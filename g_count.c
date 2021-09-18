
#include "g_count.h"
#include "str_util.h"

gc_node *init_gc_node(){
    gc_node *n = calloc(1, sizeof(gc_node));

    if (n == NULL){
        err_msg(-1, 0, "init_gc_node: %s", strerror(errno));
        return NULL;
    }

    gc_node_set_zero(n);

    return n;
}

gc_node *destroy_gc_node(gc_node *n){
    gc_node *next = n->next;
    free(n);
    return next;
}

void gc_node_set_zero(gc_node *n){
    n->ix = -1;
    int i;
    for (i = 0; i < N_SPL; ++i) n->counts[i] = 0;
    n->next = NULL;
}

bc_gc *init_bc_gc(){
    bc_gc *a = calloc(1, sizeof(bc_gc));

    if (a == NULL){
        err_msg(-1, 0, "init_bc_gc:  %s", strerror(errno));
        return NULL;
    }

    a->gene_ix = init_str_map();
    a->bc_ix = init_str_map();

    if (a->gene_ix == NULL || a->bc_ix == NULL) return NULL;

    a->gcs = kh_init(gc);
    a->n_nz = 0;
    return a;
}

void destroy_bc_gc(bc_gc *a){
    if (a == NULL) return;
    khint_t k_bc;
    for (k_bc = kh_begin(a->gcs); k_bc != kh_end(a->gcs); ++k_bc){
        if (!kh_exist(a->gcs, k_bc)) continue;
        gc_node *n = (kh_val(a->gcs, k_bc)).next;
        while (n != NULL) n = destroy_gc_node(n);
    }
    kh_destroy(gc, a->gcs);
    destroy_str_map(a->gene_ix);
    destroy_str_map(a->bc_ix);
    free(a);
}

int bc_gc_add_gene_map(bc_gc *a, str_map *sm){
    destroy_str_map(a->gene_ix);
    a->gene_ix = str_map_copy(sm);
    if (a->gene_ix == NULL) return -1;
    return 0;
}

int bc_gc_add_bc_map(bc_gc *a, str_map *sm){
    destroy_str_map(a->bc_ix);
    a->bc_ix = str_map_copy(sm);
    if (a->bc_ix == NULL) return -1;
    return 0;
}

int bc_gc_count(bc_gc *a, Records *recs){

    khint_t k_bc;
    for (k_bc = kh_begin(recs->bc); k_bc != kh_end(recs->bc); ++k_bc){
        if (!kh_exist(recs->bc, k_bc)) continue;
        char *bc_key = kh_key(recs->bc, k_bc);
        Barcode *bc = kh_val(recs->bc, k_bc);

        // add barcode to a->gcs hash table
        int ret;
        khint_t k_gcs = kh_put(gc, a->gcs, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "bc_gc_count: failed to add barcode %s to hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "bc_gc_count: barcode %s found twice\n", bc_key);
        }
        gc_node *hn = &(kh_val(a->gcs, k_gcs));
        gc_node_set_zero(hn); // dummy head node

        // add counts
        khint_t k_umi;
        for (k_umi = kh_begin(bc->umis); k_umi != kh_end(bc->umis); ++k_umi){
            if (!kh_exist(bc->umis, k_umi)) continue;
            UMI *umi = kh_val(bc->umis, k_umi);
            if (umi->best_rec == NULL) continue; // look at best rec call
            Rec *rec = umi->best_rec;
            int f, n_feat = rec->n_feat;
            for (f = 0; f < n_feat; ++f){
                int32_t fix = rec->feat[f];
                char *fid = str_map_str(recs->gene_ix, fix);
                fix = str_map_ix(a->gene_ix, fid);
                if (fix == -1){
                    fprintf(stderr, "error: bc_gc_count: could not find gene %s. Initialize bc_gc properly\n", fid);
                    return -1;
                }

                gc_node *n = hn;
                while (n->next != NULL && fix >= n->next->ix)
                    n = n->next;
                // n->ix now always <= fix

                if (n->ix < fix){ // if gene not found
                    gc_node *nn = init_gc_node();
                    nn->ix = fix;
                    nn->next = n->next;
                    n->next = nn;
                    n = nn;
                    a->n_nz += 1;
                }
                n->counts[rec->splice[f]] += 1;
            }
        }
    }

    return 0;
}
    
int bc_gc_write(bc_gc *a, Annotation *anno, char *fn){
    char delim[] = " ";
    if (mkpath(fn, 0755) == -1){
        fprintf(stderr, "error: bc_gc_write: failed to create output directory for %s", fn);
        return -1;
    }

    BGZF *fp;
    char *ofn = NULL;
    int ret;

    // gc matrix file
    char mtx_fn[] = "gc.mtx.gz";
    ofn = strcat2((const char*)fn, (const char*)mtx_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL){
        fprintf(stderr, "error: bc_gc_write: failed to open file %s\n", ofn);
        return -1;
    }
    free(ofn);

    char mtx_hdr[] = "%%MatrixMarket matrix coordinate integer general\n";
    ret = bgzf_write(fp, mtx_hdr, strlen(mtx_hdr));
    ret = bgzf_write(fp, "%\n", 2);
    int nstrs = 3;
    char *strs[nstrs];
    size_t lens[nstrs];
    int rets[nstrs];

    int i;
    for (i = 0; i < nstrs; ++i){
        lens[i] = 2;
        strs[i] = malloc(lens[i] * sizeof(char));
    }

    int nrow = (int)a->gene_ix->n * N_SPL;
    int ncol = (int)a->bc_ix->n;
    int n_nz = (int)a->n_nz * N_SPL;

    rets[0] = int2strp(nrow, strs + 0, lens + 0);
    rets[1] = int2strp(ncol, strs + 1, lens + 1);
    rets[2] = int2strp(n_nz, strs + 2, lens + 2);

    for (i = 0; i < 3; ++i){
        if (i) ret = bgzf_write(fp, delim, 1);
        ret = bgzf_write(fp, strs[i], rets[i]);
    }
    ret = bgzf_write(fp, "\n", 1);
    // for (i = 0; i < 3; ++i) free(strs[i]);
    if (ret < 0){
        fprintf(stderr, "error: bc_gc_write: failed to write to file %s\n", fn);
        return -1;
    }

    // write counts
    int k, s;
    for (s = SPLICE; s <= AMBIG; ++s){
        for (k = 0; k < a->bc_ix->n; ++k){
            char *bc_key = str_map_str(a->bc_ix, k);
            if (bc_key == NULL) continue;

            khint_t k_bc = kh_get(gc, a->gcs, bc_key);
            if (k_bc == kh_end(a->gcs)){
                // a barcode can be present but have no counts
                continue;
                // fprintf(stderr, "error: bc_gc_write: could not find barcode %s. Initialize bc_gc properly\n", bc_key);
                // return -1;
            }
            gc_node *n = &(kh_val(a->gcs, k_bc));
            for (n = n->next; n; n = n->next){
                int fix = n->ix + 1 + (s * a->gene_ix->n);

                //get total counts
                int i, total = 0;
                for (i = 0; i < N_SPL; ++i) total += (int)n->counts[i];

                rets[0] = int2strp(fix       , strs + 0, lens + 0);
                rets[1] = int2strp(k+1       , strs + 1, lens + 1);
                rets[2] = int2strp((int)n->counts[s], strs + 2, lens + 2);
                // rets[2] = int2strp(total     , strs + 2, lens + 2);
                // rets[3] = int2strp((int)n->counts[SPLICE]  , strs + 3, lens + 3);
                // rets[4] = int2strp((int)n->counts[UNSPLICE], strs + 4, lens + 4);
                // rets[5] = int2strp((int)n->counts[AMBIG]   , strs + 5, lens + 5);
                for (i = 0; i < nstrs; ++i){
                    if (i) ret = bgzf_write(fp, delim, 1);
                    ret = bgzf_write(fp, strs[i], rets[i]);
                }
                ret = bgzf_write(fp, "\n", 1);
            }
        }
    }

    for (i = 0; i < nstrs; ++i) free(strs[i]);
    bgzf_close(fp);

    // gene file
    char gene_fn[] = "gc.gene.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)gene_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL){
        fprintf(stderr, "error: bc_gc_write: failed to open file %s\n", ofn);
        return -1;
    }
    free(ofn);

    size_t intstrp_len = 1;
    char *intstrp = malloc(sizeof(char) * intstrp_len);

    char *spl_types[] = {"SPLICE", "UNSPLICE", "AMBIG"};
    for (s = SPLICE; s <= AMBIG; ++s){
        for (k = 0; k < a->gene_ix->n; ++k){
            char *gene_key = str_map_str(a->gene_ix, k);
            Gene *gene_obj = gene_from_name(anno, gene_key);

            if (gene_key == NULL) continue;

            char *chrm = str_map_str(anno->chrm_ix, gene_obj->chrm);
            ret = bgzf_write(fp, chrm, strlen(chrm));

            if (int2strp(gene_obj->beg, &intstrp, &intstrp_len) < 0){
                fprintf(stderr, "error: bc_gc_write: failed to convert int to str\n");
                return -1;
            }
            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, intstrp, strlen(intstrp));
            
            int2strp(gene_obj->end, &intstrp, &intstrp_len);
            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, intstrp, strlen(intstrp));

            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, &gene_obj->strand, 1);

            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, gene_obj->type, strlen(gene_obj->type));

            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, gene_obj->name, strlen(gene_obj->name));

            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, gene_key, strlen(gene_key));

            ret = bgzf_write(fp, "\t", 1);
            ret = bgzf_write(fp, spl_types[s], strlen(spl_types[s]));
            ret = bgzf_write(fp, "\n", 1);
        }
    }
    free(intstrp);
    bgzf_close(fp);

    // barcode file
    char bc_fn[] = "gc.barcodes.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)bc_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL){
        fprintf(stderr, "error: bc_gc_write: failed to open file %s\n", ofn);
        return -1;
    }
    free(ofn);

    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) {
            printf("barcode %i is NULL\n", k);
            continue;
        }
        ret = bgzf_write(fp, bc_key, strlen(bc_key));
        ret = bgzf_write(fp, "\n", 1);
    }
    bgzf_close(fp);

    return 0;
}

