
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
            if (n_feat > 1) continue; // only count UMIs overlapping one gene
            for (f = 0; f < n_feat; ++f){
                int32_t fix = rec->feat[f];
                char *fid = str_map_str(recs->gene_ix, fix);
                fix = str_map_ix(a->gene_ix, fid);
                if (fix == -1){
                    return err_msg(-1, 0, "bc_gc_count: could not find gene %s. Initialize bc_gc properly\n", fid);
                }

                gc_node *n = hn;
                while (n->next != NULL && fix >= n->next->ix)
                    n = n->next;
                // n->ix now always <= fix

                if (n->ix < fix){ // if gene not found
                    gc_node *nn = init_gc_node();
                    nn->ix = fix;
                    nn->counts[rec->splice[f]] += 1;
                    nn->next = n->next;
                    n->next = nn;
                    a->n_nz += 1;
                }
                else if (n->ix == fix){
                    n->counts[rec->splice[f]] += 1;
                }
                else{
                    return err_msg(-1, 0, "bc_gc_count: improper add to sorted list\n");
                }
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

    BGZF *fp[3];
    char *ofn[3];
    int ret;

    // gc matrix file
    char mtx_fn[3][14] = {"gc.spl.mtx.gz", "gc.uns.mtx.gz", "gc.amb.mtx.gz"};
    int i;
    for (i = 0; i < 3; ++i){
        ofn[i] = strcat2((const char*)fn, (const char*)mtx_fn[i]);
        if (ofn[i] == NULL) return -1;
    }
    for (i = 0; i < 3; ++i){
        fp[i] = bgzf_open(ofn[i], "wg1");
        if (fp[i] == NULL){
            err_msg(-1, 0, "bc_gc_write: failed to open file %s", ofn[i]);
            return -1;
        }
    }


    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    strs[0] = int2str((int)a->gene_ix->n, &len);
    strs[1] = int2str((int)a->bc_ix->n, &len);
    strs[2] = int2str((int)a->n_nz, &len);

    // write header
    char mtx_hdr[] = "%%MatrixMarket matrix coordinate integer general\n";
    for (i = 0; i < 3; ++i){
        ret = bgzf_write(fp[i], mtx_hdr, strlen(mtx_hdr));
        ret = bgzf_write(fp[i], "%\n", 2);

        int j;
        for (j = 0; j < 3; ++j){
            if (j) ret = bgzf_write(fp[i], delim, 1);
            ret = bgzf_write(fp[i], strs[j], strlen(strs[j]));
        }
        ret = bgzf_write(fp[i], "\n", 1);
        if (ret < 0) break;
    }
    for (i = 0; i < 3; ++i) free(strs[i]);
    if (ret < 0){
        for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
        for (i = 0; i < 3; ++i) free(ofn[i]);
        return err_msg(-1, 0, "bc_gc_write: failed to write to file %s", fn);
    }

    // write counts
    int il; // intstrp string length
    size_t intstrp_m = 1; // intstrp allocated size
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(gc, a->gcs, bc_key);
        if (k_bc == kh_end(a->gcs)){
            bci++;
            continue;
        }

        gc_node *n = &(kh_val(a->gcs, k_bc));
        for (n = n->next; n; n = n->next){
            int fix = n->ix + 1;

            for (i = 0; i < 3; ++i){
                if ((il = int2strp(fix, &intstrp, &intstrp_m)) < 0) return -1;
                ret = bgzf_write(fp[i], intstrp, il);

                ret = bgzf_write(fp[i], delim, 1);

                if ((il = int2strp(k+1, &intstrp, &intstrp_m)) < 0) return -1;
                ret = bgzf_write(fp[i], intstrp, il);

                ret = bgzf_write(fp[i], delim, 1);
            }

            if ((il = int2strp((int)n->counts[SPLICE], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[0], intstrp, il);

            if ((il = int2strp((int)n->counts[UNSPLICE], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[1], intstrp, il);

            if ((il = int2strp((int)n->counts[AMBIG], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[2], intstrp, il);
            
            for (i = 0; i < 3; ++i)
                ret = bgzf_write(fp[i], "\n", 1);
        }
        bci++;
    }
    if (k != (bci - 1)){
        fprintf(stdout, "k=%i bci=%i\n", k, bci);
    }
    for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
    for (i = 0; i < 3; ++i) free(ofn[i]);

    // gene file
    char gene_fn[] = "gc.gene.txt.gz";
    ofn[0] = strcat2((const char*)fn, (const char*)gene_fn);
    if (ofn[0] == NULL) return -1;
    fp[0] = bgzf_open(ofn[0], "wg1");
    if (fp[0] == NULL){
        fprintf(stderr, "error: bc_gc_write: failed to open file %s\n", ofn[0]);
        return -1;
    }

    for (k = 0; k < a->gene_ix->n; ++k){
        char *gene_key = str_map_str(a->gene_ix, k);
        if (gene_key == NULL){
            fprintf(stderr, "gene index %i not found\n", k);
            continue;
        }

        Gene *gene_obj = gene_from_name(anno, gene_key);

        char *chrm = str_map_str(anno->chrm_ix, gene_obj->chrm);
        if (chrm == NULL){
            fprintf(stderr, "chromosome index %i not found\n", gene_obj->chrm);
            continue;
        }
        ret = bgzf_write(fp[0], chrm, strlen(chrm));

        if ((il = int2strp(gene_obj->beg, &intstrp, &intstrp_m)) < 0) return -1;
        ret = bgzf_write(fp[0], "\t", 1);
        ret = bgzf_write(fp[0], intstrp, il);

        if ((il = int2strp(gene_obj->end, &intstrp, &intstrp_m)) < 0) return -1;
        ret = bgzf_write(fp[0], "\t", 1);
        ret = bgzf_write(fp[0], intstrp, il);

        ret = bgzf_write(fp[0], "\t", 1);
        ret = bgzf_write(fp[0], &gene_obj->strand, 1);

        ret = bgzf_write(fp[0], "\t", 1);
        ret = bgzf_write(fp[0], gene_obj->type, strlen(gene_obj->type));

        ret = bgzf_write(fp[0], "\t", 1);
        ret = bgzf_write(fp[0], gene_obj->name, strlen(gene_obj->name));

        ret = bgzf_write(fp[0], "\t", 1);
        ret = bgzf_write(fp[0], gene_key, strlen(gene_key));

        ret = bgzf_write(fp[0], "\n", 1);
    }
    bgzf_close(fp[0]);
    free(ofn[0]);

    free(intstrp);

    // barcode file
    char bc_fn[] = "gc.barcodes.txt.gz";
    ofn[0] = strcat2((const char*)fn, (const char*)bc_fn);
    if (ofn[0] == NULL) return -1;
    fp[0] = bgzf_open(ofn[0], "wg1");
    if (fp[0] == NULL){
        fprintf(stderr, "error: bc_gc_write: failed to open file %s\n", ofn[0]);
        return -1;
    }

    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) {
            printf("barcode %i is NULL\n", k);
            continue;
        }
        ret = bgzf_write(fp[0], bc_key, strlen(bc_key));
        ret = bgzf_write(fp[0], "\n", 1);
    }
    bgzf_close(fp[0]);
    free(ofn[0]);

    return 0;
}

