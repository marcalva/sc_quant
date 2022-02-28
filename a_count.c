
#include "a_count.h"
#include "str_util.h"
#include "variants.h"

ac_node *init_ac_node(){
    ac_node *n = calloc(1, sizeof(ac_node));

    if (n == NULL){
        err_msg(-1, 0, "init_ac_node: %s", strerror(errno));
        return NULL;
    }

    ac_node_set_zero(n);
    return n;
}

ac_node *destroy_ac_node(ac_node *n){
    ac_node *next = n->next;
    free(n);
    return next;
}

void ac_node_set_zero(ac_node *n){
    n->ix = -1;
    int i;
    for (i = 0; i < N_ALLELE; ++i) n->counts[i] = 0;
    n->next = NULL;
}

bc_ac *init_bc_ac(){
    bc_ac *a = calloc(1, sizeof(bc_ac));

    if (a == NULL){
        err_msg(-1, 0, "init_bc_ac: %s", strerror(errno));
        return NULL;
    }

    a->var_ix = init_str_map();
    a->bc_ix = init_str_map();

    if (a->var_ix == NULL || a->bc_ix == NULL) return NULL;

    a->acs = kh_init(ac);
    a->n_nz = 0;
    return a;
}

void destroy_bc_ac(bc_ac *a){
    khint_t k_bc;
    for (k_bc = kh_begin(a->acs); k_bc != kh_end(a->acs); ++k_bc){
        if (!kh_exist(a->acs, k_bc)) continue;
        ac_node *n = (kh_val(a->acs, k_bc)).next;
        while (n != NULL) n = destroy_ac_node(n);
    }
    kh_destroy(ac, a->acs);
    destroy_str_map(a->var_ix);
    destroy_str_map(a->bc_ix);
    free(a);
}

int bc_ac_add_var_map(bc_ac *a, str_map *sm){
    destroy_str_map(a->var_ix);
    a->var_ix = str_map_copy(sm);
    if (a->var_ix == NULL) return -1;
    return 0;
}

int bc_ac_add_bc_map(bc_ac *a, str_map *sm){
    destroy_str_map(a->bc_ix);
    a->bc_ix = str_map_copy(sm);
    if (a->bc_ix == NULL) return -1;
    return 0;
}

int bc_ac_count(bc_ac *a, Records *recs){

    khint_t k_bc;
    for (k_bc = kh_begin(recs->bc); k_bc != kh_end(recs->bc); ++k_bc){
        if (!kh_exist(recs->bc, k_bc)) continue;
        char *bc_key = kh_key(recs->bc, k_bc);
        Barcode *bc = kh_val(recs->bc, k_bc);

        // add barcode to acs hash
        int ret;
        khint_t k_acs = kh_put(ac, a->acs, bc_key, &ret);
        if (ret == -1){
            err_msg(-1, 0, "bc_ac_count: failed to add barcode %s to hash table", bc_key);
            return -1;
        }
        else if (ret == 0){
            err_msg(-1, 0, "bc_ac_count: barcode %s found twice", bc_key);
            return -1;
        }
        ac_node *hn = &(kh_val(a->acs, k_acs));
        ac_node_set_zero(hn);

        // add counts
        khint_t k_umi;
        for (k_umi = kh_begin(bc->umis); k_umi != kh_end(bc->umis); ++k_umi){
            if (!kh_exist(bc->umis, k_umi)) continue;
            UMI *umi = kh_val(bc->umis, k_umi);
            if (umi->best_rec == NULL) continue; // look at best rec call
            Rec *rec = umi->best_rec;
            int v, n_var = rec->n_var;
            for (v = 0; v < n_var; ++v){
                int32_t vix = rec->var[v];
                char *vid = str_map_str(recs->var_ix, vix);
                vix = str_map_ix(a->var_ix, vid);
                if (vix == -1){
                    err_msg(-1, 0, "bc_ac_count: could not find variant %s. "
                            "Make sure bc_ac is initialized properly", vid);
                    return -1;
                }

                ac_node *n = hn;
                while (n->next != NULL && vix >= n->next->ix)
                    n = n->next;
                // n->ix now always <= vix

                if (n->ix < vix){ // if variant not found
                    ac_node *nn = init_ac_node();
                    nn->ix = vix;
                    nn->counts[rec->base[v]] += 1;
                    nn->next = n->next;
                    n->next = nn;
                    a->n_nz += 1;
                } 
                else if (n->ix == vix){
                    n->counts[rec->base[v]] += 1;
                }
                else{
                    return err_msg(-1, 0, "bc_ac_count: improper add to sorted list\n");
                }
            }
        }
    }

    return 0;
}
    
int bc_ac_write(bc_ac *a, char *fn){
    char delim[] = " ";
    if (mkpath(fn, 0755) == -1){
        err_msg(-1, 0, "bc_ac_write: failed to create output directory for %s", fn);
        return -1;
    }

    BGZF *fp[3];
    char *ofn[3];
    int ret;

    // ac matrix file
    char mtx_fn[3][14] = {"ac.ref.mtx.gz", "ac.alt.mtx.gz", "ac.oth.mtx.gz"};
    int i;
    for (i = 0; i < 3; ++i){
        ofn[i] = strcat2((const char*)fn, (const char*)mtx_fn[i]);
        if (ofn[i] == NULL) return -1;
    }
    for (i = 0; i < 3; ++i){
        fp[i] = bgzf_open(ofn[i], "wg1");
        if (fp[i] == NULL){
            err_msg(-1, 0, "bc_ac_write: failed to open file %s", ofn[i]);
            return -1;
        }
    }
    for (i = 0; i < 3; ++i) free(ofn[i]);

    // get strings for header
    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    strs[0] = int2str((int)a->var_ix->n, &len);
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
        return err_msg(-1, 0, "bc_ac_write: failed to write to file %s", fn);
    }

    // write counts
    int il; // intstrp string length
    size_t intstrp_m = 1; // intstrp allocated size
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(ac, a->acs, bc_key);
        if (k_bc == kh_end(a->acs)){
            bci++;
            continue;
        }
        ac_node *n = &(kh_val(a->acs, k_bc));
        for (n = n->next; n; n = n->next){
            int vix = n->ix + 1;

            for (i = 0; i < 3; ++i){
                if ((il = int2strp(vix, &intstrp, &intstrp_m)) < 0) return -1;
                ret = bgzf_write(fp[i], intstrp, il);

                ret = bgzf_write(fp[i], delim, 1);

                if ((il = int2strp(k+1, &intstrp, &intstrp_m)) < 0) return -1;
                ret = bgzf_write(fp[i], intstrp, il);

                ret = bgzf_write(fp[i], delim, 1);
            }

            if ((il = int2strp((int)n->counts[REF], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[0], intstrp, il);

            if ((il = int2strp((int)n->counts[ALT], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[1], intstrp, il);

            if ((il = int2strp((int)n->counts[OTHER], &intstrp, &intstrp_m)) < 0) return -1;
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

    free(intstrp);

    // variant file
    char var_fn[] = "ac.var.txt.gz";
    ofn[0] = strcat2((const char*)fn, (const char*)var_fn);
    if (ofn[0] == NULL) return -1;
    fp[0] = bgzf_open(ofn[0], "wg1");
    if (fp[0] == NULL){
        return err_msg(-1, 0, "bc_ac_write: failed to open file %s", ofn[0]);
    }
    free(ofn[0]);

    for (k = 0; k < a->var_ix->n; ++k){
        char *var_key = str_map_str(a->var_ix, k);
        if (var_key == NULL) continue;
        ret = bgzf_write(fp[0], var_key, strlen(var_key));
        ret = bgzf_write(fp[0], "\n", 1);
    }
    bgzf_close(fp[0]);

    // barcode file
    char bc_fn[] = "ac.barcodes.txt.gz";
    ofn[0] = strcat2((const char*)fn, (const char*)bc_fn);
    if (ofn[0] == NULL) return -1;
    fp[0] = bgzf_open(ofn[0], "wg1");
    if (fp[0] == NULL){
        bgzf_close(fp[0]);
        return err_msg(-1, 0, "bc_ac_write: failed to open file %s", ofn[0]);
    }
    free(ofn[0]);

    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) continue;
        ret = bgzf_write(fp[0], bc_key, strlen(bc_key));
        ret = bgzf_write(fp[0], "\n", 1);
    }
    bgzf_close(fp[0]);

    return 0;
}

