
#include "a_count.h"
#include "str_util.h"
#include "variants.h"

ac_node *init_ac_node(){
    ac_node *n = calloc(1, sizeof(ac_node));
    n->ix = -1;
    int i;
    for (i = 0; i < N_ALLELE; ++i) n->counts[i] = 0;
    n->next = NULL;
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
    a->var_ix = init_str_map();
    a->bc_ix = init_str_map();
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
                vid = str_map_str(a->var_ix, vix);

                ac_node *n = hn;
                while (n->next != NULL && vix >= n->next->ix)
                    n = n->next;
                // n->ix now always <= vix

                if (n->ix < vix){ // if variant not found
                    ac_node *nn = init_ac_node();
                    nn->ix = vix;
                    nn->next = n->next;
                    n->next = nn;
                    n = nn;
                    a->n_nz += 1;
                }
                n->counts[rec->base[v]] += 1;
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

    BGZF *fp;
    char *ofn = NULL;
    int ret;

    // ac matrix file
    char mtx_fn[] = "ac.mtx.gz";
    ofn = strcat2((const char*)fn, (const char*)mtx_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL){
        err_msg(-1, 0, "bc_ac_write: failed to open file %s", ofn);
        return -1;
    }
    free(ofn);

    char mtx_hdr[] = "%%MatrixMarket matrix coordinate integer general\n";
    ret = bgzf_write(fp, mtx_hdr, strlen(mtx_hdr));
    ret = bgzf_write(fp, "%\n", 2);
    size_t len;
    char *strs[5];

    // var bc non-zero lengths
    strs[0] = int2str((int)a->var_ix->n, &len);
    strs[1] = int2str((int)a->bc_ix->n, &len);
    strs[2] = int2str((int)a->n_nz, &len);

    int i;
    for (i = 0; i < 3; ++i){
        if (i) ret = bgzf_write(fp, delim, 1);
        ret = bgzf_write(fp, strs[i], strlen(strs[i]));
    }
    ret = bgzf_write(fp, "\n", 1);
    for (i = 0; i < 3; ++i) free(strs[i]);
    if (ret < 0){
        err_msg(-1, 0, "bc_ac_write: failed to write to file %s", ofn);
        return -1;
    }

    // write counts
    int k, bci = 1;
    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(ac, a->acs, bc_key);
        if (k_bc == kh_end(a->acs)){
            err_msg(-1, 0, "bc_ac_write: could not find barcode %s. "
                    "Make sure bc_ac is initialized properly", bc_key);
            return -1;
        }
        ac_node *n = &(kh_val(a->acs, k_bc));
        for (n = n->next; n; n = n->next){
            int vix = n->ix + 1;
            strs[0] = int2str(vix, &len);
            strs[1] = int2str(bci, &len);
            strs[2] = int2str((int)n->counts[REF], &len);
            strs[3] = int2str((int)n->counts[ALT], &len);
            strs[4] = int2str((int)n->counts[OTHER], &len);
            for (i = 0; i < 5; ++i){
                if (i) ret = bgzf_write(fp, delim, 1);
                ret = bgzf_write(fp, strs[i], strlen(strs[i]));
            }
            ret = bgzf_write(fp, "\n", 1);
            for (i = 0; i < 5; ++i) free(strs[i]);
        }
        bci++;
    }
    bgzf_close(fp);

    // variant file
    char var_fn[] = "ac.var.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)var_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL){
        err_msg(-1, 0, "bc_ac_write: failed to open file %s", ofn);
        return -1;
    }
    free(ofn);

    for (k = 0; k < a->var_ix->n; ++k){
        char *var_key = str_map_str(a->var_ix, k);
        if (var_key == NULL) continue;
        ret = bgzf_write(fp, var_key, strlen(var_key));
        ret = bgzf_write(fp, "\n", 1);
    }
    bgzf_close(fp);

    // barcode file
    char bc_fn[] = "ac.barcodes.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)bc_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL){
        err_msg(-1, 0, "bc_ac_write: failed to open file %s", ofn);
        return -1;
    }
    free(ofn);

    for (k = 0; k < a->bc_ix->n; ++k){
        char *bc_key = str_map_str(a->bc_ix, k);
        if (bc_key == NULL) continue;
        ret = bgzf_write(fp, bc_key, strlen(bc_key));
        ret = bgzf_write(fp, "\n", 1);
    }
    bgzf_close(fp);

    return 0;
}

