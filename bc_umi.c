
#include "bc_umi.h"
#include "htslib/hts_endian.h"
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <inttypes.h>
#include <errno.h>

Rec *init_rec(){
    Rec *r = (Rec *)calloc(1, sizeof(Rec));

    if (r == NULL){
        err_msg(-1, 0, "init_rec: %s", strerror(errno));
        return NULL;
    }

    r->n_feat = 0;
    r->m_feat = 1;
    r->feat = (int32_t *)malloc(r->m_feat * sizeof(int32_t));
    r->splice = (uint8_t *)malloc(r->m_feat * sizeof(uint8_t));
    r->n_var = 0;
    r->m_var = 1;
    r->var = (int32_t *)malloc(r->m_var * sizeof(int32_t));
    r->base = (uint8_t *)malloc(r->m_var * sizeof(uint8_t));
    r->qual = (uint8_t *)malloc(r->m_var * sizeof(uint8_t));

    if (r->qual == NULL){
        err_msg(-1, 0, "init_rec: %s", strerror(errno));
        return NULL;
    }

    r->n_read = 0;
    r->next = NULL;
    return(r);
}

Rec *init_rec_s(uint8_t n_feat, uint8_t n_var){
    Rec *r = (Rec *)calloc(1, sizeof(Rec));

    if (r == NULL){
        err_msg(-1, 0, "init_rec_s: %s", strerror(errno));
        return NULL;
    }

    r->n_feat = n_feat;
    r->m_feat = n_feat;
    r->feat = (int32_t *)malloc(r->m_feat * sizeof(int32_t));
    r->splice = (uint8_t *)malloc(r->m_feat * sizeof(uint8_t));
    r->n_var = n_var;
    r->m_var = n_var;
    r->var = (int32_t *)malloc(r->m_var * sizeof(int32_t));
    r->base = (uint8_t *)malloc(r->m_var * sizeof(uint8_t));
    r->qual = (uint8_t *)malloc(r->m_var * sizeof(uint8_t));

    if (r->qual == NULL){
        err_msg(-1, 0, "init_rec_s: %s", strerror(errno));
        return NULL;
    }

    r->n_read = 0;
    r->next = NULL;
    return(r);
}

Rec *destroy_rec(Rec *r){

    if (r == NULL) return NULL;

    free(r->feat);
    free(r->splice);
    free(r->var);
    free(r->base);
    free(r->qual);
    Rec *next = r->next;
    free(r);
    return(next);
}

// also copies next member
Rec *copy_rec(Rec *r){
    if (r == NULL) return NULL;
    Rec *r_cpy = (Rec *)calloc(1, sizeof(Rec));

    if (r_cpy == NULL){
        err_msg(-1, 0, "copy_rec: %s", strerror(errno));
        return NULL;
    }

    r_cpy->n_feat = r->n_feat;
    r_cpy->m_feat = r_cpy->n_feat;
    r_cpy->feat = (int32_t *)malloc(r_cpy->m_feat * sizeof(int32_t));
    r_cpy->splice = (uint8_t *)malloc(r_cpy->m_feat * sizeof(uint8_t));

    if (r_cpy->splice == NULL){
        err_msg(-1, 0, "copy_rec: %s", strerror(errno));
        return NULL;
    }

    int i;
    for (i = 0; i < r_cpy->n_feat; ++i){
        r_cpy->feat[i] = r->feat[i];
        r_cpy->splice[i] = r->splice[i];
    }
    r_cpy->n_var = r->n_var;
    r_cpy->m_var = r_cpy->n_var;
    r_cpy->var = (int32_t *)malloc(r_cpy->m_var * sizeof(int32_t));
    r_cpy->base = (uint8_t *)malloc(r_cpy->m_var * sizeof(uint8_t));
    r_cpy->qual = (uint8_t *)malloc(r_cpy->m_var * sizeof(uint8_t));
    
    if (r_cpy->qual == NULL){
        err_msg(-1, 0, "copy_rec: %s", strerror(errno));
        return NULL;
    }

    for (i = 0; i < r_cpy->n_var; ++i){
        r_cpy->var[i] = r->var[i];
        r_cpy->base[i] = r->base[i];
        r_cpy->qual[i] = r->qual[i];
    }

    r_cpy->n_read = r->n_read;
    r_cpy->next = r->next;
    return r_cpy;
}

UMI *init_umi(){
    UMI *u = (UMI *)calloc(1, sizeof(UMI));

    if (u == NULL){
        err_msg(-1, 0, "init_umi: %s", strerror(errno));
        return NULL;
    }

    u->recs = NULL;
    u->best_rec = NULL;
    return(u);
}

void destroy_umi(UMI *u){

    if (u == NULL) return;

    Rec *r = u->recs;
    while (r) r = destroy_rec(r);
    if (u->best_rec) destroy_rec(u->best_rec);
    free(u);
}

Barcode *init_barcode(){
    Barcode *b = (Barcode *)calloc(1, sizeof(Barcode));

    if (b == NULL){
        err_msg(-1, 0, "init_barcode: %s", strerror(errno));
        return NULL;
    }

    b->umis = kh_init(umi);
    b->umis_a = NULL;
    return(b);
}

void destroy_barcode(Barcode *b){

    if (b == NULL) return;

    khint_t k;
    for (k = kh_begin(b->umis); k != kh_end(b->umis); ++k){
        if (!kh_exist(b->umis, k)) continue;
        UMI *u = kh_val(b->umis, k);
        if (u) destroy_umi(u);
    }
    kh_destroy(umi, b->umis);
    if (b->umis_a != NULL) free(b->umis_a);
    free(b);
}

Records *init_records(){
    Records *recs = (Records *)calloc(1, sizeof(Records));

    if (recs == NULL){
        err_msg(-1, 0, "init_records: %s", strerror(errno));
        return NULL;
    }

    recs->bc = kh_init(bc);
    if ((recs->bc_ix = init_str_map()) == NULL) return NULL;
    if ((recs->umi_ix = init_str_map()) == NULL) return NULL;
    if ((recs->gene_ix = init_str_map()) == NULL) return NULL;
    if ((recs->var_ix = init_str_map()) == NULL) return NULL;

    return(recs);
}

void destroy_records(Records *recs){

    if (recs == NULL) return;

    khint_t k;

    // barcode objects
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *p = kh_val(recs->bc, k);
        if (p) destroy_barcode(p);
    }
    kh_destroy(bc, recs->bc);

    destroy_str_map(recs->bc_ix);
    destroy_str_map(recs->umi_ix);
    destroy_str_map(recs->gene_ix);
    destroy_str_map(recs->var_ix);

    free(recs);
}

int del_var(Rec *rec, int i){
    if (i < 0 || i >= rec->n_var){
        return err_msg(-1, 0, "del_var: variant index %i must be between 0 and %hu", 
                i, rec->n_var);
    }

    while (i < (rec->n_var - 1)){
        rec->var[i] = rec->var[i+1];
        rec->base[i] = rec->base[i+1];
        rec->qual[i] = rec->qual[i+1];
        ++i;
    }

    rec->n_var -= 1;

    return 0;
}

int add_bc2ix_file(Records *recs, const char *path){
    const char *modes[3] = {"ru\0", "r1\0", "r1\0"};
    BGZF* fp =  bgzf_open(path, "ru");
    if (fp == 0)
        return err_msg(-1, 0, "add_bc2ix_file: failed to open file %s", path);

    int c = bgzf_compression(fp);
    bgzf_close(fp);
    fp = bgzf_open(path, modes[c]);
    if (fp == 0)
        return err_msg(-1, 0, "add_bc2ix_file: failed to open file %s", path);

    kstring_t line = KS_INITIALIZE;
    while ( bgzf_getline(fp, '\n', &line) > 0 ){
        int found = 0;
        if (add2str_map(recs->bc_ix, line.s, &found) < 0) return -1;
        ks_free(&line);
    }
    bgzf_close(fp);
    return(recs->bc_ix->n);
}

int add_genes2recs(Records *recs, Annotation *a){
    int i, found, G = a->gene_ix->n;
    for (i = 0; i < G; ++i){
        char *gene_id = str_map_str(a->gene_ix, i);
        if (add2str_map(recs->gene_ix, gene_id, &found) < 0) return -1;
    }
    return 0;
}

int recs_cmp(Rec *r1, Rec *r2){
    if (r1 == NULL || r2 == NULL) return -1;
    if (r1->n_feat != r2->n_feat) return -1;

    int i;
    for (i = 0; i < r1->n_feat; ++i){
        int j;
        for (j = 0; j < r2->n_feat; ++j){
            int eq = (r1->feat[i] == r2->feat[j]) && (r1->splice[i] == r2->splice[j]);
            if (eq) break;
        }
        if (j == r2->n_feat) return -1;
    }

    if (r1->n_var != r2->n_var) return -1;

    for (i = 0; i < r1->n_var; ++i){
        int j;
        for (j = 0; j < r2->n_var; ++j){
            int eq = (r1->var[i] == r2->var[j]) && (r1->base[i] == r2->base[j]);
            if (eq) break;
        }
        if (j == r2->n_var) return -1;
    }

    return 0;
}

// return 1 if r found in list, 0 if not (and added to list).
int add_rec2list(Rec **list, Rec *r){
    int i, found = 0;
    if (*list == NULL){
        r->n_read = 1;
        *list = r;
    }
    else {
        Rec *last = NULL;
        Rec *trec = *list;
        while (trec){
            if (recs_cmp(trec, r) == 0){
                // update to highest quality score
                for (i = 0; i < r->n_var; ++i){
                    if (r->qual[i] > trec->qual[i]) trec->qual[i] = r->qual[i];
                }
                (trec->n_read)++;
                found = 1;
                break;
            }
            if (trec->next == NULL) last = trec;
            trec = trec->next;
        }
        if (last) last->next = r;
    }
    return found;
}

int add_rec(Records *recs, const char *bc, const char *umi, 
        uint8_t n_feat, const char **feat, const uint8_t *splice, 
        uint8_t n_var, const char **var, const uint8_t *base, const uint8_t *qual){ 
    int ret, found = 0;

    int bc_i = add2str_map(recs->bc_ix, bc, &found);
    int umi_i = add2str_map(recs->umi_ix, umi, &found);
    if (bc_i < 0 || umi_i < 0) return -1;

    khint_t k;

    Barcode *tbc = NULL;
    k = kh_get(bc, recs->bc, (char *)bc);
    if ( k == kh_end(recs->bc) ) {
        tbc = init_barcode();
        
        if (tbc == NULL) return -1;

        char *bc_cpy = str_map_str(recs->bc_ix, bc_i);

        k = kh_put(bc, recs->bc, bc_cpy, &ret);
        if (ret == -1)
            return err_msg(-1, 0, "add_rec: failed to add %s to hash table", bc_cpy);

        kh_val(recs->bc, k) = tbc;
    } else {
        tbc = kh_val(recs->bc, k);
    }
    
    UMI *tumi;
    k = kh_get(umi, tbc->umis, (char *)umi);
    if ( k == kh_end(tbc->umis) ) {
        tumi = init_umi();

        if (tumi == NULL) return -1;

        char *umi_cpy = str_map_str(recs->umi_ix, umi_i);

        k = kh_put(umi, tbc->umis, umi_cpy, &ret);
        if (ret == -1)
            return err_msg(-1, 0, "add_rec: failed to add %s to hash table", umi_cpy);

        kh_val(tbc->umis, k) = tumi;
    } else {
        tumi = kh_val(tbc->umis, k);
    }

    /* populate Rec object */
    Rec *r = init_rec_s(n_feat, n_var);
    int i;
    for (i = 0; i < n_feat; ++i){
        int32_t fi = (int32_t)add2str_map(recs->gene_ix, feat[i], &found);
        if (fi < 0) {
            destroy_rec(r);
            return -1;
        }
        r->feat[i] = fi;
        r->splice[i] = splice[i];
    }
    for (i = 0; i < n_var; ++i){
        int32_t vi = (int32_t)add2str_map(recs->var_ix, var[i], &found);
        if (vi < 0){
            destroy_rec(r);
            return -1;
        }
        r->var[i] = vi;
        r->base[i] = base[i];
        r->qual[i] = qual[i];
    }
    /* */

    found = add_rec2list(&(tumi->recs), r);
    if (found) destroy_rec(r);
    return found;
}

char **bc_umi_cov(Records *recs, int u1, int u2, int *n_bc){
    int bci;
    *n_bc = 0;
    char **bcs = malloc(recs->bc_ix->n * sizeof(char *));

    if (bcs == NULL){
        err_msg(-1, 0, "bc_umi_cov: %s", strerror(errno));
        return NULL;
    }

    for (bci = 0; bci < recs->bc_ix->n; ++bci){
        char *bc_key = str_map_str(recs->bc_ix, bci);
        if (bc_key == NULL) continue;
        khint_t k_bc = kh_get(bc, recs->bc, bc_key);
        if (k_bc == kh_end(recs->bc)) return NULL;
        Barcode *bc = kh_val(recs->bc, k_bc);

        int n_umi = 0;
        khint_t k_umi;
        for (k_umi = kh_begin(bc->umis); k_umi != kh_end(bc->umis); ++k_umi){
            if (!kh_exist(bc->umis, k_umi)) continue;
            UMI *umi = kh_val(bc->umis, k_umi);
            if (umi->best_rec == NULL) continue;
            n_umi++;
        }
        if (n_umi >= u1 && n_umi <= u2) bcs[(*n_bc)++] = bc_key;
    }
    bcs = (char **)realloc(bcs, (*n_bc)*sizeof(char *));

    if (bcs == NULL){
        err_msg(-1, 0, "bc_umi_cov: %s", strerror(errno));
        return NULL;
    }

    return bcs;
}


/* */
/* Records filtering */
/* */

/* Call UMIs.
 * Fill in gene and variant data from reads.
 * Read filtering algorithm:
 *   If a read contains multiple records, pick the one with the best alignment score.
 *   If there are more than one records with the same highest score, remove the ambiguous read.
 * UMI filtering algorithm
 *   Choose the record with the highest number of supporting reads, 
 *     gene-splice-variant-base configuration. If ties, remove the UMI.
 *
 * Each UMI has a best record associated with the object.
 * If a UMI doesn't have a best record, it is removed.
 * If a barcode ends up with 0 UMIs, it is removed
 */
int call_umis(Records *recs){

    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        char *bc_key = kh_key(recs->bc, k);
        khint_t l;
        int n_umi = 0;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            if (u->recs == NULL) continue;
            u->best_rec = NULL;

            /* Get most prevalent read for UMI */
            // Rec *rec = umi_recs;
            Rec *rec = u->recs;
            Rec *mf_rec = NULL;
            int most_n = -1;
            int most_n_n = 0;
            while (rec){
                if (rec->n_read > most_n){
                    mf_rec = rec;
                    most_n = rec->n_read;
                    most_n_n = 1;
                }
                else if (rec->n_read == most_n){
                    most_n_n++;
                }
                rec = rec->next;
            }
            if (most_n_n == 1){
                u->best_rec = copy_rec(mf_rec);
                u->best_rec->next = NULL;
                n_umi++;
            }

            // remove if no variants or features
            if (u->best_rec && u->best_rec->n_feat == 0 && u->best_rec->n_var == 0)
                u->best_rec = destroy_rec(u->best_rec);

            rec = u->recs;
            while (rec) rec = destroy_rec(rec);
            u->recs = NULL;
            // remove UMI if no best_rec
            if (u->best_rec == NULL) {
                destroy_umi(u);
                kh_del(umi, bc->umis, l);
            }
        }
        // remove if barcode has 0 UMIs
        if (n_umi == 0){
            str_map_del(recs->bc_ix, bc_key);
            destroy_barcode(bc);
            kh_del(bc, recs->bc, k);
        }
    }
    return 0;
}

int fltr_bq(Records *recs, uint8_t min_qual){
    fprintf(stdout, "filtering variants base quality\n");
    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        // char *bc_key = kh_key(recs->bc, k);
        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            Rec *rec = u->best_rec;
            if (rec == NULL) continue;

            int i;
            for (i = 0; i < rec->n_var; ++i){
                if (rec->qual[i] < min_qual){
                    if (del_var(rec, i) < 0) return -1;
                    --i;
                }
            }

            // remove if no variants or features
            if (u->best_rec && u->best_rec->n_feat == 0 && u->best_rec->n_var == 0){
                u->best_rec = destroy_rec(u->best_rec);
                destroy_umi(u);
                kh_del(umi, bc->umis, l);
            }
        }
    }
    return 0;
}

int subset_vars_rec(Records *recs, Rec *rec, str_map *vars){
    int i;
    for (i = 0; i < rec->n_var; ++i){
        int vix = rec->var[i];
        char *vid = str_map_str(recs->var_ix, vix);
        int mapix = str_map_ix(vars, vid);
        if (mapix == -1){
            if (del_var(rec, i) < 0) return -1;
            --i;
        } else {
            rec->var[i] = mapix; // switch to new map
        }
    }
    return 0;
}

int subset_vars_recs(Records *recs, str_map *vars){
    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);

            Rec *rec = u->recs;
            while (rec != NULL){
                if (subset_vars_rec(recs, rec, vars) < 0) return -1;
                rec = rec->next;
            }

            rec = u->best_rec;
            if (rec != NULL){
                if (subset_vars_rec(recs, rec, vars) < 0) return -1;
                rec = rec->next;
            }
        }
    }

    // switch to new map
    destroy_str_map(recs->var_ix);
    recs->var_ix = str_map_copy(vars);
    return 0;
}

int fill_umi_array(Records *recs){
    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        int n_umi = kh_size(bc->umis), umi_ix = 0;
        if (n_umi == 0) continue;
        bc->umis_a = (UMI **)calloc(n_umi, sizeof(UMI *));

        if (bc->umis_a == NULL)
            return err_msg(-1, 0, "fill_umi_array: %s", strerror(errno));

        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            bc->umis_a[umi_ix++] = u;
        }
    }
    return 0;
}

/* */
/* Writing and reading to file */
/* */

/* Storage format for records
 * 
 * Records are stored as little-endian bytes in the following order
 *
 * block size uint32_t
 * barcode index uint32_t
 * umi index uint32_t
 * n_feat uint8_t
 * feat  int32_t * n_feat
 * splice uint8_t * n_feat
 * n_var uint8_t
 * var int32_t * n_var
 * base uint8_t * n_var
 * qual uint8_t * n_var
 * n_read uint16_t
 *
 * the block size contains the length of the record (not including the block size part)
 *
 * the file begins with char[5] "recs\0", and continues with records
 */

int write_rec_bn(Rec *rec, uint32_t bc_ix, uint32_t umi_ix, FILE *fp){
    if (rec == NULL) return 0;

    uint32_t bs = 0;
    bs += 2 * sizeof(uint32_t); // bc, umi indicies
    bs += sizeof(uint8_t); // n_feat
    bs += rec->n_feat * sizeof(int32_t); // feat
    bs += rec->n_feat * sizeof(uint8_t); // splice
    bs += sizeof(uint8_t); // n_var
    bs += rec->n_var * sizeof(int32_t); // var
    bs += rec->n_var * sizeof(uint8_t); // base
    bs += rec->n_var * sizeof(uint8_t); // qual
    bs += sizeof(uint16_t); // n_read

    uint8_t buf[8];
    int i;

    // block size
    u32_to_le(bs, buf);
    fwrite(buf, 1, 4, fp);
    // fprintf(stdout, "%u\n", *((unsigned int *)buf));

    // bc_ix
    u32_to_le(bc_ix, buf);
    fwrite(buf, 1, 4, fp);
    // fprintf(stdout, "%u\n", *((unsigned int *)buf));

    // umi_ix
    u32_to_le(umi_ix, buf);
    fwrite(buf, 1, 4, fp);

    // n_feat
    buf[0] = rec->n_feat;
    fwrite(buf, 1, 1, fp);

    // feat
    for (i = 0; i < rec->n_feat; ++i){
        i32_to_le(rec->feat[i], buf);
        fwrite(buf, 1, 4, fp);
    }

    // splice
    for (i = 0; i < rec->n_feat; ++i){
        buf[0] = rec->splice[i];
        fwrite(buf, 1, 1, fp);
    }

    // n_var
    buf[0] = rec->n_var;
    fwrite(buf, 1, 1, fp);
    // fprintf(stdout, "%u\n", buf[0]);

    // var
    for (i = 0; i < rec->n_var; ++i){
        i32_to_le(rec->var[i], buf);
        fwrite(buf, 1, 4, fp);
        // fprintf(stdout, "\t%i", rec->var[i]);
    }
    // fprintf(stdout, "\n");

    // base
    for (i = 0; i < rec->n_var; ++i){
        buf[0] = rec->base[i];
        fwrite(buf, 1, 1, fp);
    }

    // qual
    for (i = 0; i < rec->n_var; ++i){
        buf[0] = rec->qual[i];
        fwrite(buf, 1, 1, fp);
    }

    // n_read
    u16_to_le(rec->n_read, buf);
    int len = fwrite(buf, 1, 2, fp);
    if (len != 2) // check write
        return err_msg(-1, 0, "write_rec_bn: %s", strerror(errno));

    return bs;
}

void read_rec_bn(uint8_t *buf, Rec *rec, uint32_t *bc_ix, uint32_t *umi_ix){
    int i;

    // bc_ix
    *bc_ix = le_to_u32(buf);
    buf = buf + 4;

    // umi_ix
    *umi_ix = le_to_u32(buf);
    buf = buf + 4;

    // n_feat
    rec->n_feat = *buf;
    buf++;

    // m_feat
    if (rec->n_feat > rec->m_feat){
        rec->m_feat = rec->n_feat;
        rec->feat = (int32_t *)realloc(rec->feat, rec->m_feat * sizeof(int32_t));
    }

    // feat
    for (i = 0; i < rec->n_feat; ++i){
        (rec->feat)[i] = le_to_i32(buf);
        buf = buf + 4;
    }

    // splice
    rec->splice = (uint8_t *)realloc(rec->splice, rec->m_feat * sizeof(uint8_t));
    for (i = 0; i < rec->n_feat; ++i){
        (rec->splice)[i] = *buf;
        buf++;
    }

    // n_var
    rec->n_var = *buf;
    buf++;

    // m_var
    if (rec->n_var > rec->m_var){
        rec->m_var = rec->n_var;
        rec->var = (int32_t *)realloc(rec->var, rec->m_var * sizeof(int32_t));
    }   

    // var
    for (i = 0; i < rec->n_var; ++i){
        (rec->var)[i] = le_to_i32(buf);
        buf = buf + 4;
    }

    // base
    rec->base = (uint8_t *)realloc(rec->base, rec->m_var * sizeof(uint8_t));
    for (i = 0; i < rec->n_var; ++i){
        (rec->base)[i] = *buf;
        buf++;
    }

    // qual
    rec->qual = (uint8_t *)realloc(rec->qual, rec->m_var * sizeof(uint8_t));
    for (i = 0; i < rec->n_var; ++i){
        (rec->qual)[i] = *buf;
        buf++;
    }

    // n_read
    rec->n_read = le_to_u16(buf);
}

int write_ix_aux(const char *fn1, const char *fn2, str_map *sm){
    char *ofn = strcat2((const char*)fn1, (const char*)fn2);
    if (ofn == NULL) return -1;

    int ret = 0;
    ret = write_str_map(sm, ofn, '\t', '\n');
    free(ofn);
    if (ret < 0) return -1;
    else return 0;
    return 0;
}

int write_records_bn(Records *recs, const char *fn){

    // create directory
    char *fn_cpy = calloc(strlen(fn) + 1, sizeof(char));
    strcpy(fn_cpy, fn);
    if (mkpath(fn_cpy, 0755) == -1){
        free(fn_cpy);
        return err_msg(-1, 0, "write_records_bn: failed to create output directory for %s", fn);
    }
    free(fn_cpy);

    int ret = 0;
    char recsfn[] = "recs.bin";
    char bcixfn[] = "bc_ix.tsv.gz";
    char umiixfn[] = "umi_ix.tsv.gz";
    char varixfn[] = "var_ix.tsv.gz";
    char geneixfn[] = "gene_ix.tsv.gz";

    char *outrecsfn = strcat2((const char*)fn, (const char*)recsfn);
    if (outrecsfn == NULL) return -1;

    /* barcode index */
    ret = write_ix_aux(fn, bcixfn, recs->bc_ix);
    if (ret == -1) goto cleanup;
    /* */

    /* UMI ix */
    ret = write_ix_aux(fn, umiixfn, recs->umi_ix);
    if (ret == -1) goto cleanup;
    /* */

    /* var ix */
    ret = write_ix_aux(fn, varixfn, recs->var_ix);
    if (ret == -1) goto cleanup;
    /* */

    /* gene ix */
    ret = write_ix_aux(fn, geneixfn, recs->gene_ix);
    if (ret == -1) goto cleanup;
    /* */


    FILE *fp = fopen(outrecsfn, "wb");
    if (fp == NULL){
        err_msg(-1, 0, "write_records_bn: failed to open output file %s: %s", fn, strerror(errno));
        ret = -1;
        goto cleanup;
    }
    
    char beg[5] = "recs\0";
    int len = fwrite((uint8_t *)beg, 1, 5, fp);
    if (len < 5){
        err_msg(-1, 0, "write_records_bn: failed to write to output file %s: %s", fn, strerror(errno));
        ret = -1;
        goto cleanup;
    }

    uint32_t bc_ix, umi_ix;
    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        char *bc_key = kh_key(recs->bc, k);
        bc_ix = (uint32_t)str_map_ix(recs->bc_ix, bc_key);
        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            char *u_key = kh_key(bc->umis, l);
            umi_ix = (uint32_t)str_map_ix(recs->umi_ix, u_key);
            Rec *rec = u->recs;
            while (rec) {
                int wret = write_rec_bn(rec, bc_ix, umi_ix, fp);
                if (wret == -1){
                    ret = -1;
                    goto cleanup;
                }
                rec = rec->next;
            }
        }
    }

    fclose(fp);

cleanup:
    free(outrecsfn);
    return ret;
}

// TODO implement safe alloc and realloc w error checking
int read_records_bn(const char *fn, Records *recs){

    if (recs == NULL) return -1;

    int ret = 0, i;
    
    char recsfn[] = "recs.bin";
    char bcixfn[] = "bc_ix.tsv.gz";
    char umiixfn[] = "umi_ix.tsv.gz";
    char varixfn[] = "var_ix.tsv.gz";
    char geneixfn[] = "gene_ix.tsv.gz";

    char *outrecsfn = strcat2((const char*)fn, (const char*)recsfn);
    char *outbcixfn = strcat2((const char*)fn, (const char*)bcixfn);
    char *outumiixfn = strcat2((const char*)fn, (const char*)umiixfn);
    char *outvarixfn = strcat2((const char*)fn, (const char*)varixfn);
    char *outgeneixfn = strcat2((const char*)fn, (const char*)geneixfn);

    if (outgeneixfn == NULL) return -1;

    FILE *fp = NULL;
    uint8_t *buf = NULL;
    uint8_t feat_m = 0;
    char **feat = NULL;
    uint8_t var_m = 0;
    char **var = NULL;
    str_map *bc_map = NULL, *umi_map = NULL, *gene_map = NULL, *var_map = NULL;

    Rec *r = init_rec();

    if (outrecsfn == NULL || outbcixfn == NULL || outumiixfn == NULL || 
            outvarixfn == NULL || outgeneixfn == NULL){
        ret = -1;
        goto cleanup;
    }

    /* get arrays for BC, UMI, and read */
    bc_map = read_str_map((const char *)outbcixfn);
    umi_map = read_str_map((const char *)outumiixfn);
    gene_map = read_str_map((const char *)outgeneixfn);
    var_map = read_str_map((const char *)outvarixfn);
    if ( bc_map == NULL || umi_map == NULL || gene_map == NULL || var_map == NULL ){
        ret = -1;
        goto cleanup;
    }
    /* */
    
    // add barcodes
    if (add_from_str_map(recs->bc_ix, bc_map) < 0){
        ret = -1;
        goto cleanup;
    }

    // add UMIs
    if (add_from_str_map(recs->umi_ix, umi_map) < 0){
        ret = -1;
        goto cleanup;
    }

    // add genes
    if (add_from_str_map(recs->gene_ix, gene_map) < 0){
        ret = -1;
        goto cleanup;
    }

    // add variants
    if (add_from_str_map(recs->var_ix, var_map) < 0){
        ret = -1;
        goto cleanup;
    }

    fp = fopen(outrecsfn, "r");
    if (fp == NULL){
        err_msg(-1, 0, "read_records_bn: failed to open input file %s: %s", outrecsfn, strerror(errno));
        ret = -1;
        goto cleanup;
    }

    size_t len;
    size_t buf_size = 5;
    char beg[5] = "recs\0";
    buf = malloc(buf_size);

    if (buf == NULL){
        err_msg(-1, 0, "read_records_bn: cannot read header %s", strerror(errno));
        ret = -1;
        goto cleanup;
    }

    len = fread(buf, 1, 5, fp);

    if (strcmp(beg, (char*)buf) != 0){
        err_msg(-1, 0, "read_records_bn: input file %s does not have correct format", outrecsfn);
        ret = -1;
        goto cleanup;
    }

    while ( !feof(fp) && !ferror(fp)){
        if (feof(fp)) break;
        
        /* block size */
        uint32_t block_size;
        len = fread(buf, 1, 4, fp);
        if (feof(fp)) break;
        block_size = le_to_u32((uint8_t *)buf);
        if (buf_size < block_size){
            buf_size = block_size;
            if (buf_size == 0) fprintf(stderr, "bufsize for block = 0\n");
            buf = (uint8_t *)realloc(buf, buf_size);
            if (buf_size > 0 && buf == NULL){
                err_msg(-1, 0, "read_records_bn: %s", strerror(errno));
                ret = -1;
                goto cleanup;
            }
        }
        /* */

        /* record */
        len = fread(buf, 1, block_size, fp);
        if (len != block_size){
            err_msg(-1, 0, "read_records_bn: error reading records from %s: %s\n", outrecsfn, strerror(errno));
            ret = -1;
            goto cleanup;
        }

        uint32_t bc_ix=-1, umi_ix=-1;
        read_rec_bn(buf, r, &bc_ix, &umi_ix);

        char *bc = str_map_str(bc_map, bc_ix);
        char *umi = str_map_str(umi_map, umi_ix);
        if (r->n_feat > feat_m){
            feat_m = r->n_feat;
            if (feat_m == 0) fprintf(stderr, "feat_m = 0\n");
            feat = realloc(feat, feat_m * sizeof(char *));
            if (feat_m > 0 && feat == NULL){
                err_msg(-1, 0, "read_records_bn: %s", strerror(errno));
                ret = -1;
                goto cleanup;
            }
        }
        for (i = 0; i < r->n_feat; ++i)
            feat[i] = str_map_str(gene_map, r->feat[i]);
        if (r->n_var > var_m){
            var_m = r->n_var;
            if (var_m == 0) fprintf(stderr, "var_m = 0\n");
            var = realloc(var, var_m * sizeof(char *));
            if (var_m > 0 && var == NULL){
                err_msg(-1, 0, "read_records_bn: %s", strerror(errno));
                ret = -1;
                goto cleanup;
            }
        }
        for (i = 0; i < r->n_var; ++i)
            var[i] = str_map_str(var_map, r->var[i]);

        ret = add_rec(recs, (const char *)bc, (const char *)umi, 
                r->n_feat, (const char **)feat, (const uint8_t *)r->splice, 
                r->n_var, (const char **)var, (const uint8_t *)r->base, (const uint8_t *)r->qual);
        if (ret < 0){
            err_msg(-1, 0, "read_records_bn: failed to add record");
            ret = -1;
            break;
        }
        /* */
    }
    str_map_resize(recs->bc_ix);
    str_map_resize(recs->umi_ix);
    str_map_resize(recs->gene_ix);
    str_map_resize(recs->var_ix);

cleanup:
    if (fp) fclose(fp);
    if (bc_map) destroy_str_map(bc_map);
    if (umi_map) destroy_str_map(umi_map);
    if (gene_map) destroy_str_map(gene_map);
    if (var_map) destroy_str_map(var_map);
    if (feat) free(feat);
    if (var) free(var);
    free(outrecsfn);
    free(outbcixfn); free(outumiixfn); free(outgeneixfn); free(outvarixfn);
    if (buf) free(buf);
    destroy_rec(r);
    return ret;
}

Records *read_records_file(char *fn, int file_list, int verbose){
    Records *records = init_records();

    if (file_list == 0){
        if (verbose) fprintf(stdout, "reading records from %s\n", fn);
        if ((read_records_bn(fn, records)) < 0) return NULL;
    } else {
        int nfn = 0;
        char **recsfns = read_lines(fn, &nfn);
        if (recsfns == NULL) return NULL;
        int k;
        for (k = 0; k < nfn; ++k){
            if (verbose) fprintf(stdout, "samba: reading records from %s\n", recsfns[k]);
            if ((read_records_bn(recsfns[k], records)) < 0) return NULL;
        }
        for (k = 0; k < nfn; ++k) free(recsfns[k]);
        free(recsfns);
    }
    return records;
}

int n_recs(Records *recs){
    int nrec = 0;

    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            Rec *rec = u->recs;
            while (rec){
                nrec++;
                rec = rec->next;
            }
        }
    }
    return nrec;
}

void print_summary(Records *recs){
    uint32_t nrec = 0, nspl = 0, nunspl = 0, namb = 0, nref = 0, nalt = 0;

    uint32_t nfeat = (uint32_t)recs->gene_ix->n;
    uint32_t nsnp = (uint32_t)recs->var_ix->n;


    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            Rec *rec = u->recs;
            while (rec){
                nrec++;
                int i;
                for (i = 0; i < rec->n_feat; ++i){
                    if ((rec->splice)[i] == SPLICE) nspl++;
                    else if ((rec->splice)[i] == UNSPLICE) nunspl++;
                    else namb++;
                }
                for (i = 0; i < rec->n_var; ++i){
                    if ( (rec->base)[i] == REF ) nref++;
                    else nalt++;
                }
                rec = rec->next;
            }
        }
    }
    fprintf(stdout, "Num. records=%u\n", nrec);
    fprintf(stdout, "Num. genes=%u\n", nfeat);
    fprintf(stdout, "Num. spliced=%u\n", nspl);
    fprintf(stdout, "Num. unspliced=%u\n", nunspl);
    fprintf(stdout, "Num. ambiguous=%u\n", namb);
    fprintf(stdout, "Num. SNPs=%u\n", nsnp);
    fprintf(stdout, "Num. REF=%u\n", nref);
    fprintf(stdout, "Num. ALT=%u\n", nalt);
}
    
void print_best_umi(Records *recs){
    uint32_t numi = 0, nspl = 0, nunspl = 0, namb = 0, nfeat = 0, nsnp = 0, nref = 0, nalt = 0;

    int i;
    khint_t k;
    for (k = kh_begin(recs->bc); k != kh_end(recs->bc); ++k){
        if (!kh_exist(recs->bc, k)) continue;
        Barcode *bc = kh_val(recs->bc, k);
        khint_t l;
        for (l = kh_begin(bc->umis); l != kh_end(bc->umis); ++l){
            if (!kh_exist(bc->umis, l)) continue;
            UMI *u = kh_val(bc->umis, l);
            char *key = kh_key(bc->umis, l);
            Rec *r = u->best_rec;
            if (r == NULL) continue;
            numi++;
            fprintf(stdout, "%s: ", key);
            fprintf(stdout, "n_feat=%u", r->n_feat);
            nfeat += r->n_feat;
            fprintf(stdout, "; feat=");
            for (i = 0; i < r->n_feat; ++i){
                if (i) fprintf(stdout, ",");
                fprintf(stdout, "%i", r->feat[i]);
            }
            fprintf(stdout, "; splice=");
            for (i = 0; i < r->n_feat; ++i){
                if (i) fprintf(stdout, ",");
                fprintf(stdout, "%i", r->splice[i]);
                if (r->splice[i] == SPLICE) nspl++;
                else if (r->splice[i] == UNSPLICE) nunspl++;
                else namb++;
            }
            fprintf(stdout, "; n_var=%u", r->n_var);
            nsnp += r->n_var;
            fprintf(stdout, "; vars=");
            for (i = 0; i < r->n_var; ++i){
                if (i) fprintf(stdout, ",");
                fprintf(stdout, "%i", r->var[i]);
            }
            fprintf(stdout, "; base=");
            for (i = 0; i < r->n_var; ++i){
                if (i) fprintf(stdout, ",");
                fprintf(stdout, "%u", r->base[i]);
                if (r->base[i] == REF) nref++;
                else if (r->base[i] == ALT) nalt++;
            }
            fprintf(stdout, "; qual=");
            for (i = 0; i < r->n_var; ++i){
                if (i) fprintf(stdout, ",");
                fprintf(stdout, "%u", r->qual[i]);
            }
            fprintf(stdout, "; n_read=%i\n", r->n_read);
        }
    }
    fprintf(stdout, "Num. UMIs=%u\n", numi);
    fprintf(stdout, "Num. genes=%u\n", nfeat);
    fprintf(stdout, "Num. spliced=%u\n", nspl);
    fprintf(stdout, "Num. unspliced=%u\n", nunspl);
    fprintf(stdout, "Num. ambiguous=%u\n", namb);
    fprintf(stdout, "Num. SNPs=%u\n", nsnp);
    fprintf(stdout, "Num. REF=%u\n", nref);
    fprintf(stdout, "Num. ALT=%u\n", nalt);
}

