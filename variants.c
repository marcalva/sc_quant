
#include "variants.h"
#include "bins.h"
#include "gtf_anno.h"
#include "overlap.h"
#include "str_util.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>

// TODO it b->n_allele number of ref + alt alleles?
int get_var_len(bcf1_t *b){
    int len = 1;
    int i;
    for (i = 0; i < b->n_allele; i++){
        int ilen = strlen(b->d.allele[i]);
        if (ilen > len) 
            len = ilen;
    }
    return len;
}

// TODO is b->pos 0-length?
int get_var_bin(bcf1_t *b){
    int beg = b->pos;
    int len = get_var_len(b);
    int end = beg + len;
    int bin = reg2bin(beg, end);
    return bin;
}

char *var_id(const bcf_hdr_t *h, bcf1_t *b){

    char delim = ' ';

    // chromosome
    const char *chrm = bcf_seqname(h, b);
    size_t nchrm = strlen(chrm);

    // pos
    size_t npos = 0;
    char *pos = int2str(b->pos + 1, &npos);
    if (pos == NULL) return NULL;

    // ID
    size_t nid = strlen(b->d.id);

    // alleles
    size_t allele_len = 0;
    int *allele_lens = (int *)malloc(sizeof(int) * b->n_allele);
    int i;
    for (i = 0; i < b->n_allele; ++i){
        allele_lens[i] = strlen(b->d.allele[i]);
        if (i > 1) allele_lens[i] += 1;
        allele_len += allele_lens[i];
    }

    // total
    size_t ntotal = nchrm + npos + nid + allele_len + 5;
    char *out = (char *)calloc(ntotal, sizeof(char));
    if (out == NULL){
        err_msg(-1, 0, "var_id: %s", strerror(errno));
        return(NULL);
    }
    char *tmp = out;

    memcpy(tmp, chrm, nchrm * sizeof(char));
    tmp = tmp + nchrm;
    *(tmp++) = delim;

    memcpy(tmp, pos, npos * sizeof(char));
    tmp = tmp + npos;
    *(tmp++) = delim;

    memcpy(tmp, b->d.id, nid * sizeof(char));
    tmp = tmp + nid;
    *(tmp++) = delim;

    memcpy(tmp, b->d.allele[0], allele_lens[0] * sizeof(char));
    tmp = tmp + allele_lens[0];
    *(tmp++) = delim;

    for (i = 1; i < b->n_allele; ++i){
        if (i > 1) *(tmp++) = ';';
        memcpy(tmp, b->d.allele[i], allele_lens[0] * sizeof(char));
        tmp = tmp + allele_lens[0];
    }

    free(pos);
    free(allele_lens);

    return out;
}

uint8_t base_ref_alt(bcf1_t *b, char base){
    if ( base == b->d.allele[0][0] ) return ((uint8_t)REF);
    else if ( base == b->d.allele[1][0] ) return ((uint8_t)ALT);
    else return((uint8_t)OTHER);
}

Var *init_var(){
    Var *v = (Var*)calloc(1, sizeof(Var));
    if (v == NULL){
        err_msg(-1, 0, "var_id: %s", strerror(errno));
        return NULL;
    }
    v->next = NULL;
    return v;
}

GenomeVar *init_genomevar(){
    int init_n = 1<<8;
    GenomeVar *gv = (GenomeVar*)calloc(1, sizeof(GenomeVar));
    
    if (gv == NULL){
        err_msg(-1, 0, "init_genomevar: %s", strerror(errno));
        return NULL;
    }

    gv->chrm_ix = init_str_map();
    gv->var_ix = init_str_map();
    gv->vars = kh_init(var);
    gv->chrms = (ChrmVar**)calloc(init_n, sizeof(ChrmVar*));

    if (gv->chrms == NULL){
        err_msg(-1, 0, "init_genomevar: %s", strerror(errno));
        return NULL;
    }

    gv->chrms_m = init_n;
    return gv;
}

ChrmVar *init_ChrmVar(){
    ChrmVar *chrm = (ChrmVar*)calloc(1, sizeof(ChrmVar));

    if (chrm == NULL){
        err_msg(-1, 0, "init_ChrmVar: %s", strerror(errno));
        return NULL;
    }

    int i;
    for (i = 0; i < MAX_BIN; i++){
        chrm->bins[i] = NULL;
        chrm->vars_n[i] = 0;
    }
    return chrm;
}

int add_var(GenomeVar *gv, bcf1_t *b, const bcf_hdr_t *hdr){
    Var *tv = init_var();
    if (tv == NULL) return -1;
    tv->b = bcf_dup(b);
    bcf_unpack(tv->b, BCF_UN_INFO); // unpack up to and including info field

    int bin = get_var_bin(b);

    /* add variant index */
    khint_t k;
    int ret;
    char *vid = var_id(hdr, b);
    if (vid == NULL) return -1;

    int found;
    int vix = add2str_map(gv->var_ix, (const char *)vid, &found);
    if (vix < 0) return -1;

    if (found)
        return err_msg(-1, 0, "add_var: variant %s found twice", vid);

    tv->vid = str_map_str(gv->var_ix, vix);
    free(vid);
    /* */

    int rid = tv->b->rid;
    const char *seqname = bcf_hdr_id2name(hdr, rid);

    int chr_ix;
    // add chromosome ID
    if ( (chr_ix = add2str_map(gv->chrm_ix, (const char*)seqname, &found)) < 0 ) return -1;
    if (found == 0){
        while (gv->chrm_ix->n >= gv->chrms_m){
            gv->chrms_m = (gv->chrms_m)<<1;
            gv->chrms = realloc(gv->chrms, (gv->chrms_m)*sizeof(ChrmVar*));
            if (gv->chrms == NULL)
                return err_msg(-1, 0, "add_var: %s", strerror(errno));
        }
        gv->chrms[chr_ix] = init_ChrmVar();
        if (gv->chrms[chr_ix] == NULL) return -1;
    }

    // add to bins
    if (gv->chrms[chr_ix]->bins[bin]){
        Var *gvv = gv->chrms[chr_ix]->bins[bin];
        while (gvv->next){
            gvv = gvv->next;
        }
        gvv->next = tv;
    }
    else {
        gv->chrms[chr_ix]->bins[bin] = tv;
    }
    (gv->chrms[chr_ix]->vars_n[bin])++;

    // add to var hash
    k = kh_put(var, gv->vars, tv->vid, &ret);
    if (ret == -1){
        fprintf(stderr, "samba: failed to add variant ID %s to vars hash.\n", tv->vid);
        return(-1);
    }
    kh_val(gv->vars, k) = tv;

    return 0;
}

GenomeVar *vcf2gv(bcf_srs_t *sr, bcf_hdr_t *vcf_hdr){
    GenomeVar *gv = init_genomevar();
    if (gv == NULL) return NULL;

    while ( bcf_sr_next_line(sr) ){
        bcf1_t *vcf_r = bcf_sr_get_line(sr, 0);
        if (vcf_r == NULL) fprintf(stderr, "warning: a VCF record is NULL\n");

        int vt;
        vt = bcf_get_variant_types(vcf_r);
        if (vt != VCF_SNP)
            continue;
        int na = vcf_r->n_allele;
        if (na != 2){
            continue;
        }
        if ( add_var(gv, vcf_r, vcf_hdr) < 0 ) return NULL;
    }

    return gv;
}

void destroy_gv(GenomeVar *gv){

    if (gv == NULL) return;

    // chrms
    int i, j;
    for (i = 0; i < gv->chrm_ix->n; i++){
        for (j = 0; j < MAX_BIN; j++){
            Var *v = gv->chrms[i]->bins[j];
            while (v){
                bcf1_t *b = v->b;
                if (b) bcf_destroy(b);
                Var *vn = v->next;
                free(v);
                v = vn;
            }
        }
        free(gv->chrms[i]);
    }
    free(gv->chrms);
    destroy_str_map(gv->chrm_ix);
    destroy_str_map(gv->var_ix);
    kh_destroy(var, gv->vars);

    free(gv);
}

void free_bcf1_t(GenomeVar *gv){
    int i, j;
    for (i = 0; i < gv->chrm_ix->n; i++){
        for (j = 0; j < MAX_BIN; j++){
            Var *v = gv->chrms[i]->bins[j];
            while (v){
                bcf1_t *b = v->b;
                if (b) bcf_destroy(b);
                v->b = NULL;
                v = v->next;
            }
        }
    }
}

int is_gt_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b){
    bcf_unpack(b, BCF_UN_FMT);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return 0;

    // variant must be bi-allelic
    if (b->n_allele != 2) return 0;

    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GT", &gt_arr, &ndst, BCF_HT_INT); // GT must be given as BCF_HT_INT

    // check if tag is present
    if (ngt < 0){
        free(gt_arr);
        return 0;
    }
    int n_val = ngt / n_samples;

    // check max n alleles
    if (n_val > 2){
        free(gt_arr);
        return 0;
    }

    // check missing
    void *p = gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        int t;
        for (t = 0; t < n_val; t++){ // each allele
            int32_t *a = (int32_t *)p;
            a = a + (s * n_val) + t;

            if ( bcf_gt_is_missing(*a) ){
                n_miss++;
                break;
            }
        }
    }
    free(gt_arr);
    if (n_miss == n_samples) return -1;

    return 1;
}
        
int is_gp_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b){
    bcf_unpack(b, BCF_UN_FMT);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return 0;

    // variant must be bi-allelic
    if (b->n_allele != 2) return 0;

    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GP", &gt_arr, &ndst, BCF_HT_REAL);

    // check if tag is present
    if (ngt < 0) return 0;
    int n_val = ngt / n_samples;

    // check max n alleles
    if (n_val > 3){
        free(gt_arr);
        return 0;
    }

    // check missing
    void *p = gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        int t;
        for (t = 0; t < n_val; t++){ // each allele
            int32_t *a = (int32_t *)p;
            a = a + (s * n_val) + t;

            if ( *a == bcf_float_missing ){
                n_miss++;
                break;
            }
        }
    }
    free(gt_arr);
    if (n_miss == n_samples) return -1;

    return 1;
}
        
float *bcf1_ap_gt(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra){

    bcf_unpack(b, BCF_UN_FMT);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return NULL;

    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GT", &gt_arr, &ndst, BCF_HT_INT); // GT must be given as BCF_HT_INT
    if (ngt < 0){
        err_msg(-1, 0, "bcf1_ap_gt: failed to get genotypes");
        free(gt_arr);
        return NULL;
    }
    int n_val = ngt / n_samples;

    // allocate dose array
    float *dose = malloc((n_samples + extra) * sizeof(float));
    if (dose == NULL){
        err_msg(-1, 0, "bcf1_ap_gt: %s", strerror(errno));
        return NULL;
    }

    // fill dose array
    void *p = gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        float gt = 0, total = 0;
        int t;
        for (t = 0; t < n_val; t++){ // each allele
            int32_t *a = (int32_t *)p;
            a = a + (s * n_val) + t;

            // if end, sample has smaller ploidy, break
            if ( *a == bcf_int32_vector_end ) break;

            // if any missing
            if ( bcf_gt_is_missing(*a) ) {
                n_miss++;
                gt = -1.0;
                total = 1.0;
                break;
            }

            gt += (float)bcf_gt_allele(*a);
            total++;
        }
        
        if (total == 0) dose[s] = 0;
        else dose[s] = gt / total;
    }
    free(gt_arr);

    // if all missing, return NULL
    // if (n_miss == n_samples) return NULL;

    return dose;
}

float *bcf1_ap_gp(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra){

    bcf_unpack(b, BCF_UN_FMT);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return NULL;

    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GP", &gt_arr, &ndst, BCF_HT_REAL);
    if (ngt < 0){
        err_msg(-1, 0, "bcf1_dose_gp: failed to get genotypes %s", b->d.id);
        free(gt_arr);
        return NULL;
    }
    int n_val = ngt / n_samples;

    bcf_fmt_t *fmt = bcf_get_fmt(vcf_hdr, b, "GP");
    if (fmt->n != n_val)
        fprintf(stderr, "fmt->n %i doesn't match n_val %i for SNP %s\n", fmt->n, n_val, b->d.id);

    // allocate dose array
    float *dose = malloc((n_samples + extra) * sizeof(float));
    if (dose == NULL){
        err_msg(-1, 0, "bcf1_dose_gp: %s", strerror(errno));
        free(gt_arr);
        return NULL;
    }

    // fill dose array
    void *p = gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        float gt = 0;
        int t, total = 0;
        for (t = 0; t < n_val; t++){ // each allele
            float *a = (float *)p;
            a = a + (s * n_val) + t;

            // if end, sample has smaller ploidy, break
             if ( *a == bcf_float_vector_end ) break;

            // if any are missing, set to -1
            if ( *a == bcf_float_missing ){
                fprintf(stdout, "missing\n");
                n_miss++;
                gt = -1.0;
                total = 1;
                break;
            }

            gt += t * *a;
            total++;
        }

        switch ( total ) {
            case 0: 
                dose[s] = -1;
                break;
            case 1:
                dose[s] = gt;
                break;
            case 2:
                dose[s] = -1;
                break;
            case 3:
                dose[s] = gt / 2;
                break;
            default:
                err_msg(-1, 0, "bcf1_dose_gp: total number of genotypes for %s is inconsistent: %i", 
                        b->d.id, total);
                return NULL;
        }
    }
    free(gt_arr);

    // if all missing, return NULL
    if (n_miss == n_samples) return NULL;

    return dose;
}

float **ap_array_gt(GenomeVar *gv, bcf_hdr_t *vcf_hdr, char **ids, int ni, char *field){

    // check input
    if ( strcmp(field, "GP") != 0 && strcmp(field, "GT") != 0){
        err_msg(-1, 0, "ap_array_gt: %s: genotype field must be one of GT or GP\n", field);
        return NULL;
    }

    // check if field is present in header
    int gp_tag_id = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, "GP");
    int gt_tag_id = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, "GT");

    if (gp_tag_id < 0 && gt_tag_id < 0){
        err_msg(-1, 0, "ap_array_gt: neither GP nor GT is present in the VCF header\n");
        return NULL;
    }

    float **ap = malloc(ni * sizeof(float *));
    if (ap == NULL){
        err_msg(-1, 0, "ap_array_gt: %s", strerror(errno));
        return NULL;
    }
    int i;
    for (i = 0; i < ni; ++i) ap[i] = NULL;


    for (i = 0; i < ni; ++i){ // for each variant
        char *vid = ids[i];
        if (vid == NULL) continue;

        // get Var
        khint_t k = kh_get(var, gv->vars, vid);
        if (k == kh_end(gv->vars)){
            err_msg(-1, 0, "ap_array_gt: variant %s not found", vid);
            int j;
            for (j = 0; j < i; ++j) free(ap[j]);
            free(ap);
            return NULL;
        }
        Var *v = kh_val(gv->vars, k);

        // check if genotype fields are valid
        int gp_valid, gt_valid;
        if ( (gp_valid = is_gp_valid(vcf_hdr, v->b)) == -2) return NULL;
        if ( (gt_valid = is_gt_valid(vcf_hdr, v->b)) == -2) return NULL;

        if ( strcmp(field, "GP") == 0 && gp_valid > 0 )
            ap[i] = bcf1_ap_gp(vcf_hdr, v->b, 1);
        else if ( gt_valid > 0 )
            ap[i] = bcf1_ap_gt(vcf_hdr, v->b, 1);
        else{
            err_msg(-1, 0, "ap_array_gt: error getting genotypes from variant %s", vid);
            int j;
            for (j = 0; j < i; ++j) free(ap[j]);
            free(ap);
            return NULL;
        }

        // if error getting genotypes
        if ( ap[i] == NULL ){
            err_msg(-1, 0, "ap_array_gt: error getting genotypes from variant %s", vid);
            int j;
            for (j = 0; j < i; ++j) free(ap[j]);
            free(ap);
            return NULL;
        }
    }
    return ap;
}

// return pointers to Var objects
int vars_from_region(GenomeVar *gv, const char* ref, hts_pos_t beg, 
        hts_pos_t end, Var ***vars, int *vars_m){
    int tid = str_map_ix(gv->chrm_ix, (char *)ref);
    if (tid < 0)
        return err_msg(-1, 0, "vars_from_region: chromosome %s not found", ref);

    double reg_len = (double)end - (double)beg;
    if (reg_len < 0)
        return err_msg(-1, 0, "vars_from_region: end (%i) < beg (%i)", beg, end);

    int nvars = 0;
    if (*vars_m <= 0){
        *vars_m = 1;
        *vars = (Var **)realloc(*vars, *vars_m * sizeof(Var *));
        if (*vars == NULL)
            return err_msg(-1, 0, "vars_from_region: %s", strerror(errno));
    }

    uint16_t list[MAX_BIN];
    int n_bin = reg2bins((int)beg, (int)end, list);
    int i;
    for (i = 0; i < n_bin; ++i){
        uint16_t bin = list[i];
        Var *v = gv->chrms[tid]->bins[bin]; 
        for (; v; v = v->next){
            bcf1_t *b = v->b;
            int ovrlp = bp_overlap((int)beg, (int)end, '.', b->pos, b->pos + b->rlen, '.');
            if (ovrlp < 0)
                return err_msg(-1, 0, "vars_from_region: failed to get bp_overlap");

            double frac_ovrlp = (double)ovrlp / reg_len;
            if (frac_ovrlp == 0.0) continue;
            if (nvars >= *vars_m){
                *vars_m = (*vars_m)<<1;
                *vars = (Var **)realloc(*vars, (*vars_m) * sizeof(Var *));
                if (*vars == NULL)
                    return err_msg(-1, 0, "vars_from_region: %s", strerror(errno));
            }
            (*vars)[nvars++] = v;
        }
    }
    return(nvars);
}

int n_snp(GenomeVar *gv, int *n_snp){
    *n_snp = 0;
    int nc = gv->chrm_ix->n;
    int i, j;
    for (i = 0; i < nc; i++){
        for (j = 0; j < MAX_BIN; j++){
            *n_snp = *n_snp + gv->chrms[i]->vars_n[j];
        }
    }
    return nc;
}

void print_dose(double **dose, char **var_ids, int dose_len, int nsamples){
    int v;
    for (v = 0; v < dose_len; ++v){
        int s;
        fprintf(stdout, "%s", var_ids[v]);
        for (s = 0; s < nsamples; ++s) fprintf(stdout, "\t%f", dose[v][s]);
        fprintf(stdout, "\n");
    }
}


void print_bcf_hdr(bcf_hdr_t *h){

    int igt, a, b;
    igt = 0;
    bcf_gt2alleles(igt, &a, &b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);
    igt = 1;
    bcf_gt2alleles(igt, &a, &b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);
    igt = 2;
    bcf_gt2alleles(igt, &a, &b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);

    a = 0; b =0;
    igt = bcf_alleles2gt(a, b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);
    a = 0; b = 1;
    igt = bcf_alleles2gt(a, b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);
    a = 1; b =0;
    igt = bcf_alleles2gt(a, b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);
    a = 1; b =1;
    igt = bcf_alleles2gt(a, b);
    fprintf(stdout, "igt=%i, a=%i, b=%i\n", igt, a, b);
    
    fprintf(stdout, "\n============= header info =============\n\n");

    // n object
    fprintf(stdout, "n[BCF_DT_ID] = %i; n[BCF_DT_CTG] = %i; n[BCF_DT_SAMPLE] = %i\n", 
            h->n[BCF_DT_ID], h->n[BCF_DT_CTG], h->n[BCF_DT_SAMPLE]); 

    // id (BCF_DT_ID)
    fprintf(stdout, "The BCF_DT_ID keys are: \n");
    int i, j;
    for (i = 0; i < h->n[BCF_DT_ID]; i++){
        // unsigned int infoN0 = (unsigned int)bcf_hdr_id2type(h, 0, i);
        // unsigned int infoN1 = (unsigned int)bcf_hdr_id2type(h, 1, i);
        // unsigned int infoN2 = (unsigned int)bcf_hdr_id2type(h, 2, i);
        fprintf(stdout, "key: %s", h->id[BCF_DT_ID][i].key);
        fprintf(stdout, ", id: %i", h->id[BCF_DT_ID][i].val->id);
        // fprintf(stdout, ", info[0] num = %u; info[1] num = %u; info[0] num = %u", infoN0, infoN1, infoN2);
        // The info keys are not printable.
        // fprintf(stdout, ", info[0]: %"PRIu64", info[1]: %"PRIu64", info[2]: %"PRIu64, 
        //         h->id[BCF_DT_ID][i].key, h->id[BCF_DT_ID][i].val->info[0], h->id[BCF_DT_ID][i].val->info[1], h->id[BCF_DT_ID][i].val->info[2]);
        // fprintf(stdout, ", hrec type[0]: %i, hrec type[1]: %i, hrec type[2]: %i", 
        //         h->id[BCF_DT_ID][i].val->hrec[0]->type, h->id[BCF_DT_ID][i].val->hrec[2]->type, h->id[BCF_DT_ID][i].val->hrec[2]->type);
        fprintf(stdout, "\n");
    }
    // fprintf(stdout, "\n");

    // dict
    // don't think I need to use this member


    // samples
    fprintf(stdout, "The %i samples are: ", h->n[BCF_DT_SAMPLE]);
    for (i = 0; i < h->n[BCF_DT_SAMPLE]; i++){
        fprintf(stdout, "\t%s", h->samples[i]);
    }
    fprintf(stdout, "\n");
    

    // hrec object
    fprintf(stdout, "The header has %i hrecs:\n", h->nhrec);
    for (i = 0; i < h->nhrec; i++){
        fprintf(stdout, "type = %i; nkeys = %i; key = %s; value = %s", 
                h->hrec[i]->type, h->hrec[i]->nkeys, h->hrec[i]->key, h->hrec[i]->value);
        for (j = 0; j < h->hrec[i]->nkeys; j++){
            fprintf(stdout, "; %s=%s", h->hrec[i]->keys[j], h->hrec[i]->vals[j]);
        }
        fprintf(stdout, "\n");
        
        //fprintf("type = %i\n", h->hrec[i]->type);
    }

    // 

    fprintf(stdout, "\n=============             =============\n\n");
}

void print_var(GenomeVar *gv, const bcf_hdr_t *vcf_hdr){
    int nc = gv->chrm_ix->n;
    int i, j, k;
    for (i = 0; i < nc; i++){
        for (j = 0; j < MAX_BIN; j++){
            Var *v = gv->chrms[i]->bins[j];
            while (v){
                bcf1_t *b = v->b;
                bcf_unpack(b, BCF_UN_ALL);
                /* b members */
                fprintf(stdout, "\n***** %s *****\n", b->d.id);
                char *vid = var_id(vcf_hdr, b);
                fprintf(stdout, "vid=%s\n", vid);
                free(vid);
                fprintf(stdout, "pos=%li; rlen=%li; rid=%i; qual=%f; n_info=%i; n_allele=%i; n_fmt=%i; n_sample=%i", 
                        b->pos, b->rlen, b->rid, b->qual, b->n_info, b->n_allele, b->n_fmt, b->n_sample);
                
                fprintf(stdout, "; max_unpack=%i; unpacked=%i; errcode=%i", 
                        b->max_unpack, b->unpacked, b->errcode);

                /* b->unpack_size not sure what this is */
                fprintf(stdout, "; unpack_size_1=%i; unpack_size_2=%i; unpack_size_3=%i\n", 
                        b->unpack_size[0], b->unpack_size[1], b->unpack_size[2]);

                /* b->shared cannot just print out */
                // fprintf(stdout, "\tshared %s\n", b->shared.s);

                /* b->indiv cannot just print out */
                // fprintf(stdout, "\tindiv %s\n", b->indiv.s);

                /* b->d */
                fprintf(stdout, "\td->m_fmt=%i; d->m_info=%i; d->m_id=%i; d->m_als=%i; d->m_allele=%i; d->m_flt: %i", 
                        b->d.m_fmt, b->d.m_info, b->d.m_id, b->d.m_als, b->d.m_allele, b->d.m_flt);
                for (k = 0; k < b->d.n_flt; k++){
                    fprintf(stdout, "; d->flt[%i]=%i", k, b->d.flt[k]);
                }
                fprintf(stdout, "; d->als=%s; d->n_flt=%i; d->n_var=%i; d->var_type=%i; d->shared_dirty=%i, d->indiv_dirty=%i\n", 
                        b->d.als, b->d.n_flt, b->d.n_var, b->d.var_type, b->d.shared_dirty, b->d.indiv_dirty);

                /*  */
                for (k = 0; k < b->d.n_flt; k++){
                    fprintf(stdout, "\tflt %i: %i", k, b->d.flt[k]);
                }
                fprintf(stdout, "\n");

                /* b->d.alleles */
                for (k = 0; k < b->n_allele; k++){
                    if (k) fprintf(stdout, "; ");
                    else fprintf(stdout, "Alleles: ");
                    fprintf(stdout, "allele[%i]=%s", k, b->d.allele[k]);
                }
                fprintf(stdout, "\n");

                /* b->d.info */
                for (k = 0; k < b->n_info; k++){
                    int infokey = b->d.info[k].key;
                    int infotype = b->d.info[k].type;
                    // uint8_t *vptr = b->d.info[k].vptr;
                    uint32_t vptr_len = b->d.info[k].vptr_len;
                    int len = b->d.info[k].len;
                    // char *strkey = vcf_hdr->id[BCF_DT_ID][infokey].key;
                    // same as:
                    const char *strkey = bcf_hdr_int2id(vcf_hdr, BCF_DT_ID, infokey);
                    uint32_t hdrtype = bcf_hdr_id2type(vcf_hdr, BCF_HL_INFO, infokey); // returns one of BCF_HT_*
                    // uint32_t hdrtype = bcf_hdr_id2coltype(vcf_hdr, BCF_HL_INFO, infokey); // returns one of BCF_HL_* Should always be 1

                    fprintf(stdout, "INFO[%i]={", k);
                    fprintf(stdout, "key_ix=%i; strkey=%s; type=%i; vptr_len=%u; len=%i; hdrtype=%u", 
                            infokey, strkey, infotype, vptr_len, len, hdrtype);
                    // get INFO value
                    int si = 0;
                    switch (infotype) {
                        case BCF_BT_NULL: si = 1 * len; break; 
                        case BCF_BT_INT8: si = 1 * len; break; 
                        case BCF_BT_INT16: si = 2 * len; break; 
                        case BCF_BT_INT32: si = 4 * len; break; 
                        case BCF_BT_INT64: si = 8 * len; break; 
                        case BCF_BT_FLOAT: si = 8 * len; break; 
                    }

                    // void *val = malloc(si);
                    void *val = NULL; // will realloc
                    si = 0;
                    int ret;
                    ret = bcf_get_info_values(vcf_hdr, b, (const char*)strkey, &val, &si, hdrtype);
                    if (ret < 0){
                        fprintf(stdout, "RETURN %i < 0\n", ret);
                        exit(1);
                    }

                    switch (infotype) {
                        case BCF_BT_NULL: fprintf(stdout, "; value=None"); break;
                        case BCF_BT_INT8: fprintf(stdout, "; value=%"PRIu8, *(int8_t*)val); break;
                        case BCF_BT_INT16: fprintf(stdout, "; value=%"PRIu16, *(int16_t*)val); break;
                        case BCF_BT_INT32: fprintf(stdout, "; value=%"PRIu32, *(int32_t*)val); break;
                        case BCF_BT_INT64: fprintf(stdout, "; value=%"PRIu64, *(int64_t*)val); break;
                        case BCF_BT_FLOAT: fprintf(stdout, "; value=%.6f", *(float*)val); break;
                        case BCF_BT_CHAR: fprintf(stdout, "; value=%s", (char*)val); break;
                    }
                    fprintf(stdout, "}");
                    free(val);

                }
                fprintf(stdout, "\n");

                /* FMT */
                for (k = 0; k < b->n_fmt; k++){
                    int fmt_id = b->d.fmt[k].id; // numeric TAG id
                    int fmt_n = b->d.fmt[k].n; // number values in each sample
                    int fmt_size = b->d.fmt[k].size; // number of bytes
                    int fmt_type = b->d.fmt[k].type; // one of BCF_BT_*
                    uint32_t fmt_plen = b->d.fmt[k].p_len;
                    uint32_t fmt_poff = b->d.fmt[k].p_off;
                    uint32_t fmt_pfree = b->d.fmt[k].p_free;
                    uint32_t nsamples = bcf_hdr_nsamples(vcf_hdr);

                    const char *strkey = bcf_hdr_int2id(vcf_hdr, BCF_DT_ID, fmt_id);
                    uint32_t hdrtype = bcf_hdr_id2type(vcf_hdr, BCF_HL_FMT, fmt_id); // returns one of BCF_HT_*

                    fprintf(stdout, "FMT[%i]={", k);
                    fprintf(stdout, "ID=%i; key=%s; n=%i; size=%i; type=%u; plen=%u; poff=%u; pfree=%u; nsamples=%u; hdrtype=%u}", 
                            fmt_id, strkey, fmt_n, fmt_size, fmt_type, fmt_plen, fmt_poff, fmt_pfree, nsamples, hdrtype);


                    // get values 
                    // void *val = malloc(1); // will realloc
                    void *val = NULL;
                    int ndst = 0;
                    int ret;
                    if (strcmp(strkey, "GT") == 0){
                        ret = bcf_get_format_values(vcf_hdr, b, strkey, &val, &ndst, BCF_HT_INT);
                    }
                    else{
                        ret = bcf_get_format_values(vcf_hdr, b, strkey, &val, &ndst, hdrtype);
                    }
                    if (ret < 0){
                        fprintf(stdout, "RETURN %i < 0\n", ret);
                        exit(1);
                    }

                    fprintf(stdout, "%i values: {", ret);
                    int s;
                    void *p = val;
                    for (s = 0; s < nsamples; s++){
                        int t;
                        for (t = 0; t < fmt_n; t++){
#define BRANCH(type_t, print_m){ uint32_t *a = (uint32_t *)p; a = a + (s*fmt_n) + t; type_t *b = (type_t *)a; fprintf(stdout, "; fmt[%u][%i]=%"print_m, s, t, *b);}
#define BRANCH_GT(type_t, print_m){ uint32_t *a = (uint32_t *)p; a = a + (s*fmt_n) + t; type_t *b = (type_t *)a; fprintf(stdout, "; fmt[%u][%i]=%"print_m, s, t, bcf_gt_allele(*b));}
                            if (strcmp(strkey, "GT") == 0){
                                BRANCH_GT(int32_t, PRIi32);
                            }else {
                                switch (fmt_type) {
                                    case BCF_BT_NULL: fprintf(stdout, "; fmt[%u][%i]=None", s, t); break;
                                    case BCF_BT_INT8: ; BRANCH(int8_t, PRIi8); break;
                                    case BCF_BT_INT16: ; BRANCH(int16_t, PRIi16); break;
                                    case BCF_BT_INT32: ; BRANCH(int32_t, PRIi32); break;
                                    case BCF_BT_FLOAT: ; float *a = (float *)p; a = a + (s*fmt_n) + t; fprintf(stdout, "; fmt[%u][%i]==%f", s, t, *a); break;
                                                       // case BCF_BT_CHAR: ; char *tc = (char*)p; fprintf(stdout, "; fmt[%u][%i]==%s", s, t, tc); break;
                                }
                            }
                        }
                        if (fmt_type == BCF_BT_CHAR){
                            char *sp = (char *)val;
                            int t;
                            for (t = 0; t < fmt_n; t++){
                                sp = sp + (fmt_n*s) + t;
                                if (*sp == bcf_str_missing) fprintf(stdout, "; fmt[%u][%i]=%s", s, t, "MISS");
                                else if (*sp == bcf_str_vector_end) fprintf(stdout, "; fmt[%u][%i]=%s", s, t, "END");
                                else fprintf(stdout, "; fmt[%u][%i]=%s", s, t, sp);
                            }
                        }
                    }
                    fprintf(stdout, "}\n");

                    /* an equivalent way to place FMT data into string
                    kstring_t stmp = {0,0,0};
                    bcf_fmt_array(&stmp, fmt_n*nsamples, fmt_type, b->d.fmt[k].p);
                    fprintf(stdout, "FORMATTED: %s\n", stmp.s);
                    ks_free(&stmp);
                    */
                    

                    free(val);
                }
                v = v->next;
                fprintf(stdout, "\n***** *****\n");
            }
        }
    }
}


int test_vcf(char *vcffn){

    int region_set = 0;
    char *region = NULL;

    /* Read in VCF file */
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    // set region
    if ( region_set ){
        int tret = bcf_sr_set_regions(sr, region, 0);
        // int tret = bcf_sr_set_regions(sr, region, 0);
        if ( tret == -1 ){
            fprintf(stderr, "Could not set region %s\n", region);
            exit(EXIT_FAILURE);
        }
    }
    if ( bcf_sr_add_reader(sr, vcffn) == 0 ){
        fprintf(stderr, "Could not add VCF file %s\n", vcffn);
        return 1;
    }
    bcf_hdr_t *vcf_hdr = sr->readers[0].header;
    /* */

    /* Store VCF variants to GenomeVar object */
    GenomeVar *gv = vcf2gv(sr, vcf_hdr);

    print_bcf_hdr(vcf_hdr);
    
    print_var(gv, vcf_hdr);

    int nsnp;
    int nchr = n_snp(gv, &nsnp);
    fprintf(stdout, "nchr = %i and nsnp = %i\n", nchr, nsnp);

    destroy_gv(gv);

    bcf_sr_destroy(sr);

    return 0;

}


