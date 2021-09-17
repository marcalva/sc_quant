
#include "overlap.h"
#include "sam_read.h"
#include "bc_umi.h"
#include "str_util.h"
#include "bins.h"
#include "htslib/khash.h"
#include <string.h>

int overlap_bam1_feats(const sam_hdr_t *h, bam1_t *b, const Annotation *a, uint8_t *n_feat, uint8_t *m_feat, 
        char ***feat, uint8_t **splice){

    int32_t tid = b->core.tid;
    const char *ref = sam_hdr_tid2name(h, (int)tid);
    
    hts_pos_t b_beg = b->core.pos;
    hts_pos_t b_end = bam_endpos(b);
    if (b_end == b_beg + 1) return(-1); // unmapped;

    char strand = '+';
    if (bam_is_rev(b)) strand = '-';

    // get features that contain bam1_t read b.
    int genes_len = 1;
    Gene **genes = malloc(genes_len * sizeof(Gene *)); // must be freed.

    if (genes == NULL)
        return err_msg(-1, 0, "overlap_bam1_feats: %s", strerror(errno));

    int ngenes = feats_from_region(a, ref, b_beg, b_end, 1, strand, &genes, &genes_len);
    if (ngenes < 0) return -1;

    // reallocate memory for feat and splice
    if (ngenes > *m_feat){
        *m_feat = ngenes;
        *feat = realloc(*feat, *m_feat * sizeof(char *));
        if (*feat == NULL)
            return err_msg(-1, 0, "overlap_bam1_feats: %s", strerror(errno));
        *splice = realloc(*splice, *m_feat * sizeof(uint8_t));
        if (*splice == NULL)
            return err_msg(-1, 0, "overlap_bam1_feats: %s", strerror(errno));
    }

    // add gene ID and splice status
    int i, j;
    for (i = 0; i < ngenes; ++i){
        if (genes[i] == NULL)
            return err_msg(-1, 0, "overlap_bam1_feats: genes cannot be NULL");

        char *id = genes[i]->id;
        if (id == NULL)
            return err_msg(-1, 0, "overlap_bam1_feats: gene name not found");

        // sorted insert (by string)
        for (j = 0; j < i; j++){
            if (strcmp(id, (*feat)[j]) < 0){ // if gene[i] < gene[j]
                int k;
                for (k = i; k > j;--k){
                    (*feat)[k] = (*feat)[k - 1];
                    (*splice)[k] = (*splice)[k - 1];
                }
                break;
            }
        }

        (*feat)[j] = id;
        int sret = 0;
        (*splice)[j] = bam1_spliced(b, genes[i], &sret);
        if (sret < 0)
            return err_msg(-1, 0, "overlap_bam1_feats: splice failed");
    }
    *n_feat = (uint32_t)ngenes;

    free(genes);

    return(0);
}

/* 
 * Detect whether a bam record is spliced. 
 * The record must overlap a splice junction and 
 * only overlap exons.
 *
 * A record overlaps an annotated 
 * splice junction if the beginning and end of any 
 * 'N' CIGAR tag is equal to the base after the end 
 * of one exon and the base before the start of the 
 * next exon, respectively.
 *
 * If a record only overlaps exons but doesn't 
 * overlap a splice junction, it is ambiguous.
 *
 * If a record has any overlap outside of exons, 
 * it is considered intronic.
 * 
 */
uint8_t bam1_spliced(bam1_t *b, Gene *g, int *ret){
    *ret = 0;

    if (g == NULL){
        *ret = -1;
        return 0;
    }

    int n_iso = g->isoforms_n;

    if (n_iso <= 0 || g->isoforms == NULL){
        err_msg(-1, 0, "bam1_spliced: %s has no isoforms", g->id);
        *ret = -1;
        return 0;
    }

    double *ro = calloc(n_iso, sizeof(double)); // length of read overlap
    double *rl = calloc(n_iso, sizeof(double)); // length of read
    double *pe = calloc(n_iso, sizeof(double)); // ro / rl
    uint8_t *sj = calloc(n_iso, sizeof(uint8_t)); // splice junctions

    if (ro == NULL || rl == NULL || pe == NULL || sj == NULL){
        return err_msg(-1, 0, "bam1_spliced: %s", strerror(errno));
        *ret = -1;
        return 0;
    }

    int i;
    for (i = 0; i < n_iso; ++i){
        ro[i] = 0;
        rl[i] = 0;
        sj[i] = 0;
    }

    // loop through isoforms. Update struct to make this faster
    khint_t k;
    int iso_i = 0; // isoform index
    for (k = kh_begin(g->isoforms); k != kh_end(g->isoforms); ++k){
        if (!kh_exist(g->isoforms, k)) continue;
        Isoform *iso = kh_val(g->isoforms, k);

        if (iso == NULL){
            err_msg(-1, 0, "bam1_spliced: isoform not stored properly");
            *ret = -1;
            return 0;
        }

        Exon *exon = iso->exons;
        if (exon == NULL) continue;

        // pos is leftmost base of first CIGAR op that consumes reference.
        char strand = '+';
        if (bam_is_rev(b)) strand = '-';

        hts_pos_t pos = b->core.pos;
        hts_pos_t beg = pos, end;

        uint32_t *cigar_raw = bam_get_cigar(b);
        uint32_t n_cigar = b->core.n_cigar;

        int i;
        for (i = 0; i < n_cigar; ++i){ // for each CIGAR op
            uint32_t op = bam_cigar_op(cigar_raw[i]); // lower 4 bits is cigar length
            uint32_t len = bam_cigar_oplen(cigar_raw[i]); // higher 24 bits is length
            end = beg + len; // b and e are coordinates of CIGAR. Not inclusive [beg,end)

            int cr = bam_cigar_type(op)&2; // consumes reference
            int cq = bam_cigar_type(op)&1; // consumes query

            /*********************************************************
             * check splice junction N
             *********************************************************/

            if (op == BAM_CREF_SKIP){
                exon = iso->exons;
                while (exon != NULL && exon->next != NULL && exon->end <= beg){
                    if (exon->end == beg && exon->next->beg == end){
                        sj[iso_i] = 1;
                        break;
                    }
                    exon = exon->next;
                }
            }

            /*********************************************************
             * check exon overlap
             * only check exon-query overlap with query 
             *   segments that consume reference and query
             *********************************************************/

            if ( cr == 0 || cq == 0 ){
                beg = end;
                continue;
            }

            rl[iso_i] += (double)len;

            // skip to overlapping exons
            exon = iso->exons;
            while (exon != NULL && exon->end <= beg)
                exon = exon->next;

            // get overlap of exons
            int ovrlp = 0;
            while (exon != NULL && exon->beg < end){
                int tmp_ovrlp = bp_overlap(beg, end, strand, exon->beg, exon->end, g->strand);
                if (tmp_ovrlp < 0){
                    *ret = -1;
                    return 0;
                }
                ovrlp += tmp_ovrlp;
                exon = exon->next;
            }
            
            // save overlap of CIGAR segment
            ro[iso_i] += (double)ovrlp;

            beg = end;
        }
        pe[iso_i] = ro[iso_i] / rl[iso_i];
        iso_i++;
    }

    if (iso_i != n_iso){
        return err_msg(-1, 0, "bam1_spliced: number of isoforms are inconsistent");
        *ret = -1;
        return 0;
    }

    /*********************************************************
     * select isoforms by precedence
     *   sort by sj=1 > sj=0, then by pe
     *********************************************************/

    int max_i = 0;
    for (i = 1; i < n_iso; ++i){
        if (sj[i] > sj[max_i]){
            max_i = i;
            continue;
        }
        if (pe[i] > pe[max_i]){
            max_i = i;
            continue;
        }
    }

    /*********************************************************
     * if sj = 1 and pe = 1 SPLICE
     * if sj = 0 and pe = 1 UNSPLICE
     * else AMBIG
     *********************************************************/

    uint8_t spl_stat;
    if (sj[max_i] == 1 && pe[max_i] == 1) spl_stat = SPLICE;
    else if (sj[max_i] == 0 && pe[max_i] == 1) spl_stat = AMBIG;
    else spl_stat = UNSPLICE;

    free(pe);
    free(ro);
    free(rl);
    free(sj);

    return(spl_stat);
}

int overlap_bam1_vars(const sam_hdr_t *h, bam1_t *b, GenomeVar *gv, uint8_t *n_var, uint8_t *m_var,
        char ***var, uint8_t **base, uint8_t **qual, uint8_t min_qual){

    int32_t tid = b->core.tid;
    const char *ref = sam_hdr_tid2name(h, (int)tid);
    
    hts_pos_t b_beg = b->core.pos;
    hts_pos_t b_end = bam_endpos(b);
    if (b_end == b_beg + 1) return(-1); // unmapped;

    int var_q_n = 0;
    int var_q_m = 4;
    Var **var_q = calloc(var_q_m, sizeof(Var *));
    int var_cig_n = 0;
    int var_cig_m = 4;
    Var **var_cig = calloc(var_cig_m, sizeof(Var *));
    if (var_q == NULL || var_cig == NULL){
        err_msg(-1, 0, "overlap_bam1_vars: %s", strerror(errno));
        return -1;
    }

    uint32_t *cigar_raw = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;
    int i;
    for (i = 0; i < n_cigar; ++i){
        // get CIGAR op
        uint32_t op = bam_cigar_op(cigar_raw[i]);
        uint32_t len = bam_cigar_oplen(cigar_raw[i]);
        int cr = bam_cigar_type(op)&2; // consumes reference
        int cq = bam_cigar_type(op)&1; // consumes query
        if (cr) b_end = b_beg + len;
        else continue;
        if ( !cq ) {
            b_beg = b_end;
            continue;
        }

        // get overlapping variants in CIGAR
        var_cig_n = vars_from_region(gv, ref, b_beg, b_end, &var_cig, &var_cig_m);
        if (var_cig_n < 0) return(-1);

        var_q_n += var_cig_n;
        if (var_q_n > var_q_m){
            var_q_m = var_q_n << 1;
            var_q = realloc(var_q, var_q_m * sizeof(Var *));
        }
        Var **tmp = var_q + (var_q_n - var_cig_n);
        memcpy(tmp, var_cig, var_cig_n * sizeof(Var *));

        b_beg = b_end;
    }
    free(var_cig);

    // fprintf(stdout, "Found %i feats from region. genes_len=%i\n", ngenes, genes_len);

    // reallocate memory for var, base, and qual
    if (var_q_n > *m_var){
        *m_var = var_q_n;
        *var = realloc(*var, *m_var * sizeof(char *));
        *base = realloc(*base, *m_var * sizeof(uint8_t));
        *qual = realloc(*qual, *m_var * sizeof(uint8_t));
        if (*var == NULL || *base == NULL || *qual == NULL){
            err_msg(-1, 0, "overlap_bam1_vars: %s", strerror(errno));
            return(-1);
        }
    }

    // fill var, base, qual
    for (i = 0; i < var_q_n; ++i){
        char t_base = 'N';
        uint8_t t_qual = 0;
        if ( get_ovrlp_base(b, tid, var_q[i]->b->pos, &t_base, &t_qual) < 0 )
            return -1;
        
        if (t_qual < min_qual){
            int k;
            for (k = i; k < var_q_n - 1; ++k){
                var_q[k] = var_q[k+1];
            }
            var_q[k] = NULL;
            i--;
            var_q_n--;
            continue;
        }
        char *id = var_q[i]->vid;
        int j = 0;
        for (; j < i; j++){
            if (strcmp(id, (*var)[j]) < 0){
                int k;
                for (k = i; k > j;--k){
                    (*var)[k] = (*var)[k-1];
                    (*base)[k] = (*base)[k-1];
                    (*qual)[k] = (*qual)[k-1];
                }
                break;
            }
        }
        (*var)[j] = id;
        (*base)[j] = base_ref_alt(var_q[i]->b, t_base);
        (*qual)[j] = t_qual;
    }

    // debugging
    int printout = 0;
    if (var_q_n > 1 && printout){
        fprintf(stdout, "var_q_n: %i\n", var_q_n);
        int rstart = b->core.pos;
        int rend = rstart;
        for (i = 0; i < n_cigar; ++i){
            // get CIGAR op
            uint32_t op = bam_cigar_op(cigar_raw[i]);
            uint32_t len = bam_cigar_oplen(cigar_raw[i]);
            char cig = bam_cigar_opchr(op);
            int cr = bam_cigar_type(op)&2; // consumes reference
            int cq = bam_cigar_type(op)&1; // consumes query
            if ( cr ) rend = rstart + len; // b and e are coordinates of CIGAR
            fprintf(stdout, "\t%c(%i:CR%i:CQ%i):%i-%i", cig, len, cr, cq, rstart, rend);
            rstart = rend;
        }
        fprintf(stdout, "\n");

        for (i = 0; i < var_q_n; ++i){
            fprintf(stdout, "\t%s\tb=%i\tbq=%u\n", var_q[i]->vid, (*base)[i], (*qual)[i]);
        }
    }

    *n_var = (uint8_t)var_q_n;

    free(var_q);
    
    return 0;
}

int bp_overlap(int64_t a1, int64_t a2, char a_strand, int64_t b1, int64_t b2, char b_strand){
    if ( (a_strand == '+') && (b_strand == '-') )
        return 0;
    if ( (a_strand == '-') && (b_strand == '+') )
        return 0;

    if (a2 <= a1 || b2 <= b1){
        return err_msg(-1, 0, "bp_overlap: incorrect parameters A=[%i, %i) B=[%i, %i)", 
                a1, a2, b1, b2);
    }

    if (a2 <= b1)
        return 0;
    if (b2 <= a1)
        return 0;

    int c1, c2;
    if (a1 <= b1)
        c1 = b1;
    else
        c1 = a1;

    if (a2 <= b2)
        c2 = a2;
    else
        c2 = b2;

    return(c2 - c1);
}

