
#ifndef VARIANTS_H
#define VARIANTS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "str_util.h"
#include "bins.h"

/* 
#define REF 0
#define ALT 1
#define OTHER 2
#define NA_ALLELE 3
*/
enum alleles {REF, ALT, OTHER, NA_ALLELE};
#define N_ALLELE 4

// #define MAX_BIN (((1<<18)-1)/7)

typedef struct Var { 
    bcf1_t *b;
    char *vid; // points to same array as variant in hash GenomeVar.var_ix
               // doesn't need to be freed, freed with destroy_gv call.
    struct Var *next;
} Var;

KHASH_INIT(var, char*, Var*, 1, kh_str_hash_func, kh_str_hash_equal);

typedef struct chrmVar {
    Var *bins[MAX_BIN];
    uint16_t vars_n[MAX_BIN];
} ChrmVar;

typedef struct {
    str_map *chrm_ix; // chromosome index
    str_map *var_ix; // variant index
    str_map *var_bin; // variant bin
    khash_t(var) *vars; // ID to Var hash table. Access Var object through ID. 
    ChrmVar **chrms; // array to array of chromosomes
    int chrms_m; // max number of elements
} GenomeVar;

// Functions

/* Return a variant ID for a BCF line.
 *
 * Since the ID column may be missing, a variant can be identified using 
 * CHR, POS, ID, REF, and ALT. The variant ID is tab-delimited.
 *
 * This returns a pointer to char that contains the NULL terminated string. 
 * The function allocates the memory, and it is the caller's job to free 
 * the memory.
 *
 * @param h pointer to VCF header object.
 * @param b pointer to bcf1_t object.
 * @param delim delimiter between fields
 *
 * @return pointer to char that contains the NULL terminated variant ID string.
 * variant ID is of the format CHRM\tPOS\tID\tREF\tALT. If there are multiple 
 * ALT alleles, they are separated by commas.
 */
char *var_id(const bcf_hdr_t *h, bcf1_t *b, char delim);

/* get SNP allele
 *
 * return whether the given base is REF, ALT, or OTHER in the VCF/BCF record.
 * @param b bcf record with ref and alt allele to test against
 * @param base base allele to test
 * @return REF if base == b->ref, ALT if base == b->alt, OTHER if neither
 */
uint8_t base_ref_alt(bcf1_t *b, char base);

int get_var_len(bcf1_t *b);

int get_var_bin(bcf1_t *b);

Var *init_var();

GenomeVar *init_genomevar();

/* initialize chrmVar object
 * return NULL if memory wasn't allocated
 */ 
ChrmVar *init_ChrmVar();

int add_var(GenomeVar *gv, bcf1_t *b, const bcf_hdr_t *hdr);

GenomeVar *vcf2gv(bcf_srs_t *sr, bcf_hdr_t *vcf_hdr);

void destroy_gv(GenomeVar *gv);

// destroy bcf1_t records to preserve memory
void free_bcf1_t(GenomeVar *gv);

/* check if fmt tag is valid for getting allele probabilities from a VCF line
 *
 * to be valid, the variant must be 
 *  bi-allelic
 *  fmt must be present in the header
 *  fmt must be present in the VCF line
 *  is monoploid or diploid
 *  has at least one non-missing sample
 *
 * @return 
 *   -2 if error 
 *   -1 if all samples are missing 
 *   0 if invalid 
 *   1 if valid
 */
int is_gt_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b);
int is_gp_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b);

/* get dosage for bi-allelic SNPs 
 *
 * This returns the dose (fraction of alternate alleles) for all samples 
 * in the bcf line. If missing, the dose is equal to -1 for that sample.
 * If the samples were subsetted with bcf_hdr_set_samples, 
 * then only those samples will be returned.
 *
 * 
 * @param vcf_hdr VCF header
 * @param b VCF line
 * @param extra add this many elements to the allocated array
 * @param tag_id the numeric ID of the genotype tag (GT or GP)
 * @return double array of dose, or NULL on failure.
 * 
 * The returned array must be freed by caller.
 */
float *bcf1_ap_gt(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra);
float *bcf1_ap_gp(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra);

/* Return dose in a two-dimensional array.
 *
 * @param gv pointer to GenomeVar object.
 * @param vcf_hdr pointer to VCF header of the VCF lines in @p gv.
 * @param ids array of character arrays containing the variant IDs to return dose for.
 *   The order of these variants will be preserved in @p dose.
 * @param ni length of @p ids array.
 * @param field one of "GT" or "GP".
 * @return float matrix, NULL on error.
 */
float **ap_array_gt(GenomeVar *gv, bcf_hdr_t *vcf_hdr, char **ids, int ni, char *field);

/* Get variants that overlap given region [beg, end).
 *
 * @param gv GenomeVar object to retrieve variants from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param vars pointer to array of Var pointers.
 * @param vars_m Pointer to integer that contains the length of the vars array.
 *
 * @return Number of overlapping variants in @p vars.
 *  -1 if the reference is not found in GenomeVar. -2 on error
 *
 * The @p vars array must be an allocated array, or NULL. The @p vars_len 
 * gives the length of the @p vars array. The function can increase the size of 
 * @p vars by calling realloc, and will change @p vars_len to reflect the udpated 
 * size.
 */
int vars_from_region(GenomeVar *gv, const char* ref, hts_pos_t beg, 
        hts_pos_t end, Var ***vars, int *vars_m);

int n_snp(GenomeVar *gv, int *n_snp);

#endif // VARIANTS_H

