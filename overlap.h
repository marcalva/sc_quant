
#ifndef OVERLAP_H
#define OVERLAP_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "gtf_anno.h"
#include "variants.h"
#include "bc_umi.h"

/* Overlap a bam1_t read with features given Annotation Find the features that
 * the read fully intersects, i.e. is fully contained.  For each feature, give
 * whether the read is spliced, unspliced, or ambiguous.
 *
 * The features and splice status are returned as integer indices in an array.
 * The number of elements in the array are stored in @p n_feat, while the
 * allocated size of the array is stored in @p m_feat.
 *
 * @param h pointer to header for bam record @p b.  @param b pointer to bam1_t
 * object for overlap.  @param a pointer to Annotation object with features.
 * @param n_feat pointer to integer giving the number of features that overlap.
 * @param m_feat pointer to integer giving the allocated size of the updated
 * arrays @p feat and @p splice.  @param feat pointer to array of char arrays
 * from @p a containing variant IDs.  @param splice integer array containing
 * splice status for each overlapping feature in @p feat.  @return 0 if
 * success, -1 if fail
 */
int overlap_bam1_feats(const sam_hdr_t *h, bam1_t *b, const Annotation *a,
        uint8_t *n_feat, uint8_t *m_feat, char ***feat, uint8_t **splice);

/* Calculate whether a read is spliced, unspliced, or ambiguous.
 * 
 * @param b pointer to bam record.
 * @param gene pointer to gene.
 * @param ret set to 0 on success, -1 on error
 * @return One of SPLICE UNSPLICE AMBIG
 */
uint8_t bam1_spliced(bam1_t *b, Gene *g, int *ret);

/* Overlap a bam1_t read with variants.
 *
 * @param h pointer to header for bam record @p b.
 * @param b pointer to bam1_t object for overlap.
 * @param gv pointer to GenomeVar object with variants for overlap.
 * @param n_var pointer to integer giving the number of variants that overlap.
 * @param m_var pointer to integer giving the allocated size of the updated @p var @p base and @p qual arrays.
 * @param var pointer to array of char arrays containing the IDs of the overlapping variants.
 * @param base pointer to integer array containing the base calls for each variant.
 * @param qual pointer to integer array containing the quality of the base calls for each variant.
 * @param min_qual only return overlapping variants that have a base quality score of at least this value.
 * 
 * @return 0 if success, -1 if fail.
 *
 */
int overlap_bam1_vars(const sam_hdr_t *h, bam1_t *b, GenomeVar *gv, uint8_t *n_var, uint8_t *m_var, 
        char ***var, uint8_t **base, uint8_t **qual, uint8_t min_qual);

/* Return the number of bases that overlap between [a1,a2) and [b1,b2) features. 
 * a1 must be less than a2, and b1 must be less than b2.
 * If a_strand and b_strand are opposite strands ('+' or '-'), return 0 overlap.
 * Otherwise, the strand is ignored.
 * the a and b regions are half open, with a2 and b2 open.
 *
 * @return Number of base pairs that overlap, or -1 on error.
 * */
int bp_overlap(int64_t a1, int64_t a2, char a_strand, int64_t b1, int64_t b2, char b_strand);

#endif // OVERLAP_H
