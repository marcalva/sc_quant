
#ifndef SAM_READ_H
#define SAM_READ_H

#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/hts.h"

// Overlap tid
int ovrlp_tid(sam_hdr_t *sam_hdr, bcf_hdr_t *bcf_hdr, int **t1, int **t2);

// Return the base that overlaps the BAM record. If no overlap, return 'N'
// base_ref is the tid chromosome that corresponds to the BAM record
int get_ovrlp_base(const bam1_t *b, int32_t base_ref, hts_pos_t base_pos, 
        char *base, uint8_t *qual);

char *get_tag(const bam1_t *b, char tag[2]);

int64_t get_itag(const bam1_t *b, char tag[2]);

void print_bam1_t(const bam1_t *b);

#endif // SAM_READ_H
