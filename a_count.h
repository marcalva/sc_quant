
#ifndef A_COUNT_H
#define A_COUNT_H

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "str_util.h"
#include "variants.h"
#include "bc_umi.h"

/* hold allele counts for a variant */
typedef struct ac_node {
    int ix; // variant index from var_ix
    uint32_t counts[N_ALLELE]; // ref, alt, other, na counts
    struct ac_node *next; // pointer to next ac_node, NULL if none
} ac_node;

// value in hash is a dummy head node with index=-1
KHASH_INIT(ac, char*, ac_node, 1, kh_str_hash_func, kh_str_hash_equal);

/* hold allele counts for barcodes */
typedef struct {
    str_map *var_ix; // variants
    str_map *bc_ix; // barcodes
    khash_t(ac) *acs; // allele counts per SNP
                      // key is barcode ID string. Key is same memory as bc_ix. 
                      // val is ac_node
    int n_nz; // store number of non-zero elements
} bc_ac;

ac_node *init_ac_node();

// return next
ac_node *destroy_ac_node(ac_node *n);

void ac_node_set_zero(ac_node *n);

bc_ac *init_bc_ac();

void destroy_bc_ac(bc_ac *a);

int bc_ac_add_var_map(bc_ac *a, str_map *sm);
int bc_ac_add_bc_map(bc_ac *a, str_map *sm);

int bc_ac_count(bc_ac *a, Records *recs);

int bc_ac_write(bc_ac *a, char *fn);

#endif // A_COUNT_H
