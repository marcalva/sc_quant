
#ifndef G_COUNT_H
#define G_COUNT_H

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "str_util.h"
#include "bc_umi.h"

/* hold counts for a gene */
typedef struct gc_node {
    int ix; // gene index from gene_ix
    uint32_t counts[N_SPL];
    struct gc_node *next; // pointer to next gc_node, NULL if none
} gc_node;

// value in hash is a dummy head node with index=-1
KHASH_INIT(gc, char*, gc_node, 1, kh_str_hash_func, kh_str_hash_equal);

/* hold allele counts for barcodes */
typedef struct {
    str_map *gene_ix; // variants
    str_map *bc_ix; // barcodes
    khash_t(gc) *gcs; // counts per gene
                      // key is barcode ID string. Key is same memory as bc_ix. 
                      // val is gc_node
    int n_nz; // store number of non-zero elements
} bc_gc;

gc_node *init_gc_node();

// return next
gc_node *destroy_gc_node(gc_node *n);

void gc_node_set_zero(gc_node *n);

bc_gc *init_bc_gc();

void destroy_bc_gc(bc_gc *a);

int bc_gc_add_gene_map(bc_gc *a, str_map *sm);
int bc_gc_add_bc_map(bc_gc *a, str_map *sm);

int bc_gc_count(bc_gc *a, Records *recs);

int bc_gc_write(bc_gc *a, Annotation *anno, char *fn);

#endif // G_COUNT_H
