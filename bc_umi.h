
#ifndef BC_UMI_H
#define BC_UMI_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "str_util.h"
#include "variants.h"
#include "gtf_anno.h"

enum spl {SPLICE, UNSPLICE, AMBIG, NA};
#define N_SPL 3

/* Record struct
 * A record holds information on a UMI alignment.
 */
typedef struct Rec {
    uint8_t n_feat; // number of elements in feat array
    uint8_t m_feat; // allocated size of feat array
    int32_t *feat; // array of feature indices, IDs in records->gene_ix hash table
    uint8_t *splice; // array of splice values (from enum spl)
    uint8_t n_var; // number of elements in var, base, qual arrays
    uint8_t m_var; // allocated size of var, base, qual arrays
    int32_t *var; // array of variant indices, IDs in records->var_ix hash
    uint8_t *base; // base calls: ref, alt, other (from enum alleles)
    uint8_t *qual; // quality scores: phred quality 0-40
    uint16_t n_read; // number of reads that support this rec
    struct Rec *next; // Next record in list in UMI. NULL if none.
} Rec;

/* UMI struct
 * @field recs Linked list of records
 * @field best_rec Best read/alignment. This pointer is initialized to NULL, and filled after 
 *   calling the function call_umis. If no best rec can be called, it remains NULL.
 */
typedef struct {
    struct Rec *recs;
    struct Rec *best_rec;
} UMI;

KHASH_INIT(umi, char*, UMI*, 1, kh_str_hash_func, kh_str_hash_equal);

/* Barcode struct */
typedef struct {
    khash_t(umi) *umis; // key is same memory as Records->umi_ix
} Barcode;

KHASH_INIT(bc, char*, Barcode*, 1, kh_str_hash_func, kh_str_hash_equal);

/* Records stores the BC-UMI record data */
typedef struct {
    khash_t(bc) *bc; // key is barcode string, same memory as Records->bc_ix

    /* ix objects are str_map. 
     * use function str_map_str to get ID from index
     * use function str_map_ix to get index from ID */
    str_map *bc_ix;
    str_map *umi_ix;
    str_map *gene_ix;
    str_map *var_ix;
} Records;

/*****************************
 * Records
 *****************************/

/* Initialize an empty Rec object.
 * Return pointer to object.
 * Data must be freed with destroy_rec()
 */
Rec *init_rec();

/* Initialize an empty Rec object.
 * Allocates storage for feat and var arrays
 * @param n_feat number of genes to allocate space for
 * @param n_var number of variants to allocate space for
 * @return pointer to initialized Rec object.
 *
 * returned object must be freed with destroy_rec()
 */
Rec *init_rec_s(uint8_t n_feat, uint8_t n_var);

/* Destroy rec object.
 * Returns pointer to next object. NULL if no other Rec objects are left.
 */
Rec *destroy_rec(Rec *r);

/* Return a copy of rec
 * Returns NULL if failed.
 */
Rec *copy_rec(Rec *r);

/* Initilize UMI object
 * Return pointer to object.
 * Data must be freed with destroy_umi
 */
UMI *init_umi();

/* Destroy UMI object.
 * Free allocated memory.
 */
void destroy_umi(UMI *u);

/* Initialize empty Barcode object.
 * Return pointer to object.
 * Data must be freed with destroy_bc
 */
Barcode *init_barcode();

/* Destroy Barcode object
 * Free allocated memory.
 */
void destroy_barcode(Barcode *b);

/* Initialize an empty Records object. 
 * Return pointer to object.
 * Data must be freed with destroy_records() 
 */
Records *init_records();

/* Destroy Records object, freeing allocated memory.
*/
void destroy_records(Records *recs);

/* Delete variant from a Rec object
 * @return 0 on success, -1 on error
 */
int del_var(Rec *rec, int i);

/* Add a barcode to str_map in Records recs object.
*/
int add_bc2ix(Records *recs, char *str);

/* Add UMI to str_map in Records recs object.
*/
int add_umi2ix(Records *recs, char *str);

/* Add barcodes from path to str_map in Records recs object.
*/
int add_bc2ix_file(Records *recs, const char *path);

/* Add genes from annotation to str_map in Records recs object.
 */
int add_genes2recs(Records *recs, Annotation *a);


/* compare two rec objects for equality
 *
 * compares n_feat, feat, splice, n_var, var, base, qual
 * in @p r1 and @p r2
 * @return 0 if equal, -1 if not equal
 */
int recs_cmp(Rec *r1, Rec *r2);

int add_rec2list(Rec **list, Rec *r);

/* add rec to Records
 *
 * add a record to Records object. If identical record is 
 * present, update rec->n_read and return 1. If no identical record was 
 * present, return 0. Return -1 on error.
 * Records are freed with function destroy_records.
 */
int add_rec(Records *recs, const char *bc, const char *umi, 
        uint8_t n_feat, const char **feat, const uint8_t *splice, 
        uint8_t n_var, const char **var, const uint8_t *base, const uint8_t *qual); 

/* Get barcodes with UMI given in a range [u1, u2]
 * @param recs pointer to Records object
 * @param u1 lower bound
 * @param u2 upper bound
 * @param n_bc pointer to int that is updated to number of elements in returned array.
 * @return array of char arrays. This top-level array must be freed, but not the 
 *   individual elements.
 */
char **bc_umi_cov(Records *recs, int u1, int u2, int *n_bc);

/* Call UMIs from records */
int call_umis(Records *recs);

/* filter out variants by base quality 
 * applies to best_rec, so call_umis must be called before
 *
 * @param recs pointer to Records object
 * @param min_qual phred based quality score (>=0).
 * @return 0 on success, -1 on error
 */
int fltr_bq(Records *recs, uint8_t min_qual);

/* subset to variants in vars only */
int subset_vars_rec(Records *recs, Rec *rec, str_map *vars);
int subset_vars_recs(Records *recs, str_map *vars);




/*****************************
 * I/O
 *****************************/

/* helper function to write out indexes */
int write_ix_aux(const char *fn1, const char *fn2, str_map *sm);

/* Write data in rec and *ix into fp.
 * @return Number of bytes written, or -1 on error.
 */
int write_rec_bn(Rec *rec, uint32_t bc_ix, uint32_t umi_ix, FILE *fp);

/* write records to file in binary format.
 *
 * Writes the list of records to a binary file with suffix "recs.bin".
 * It writes the UMI, barcode, and read indices to bgzipped 
 * tab-separated text files. 
 * The first column has the barcode, UMI, or read, and the second 
 * column has the index. These files have the suffix 
 * "bc_ix.tsv.gz", "umi_ix.tsv.gz", and "rd_x.tsv.gz"
 *
 * @param recs Records to write out
 * @param fn Prefix of output file names. If it is in a directory 
 *   that doesn't exist, it will be created.
 *
 * @return 0 on success, or -1 on error.
 */
int write_records_bn(Records *recs, const char *fn);

/* Read from buffer into a rec object. Decode the data in 
 * @p buf and fill into @p rec and @p *ix objects.
 *
 * @param buf pointer to uint8_t array containing rec data
 * @param rec pointer to initialized Rec object
 * @param bc_ix pointer to uint32_t to store index of barcode
 * @param umi_ix pointer to uint32_t to store index of UMI
 * @param rd_ix pointer to uint32_t to store index of read
 *
 * */
void read_rec_bn(uint8_t *buf, Rec *rec, uint32_t *bc_ix, uint32_t *umi_ix);

/* read records from file in binary format
 * Read in records from file @p fn and add to existing @p recs object.
 * The @p recs object can be empty, but it has to be initialized with 
 * function init_records(). If @p recs is NULL, return -1.
 *
 * @param fn Prefix to input file names. The full filenames are generated by 
 *   adding the appropriate suffix, such as recs.bin or bc_ix.tsv.gz.
 * @param recs pointer to initialized Record object.
 *
 * @return 0 on success, -1 on error.
 */
int read_records_bn(const char *fn, Records *recs);

Records *read_records_file(char *fn, int file_list, int verbose);

int n_recs(Records *recs);


/*****************************
 * test
 *****************************/



/* print summary of Records object
*/
void print_summary(Records *recs);
void print_best_umi(Records *recs);

#endif // BC_UMI_H

