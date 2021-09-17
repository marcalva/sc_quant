
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <errno.h>
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "str_util.h"
#include "variants.h"
#include "bc_umi.h"
#include "a_count.h"

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "samba v0.1.0: single-cell demultiplexing\n"
            "Usage:    sc_quant ac --recs records --out results\n"
            "\n"
            "Options:\n"
            "\n"
            "Required options:\n"
            "\n"
            "  -r, --rec           Prefix of records output from samba plp.\n"
            "  -R, --recs          Text file that lists multiple plp outputs, one per line.\n"
            "  -o, --out           Prefix for results output file.\n"
            "  -v, --vcf           An optional VCF file. Will output counts for the SNPs in this file.\n"
            "  -V, --verbose       Write status on output.\n"
            "\n");
    exit(exit_status);
}

int allele_count(int argc, char *argv[]){

    if (argc == 1) usage(stderr, EXIT_FAILURE);

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"rec", required_argument, NULL, 'r'}, 
        {"recs", required_argument, NULL, 'R'}, 
        {"vcf", required_argument, NULL, 'v'}, 
        {"out", required_argument, NULL, 'o'}, 
        {"verbose", no_argument, NULL, 'V'},
        {NULL, 0, NULL, 0}
    };

    /* parameters */
    int ret = EXIT_SUCCESS;
    // int dlen = 20;
    // char dt[dlen];
    int hs = 0;
    char *recfn = NULL;
    char *recsfn = NULL;
    char *vcffn = NULL;
    char *out_dltf = "samba.";
    char *out = strdup(out_dltf);
    int verbose = 0;
    GenomeVar *gv = NULL;
    /* */

    /* local variables */
    bc_ac *a = NULL;
    Records *records = NULL;
    bcf_srs_t *sr = NULL;
    /* */

    int option_index = 0;
    int cm;
    while ((cm = getopt_long_only(argc, argv, "hr:R:v:o:V", loptions, &option_index)) != -1){
        switch(cm){
            case 'h': hs = 1; break;
            case 'r': recfn = strdup(optarg); break;
            case 'R': recsfn = strdup(optarg); break; 
            case 'v': vcffn = strdup(optarg); break; 
            case 'o': free(out); out = strdup(optarg); break; 
            case 'V': verbose = 1; break;
            default: exit(EXIT_FAILURE);
        }
    }

    if ( hs ) usage(stdout, EXIT_SUCCESS);

    /* Check files */
    if ( recfn == NULL && recsfn == NULL ){
        err_msg(1, 0, "allele_count: one of --rec or --recs file must be given");
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if ( recfn != NULL && recsfn != NULL ){
        err_msg(1, 0, "allele_count: only one of '--rec' or '--recs' must be given");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* print options set */
    if ( verbose ){
                   fprintf(stdout, "options set:\n");
        if (recfn){fprintf(stdout, "    records:         %s\n", recfn);}
              else{fprintf(stdout, "    records:         %s\n", recsfn);}
        if (vcffn){fprintf(stdout, "    VCF:             %s\n", vcffn);}
                   fprintf(stdout, "    out:             %s\n", out);
        fflush(stdout);
    }
    /* */


    /* Instantiate and read in records */
    if (recfn) records = read_records_file(recfn, 0, verbose);
    else records = read_records_file(recsfn, 1, verbose);
    if (records == NULL){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    // call UMIs
    if (verbose) fprintf(stdout, "filtering UMI reads\n");
    fflush(stdout);
    call_umis(records);

    /* Store VCF variants  */
    if (vcffn != NULL){
        if (verbose) fprintf(stdout, "reading VCF file %s\n", vcffn);
        fflush(stdout);
        sr = bcf_sr_init();
        bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
        if ( bcf_sr_add_reader(sr, vcffn) == 0 ){
            err_msg(1, 0, "allele_count: could not read VCF file %s", vcffn);
            ret = EXIT_FAILURE;
            goto cleanup;
        }
        bcf_hdr_t *vcf_hdr = sr->readers[0].header;
        gv = vcf2gv(sr, vcf_hdr);
    }
    /* */

    a = init_bc_ac();

    str_map *vars = NULL;
    if (gv) vars = gv->var_ix;
    else vars = records->var_ix;
    str_map *bcs = records->bc_ix;

    if (bc_ac_add_var_map(a, vars) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (bc_ac_add_bc_map(a, bcs) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (bc_ac_count(a, records) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (bc_ac_write(a, out) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose) fprintf(stdout, "finished\n");

cleanup:
    if (recfn) free(recfn);
    if (recsfn) free(recsfn);
    if (vcffn) free(vcffn);
    if (out) free(out);
    if (gv) destroy_gv(gv);
    if (records) destroy_records(records);
    if (sr) bcf_sr_destroy(sr);
    if (a) destroy_bc_ac(a);

    return(ret);
}


