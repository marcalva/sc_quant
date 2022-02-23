
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
#include "sam_read.h"
#include "overlap.h"
#include "a_count.h"

void print_status_ac(const char* k, int i, const char *chr, int pos){
    time_t now;
    time(&now);
    char dt[20];
    get_time(dt, 20);
    fprintf(stdout, "%s: processed %i %s through %s:%i\n", dt, i, k, chr, pos);
    fflush(stdout);
}

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "sc_quant v0.1.0: Generate counts from single-cell experiments\n"
            "Usage:    sc_quant ac [options] --bam bamfile --vcf vcffile --out outfile \n"
            "\n"
            "Options:\n"
            "\n"
            "Required options:\n"
            "\n"
            "  -b, --bam           Indexed BAM file.\n"
            "  -v, --vcf           VCF file. Will output counts for the SNPs in this file.\n"
            "  -o, --out           Output file prefix [ac.].\n"
            "  -B, --bc-tag        BAM tag for barcode [CB].\n"
            "  -U, --umi-tag       BAM tag for UMI [UB].\n"
            "  -H, --nh-tag        BAM tag for the number of alignments of a read [NH].\n"
            "  -m, --max-nh        Only process reads with a maximum of this many alignments [10]. Set to 0 to ignore.\n"
            "  -p, --min-phredq    Minimum base phred quality score in read [30].\n"
            "  -c, --barcodes      File containing list of barcode IDs, one per line.\n"
            "  -r, --region        Region (hts format), for example 'chr21,chr21:10-,chr21-10-20'.\n"
            "  -V, --verbose       Write status on output.\n"
            "\n");
    exit(exit_status);
}

int allele_count(int argc, char *argv[]){

    if (argc == 1) usage(stderr, EXIT_FAILURE);

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"bam", required_argument, NULL, 'b'}, 
        {"vcf", required_argument, NULL, 'v'}, 
        {"out", required_argument, NULL, 'o'}, 
        {"bc-tag", required_argument, NULL, 'B'},
        {"umi-tag", required_argument, NULL, 'U'},
        {"nh-tag", required_argument, NULL, 'H'}, 
        {"max-nh", required_argument, NULL, 'm'}, 
        {"min-phredq", required_argument, NULL, 'p'},
        {"barcodes", required_argument, NULL, 'c'},
        {"region", required_argument, NULL, 'r'},
        {"verbose", no_argument, NULL, 'V'},
        {NULL, 0, NULL, 0}
    };

    /* parameters */
    int ret = EXIT_SUCCESS;
    int hs = 0;
    int verbose = 0;

    /* output prefix */
    char *out_dflt = "ac.";
    char *out = strdup(out_dflt);

    /* variant data */
    char *vcffn = NULL;
    bcf_srs_t *sr = NULL;
    GenomeVar *gv = NULL;

    /* bam file */
    char *bamfn = NULL;
    samFile *bam = NULL;
    sam_hdr_t *bam_hdr = NULL;
    hts_idx_t *bam_idx = NULL;
    hts_itr_t *bam_itr = NULL;
    bam1_t *bam_r = NULL;

    /* bam tags */
    char bc_tag[] = "CB";
    char umi_tag[] = "UB";
    char nh_tag[] = "NH";

    /* read filters */
    int max_nh = 10; // max number mappings per read
    int min_phred = 30; // min phred quality score
    char *bc_fn = NULL;

    /* genomic region to iterate over */
    char *region = ".";
    int region_set = 0;

    /* local variables */
    Records *records = NULL;
    bc_ac *a = NULL;

    uint8_t n_feat = 0, m_feat = 4;
    char **feat = (char **)calloc(m_feat, sizeof(char *));
    uint8_t *splice = (uint8_t *)calloc(m_feat, sizeof(uint8_t));;
    uint8_t n_var = 0, m_var = 4;
    char **var = (char **)calloc(m_var, sizeof(char *));
    uint8_t *base = (uint8_t *)calloc(m_var, sizeof(uint8_t));
    uint8_t *qual = (uint8_t *)calloc(m_var, sizeof(uint8_t));
    char *p_end = NULL;
    int option_index = 0;
    int cm;
    while ((cm = getopt_long_only(argc, argv, "hb:v:o:B:U:H:m:p:c:r:V", 
                    loptions, &option_index)) != -1){
        switch(cm){
            case 'h': hs = 1; break;
            case 'b': bamfn = strdup(optarg); break;
            case 'v': vcffn = strdup(optarg); break; 
            case 'o': free(out); out = strdup(optarg); break; 
            case 'B': if (strlen(optarg) != 2){
                          ret = err_msg(1, 0, "ac: bc-tag must have a length of 2");
                          goto cleanup;
                      }
                      strncpy(bc_tag, optarg, 2);
                      break;
            case 'U': if (strlen(optarg) != 2){
                          ret = err_msg(1, 0, "ac: umi-tag must have a length of 2");
                          goto cleanup;
                      }
                      strncpy(umi_tag, optarg, 2);
                      break;
            case 'H': if (strlen(optarg) != 2){
                          ret = err_msg(1, 0, "ac: nh-tag must have a length of 2");
                          goto cleanup;
                      }
                      strncpy(nh_tag, optarg, 2);
                      break;

            case 'm': errno = 0;
                      max_nh = (int)strtol(optarg, &p_end, 10);
                      if (max_nh == 0 && errno > 0){
                          return err_msg(EXIT_FAILURE, 0, 
                                  "ac: could not convert --max-nh %s to int: %s", 
                                  optarg, strerror(errno));
                      }
                      break;
            case 'p':
                      errno = 0;
                      min_phred = (int)strtol(optarg, &p_end, 10);
                      if (min_phred == 0 && errno > 0){
                          return err_msg(EXIT_FAILURE, 0, 
                                  "ac: could not convert --min-phred %s to int: %s", 
                                  optarg, strerror(errno));
                      }
                      break;
            case 'c': bc_fn = strdup(optarg); break;
            case 'r': region = strdup(optarg); region_set = 1; break;
            case 'V': verbose = 1; break;
            default: exit(EXIT_FAILURE);
        }
    }

    if ( hs ) usage(stdout, EXIT_SUCCESS);

    /* Check files */
    if (bamfn == NULL){
        ret = err_msg(1, 0, "ac: BAM file must be given with --bam option");
        goto cleanup;
    }
    if (vcffn == NULL){
        ret = err_msg(1, 0, "ac: VCF file must be given with --vcf option");
        goto cleanup;
    }

    /* print options set */
    if ( verbose ){
                   fprintf(stdout, "options set:\n");
                   fprintf(stdout, "    BAM:             %s\n", bamfn);
                   fprintf(stdout, "    VCF:             %s\n", vcffn);
                   fprintf(stdout, "    out:             %s\n", out);
                   fprintf(stdout, "    bc-tag:          %s\n", bc_tag);
                   fprintf(stdout, "    umi-tag:         %s\n", umi_tag);
                   fprintf(stdout, "    nh-tag:          %s\n", nh_tag);
                   fprintf(stdout, "    max-nh:          %i\n", max_nh);
                   fprintf(stdout, "    min-phredq:      %i\n", min_phred);
                   fprintf(stdout, "    barcodes:        %s\n", bc_fn);
   if (region_set){fprintf(stdout, "    region:          %s\n", region);}
        fflush(stdout);
    }
    /* */

    /* Instantiate records */
    records = init_records();

    int n_reads = 0;
    int n_bc = 0;
    /* */

    /* set barcodes */
    if ( bc_fn != NULL ){
        n_bc = add_bc2ix_file(records, bc_fn);
        if (n_bc < 0){
            ret = err_msg(1, 0, "could not read barcodes from %s\n", bc_fn);
            goto cleanup;
        }
        if ( verbose ){
            fprintf(stdout, "read %i barcodes from %s\n", n_bc, bc_fn);
            fflush(stdout);
        }
    }
    /* */

    /* Store VCF variants  */
    if (verbose){
        fprintf(stdout, "reading VCF file %s\n", vcffn);
        fflush(stdout);
    }
    sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    // set region
    if ( region_set ){
        int tret = bcf_sr_set_regions(sr, region, 0);
        if ( tret == -1 ){
            ret = err_msg(EXIT_FAILURE, 0, "could not set region %s in VCF file", region);
            goto cleanup;
        }
    }
    if ( bcf_sr_add_reader(sr, vcffn) == 0 ){
        ret = err_msg(EXIT_FAILURE, 0, "could not read VCF file %s", vcffn);
        goto cleanup;
    }
    bcf_hdr_t *vcf_hdr = sr->readers[0].header;

    gv = vcf2gv(sr, vcf_hdr);
    /* */

    /* Open BAM file */
    bam = sam_open(bamfn, "r");
    if (bam == NULL){
        ret = err_msg(1, 0, "could not open BAM file %s", bamfn);
        goto cleanup;
    }
    bam_hdr = sam_hdr_read(bam);
    if (bam_hdr == NULL){
        ret = err_msg(1, 0, "could not read header for BAM file %s", bamfn);
        goto cleanup;
    }
    bam_idx = sam_index_load(bam, bamfn);
    if (bam_idx == NULL){
        ret = err_msg(1, 0, "could not load index for BAM file %s", bamfn);
        goto cleanup;
    }
    /* */

    /* Set BAM iterator to region */
    bam_itr = sam_itr_querys(bam_idx, bam_hdr, region);
    if (bam_itr == NULL){
        ret = err_msg(EXIT_FAILURE, 0, "could not set region %s in BAM file", region);
        goto cleanup;
    }
    /* */
            
    // Loop over BAM records/
    if (verbose){
        fprintf(stdout, "processing alignments in %s\n", bamfn); fflush(stdout);
    }
    int iter_ret;
    bam_r = bam_init1();
    while ( (iter_ret = sam_itr_next(bam, bam_itr, bam_r)) >= 0){
        if ( (verbose) && (n_reads % (int)1e6 == 0)){
            const char *chr = sam_hdr_tid2name(bam_hdr, bam_r->core.tid);
            print_status_ac("alignments", n_reads, chr, (int)(bam_r->core.pos+1));
        }
        n_reads++;

        // print_bam1_t(bam_r);
        const char *qname = bam_get_qname(bam_r);

        uint8_t *tag_ptr;
        const char *b_umi = get_tag(bam_r, umi_tag);
        if (b_umi == NULL) continue;
        const char *b_cb = get_tag(bam_r, bc_tag);
        if (b_cb == NULL) continue;
        int scnd_aln = ((bam_r->core.flag)&(BAM_FSECONDARY)) != 0; // secondary alignment
        if (scnd_aln) continue;

        if (max_nh > 0){
            tag_ptr = bam_aux_get(bam_r, nh_tag);
            if (tag_ptr == NULL){
                ret = err_msg(EXIT_FAILURE, 0, 
                        "NH set but could not find tag %s in read %s\n", 
                        nh_tag, qname);
                goto cleanup;
            }
            int nh = (int)bam_aux2i(tag_ptr);
            if (nh > max_nh) continue;
        }

        // Get overlapping variants, base calls, and base qualities.
        int vret = overlap_bam1_vars(bam_hdr, bam_r, gv, &n_var, &m_var, 
                &var, &base, &qual, min_phred);
        if (vret < 0){
            ret = err_msg(EXIT_FAILURE, 0, "ac: could not overlap variants with BAM record %s\n", 
                    qname);
            goto cleanup;
        }

        if (n_var > 0 || n_feat > 0){
            int rret = add_rec(records, b_cb, b_umi,  
                    n_feat, (const char **)feat, (const uint8_t *)splice, 
                    n_var, (const char **)var, (const uint8_t *)base, (const uint8_t *)qual); 
            if (rret < 0){
                ret = EXIT_FAILURE;
                goto cleanup;
            }
        }
    }

    // if sam_itr_next returns an error
    if (iter_ret < -1){
        ret = err_msg(EXIT_FAILURE, 0, "ac: could not get next read from BAM file\n");
        goto cleanup;
    }

    /* print out record stats */
    if (verbose){
        int nrec = n_recs(records);
        fprintf(stdout, "stored %i records from %i UMIs in %i barcodes\n", 
                nrec, records->umi_ix->n, records->bc_ix->n);
        fflush(stdout);
    }

    // if (verbose) print_summary(records);

    // call UMIs
    if (verbose){
        fprintf(stdout, "filtering UMI reads\n");
        fflush(stdout);
    }
    call_umis(records);

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

    if (verbose){
        fprintf(stdout, "finished\n");
        fflush(stdout);
    }

cleanup:
    if (vcffn) free(vcffn);
    if (out) free(out);
    if (gv) destroy_gv(gv);
    if (records) destroy_records(records);
    if (sr) bcf_sr_destroy(sr);
    if (a) destroy_bc_ac(a);

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);
    sam_hdr_destroy(bam_hdr);
    if (bam) sam_close(bam);
    if (bam_idx) hts_idx_destroy(bam_idx);

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);
    if (sr) bcf_sr_destroy(sr);
    sam_hdr_destroy(bam_hdr);
    if (bam) sam_close(bam);
    if (bam_idx) hts_idx_destroy(bam_idx);

    free(feat); 
    free(splice);
    free(var); 
    free(base); 
    free(qual);

    return(ret);
}


