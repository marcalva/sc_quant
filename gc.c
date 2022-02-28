
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <errno.h>
#include "htslib/sam.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "sam_read.h"
#include "str_util.h"
#include "variants.h"
#include "gtf_anno.h"
#include "g_count.h"
#include "bc_umi.h"
#include "overlap.h"

static void print_status_gc(const char* k, int i, const char *chr, int pos){
    log_msg("processed %i %s through %s:%i", i, k, chr, pos);
}

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "sc_quant v0.1.0: Generate counts from single-cell experiments\n"
            "Usage:    sc_quant gc [options] --bam bamfile --gtf gtffile --out outfile \n"
            "\n"
            "Options:\n"
            "  -b, --bam           Indexed BAM file.\n"
            "  -g, --gtf           GTF file.\n"
            "  -o, --out           Output file prefix [gc.].\n"
            "  -B, --bc-tag        BAM tag for barcode [CB].\n"
            "  -U, --umi-tag       BAM tag for UMI [UB].\n"
            "  -H, --nh-tag        BAM tag for the number of alignments of a read [NH].\n"
            "  -m, --max-nh        Only process reads with a maximum of this many alignments [1]. Set to 0 to ignore.\n"
            "  -t, --tx-basic      Read only transcripts tagged as 'basic' in the GTF file.\n"
            "  -c, --barcodes      File containing list of barcode IDs, one per line.\n"
            "  -r, --region        Region (hts format), for example 'chr21,chr21:10-,chr21-10-20'.\n"
            "  -V, --verbose       Write status on output.\n"
            "\n");
    exit(exit_status);
}


int gene_count(int argc, char *argv[]){

    if (argc == 1) usage(stderr, EXIT_FAILURE);

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"bam", required_argument, NULL, 'b'}, 
        {"gtf", required_argument, NULL, 'g'},
        {"out", required_argument, NULL, 'o'}, 
        {"bc-tag", required_argument, NULL, 'B'},
        {"umi-tag", required_argument, NULL, 'U'},
        {"nh-tag", required_argument, NULL, 'H'}, 
        {"max-nh", required_argument, NULL, 'm'}, 
        {"tx-basic", no_argument, NULL, 't'}, 
        {"barcodes", required_argument, NULL, 'c'},
        {"region", required_argument, NULL, 'r'},
        {"verbose", no_argument, NULL, 'V'},
        {NULL, 0, NULL, 0}
    };

    /* parameters */
    int ret = 0;
    int hs = 0;
    char bc_tag[] = "CB";
    char umi_tag[] = "UB";
    char nh_tag[] = "NH";
    char *bamfn = NULL;
    char *gtffn = NULL;
    char *outfn = "gc.";
    int max_nh = 1;
    int tx_basic = 0;
    char *p_end = NULL;
    int bc_set = 0;
    char *bc_fn = NULL;
    char *region = ".";
    int region_set = 0;
    int verbose = 0;
    /* */

    /* local variables */
    Records *records = NULL;
    Annotation *a = NULL;
    bc_gc *g_counts = NULL;
    GenomeVar *gv = NULL;
    samFile *bam = NULL;
    sam_hdr_t *bam_hdr = NULL;
    hts_idx_t *bam_idx = NULL;
    hts_itr_t *bam_itr = NULL;
    bam1_t *bam_r = NULL;

    uint8_t n_feat = 0, m_feat = 4;
    char **feat = (char **)calloc(m_feat, sizeof(char *));
    uint8_t *splice = (uint8_t *)calloc(m_feat, sizeof(uint8_t));;
    uint8_t n_var = 0;
    char **var = NULL;
    uint8_t *base = NULL;
    uint8_t *qual = NULL;
    /* */

    int option_index = 0;
    int cm;
    while ((cm = getopt_long_only(argc, argv, "hb:g:o:B:U:H:m:p:tc:r:V", loptions, &option_index)) != -1){
        switch(cm){
            case 'h': hs = 1; break;
            case 'b': bamfn = optarg; break; 
            case 'g': gtffn = optarg; break;
            case 'o': outfn = optarg; break; 
            case 'B':
                      if (strlen(optarg) != 2){
                          ret = err_msg(EXIT_FAILURE, 0, 
                                  "bc-tag must be a 2 character string: %s given", optarg);
                          goto cleanup;
                      }
                      strncpy(bc_tag, optarg, 2);
                      break;
            case 'U':
                      if (strlen(optarg) != 2){
                          ret = err_msg(EXIT_FAILURE, 0, 
                                  "umi-tag must be a 2 character string: %s given", optarg);
                          goto cleanup;
                      }
                      strncpy(umi_tag, optarg, 2);
                      break;
            case 'H':
                      if (strlen(optarg) != 2){
                          ret = err_msg(EXIT_FAILURE, 0, 
                                  "nh-tag must be a 2 character string: %s given", optarg);
                          goto cleanup;
                      }
                      strncpy(nh_tag, optarg, 2);
                      break;
            case 'm':
                      errno = 0;
                      max_nh = (int)strtol(optarg, &p_end, 10);
                      if (max_nh == 0 && errno > 0){
                          return err_msg(EXIT_FAILURE, 0, 
                                  "plp: could not convert --max-nh %s to int: %s", 
                                  optarg, strerror(errno));
                      }
                      break;
            case 't':
                      tx_basic = 1;
                      break;
            case 'c':
                      bc_fn = optarg;
                      bc_set = 1;
                      break;
            case 'r':
                      region = optarg;
                      region_set = 1;
                      break;
            case 'V':
                      verbose = 1;
                      break;
            default: 
                      err_msg(EXIT_FAILURE, 0, 
                              "unrecognized option: %s", loptions[option_index].name);
                      usage(stdout, EXIT_FAILURE);
        }
    }

    if ( hs ) usage(stdout, EXIT_SUCCESS);

    /* Check files */
    if ( !bamfn ){
        err_msg(1, 0, "BAM file must be given");
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if ( !gtffn ){
        err_msg(1, 0, "GTF file must be given");
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if ( !outfn ){
        err_msg(1, 0, "Out file prefix must be given");
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    /* */

    /* print options set */
    if ( verbose ){
        fprintf(stdout, "options set:\n");
        fprintf(stdout, "    BAM:             %s\n", bamfn);
        fprintf(stdout, "    GTF:             %s\n", gtffn);
        fprintf(stdout, "    Out:             %s\n", outfn);
        fprintf(stdout, "    BC tag:          %s\n", bc_tag);
        fprintf(stdout, "    UMI tag:         %s\n", umi_tag);
        fprintf(stdout, "    NH tag:          %s\n", nh_tag);
        fprintf(stdout, "    Max NH:          %i\n", max_nh);
        if ( bc_set ){
            fprintf(stdout, "    Barcode file:    %s\n", bc_fn);
        }
        if ( region_set ){
            fprintf(stdout, "    Region:          %s\n", region);
        }
        fflush(stdout);
    }
    /* */

    /* Instantiate records */
    records = init_records();

    /* Read in GTF file */
    if (verbose)
        log_msg("reading annotations from %s", gtffn);
    a = read_from_gtf(gtffn, tx_basic);
    if (a == NULL){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (verbose)
        log_msg("read %i genes", a->gene_ix->n); 


    /* Read in BAM file */
    bam = sam_open(bamfn, "r");
    if (bam == NULL){
        err_msg(1, 0, "could not open BAM file %s", bamfn);
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    bam_hdr = sam_hdr_read(bam);
    if (bam_hdr == NULL){
        err_msg(1, 0, "could not read header for BAM file %s", bamfn);
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    // may remove this
    bam_idx = sam_index_load(bam, bamfn);
    if (bam_idx == NULL){
        err_msg(1, 0, "could not load index for BAM file %s", bamfn);
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    /* */

    int n_reads = 0;
    int n_bc = 0;
    /* */

    /* set barcodes */
    if ( bc_set ){
        n_bc = add_bc2ix_file(records, bc_fn);
        if (n_bc < 0){
            ret = err_msg(1, 0, "could not read barcodes from %s", bc_fn);
            goto cleanup;
        }
        if ( verbose )
            log_msg("read %i barcodes from %s", n_bc, bc_fn);
    }
    /* */

    /* set genes */
    if (add_from_str_map(records->gene_ix, a->gene_ix) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* Set BAM iterator to region */
    bam_itr = sam_itr_querys(bam_idx, bam_hdr, region);
    if (bam_itr == NULL){
        err_msg(1, 0, "could not set region %s in BAM file", region);
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    /* */
            
    // Loop over BAM records/

    log_msg("processing alignments in %s", bamfn);
    int iter_ret;
    bam_r = bam_init1();
    while ( (iter_ret = sam_itr_next(bam, bam_itr, bam_r)) >= 0){
        if ( (verbose) && (n_reads % (int)10e6 == 0)){
            const char *chr = sam_hdr_tid2name(bam_hdr, bam_r->core.tid);
            print_status_gc("million alignments", n_reads/1e6, chr, (int)(bam_r->core.pos+1));
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
                err_msg(1, 0, "NH set but could not find tag %s in read %s\n", nh_tag, qname);
                ret = EXIT_FAILURE;
                goto cleanup;
            }
            int nh = (int)bam_aux2i(tag_ptr);
            if (nh > max_nh) continue;
        }

        // Get overlapping features
        int fret = overlap_bam1_feats(bam_hdr, bam_r, a, &n_feat, &m_feat, &feat, &splice);
        if (fret < 0){
            ret = EXIT_FAILURE;
            goto cleanup;
        }

        if (n_feat > 0){
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
        err_msg(1, 0, "could not get next read from BAM file\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (verbose){
        int nrec = n_recs(records);
        log_msg("stored %i records from %i UMIs in %i barcodes", 
                nrec, records->umi_ix->n, records->bc_ix->n);
    }

    // if (verbose) print_summary(records);

    // call UMIs
    if (verbose) log_msg("filtering UMI reads");
    fflush(stdout);
    call_umis(records);

    // gene counts
    g_counts = init_bc_gc();

    if (verbose) log_msg("counting UMIs");
    if (bc_gc_add_gene_map(g_counts, records->gene_ix) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (bc_gc_add_bc_map(g_counts, records->bc_ix) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (bc_gc_count(g_counts, records) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    if (bc_gc_write(g_counts, a, outfn) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }


cleanup:

    destroy_gv(gv);
    destroy_anno(a);
    destroy_records(records);
    destroy_bc_gc(g_counts);

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);
    sam_hdr_destroy(bam_hdr);
    if (bam) sam_close(bam);
    if (bam_idx) hts_idx_destroy(bam_idx);

    free(feat); 
    free(splice);
    free(var); 
    free(base); 
    free(qual);

    if (verbose) log_msg("finished counting UMIs");

    return(ret);
}

