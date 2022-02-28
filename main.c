
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <errno.h>
#include "scq_version.h"

int gene_count(int argc, char *argv[]);
int allele_count(int argc, char *argv[]);

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "sc_quant: Generate counts from single-cell experiments\n"
            "Version:  %s\n"
            "Usage:    sc_quant <command> [options]\n"
            "\n"
            "Commands:\n"
            "  gc           Generate gene-barcode counts.\n"
            "  ac           Generate variant allele-barcode counts.\n"
            "\n"
            , SCQ_VERSION);
    exit(exit_status);
}


int main(int argc, char *argv[]){
    errno = 0;

    if (argc < 2) usage(stderr, EXIT_FAILURE);

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0){
        usage(stderr, EXIT_SUCCESS);
    }

    int ret = 0;
    if (strcmp(argv[1], "gc") == 0)
        ret = gene_count(argc - 1, argv + 1);
    else if (strcmp(argv[1], "ac") == 0)
        ret = allele_count(argc - 1, argv + 1);
    else {
        fprintf(stderr, "command %s not recognized\n", argv[1]);
        return EXIT_FAILURE;
    }

    return(ret);
}

