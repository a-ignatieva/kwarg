/***************************************************************************
 * 
 *    flip.c: Flip specified nucleotides within a dataset.
 * 
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>

#include "gene.h"
#include "common.h"
#include "bounds.h"

static void _print_usage(FILE *f, char *name)
{
    fprintf(f, "Usage: %s [options] < [input]\n", name);
    pretty_print(f, "The program reads data from the input file specified, flips the specified nucleotide (given by sequence and site number, starting at 1) and output the result to a specified file.", 80, 0);
    fprintf(f, "Legal options are:\n");
    print_option(f, "-b[]", "Output file name to write the amended sequence data.", 80, -1);
    print_option(f, "-q[]", "Sequence id to alter.", 80, -1);
    print_option(f, "-s[]", "Site number to alter.", 80, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 80, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 80, -1);
    print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site).", 80, -1);
    print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site).", 80, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 80, -1);
}

int main(int argc, char **argv)
{
    Genes *g;
    AnnotatedGenes *a;
    int i;
    int q;
    int s;
    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;
    FILE *fp;
    fp = stdout;
    LList *sequences = MakeLList();
    LList *sites = MakeLList();
    char *token;
    char *endptr;
    errno = 0;
    
    /* Analyse command line options */
    #define SIMPLIFY_OPTIONS "b::q:s:ofanhH?"
    
    /* Parse command line options */
    while ((i = getopt(argc, argv, SIMPLIFY_OPTIONS)) >= 0){
        switch(i){
            case 'b':
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -b requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        _print_usage(stderr, argv[0]);
                        exit(1);
                    }
                }
                
                break;
            case 'q':
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        q = strtol(token, &endptr, 10);
                        if(errno != 0 || *endptr != '\0' || q <= 0) {
                            fprintf(stderr, "Sequence number should be a positive integer.\n");
                            exit(1);
                        }
                        Enqueue(sequences, (void *)q);
                        token = strtok(NULL, ",");
                    }
                }
                break;
            case 's':
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        q = strtol(token, &endptr, 10);
                        if(errno != 0 || *endptr != '\0' || q <= 0) {
                            fprintf(stderr, "Site number should be a positive integer.\n");
                            exit(1);
                        }
                        Enqueue(sites, (void *)q);
                        token = strtok(NULL, ",");
                    }
                }
                break;
            case 'o':
                format = GENE_BEAGLE;
                break;
            case 'f':
                format = GENE_FASTA;
                break;
            case 'a':
                seqtype = GENE_AMINO;
                break;
            case 'n':
                seqtype = GENE_NUCLEIC;
                break;
            case 'h':
            case 'H':
            case '?':
                _print_usage(stdout, argv[0]);
                exit(0);
            case ':':
                _print_usage(stderr, argv[0]);
                exit(1);
        }
    }
    
    
    /* Read data */
    if (argc > optind){
            fprintf(stderr, "Not a valid option: %s\n", argv[optind]);
            exit(1);
    }
    else{
        if ((a = read_genes(NULL, format, seqtype)) == NULL){
            fprintf(stderr, "Could not parse input as valid data\n");
            exit(1);
        }
    }
    
    g = a->g;
    
    char c;
    while (Length(sequences) !=0) {
        q = (int)Pop(sequences);
        s = (int)Pop(sites);
        c = get_genes_character(g, q-1, s-1);
        switch(c) {
            case 0:
                set_genes_character(g, q-1, s-1, 1);
                break;
            case 1:
                set_genes_character(g, q-1, s-1, 0);
                break;
            case 2:
                break;
        }
    }
    
    output_genes(g, fp, NULL);
    
    // Clean up
    free_annotatedgenes(a);
    DestroyLList(sequences);
    DestroyLList(sites);
    
    return 0;
}
