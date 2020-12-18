/***************************************************************************
 * 
 *    simplify.c: Import a dataset and run the Clean algorithm
 * 
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "gene.h"
#include "common.h"
#include "bounds.h"

static void _print_usage(FILE *f, char *name)
{
    fprintf(f, "Usage: %s [options] < [input]\n", name);
    pretty_print(f, "The program reads data from the input file specified, deletes all sequences that can coalesce, and removes all singleton mutations.", 80, 0);
    fprintf(f, "Legal options are:\n");
    print_option(f, "-b[name]", "Output to file name.", 80, -1);
    print_option(f, "-k", "Whether common ancestor is known. If data is in binary format, the all zero sequence is assumed to be ancestral (whether or not this is present in the data). Otherwise, see the '-a' and '-n' options.", 80, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 80, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 80, -1);
    print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the '-k' option is used, the first sequence is assumed to be ancestral.", 80, -1);
    print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the '-k' option is used, the first sequence is assumed to be ancestral.", 80, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 80, -1);
}

int main(int argc, char **argv)
{
    Genes *g;
    AnnotatedGenes *a;
    int i, j, quiet = 0;
    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;
    FILE *fp;
    fp = stdout;
    
    eventlist = MakeLList();
    
    /* Analyse command line options */
    #define SIMPLIFY_OPTIONS "b::kofanQhH?"
    
    /* Parse command line options */
    while ((i = getopt(argc, argv, SIMPLIFY_OPTIONS)) >= 0){
        switch(i){
            case 'b':
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -b requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating compuation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        _print_usage(stderr, argv[0]);
                        exit(1);
                    }
                }
                break;
            case 'k':
                gene_knownancestor = 1;
                break;
            case 'Q':
                quiet = 1;
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
    
    if ((gene_knownancestor) && (seqtype != GENE_BINARY)) {
        /* First sequence only included to specify known common ancestor */
        remove_annotatedgene(a, 0);
    }
    g = a->g;
    

    //Initialise the elements array (this will track the number of recombinations which each of the sequences has undergone)
    elements = elist_make();
    sites = elist_make();
    if ((gene_knownancestor) && (seqtype != GENE_BINARY)) {
        for(i=0; i < g->n; i++) {
            elist_append(elements, (void *)(i+1));
        }
    } else {
        for(i=0; i < g->n; i++) {
            elist_append(elements, (void *)i);
        }
    }
    // Initialise the list of sites
    for(i=0; i < g->length; i++) {
        elist_append(sites, (void *)i);
    }
    
    // Print stats for input dataset
//         printf("Input dataset has %d sequences and %d sites\n", g->n, g->length);
    printf("Input dataset: %d sequences, %d sites\n", g->n, g->length);
    
    implode_genes(g);
    
    // Print stats for reduced dataset
//         printf("Reduced dataset has %d sequences and %d sites\n", g->n, g->length);
    printf("Reduced dataset: %d sequences, %d sites\n", g->n, g->length);
    
    output_genes(g, fp, NULL);
   
    printf("Sequences:\n");
    print_elist(elements, NULL);
    
    printf("Sites:\n");
    print_elist(sites, NULL);

    // Tidying
    if (eventlist != NULL){
        while (Length(eventlist) > 0)
            free(Pop(eventlist));
        DestroyLList(eventlist);
    }
    free_annotatedgenes(a);

    
    return 0;
}
