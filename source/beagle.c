/**************************************************************************
 * 
 * beagle.c: Implementation of a branch and bound algorithm to find
 * the minimum number of recombinations in an SNP data set. 
 * 
 * The name pays homage, not to a certain ship carrying yet another Englishman
 * on a cruise, but to the infamous Beagle Boys that go by the same
 * initials as Branch & Bound.
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include "gene.h"
#include "bounds.h"
#include "exact.h"
#include "common.h"
#include "backtrack.h"

static void print_usage(FILE *f, char *name)
{
    fprintf(f, "Usage: %s [options] < [input]\n", name);
    pretty_print(f, "The program reads data from the input file specified (or from standard input if no input is specified) and determines the minimum number of recombinations required to explain the data set under the infinite sites assumption.", 80, 0);
    fprintf(f, "Legal options are:\n");
    print_option(f, "-p", " Print lower bounds on minimum number of recombinations as they become established.", 80, -1);
    print_option(f, "-s", " Suppress output of minimum number of recombinations required", 80, -1);
    print_option(f, "-b[name]", "Output a minimum recombination history to file name.", 80, -1);
    print_option(f, "-d[name]", "Output ancestral recombination graph of minimum recombination history in dot format to file name.", 80, -1);
    print_option(f, "-g[name]", "Output ancestral recombination graph of minimum recombination history in GDL format to file name.", 80, -1);
    print_option(f, "-j[name]", "Output ancestral recombination graph of minimum recombination history in GML format to file name.", 80, -1);
    print_option(f, "-t[name]", "Output list of marginal phylogenies for each site in Newick's 8:45 format to file name.", 80, -1);
    print_option(f, "-D[name]", "Output list of marginal phylogenies for each site in dot format to file name.", 80, -1);
    print_option(f, "-G[name]", "Output list of marginal phylogenies for each site in GDL format to file name.", 80, -1);
    print_option(f, "-J[name]", "Output list of marginal phylogenies for each site in GML format to file name.", 80, -1);
    print_option(f, "-I", "Marginal trees are only output one for each of the intervals between two recombination points, instead of one for each site.", 80, -1);
    print_option(f, "-vnodelabel", "Use nodelabel convention for labelling nodes in ancestral recombination graphs. The possible conventions are\nnone - do not label nodes\nid - only label nodes representing sampled sequences, using their sequence ids from the data file, and nodes representing recombinations, indicating the recombination point\nsequence - label nodes with the inferred sequences; these sequences will be in binary format, with 0 representing wild type and 1 representing mutant type, even if the original data is not in binary format\nboth - label nodes with both id and inferred sequence\none - use only one label for a node, sequence id or recombination point if avaialable and otherwise the inferred sequence\nDefault convention is id. The colour coding scheme used for the nodes is red for sequences in the input data, blue for recombination nodes, green for standard coalescent nodes, and yellow for the final coalescence into the most recent common ancestor.", 80, -1);
    print_option(f, "-i", "Sequences not having a sequence id in the data file are assigned their index in the data file as id, e.g. the first sequence in the data file would be assigned `1' as id.", 80, -1);
    print_option(f, "-e", "Label edges in ancestral recombination graphs with the sites undergoing mutation along the edge.", 80, -1);
    #ifdef HAPLOTYPE_BLOCKS
    print_option(f, "-m[name]", "Output matrix of minimum number of recombinations required in any block, with the entry of column i and row j being the minimum number of recombinations required in the region from site i to site j + 1; the matrix is output to file name.", 80, -1);
    #endif
    print_option(f, "-k", "Assume that the common ancestral sequence is known, i.e. that we know which is the wild type and which is the mutant in each site. In binary format the all-0 sequence is assumed to be the common ancestral sequence, while the common ancestral sequence has to be specified directly for amino acid and nucleic sequences (see options -a and -n)", 80, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 80, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 80, -1);
    print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known, the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set." , 80, -1);
    print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known, the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 80, -1);
    print_option(f, "-ln", "Start search from n recombinations. If n is larger than the minimum number of recombinations required, then the inferred evolutionary history may not be minimal.", 80, -1);
    print_option(f, "-un", "Only look for histories with at most n recombinations. If n is less than the minimum number of recombinations required, then no evolutionary history is inferred; it is simply stated that the data requires at least n + 1 recombinations.", 80, -1);
    print_option(f, "-r", "Usually the program tries to explore search paths in order of decreasing likelihood of success; this option changes the behaviour such that search paths are explored in random order, i.e. two consecutive runs may give different evolutionary histories.", 80, -1);
    print_option(f, "-c[N]", "This option allows a more comprehensive (although also more time consuming) random exploration of search paths than the -r option. With this option also search paths ignored a priori in the branch & bound algorithm of beagle are possible. N more than the minimum number of recombinations are allowed.", 80, -1);
    print_option(f, "-h, -H, -?", "Print this information and stop.", 80, -1);
}

int main(int argc, char **argv)
{
    Genes *g;
    AnnotatedGenes *a;
    int i, j,
    print_progress = 0,
    lower = -1,
    upper = -1,
    comprehensive_bound = -1,
    silent = 0,
    intervals = 0;
    char *tmp;
    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;
    FILE *fp;
    LList *history_files = MakeLList(),
    *dot_files = MakeLList(),
    *gml_files = MakeLList(),
    *gdl_files = MakeLList(),
    *tree_files = MakeLList(),
    *dottree_files = MakeLList(),
    *gmltree_files = MakeLList(),
    *gdltree_files = MakeLList();
    ARG *arg = NULL;
    ARGLabels nodelabel = ARGLABEL;
    int edgelabel = 0;
    int generate_id = 0;
    HashTable *t = NULL;
    #ifdef HAPLOTYPE_BLOCKS
    FILE *haploblock_file = NULL;
    #endif
    gene_knownancestor = 0;
    
    /* Initialise random number generator */
    initialise_xrandom();
    
    #ifdef ENABLE_VERBOSE
    set_verbose(1);
    #endif
    
    /* Analyse command line options */
    #define BEAGLE_OPTIONS1 "phH?l:u:rfanosb::d::g::j::t::D::G::J::Iv:iec::k"
    #ifdef HAPLOTYPE_BLOCKS
    #define BEAGLE_OPTIONS BEAGLE_OPTIONS1 "m::"
    #else
    #define BEAGLE_OPTIONS BEAGLE_OPTIONS1
    #endif
    
    /* Parse command line options */
    while ((i = getopt(argc, argv, BEAGLE_OPTIONS)) >= 0){
        switch(i){
            case 'h':
            case 'H':
            case '?':
                print_usage(stdout, argv[0]);
                exit(0);
            case ':':
                print_usage(stderr, argv[0]);
                exit(1);
            case 'p':
                print_progress = 1;
                break;
            case 's':
                silent = 1;
                break;
            case 'b':
                /* Backtrack history leading to minimum number of recombinations */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -b requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating compuation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(history_files, (void *)optarg);
                }
                else
                    Enqueue(history_files, stdout);
                break;
            case 'd':
                /* Output ancestral recombination graph of history leading to
                 * minimum number of recombinations in dot format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -d requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(dot_files, (void *)optarg);
                }
                else
                    Enqueue(dot_files, stdout);
                break;
            case 'j':
                /* Output ancestral recombination graph of history leading to
                 * minimum number of recombinations in gml format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -j requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gml_files, (void *)optarg);
                }
                else
                    Enqueue(gml_files, stdout);
                break;
            case 'g':
                /* Output ancestral recombination graph of history leading to
                 * minimum number of recombinations in gdl format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -g requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gdl_files, (void *)optarg);
                }
                else
                    Enqueue(gdl_files, stdout);
                break;
            case 't':
                /* Output marginal trees in ancestral recombination graph of
                 * history leading to minimum number of recombinations in
                 * Newick's 8:45 format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -t requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(tree_files, (void *)optarg);
                }
                else
                    Enqueue(tree_files, stdout);
                break;
            case 'D':
                /* Output marginal trees in ancestral recombination graph of
                 * history leading to minimum number of recombinations in
                 * dot format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -D requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(dottree_files, (void *)optarg);
                }
                else
                    Enqueue(dottree_files, stdout);
                break;
            case 'J':
                /* Output marginal trees in ancestral recombination graph of
                 * history leading to minimum number of recombinations in
                 * GML format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -J requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gmltree_files, (void *)optarg);
                }
                else
                    Enqueue(gmltree_files, stdout);
                break;
            case 'G':
                /* Output marginal trees in ancestral recombination graph of
                 * history leading to minimum number of recombinations in
                 * GDL format.
                 */
                /* Was a file name specified? */
                if (optarg != 0){
                    if(optarg[0] == '-') {
                        fprintf(stderr, "Option -G requires an output file.\n");
                        exit(1);
                    }
                    /* Check whether file can be written before initiating computation */
                    if ((fp = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gdltree_files, (void *)optarg);
                }
                else
                    Enqueue(gdltree_files, stdout);
                break;
            case 'I':
                intervals = 1;
                break;
            case 'v':
                if (!strcmp(optarg, "none"))
                    nodelabel = ARGNONE;
                else if (!strcmp(optarg, "id"))
                    nodelabel = ARGLABEL;
                else if (!strcmp(optarg, "sequence"))
                    nodelabel = ARGSEQUENCE;
                else if (!strcmp(optarg, "both"))
                    nodelabel = ARGBOTH;
                else if (!strcmp(optarg, "one"))
                    nodelabel = ARGLABELFIRST;
                else{
                    fprintf(stderr, "Unrecognised nodelabel convention (%s) for -v option\n", optarg);
                    print_usage(stderr, argv[0]);
                    exit(1);
                }
                break;
            case 'i':
                generate_id = 1;
                break;
            case 'e':
                edgelabel = 1;
                break;
                #ifdef HAPLOTYPE_BLOCKS
            case 'm':
                /* Compute matrix of local minimum number of recombinations required */
                /* Check whether we already opened an output file that now needs
                 * to be closed.
                 */
                if ((haploblock_file != NULL) && (haploblock_file != stdout))
                    fclose(haploblock_file);
                /* Was a file name specified? */
                if (optarg != 0){
                    if ((haploblock_file = fopen(optarg, "w")) == NULL){
                        fprintf(stderr, "Could not open file %s for output\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                }
                else
                    haploblock_file = stdout;
                break;
                #endif
            case 'k':
                gene_knownancestor = 1;
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
            case 'l':
                lower = strtol(optarg, &tmp, 10);
                if (*tmp != '\0'){
                    fprintf(stderr, "Could not recognize %s as lower bound on recombinations needed\n", optarg);
                    print_usage(stderr, argv[0]);
                    exit(1);
                }
                break;
            case 'u':
                upper = strtol(optarg, &tmp, 10);
                if (*tmp != '\0'){
                    fprintf(stderr, "Could not recognize %s as upper bound on recombinations allowed\n", optarg);
                    print_usage(stderr, argv[0]);
                    exit(1);
                }
                break;
            case 'r':
                exact_randomise = 1;
                break;
            case 'c':
                /* Was a bound specified? */
                if (optarg != 0){
                    comprehensive_bound = strtol(optarg, &tmp, 10);
                    if (*tmp != '\0'){
                        fprintf(stderr, "Could not recognize %s as number of recombinations allowed\n", optarg);
                        print_usage(stderr, argv[0]);
                        exit(1);
                    }
                }
                else
                    /* Default is histories with a minimum number of recombinations */
                    comprehensive_bound = 0;
                break;
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
    if ((gene_knownancestor) && (seqtype != GENE_BINARY))
        /* First sequnce only included to specify known common ancestor */
        remove_annotatedgene(a, 0);
    g = a->g;
    
    #ifdef ENABLE_VERBOSE
    output_annotatedgenes(a, NULL, NULL);
    #endif
    
    /* Set up structures for computation */
    if (comprehensive_bound >= 0)
        t = beagle_allocate_hashtable(g, -1);
    else if ((Length(history_files) > 0) || (Length(dot_files) > 0)
        || (Length(gml_files) > 0) || (Length(gdl_files) > 0)
        || (Length(tree_files) > 0) || (Length(dottree_files) > 0)
        || (Length(gmltree_files) > 0) || (Length(gdltree_files) > 0))
        eventlist = MakeLList();
    #ifdef HAPLOTYPE_BLOCKS
    if (haploblock_file != NULL){
        haploblocks = (int **)xmalloc((g->length - 1) * sizeof(int *));
        for (i = 0; i < g->length - 1; i++)
            haploblocks[i] = (int *)xcalloc((g->length - i - 1), sizeof(int));
    }
    #endif
    #ifdef DEBUG
    ancestral_state_trace = beagle_allocate_hashtable(NULL, -1);
    #endif
    
    /* Compute number of recombinations */
    if ((lower >= 0) || (upper >= 0)){
        if (comprehensive_bound >= 0){
            i = beagle_reusable_bounded(g, (print_progress ? stdout : NULL), lower,
                                        upper, t);
            comprehensive_bound += i;
        }
        else
            i = beagle_bounded(g, (print_progress ? stdout : NULL), lower, upper);
    }
    else{
        if (comprehensive_bound >= 0){
            i = beagle_reusable(g, (print_progress ? stdout : NULL), t);
            comprehensive_bound += i;
        }
        else
            i = beagle(g, (print_progress ? stdout : NULL));
    }
    
    if (!silent){
        if ((upper < 0) || (i <= upper))
            if (lower < i)
                printf("Minimum number of recombinations required is %d\n", i);
            else
                printf("Data can be explained with %d recombination%s\n", i, i != 1 ? "s" : "");
            else
                printf("At least %d recombination%s required\n", i, i != 1 ? "s" : "");
    }
    if ((upper >= 0) && (comprehensive_bound > upper))
        comprehensive_bound = upper;
    
    #ifdef HAPLOTYPE_BLOCKS
    if (haploblock_file != NULL){
        for (i = 1; i < g->length; i++){
            for (j = 0; j < i - 1; j++)
                fprintf(haploblock_file, "%d ", haploblocks[j][i - j - 1]);
            fprintf(haploblock_file, "%d\n", haploblocks[j][i - j - 1]);
        }
        for (i = 0; i < g->length - 1; i++)
            free(haploblocks[i]);
        free(haploblocks);
        if (haploblock_file != stdout)
            fclose(haploblock_file);
    }
    #endif
    if ((Length(history_files) > 0) || (Length(dot_files) > 0)
        || (Length(gml_files) > 0) || (Length(gdl_files) > 0)
        || (Length(tree_files) > 0) || (Length(dottree_files) > 0)
        || (Length(gmltree_files) > 0) || (Length(gdltree_files) > 0)){
        if (comprehensive_bound >= 0){
            eventlist	= beagle_randomised(g, NULL, comprehensive_bound, t);
            beagle_deallocate_hashtable(t);
        }
        while ((fp = (FILE *)Pop(history_files)) != NULL){
            if (fp != stdout)
                /* Open named file for output */
                if ((fp = fopen((char *)fp, "w")) == NULL){
                    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    continue;
                }
                /* Only remember last ARG constructed (they should all be the same) */
                if (arg != NULL)
                    arg_destroy(arg);
                arg = eventlist2history(a, fp);
        }
        if (arg == NULL)
            arg = eventlist2history(a, NULL);
        if (arg != NULL){
            /* Output ARG in dot format */
            while ((fp = (FILE *)Pop(dot_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, ARGDOT, nodelabel, edgelabel, generate_id);
                if (fp != stdout)
                    fclose(fp);
            }
            /* Output ARG in GML format */
            while ((fp = (FILE *)Pop(gml_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, ARGGML, nodelabel, edgelabel, generate_id);
                if (fp != stdout)
                    fclose(fp);
            }
            /* Output ARG in dot format */
            while ((fp = (FILE *)Pop(gdl_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, ARGGDL, nodelabel, edgelabel, generate_id);
                if (fp != stdout)
                    fclose(fp);
            }
            /* Output marginal trees of ARG in Newick's 8:45 format */
            while ((fp = (FILE *)Pop(tree_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, TREENEWICK, nodelabel, edgelabel, generate_id,
                               intervals);
                    if (fp != stdout)
                        fclose(fp);
            }
            /* Output marginal trees of ARG in dot format */
            while ((fp = (FILE *)Pop(dottree_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, TREEDOT, nodelabel, edgelabel, generate_id,
                               intervals);
                    if (fp != stdout)
                        fclose(fp);
            }
            /* Output marginal trees of ARG in GML format */
            while ((fp = (FILE *)Pop(gmltree_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, TREEGML, nodelabel, edgelabel, generate_id,
                               intervals);
                    if (fp != stdout)
                        fclose(fp);
            }
            /* Output marginal trees of ARG in GDL format */
            while ((fp = (FILE *)Pop(gdltree_files)) != NULL){
                if (fp != stdout)
                    /* Open named file for output */
                    if ((fp = fopen((char *)fp, "w")) == NULL)
                        fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
                    arg_output(arg, a, fp, TREEGDL, nodelabel, edgelabel, generate_id,
                               intervals);
                    if (fp != stdout)
                        fclose(fp);
            }
            arg_destroy(arg);
        }
        if (eventlist != NULL){
            while (Length(eventlist) > 0)
                free(Pop(eventlist));
            DestroyLList(eventlist);
        }
        }
        
        /* Clean up */
        DestroyLList(dot_files);
        DestroyLList(gml_files);
        DestroyLList(gdl_files);
        DestroyLList(tree_files);
        DestroyLList(dottree_files);
        DestroyLList(gmltree_files);
        DestroyLList(gdltree_files);
        DestroyLList(history_files);
        #ifdef DEBUG
        if (ancestral_state_trace != NULL)
            hashtable_destroy(ancestral_state_trace,
                              (void (*)(void *))free_packedgenes, NULL,
                              (void (*)(void *))free);
        #endif
            
            
        free_annotatedgenes(a);
        
        return 0;
}
