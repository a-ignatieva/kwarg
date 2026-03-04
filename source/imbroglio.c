/***************************************************************************
 *
 *    imbroglio.c
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

#include "gene.h"
#include "bounds.h"
#include "exact.h"
#include "common.h"
#include "backtrack.h"

static void _print_usage(FILE *f, char *name)
{
    fprintf(f, "Usage: %s [options] < [input]\n", name);
    pretty_print(f, "This is a silly program.", 70, 0);
    fprintf(f, "Legal options are:\n");
    print_option(f, "-V[x]", "If running a single iteration with given cost parameters, this controls the level of verbosity. \nx = 0: no extra output \nx = 1: during each neighbourhood search, output the number of neighbours explored, the move selected and its cost \nx = 2: during each neighbourhood search, output the number of neighbours explored, the resulting configuration and cost of each neighbour, the move selected and its cost.", 70, -1);
    print_option(f, "-b[name]", "Output a minimum recombination history to file name.", 70, -1);
    print_option(f, "-R[x]", "Which sequence to recombine (indexing input sequences from 0). No recombination if this is -1.", 70, -1);
    print_option(f, "-P[x]", "After which site to recombine (indexing input sites from 0). No recombination if this is -1.", 70, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 70, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 70, -1);
    print_option(f, "-k", "Assume that the common ancestral sequence is known, i.e. that we know which is the wild type and which is the mutant in each site. If the data is in binary format, the all-0 sequence is assumed to be the common ancestral sequence (this does not need to be present in the data). If the data is in amino acid or nucleotide format, the common ancestral sequence has to be specified directly and is taken to be the first sequence in the data file (see options -a and -n)", 70, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 70, -1);
}

/* Parse a floating point option argument and store it in value.
 * Return value states whether the argument could be parsed in full.
 */
static int _parse_double(char *s, double *value)
{
    int i;
    
    if (sscanf(s, "%lf%n", value, &i) != 1)
        return 0;
    
    return s[i] == '\0';
}

/* Read entire content of file name into a string */
static char *_read_file(char *name)
{
    FILE *f = stdin;
    char *s = (char *)xmalloc(8 * sizeof(char));
    int size = 0, capacity = 8;
    
    /* Open file */
    if ((name != NULL) && ((f = fopen(name, "r")) == NULL))
        return NULL;
    
    /* Read file character by character */
    while ((s[size++] = fgetc(f)) != EOF)
        if (size == capacity - 1){
            /* Ran out of buffer capacity, double its size */
            capacity *= 2;
            s = (char *)xrealloc(s, capacity * sizeof(char));
        }
        
        /* Zero-terminate string, shrink buffer to fit it, and close file */
        s[size - 1] = '\0';
    s = xrealloc(s, size);
    if (f != stdin)
        fclose(f);
    
    return s;
}

int main(int argc, char **argv)
{
    Genes *g, *h;
    AnnotatedGenes *a;
    KwargContext ctx;
    
    int i, j = 0, k = 0, l = 0, m = 0, t = 0,
    head = 1,
    silent = 0,
    intervals = 0,
    multruns = 0,
    costs_in = 0;
    int R = -1, P = -1;
    
    double n;
    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;
    FILE *print_progress = stdout;
    int (*select)(double, KwargContext *) = NULL;
    FILE *fp;
    LList *history_files = MakeLList();
    ARG *arg = NULL;
    ARGLabels nodelabel = ARGLABEL;
    int edgelabel = 0;
    int generate_id = 0;
    int ontheflyselection = 0;
    ctx.gc_enabled = 0;
    Event *e;
    LList *tmp;
    char *token;
    double timer;
    clock_t tic, toc;
    char *endptr;
    errno = 0;
    int reference = -1;
    
    ctx.eventlist = NULL; // Will be created later if needed
    ctx.elements = NULL;
    ctx.sites = NULL;
    ctx.lookup = NULL;
    ctx.seq_numbering = 0;
    ctx.gc_enabled = 0;
    ctx.rec_max = INT_MAX;
    ctx.rm_max = INT_MAX;
    ctx._greedy_functioncalls = NULL;
    ctx._greedy_beaglereusable = NULL;
    
    counter = 0;
    r_seed = 0;
    xseed = 0;
    x2seed = 0;

    // Default cost values from kwarg.c
    ctx.se_cost = -1.0;
    ctx.rm_cost = -1.0;
    ctx.r_cost = 1.0;
    ctx.rr_cost = 2.0;
    
    // Default temperature from kwarg.c
    ctx.Temp = 30.0;

    #ifdef ENABLE_VERBOSE
        ctx.howverbose = 1;
    #else
        ctx.howverbose = 0;
    #endif

    // Initialization of static selection variables from kwarg.c
    ctx._ms_w = 0;
    ctx._ms_v = DBL_MAX;
    ctx._rs_w = -DBL_MAX + 1;
    ctx._rs_n = 0;
    ctx._prs_kT = 1.0; // From _pseudoenergy_random_select logic
    ctx._prs_Z = 0;
    ctx._prs_offset = 0;

    // Initialization of static variables from exact.c
    ctx.exact_randomise = 0;
    ctx.reusable = 0;
    ctx.skip_lookup = 0;
    ctx._coalesce_compatibleandentangled_states = NULL;
    ctx._greedy_rmin = -1;
    ctx._greedy_currentstate = NULL;
    ctx._am = 0.0;
    ctx._choice_fixed = 0;
    ctx._greedy_choice = NULL;
    ctx.sc_min = DBL_MAX;
    ctx.sc_max = 0;
    ctx._lb = 0.0;
    ctx._predecessors = NULL;
    
    int T_in = 0, cost_in = 0;
    double T_array[100] = {30};
    double se_costs[100] = {0};
    double rm_costs[100] = {0};
    double r_costs[100] = {0};
    double rr_costs[100] = {0};
    
    #ifdef ENABLE_VERBOSE
    set_verbose(1);
    #endif
    
    /* Analyse command line options */
    #define KWARG_OPTIONS "V:b::R:P:kofhH?"
    
    /* Parse command line options */
    while ((i = getopt(argc, argv, KWARG_OPTIONS)) >= 0){
        switch(i){
            case 'V':
                ctx.howverbose = strtol(optarg, &endptr, 10);
                if(errno != 0 || *endptr != '\0') {
                    fprintf(stderr, "Verbosity input should be 0, 1 or 2.\n");
                    exit(1);
                }
                if(ctx.howverbose > 2 && ctx.howverbose < 0) {
                    fprintf(stderr, "Verbosity input should be 0, 1 or 2.\n");
                    exit(1);
                }
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
                        exit(1);
                    }
                }
                break;
            case 'R':
                R = strtod(optarg, &endptr);
                break;
            case 'P':
                P = strtod(optarg, &endptr);
                break;
            case 'k':
                gene_knownancestor = 1;
                break;
            case 'o':
                format = GENE_BEAGLE;
                break;
            case 'f':
                format = GENE_FASTA;
                break;
            case 'L':
                reference = strtol(optarg, &endptr, 10);
                if(errno != 0 || *endptr != '\0') {
                    fprintf(stderr, "Reference should be a positive integer.\n");
                    exit(1);
                }
                if(reference < 0) {
                    fprintf(stderr, "Reference should be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'h':
            case 'H':
            case '?':
                _print_usage(stdout, argv[0]);
                /* Clean up */
                DestroyLList(history_files);
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
    if ((gene_knownancestor) && (seqtype != GENE_BINARY))
        /* First sequnce only included to specify known common ancestor */
        remove_annotatedgene(a, 0);
    g = a->g;
    
    printf("Recombining sequence %d before site %d\n", R, P);
    
    /* Set up structures for computation */
    if ((Length(history_files) > 0)) {
        ctx.eventlist = MakeLList();
        multruns = 0;
    }
    
    T_in = 1;
    cost_in = 1;
    
    ctx.Temp = T_array[l];
    select = NULL;
        
    ctx.se_cost = se_costs[k];
    ctx.rm_cost = rm_costs[k];
    ctx.r_cost = r_costs[k];
    ctx.rr_cost = rr_costs[k];
        
    /* Initialise random number generator */
    initialise_x2random(r_seed);
    counter = 0;
    
    // Copy the data and set up the tracking lists
    h = copy_genes(g);
    ctx.seq_numbering = h->n;
    ctx.elements = elist_make();
    ctx.sites = elist_make();
    // Initialise list of sequences
    if ((gene_knownancestor) && (seqtype != GENE_BINARY)) {
        for(i=0; i < h->n; i++) {
            elist_append(ctx.elements, (void *)(i+1));
        }
    } else {
        for(i=0; i < h->n; i++) {
            elist_append(ctx.elements, (void *)i);
        }
    }
    // Initialise the list of sites
    for(i=0; i < h->length; i++) {
        elist_append(ctx.sites, (void *)i);
    }
    
    // Do the recombination event
    fprintf(fp, "input_data\n");
    output_genes(g, fp, NULL);
    if(R >= 0 & P > 0) {
        split(h, R, P, &ctx);
    }
    
    // Get a history
    tic = clock();
    n = output_coalescences(h, print_progress, &ctx, fp);
    toc = clock();
    timer = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("Time taken: %15.8f\n", timer);
    
    free_genes(h);
    ctx.elements = NULL;
    ctx.sites = NULL;
    r_seed = 0;
            
    /* Output inferred ARG */
    if ((Length(history_files) > 0)){
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
            arg = eventlist2history(a, fp, &ctx);
        }
        if (arg == NULL)
            arg = eventlist2history(a, NULL, &ctx);
        if (arg != NULL){
            arg_destroy(arg);
        }
    
    
        if (ctx.eventlist != NULL){
            while (Length(ctx.eventlist) > 0)
                free(Pop(ctx.eventlist));
            DestroyLList(ctx.eventlist);
        }
    }
        
    /* Clean up */
    if (ctx.lookup != NULL){
        elist_destroy(ctx.lookup);
    }
    
    if (ctx._greedy_beaglereusable != NULL) {
        beagle_deallocate_hashtable(ctx._greedy_beaglereusable);
        ctx._greedy_beaglereusable = NULL;
    }
    if (ctx._greedy_functioncalls != NULL) {
        hashtable_destroy(ctx._greedy_functioncalls, free, NULL, free);
        ctx._greedy_functioncalls = NULL;
    }
    
    DestroyLList(history_files);
    free_annotatedgenes(a);
    
    return 0;
}

