/***************************************************************************
 * 
 *    kwarg.c
 *  
 *    Implementation of front end for a greedy heuristic for finding plausible
 *    evolutionary histories, minimising number of recombinations and/or recurrent
 *    mutations (depending on selected cost of these events).
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
    pretty_print(f, "The program reads data from the input file specified and greedily constructs a history with a low number of recombinations and recurrent mutations (homoplasies). The history is constructed by stepping backwards in time using coalescence, mutation and recombination events. At each point in this process, all possible next events (strictly speaking, only a useful subset of all possible next events) are considered, and the resulting ancestral states are scored. The scores are used to choose an event either at random or to proceed to an ancestral state with minimum score (see option -T). This process is NOT guaranteed to lead to a history with a minimum number of recombinations or the minimum number of homoplasies.", 70, 0);
    fprintf(f, "Legal options are:\n");
    print_option(f, "-L[x]", "Provide an optional label x (should be an integer) to print at the start of each line.", 70, -1);
    print_option(f, "-S[x]", "Specify cost of a sequencing error (default: x = 0.5).", 70, -1);
    print_option(f, "-M[x]", "Specify cost of a recurrent mutation (default: x = 0.9).", 70, -1);
    print_option(f, "-R[x]", "Specify cost of a single recombination (default: x = 1.0).", 70, -1);
    print_option(f, "-C[x]", "Specify cost of two recombinations immediately following each other (default: x = 2.0).", 70, -1);
    print_option(f, "-T[F]", "Set annealing temperature 0 < F < 700 (default: F = 30). Scores are normalised to lie between 0 and 1, and the next step selected with probability proportional to exp(F * normalised_score). Setting F = 0 corresponds to random selection among possible moves (not recommended). Setting F = -1 corresponds to setting F = Infinity, and will force selection of the move with the minimum score at each step.", 70, -1);
    print_option(f, "-V[x]", "If running a single iteration with given cost parameters, this controls the level of verbosity. \nx = 0: no extra output \nx = 1: during each neighbourhood search, output the number of neighbours explored, the move selected and its cost \nx = 2: during each neighbourhood search, output the number of neighbours explored, the resulting configuration and cost of each neighbour, the move selected and its cost.", 70, -1);
    print_option(f, "-b[name]", "Output a minimum recombination history to file name.", 70, -1);
    print_option(f, "-d[name]", "Output ancestral recombination graph of minimum recombination history in dot format to file name.", 70, -1);
    print_option(f, "-g[name]", "Output ancestral recombination graph of minimum recombination history in GDL format to file name.", 70, -1);
    print_option(f, "-j[name]", "Output ancestral recombination graph of minimum recombination history in GML format to file name.", 70, -1);
    print_option(f, "-t[name]", "Output list of marginal phylogenies for each site in Newick's 8:45 format to file name.", 70, -1);
    print_option(f, "-D[name]", "Output list of marginal phylogenies for each site in dot format to file name.", 70, -1);
    print_option(f, "-G[name]", "Output list of marginal phylogenies for each site in GDL format to file name.", 70, -1);
    print_option(f, "-J[name]", "Output list of marginal phylogenies for each site in GML format to file name.", 70, -1);
    print_option(f, "-I", "Marginal trees are only output one for each of the intervals between two recombination points, instead of one for each site.", 70, -1);
    print_option(f, "-v[nodelabel]", "Use nodelabel convention for labelling nodes in ancestral recombination graphs. The possible conventions are: nodelabel = \n'none': do not label nodes\n'id': only label nodes representing sampled sequences, using their sequence ids from the data file, and nodes representing recombinations, indicating the recombination point\n'sequence': label nodes with the inferred sequences; these sequences will be in binary format, with 0 representing wild type and 1 representing mutant type, even if the original data is not in binary format\n'both': label nodes with both id and inferred sequence\n'one': use only one label for a node, sequence id or recombination point if available and otherwise the inferred sequence\nDefault convention is id. The colour coding scheme used for the nodes is red for sequences in the input data, blue for recombination nodes, green for standard coalescent nodes, and yellow for the final coalescence into the most recent common ancestor.", 70, -1);
    print_option(f, "-i", "Sequences not having a sequence id in the data file are assigned their index in the data file as id, e.g. the first sequence in the data file would be assigned '1' as id.", 70, -1);
    print_option(f, "-e", "Label edges in ancestral recombination graphs with the sites undergoing mutation along the edge.", 70, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 70, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 70, -1);
    print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known (see option -k), the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 70, -1);
    print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known (see option -k), the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 70, -1);
    print_option(f, "-k", "Assume that the common ancestral sequence is known, i.e. that we know which is the wild type and which is the mutant in each site. If the data is in binary format, the all-0 sequence is assumed to be the common ancestral sequence (this does not need to be present in the data). If the data is in amino acid or nucleotide format, the common ancestral sequence has to be specified directly and is taken to be the first sequence in the data file (see options -a and -n)", 70, -1);
    print_option(f, "-Q[x]", "Sets the number of runs.", 70, -1);
    print_option(f, "-Z[x]", "Sets the random seed z (only one run is made in this case).", 70, -1);
    print_option(f, "-X[x]", "Provide an upper bound x on the number of recombinations needed for the input dataset (solutions with more than x recombinations will be abandoned).", 70, -1);
    print_option(f, "-s", "Turns off the header row of the results table.", 70, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 70, -1);
}

/* Use the random function to draw an integer between 0 and n - 1 */
static int _unbiased_random(int n)
{
    long int l = XRAND_MAX / n;
    long int i;
    
    do{
        i = x2random() / l;
    } while (i >= n); /* i ought to always be at most n, but just to make sure */
    
    return i;
}

/* Determine whether to select a value according to a random minimal
 * value scheme. The minimum value seen so far is maintained in _ms_v
 * and the number of times this has been encountered is maintained in
 * _ms_w. Return value indicates whether the value should be selected.
 */
static int _ms_w = 0;
static double _ms_v = DBL_MAX;
static int _minimum_select(double a)
{
    if (a == _ms_v)
        /* The two values are equal - choose one at random and increment
         * number of times we have seen this value.
         */
        return (_unbiased_random(++_ms_w) == 0);
    else if (a < _ms_v){
        /* The value of a is new minimum - reset count and report this */
        _ms_v = a;
        _ms_w = 1;
        return 1;
    }
    
    /* The old value is still the minimum */
    return 0;
}

/* Determine whether to select a value according to a scheme where
 * values are selected with probability proportional to their
 * value. The sum of values seen so far is maintained in _rs_w. Values
 * less than zero or truncated to zero (as long as all values seen so
 * far are non-positive a least negative will be the
 * selection). Return value indicates whether the value should be
 * selected.
 */
static double _rs_w = -DBL_MAX + 1;
static int _rs_n = 0;
static int _random_select(double a)
{
    /* Check for non-positive values */
    if (a <= 0){
        if (_rs_w < a){
            _rs_w = a;
            _rs_n = 1;
            return 1;
        }
        else if (_rs_w == a)
            return (_unbiased_random(++_rs_n) == 0);
        else
            return 0;
    }
    
    /* Update _rs_w and perform random selection */
    if (_rs_w < 0)
        _rs_w = a;
    else
        _rs_w += a;
    
    return (a * XRAND_MAX > _rs_w * x2random());
}

/* Determine whether to select a value according to a scheme where
 * values are interpreted as energies and selected according to the
 * corresponding Boltzmann distribution at temperature _prs_kT / k.
 * The partition function so far is maintained in _srs_w and values
 * are shifted such that this is maintained to be close to 1. The
 * temperature is required to be positive. Return value indicates
 * whether the value should be selected.
 */
static double _prs_kT = 1;
static double _prs_Z = 0;
static double _prs_offset = 0;
static int _pseudoenergy_random_select(double a)
{
    if (_prs_Z == 0){
        /* First value seen */
        /* Choose offset such that partition function initially is 1 */
        _prs_offset = a;
        _prs_Z = 1;
        return 1;
    }
    else{
        /* Include contribution of a in partition function and update
         * offset if necessary.
         */
        if (a < _prs_offset){
            /* It's a good idea to change offset before proceeding */
            _prs_Z = exp((a - _prs_offset) / _prs_kT) * _prs_Z + 1;
            _prs_offset = a;
        }
        else
            /* Update partition function */
            _prs_Z += exp((_prs_offset - a) / _prs_kT);
        /* Check whether we should change the offset after the fact */
        if (_prs_Z > 2){
            /* By the precheck, inclusion of a can at most increase the
             * value of the partition function by 1 so reducing by a factor
             * of e should be sufficient.
             */
            _prs_offset -= _prs_kT;
            _prs_Z = _prs_Z / M_E;
        }
    }
    
    /* Determine whether to choose a */
    if (exp((_prs_offset - a) / _prs_kT) * XRAND_MAX < x2random() * _prs_Z)
        return 0;
    else
        return 1;
}

/* Reset variables for the various selection functions */
static void _reset_selections()
{
    _ms_w = 0;
    _ms_v = DBL_MAX;
    _rs_w = -DBL_MAX + 1;
    _rs_n = 0;
    _prs_Z = 0;
    _prs_offset = 0;
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
    int i, j = 0, k = 0, l = 0, m = 0, t = 0,
    head = 1,
    silent = 0,
    intervals = 0,
    multruns = 0,
    costs_in = 0;
    double n;
    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;
    FILE *print_progress = stdout;
    int (*select)(double) = _random_select;
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
    int ontheflyselection = 0;
    gc_enabled = 0;
    Event *e;
    LList *tmp;
    r_seed = 0;
    rec_max = INT_MAX;
    char *token;
//     int gc_ind = 0;
    double timer;
    clock_t tic, toc;
    char *endptr;
    errno = 0;
    reference = -1;
    
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
    #define KWARG_OPTIONS "S:M:R:C:T:V:X:b::d::g::j::t::D::G::J::Iv:iekofaZ:Q:sL:nhH?"
    
    /* Parse command line options */
    while ((i = getopt(argc, argv, KWARG_OPTIONS)) >= 0){
        switch(i){
            case 'S':
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        se_costs[j] = strtod(token, &endptr);
                        if(errno != 0 || *endptr != '\0') {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        if(se_costs[j] <= 0 && se_costs[j] != -1) {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        j++;
                        token = strtok(NULL, ",");
                    }
                }
                break;
            case 'M':
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        rm_costs[k] = strtod(token, &endptr);
                        if(errno != 0 || *endptr != '\0') {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        if(rm_costs[k] <= 0 && rm_costs[k] != -1) {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        k++;
                        token = strtok(NULL, ",");
                    }
                }
                break;
            case 'R':
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        r_costs[l] = strtod(token, &endptr);
                        if(errno != 0 || *endptr != '\0') {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        if(r_costs[l] <= 0 && r_costs[l] != -1) {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        l++;
                        token = strtok(NULL, ",");
                    }
                }
                break;
            case 'C':
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        rr_costs[m] = strtod(token, &endptr);
                        if(errno != 0 || *endptr != '\0') {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        if(rr_costs[m] <= 0 && rr_costs[m] != -1) {
                            fprintf(stderr, "Cost input should be a positive number.\n");
                            exit(1);
                        }
                        m++;
                        token = strtok(NULL, ",");
                    }
                }
                break;
            case 'T':
                /* Set new temperature */
                if (optarg != 0){
                    token = strtok(optarg, ",");
                    while(token != NULL) {
                        T_array[T_in] = strtod(token, &endptr);
                        if(T_array[T_in] != -1 && (errno != 0 || *endptr != '\0')) {
                            fprintf(stderr, "Annealing temperature can be -1, 0 or a positive real number.\n");
                            exit(1);
                        }
                        if(T_array[T_in] < 0 && T_array[T_in] != -1) {
                            fprintf(stderr, "Annealing temperature can be -1, 0 or a positive real number.\n");
                            exit(1);
                        }
                        if(T_array[T_in] >= log(DBL_MAX)) {
                            fprintf(stderr, "Warning: annealing temperature %.2f too high, setting to -1.\n", T_array[T_in]);
                            T_array[T_in] = -1;
                        }
                        T_in++;
                        token = strtok(NULL, ",");
                    }
                }
                break;
			case 'V':
                howverbose = strtol(optarg, &endptr, 10);
                if(errno != 0 || *endptr != '\0') {
                    fprintf(stderr, "Verbosity input should be 0, 1 or 2.\n");
                    exit(1);
                }
                if(howverbose > 2 && howverbose < 0) {
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
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(dot_files, (void *)optarg);
                }
                else
                    Enqueue(dot_files, stdout);
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
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gdl_files, (void *)optarg);
                }
                else
                    Enqueue(gdl_files, stdout);
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
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gml_files, (void *)optarg);
                }
                else
                    Enqueue(gml_files, stdout);
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
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(dottree_files, (void *)optarg);
                }
                else
                    Enqueue(dottree_files, stdout);
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
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gdltree_files, (void *)optarg);
                }
                else
                    Enqueue(gdltree_files, stdout);
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
                        exit(1);
                    }
                    fclose(fp);
                    Enqueue(gmltree_files, (void *)optarg);
                }
                else
                    Enqueue(gmltree_files, stdout);
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
                    _print_usage(stderr, argv[0]);
                    exit(1);
                }
                break;
            case 'i':
                generate_id = 1;
                break;
            case 'e':
                edgelabel = 1;
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
            case 'a':
                seqtype = GENE_AMINO;
                break;
            case 'n':
                seqtype = GENE_NUCLEIC;
                break;
            case 'Z':
                r_seed = strtod(optarg, &endptr);
                if(errno != 0 || *endptr != '\0') {
                    fprintf(stderr, "Seed input should be a positive integer.\n");
                    exit(1);
                }
                if(r_seed <= 0) {
                    fprintf(stderr, "Seed input should be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'Q':
                multruns = strtol(optarg, &endptr, 10) - 1;
                if(errno != 0 || *endptr != '\0') {
                    fprintf(stderr, "Number of iterations should be a positive integer.\n");
                    exit(1);
                }
                if(multruns < 0) {
                    fprintf(stderr, "Number of iterations should be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'X':
                rec_max = strtol(optarg, &endptr, 10);
                if(errno != 0 || *endptr != '\0') {
                    fprintf(stderr, "Upper bound on number of recombinations should be a positive integer.\n");
                    exit(1);
                }
                break;
            case 's':
                head = 0;
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
                DestroyLList(dot_files);
                DestroyLList(gml_files);
                DestroyLList(gdl_files);
                DestroyLList(tree_files);
                DestroyLList(dottree_files);
                DestroyLList(gmltree_files);
                DestroyLList(gdltree_files);
                DestroyLList(history_files);
                exit(0);
            case ':':
                _print_usage(stderr, argv[0]);
                exit(1);
        }
    }
    
    if(j != k) {
        fprintf(stderr, "Wrong number of costs specified.\n");
        exit(1);
    }
    else {
        cost_in = j;    
    }
    
    if((l > 0 || m > 0 ) && (l != m || l != j)) {
        fprintf(stderr, "Wrong number of costs specified.\n");
        exit(1);
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

    
    /* Set up structures for computation */
    if ((Length(history_files) > 0) || (Length(dot_files) > 0)
        || (Length(gml_files) > 0) || (Length(gdl_files) > 0)
        || (Length(tree_files) > 0) || (Length(dottree_files) > 0)
        || (Length(gmltree_files) > 0) || (Length(gdltree_files) > 0)) {
        eventlist = MakeLList();
        multruns = 0;
	cost_in = 1;
	T_in = 1;
    }
    
    
    initialise_xrandom();
    
    
    if(r_seed > 0) {
        multruns = 0;
        T_in = 1;
        cost_in = 1;
        se_costs[0] = (se_costs[0] !=0 ? se_costs[0] : 0.5);
        rm_costs[0] = (rm_costs[0] !=0 ? rm_costs[0] : 0.9);
        r_costs[0] = (r_costs[0] !=0 ? r_costs[0] : 1.0);
        rr_costs[0] = (rr_costs[0] != 0 ? rr_costs[0] : 2.0);
        if(howverbose > 0) {
            head = 0;
        }
    }
    else if(cost_in > 0) {
        if(multruns > 0 || cost_in > 1) {
            howverbose = 0;
        }
        for(t = 0; t < cost_in; t++) {
            se_costs[t] = (se_costs[t] !=0 ? se_costs[t] : 0.5);
            rm_costs[t] = (rm_costs[t] !=0 ? rm_costs[t] : 0.9);
            r_costs[t] = (r_costs[t] !=0 ? r_costs[t] : 1.0);
            rr_costs[t] = (rr_costs[t] != 0 ? rr_costs[t] : 2.0);
        }
    }
    else {
        howverbose = 0;
        cost_in = 13;
        double template1[13] = {-1, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01, 1};
        double template2[13] = {-1, 1.01, 0.91, 0.81, 0.71, 0.61, 0.51, 0.41, 0.31, 0.21, 0.11, 0.02, 1.1};
        double template3[13] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1};
        double template4[13] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -1};
        set_array(se_costs, template1, cost_in, 0);
        set_array(rm_costs, template2, cost_in, 0);
        set_array(r_costs, template3, cost_in, 0);
        set_array(rr_costs, template4, cost_in, 0);
    }
    
    if(T_in == 0) {
        T_in = 1;
    }
    
    if(head) {
        if(reference > 0) {
            fprintf(print_progress, "%10s %13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Ref", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
            "SE", "RM", "R", "N_states", "Time");
        }
        else {
            fprintf(print_progress, "%13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
            "SE", "RM", "R", "N_states", "Time");
        }
    }
    
    // Create the lookup array if will have multiple runs
    if(multruns > 0 || cost_in > 1){  
        if(rec_max == INT_MAX) {
            rec_max = 1000; // TODO What should this be
        }
        // Will store SE + RM in a lookup list with index = number of recombinations
        lookup = elist_make();
        for(i = 0; i <= rec_max; i++) {
            elist_append(lookup, (void *)INT_MAX);
        }
        // We need 0 se/rm for rec_max recombinations
        elist_change(lookup, rec_max, (void *)0);
    }
    
    for(l = 0; l < T_in; l++) {
        
        Temp = T_array[l];
        if(Temp == -1) {
            select = _minimum_select;
        }
        else {
            select = _random_select;
        }
        
        for(k = 0; k < cost_in; k++) {
            
            se_cost = se_costs[k];
            rm_cost = rm_costs[k];
            r_cost = r_costs[k];
            rr_cost = rr_costs[k];
            if(se_cost == -1 && rm_cost == -1 && r_cost == -1 && rr_cost == -1) {
                fprintf(stderr, "At least one type of event should be allowed (all event costs are -1).\n");
                exit(1);
            }

            /* Find a history for each iteration */
            for(j = 0; j <= multruns; j++) {
                
                /* Initialise random number generator */
                initialise_x2random(r_seed);
                counter = 0;
                
                // Copy the data and set up the tracking lists
                h = copy_genes(g);
                seq_numbering = h->n;
                elements = elist_make();
                sites = elist_make();
                // Initialise list of sequences
                if ((gene_knownancestor) && (seqtype != GENE_BINARY)) {
                    for(i=0; i < h->n; i++) {
                        elist_append(elements, (void *)(i+1));
                    }
                } else {
                    for(i=0; i < h->n; i++) {
                        elist_append(elements, (void *)i);
                    }
                }
                // Initialise the list of sites
                for(i=0; i < h->length; i++) {
                    elist_append(sites, (void *)i);
                }
                
                // Get a history
                tic = clock();
                n = ggreedy(h, print_progress, select, _reset_selections, ontheflyselection);
                toc = clock();
                timer = (double)(toc - tic) / CLOCKS_PER_SEC;
                printf("%15.8f\n", timer);
                // The ggreedy function will update rec_max and the lookup array
                
                // Tidy up for the next run
                free_genes(h);
                elist_destroy(elements);
                elements = NULL;
                elist_destroy(sites);
                sites = NULL;
                r_seed = 0;
                
            }

        }
        
    }

    
    /* Output inferred ARG */
    if ((Length(history_files) > 0) || (Length(dot_files) > 0)
        || (Length(gml_files) > 0) || (Length(gdl_files) > 0)
        || (Length(tree_files) > 0) || (Length(dottree_files) > 0)
        || (Length(gmltree_files) > 0) || (Length(gdltree_files) > 0)){
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
        if (lookup != NULL){
            elist_destroy(lookup);
        }
        
        if (_greedy_beaglereusable != NULL) {
            beagle_deallocate_hashtable(_greedy_beaglereusable);
            _greedy_beaglereusable = NULL;
        }
        if (_greedy_functioncalls != NULL) {
            hashtable_destroy(_greedy_functioncalls, free, NULL, free);
            _greedy_functioncalls = NULL;
        }
        
        DestroyLList(dot_files);
        DestroyLList(gml_files);
        DestroyLList(gdl_files);
        DestroyLList(tree_files);
        DestroyLList(dottree_files);
        DestroyLList(gmltree_files);
        DestroyLList(gdltree_files);
        DestroyLList(history_files);
        free_annotatedgenes(a);
        
        return 0;
}
