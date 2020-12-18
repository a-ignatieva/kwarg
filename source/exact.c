/*******************************************************************
 *
 *    exact.c
 *
 *    Implementation of functions to compute the exact minimum
 *    number of recombinations required under the infinite sites
 *    assumption for an SNP data set (Beagle), and a heuristic method to 
 *    obtain a (near-)minimal ARG in the presence of both recombination
 *    and recurrent mutation (KwARG).
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "gene.h"
#include "common.h"
#include "mergesort.h"
#include "bounds.h"
#include "llist.h"
#include "elist.h"
#include "exact.h"
#include "hashtable.h"
#include "bitfunctions.h"
#include "backtrack.h"
#include "common.h"

/* Select method for computing a lower bound on the number of recombinations */
#ifdef BEAGLE_HAPLOTYPE
#define LOWERBOUND(g) haplotype_bound_genes(g)
#elif defined(BEAGLE_HAPLOTYPEHEURISTIC)
#define LOWERBOUND(g) haplotype_heuristic_genes(g, INT_MAX, INT_MAX, 1)
#else
#define LOWERBOUND(g) hudson_kaplan_genes(g)
#endif

typedef struct _BeagleSplitInformation {
    HistoryFragment *g;
    int am;
    char splits;
    char representative;
} BeagleSplitInformation;

/* Variables declared globally to avoid parameter clutter */
int exact_randomise = 0;    /* Random evolutionary histories? */
static int reusable = 0;    /* Is hash table reusable? */
static int skip_lookup = 0; /* Should LOOKUP event be allowed? */

/* The smaller element is the one with the smallest number of splits,
 * or if those are equal the element with the smaller amount of
 * ancestral material left.
 */
static int compareboundandancestralmaterial(BeagleSplitInformation *a,
        BeagleSplitInformation *b)
{
    if ((a->representative < b->representative)
            || ((a->representative == b->representative)
                && ((a->splits > b->splits)
                    || ((a->splits == b->splits) && (a->am < b->am)))))
        return 1;
    else
        return 0;
}

/* Use the random function to draw an integer between 0 and n - 1 */
static int unbiased_random(int n)
{
    long int l = XRAND_MAX / n;
    long int i;

    do {
        i = xrandom() / l;
    } while (i >= n); /* i ought to always be at most n, but just to make sure */

    return i;
}

static void permute(BeagleSplitInformation *splits, int n)
{
    int i, j;
    BeagleSplitInformation tmp;

    for (i = n; i > 1;) {
        j = unbiased_random(i);
        i--;
        if (j != i) {
            memcpy(&tmp, splits + i, sizeof(BeagleSplitInformation));
            memcpy(splits + i, splits + j, sizeof(BeagleSplitInformation));
            memcpy(splits + j, &tmp, sizeof(BeagleSplitInformation));
        }
    }
}

/* Permute the elements in elist */
static void permute_elist(EList *elist)
{
    int i, j;

    for (i = elist_length(elist); i > 1;) {
        j = unbiased_random(i);
        i--;
        elist_swap(elist, i, j);
    }
}

/* Transfer the ancestral states in genes to the array splits (that
 * should be sufficiently large to hold all transfered ancestral
 * states). For each ancestral state it is checked whether it has
 * already been investigated or whether a computed lower bound on the
 * number of recombinations required plus the base number of
 * recombinations already used exceeds target. The total number of
 * ancestral states transferred is added to n.
 */
static int transfer2splitinformation(BeagleSplitInformation *splits,
                                     EList *genes, int base, int *n,
                                     HashTable *t, int target)
{
    int i, prevtarget, bound;
    void *lookup;
    Genes *g;
    PackedGenes *p;
    Event *e;

    /* Insert all feasible branches from genes into splits */
    for (i = 0; i < elist_length(genes); i++) {
        /* Check whether we have already encountered this ancestral state
         * with the same, or larger, target.
         */
        splits[*n].g = elist_get(genes, i);
        g = splits[*n].g->g;

        p = pack_genes(g);
        if (hashtable_lookup((void *)p, t, &lookup)) {
            /* This ancestral state is already present in the hash table */
            prevtarget = (int)lookup;
            splits[*n].representative = 0;
            if (reusable && (prevtarget < 0) && (base - target <= prevtarget + 1)) {
                /* We know this set of sequences needs at most as many
                 * recombinations as we have left.
                 */
                free_packedgenes(p);
                if (!skip_lookup) {
                    /* And we don't want to search histories */
                    if (eventlist != NULL) {
                        Append(eventlist, splits[*n].g->event);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = LOOKUP;
                        e->event.lookup = -prevtarget - 1;
                        Enqueue(eventlist, (void *)e);
                    }
                    free(splits[*n].g);
                    /* We found a solution - eliminate all candidates */
                    free_genes(g);
                    for (i++; i < elist_length(genes); i++) {
                        splits[*n].g = elist_get(genes, i);
                        if (eventlist != NULL) {
                            while (Length(splits[*n].g->event) > 0)
                                free(Pop(splits[*n].g->event));
                            DestroyLList(splits[*n].g->event);
                        }
                        free_genes(splits[*n].g->g);
                        free(splits[*n].g);
                    }
                    for (i = 0; i < *n; i++) {
                        free_genes(splits[i].g->g);
                        if (eventlist != NULL) {
                            while (Length(splits[i].g->event) > 0)
                                free(Pop(splits[i].g->event));
                            DestroyLList(splits[i].g->event);
                        }
                        free(splits[i].g);
                    }
                    elist_destroy(genes);
                    free(splits);
                    return 1;
                }
            }
            else if ((reusable && ((2 * (target - base) + (base > 0) < prevtarget)
                                   || (base - target > prevtarget + 1)))
                     || (!reusable && (target - base <= prevtarget))) {
                /* We know this set of sequences needs more recombinations
                 * than we have left.
                 */
                if (eventlist != NULL) {
                    while (Length(splits[*n].g->event) > 0)
                        free(Pop(splits[*n].g->event));
                    DestroyLList(splits[*n].g->event);
                }
                free(splits[*n].g);
                free_genes(g);
                free_packedgenes(p);
                continue;
            }
            else {
                /* By the structure of the method, we know that no history
                 * exists with fewer than target - base recombinations.
                 */
                hashtable_update((void *)p,
                                 (void *)((target - base) << reusable), t, NULL);
                free_packedgenes(p);
            }
        }
        else {
            splits[*n].representative = 1;

            /* Compute lower bound for this sequence set */
            bound = LOWERBOUND(g);

            /* Is there any hope for this branch */
            if (bound > target - base) {
#ifdef BEAGLE_DONOTSTORELEAVES
                free_packedgenes(p);
                free_genes(g);
#else
                if (reusable)
                    hashtable_update((void *)p, (void *)(2 * bound), t, NULL);
                else
                    hashtable_update((void *)p, (void *)(bound - 1), t, NULL);
#endif
                if (eventlist != NULL) {
                    while (Length(splits[*n].g->event) > 0)
                        free(Pop(splits[*n].g->event));
                    DestroyLList(splits[*n].g->event);
                }
                free(splits[*n].g);
                continue;
            }

            /* By the structure of the method, we know that no history
             * exists with fewer than target - base recombinations.
             */
            hashtable_update((void *)p,
                             (void *)((target - base) << reusable), t, NULL);
        }

        splits[*n].am = ancestral_material(g);
        splits[*n].splits = base;

        /* Another sequence set inserted into splits */
        (*n)++;
    }

    elist_destroy(genes);
    return 0;
}

/* Free the memory used by the sequence sets stored in elist as well
 * as by elist itself.
 */
static void free_elist_elements(EList *elist)
{
    int i;
    HistoryFragment *s;

    for (i = 0; i < elist_length(elist); i++) {
        s = (HistoryFragment *)elist_get(elist, i);
        if ((eventlist != NULL) && (s->event != NULL)) {
            while (Length(s->event) > 0)
                free(Pop(s->event));
            DestroyLList(s->event);
        }
        free_genes(s->g);
        free(s);
    }

    elist_destroy(elist);
}

/* Check whether any of the set of sequences in genes can be explained
 * without using recombinations. If so, free the memory used by the
 * sequence sets stored in genes as well as by genes itself and return
 * 1; otherwise return 0.
 */
static int check_for_bottom(EList *genes)
{
    int i;
    Genes *g;
    HistoryFragment *s;

    for (i = 0; i < elist_length(genes); i++) {
        s = (HistoryFragment *)elist_get(genes, i);
        g = s->g;
        if (no_recombinations_required(g)) {
            /* We can explain the set of sequences with no further recombinations */
#ifdef ENABLE_VERBOSE
            if (verbose()) {
                printf("Bottom of recursion:\n");
                output_genes_indexed(g, NULL);
            }
#endif
            if (eventlist != NULL) {
                Append(eventlist, s->event);
                s->event = NULL;
            }
            free_elist_elements(genes);
            return 1;
        }
    }

    return 0;
}

/* Data structure for storing pair of a node and a class */
typedef struct _NodeClass {
    int node;
    int class;
} NodeClass;
/* The recursion actually implementing the enumeration of
 * _coalesce_compatibleandentangled. The nodes stored on stack still
 * needs to be considered, the nodes stored in component are in current
 * component, E is array of edge lists for each node in the entangled
 * graph, C is array of edge lists for each node in the compatible
 * graph, components is an array containing current assignment, and
 * states is the list new ancestral states are stored in.
 */
static void _coalesce_cande_recursion(LList *stack, EList *component, Genes *g,
                                      EList *E, int **C, int *components,
                                      void (*f)(Genes *g))
{
    int i, j, n = 0, old;
    Genes *h = NULL;
    NodeClass *current, *tmp;
    Event *event;
    LList *oldevents = eventlist;
    EList *oldelements = elements;
    EList *oldsites = sites;

    /* Check whether we have reached bottom of recursion */
    if (Length(stack) == 0) {
        /* Coalesce in reverse order to avoid having to keep track of
         * where other sequences are moved when sequences disappear by
         * coalescence.
         */
        for (i = g->n - 1; i >= 0; i--){
            if (components[i] - 1 != i) {
                /* Coalesce as specified by components */
                if (h == NULL) {
                    /* First coalescence encountered; create structure for
                     * carrying out the coalescences in.
                     */
                    h = copy_genes(g);
                    if (eventlist != NULL)
                        eventlist = MakeLList();
                    if(elements != NULL) {
                        elements = elist_make();
                        elist_safeextend(elements, oldelements);
                    }
                    if(sites != NULL) {
                        sites = elist_make();
                        elist_safeextend(sites, oldsites);
                    }
                }
                coalesce(h, components[i] - 1, i);
                if (eventlist != NULL) {
                    event = (Event *)xmalloc(sizeof(Event));
                    event->type = COALESCENCE;
                    event->event.c.s1 = components[i] - 1;
                    event->event.c.s2 = i;
                    Enqueue(eventlist, event);
                }
            }
        }
        if (h != NULL) {
            /* At least one coalescence carried out, leading to a new
             * ancestral state that should be pursued.
             */
            implode_genes(h);
            f(h);
            eventlist = oldevents;
            elements = oldelements;
            sites = oldsites;
        }
    }
    else {
        /* Still nodes left in the stack to visit */
        current = Pop(stack);
        if ((components[current->node] <= 0)
                && (components[current->node] != -current->class)) {
            /* No decision has yet been taken on whether to add current node
             * to current component.
             */
            /* Save old component value for restoration at exit */
            old = components[current->node];
            /* Check whether current node is compatible with all other nodes
             * in current component.
             */
            if (current->class - 1 != current->node) {
                /* Current node does not initiate new component */
                for (i = 0; i < elist_length(component); i++)
                    if (!C[current->node][(int)elist_get(component, i)])
                        break;
            }
            else {
                /* Current node does initiate new component */
                component = elist_make();
                i = 0;
            }
            if (i == elist_length(component)) {
                /* Add current node to current component */
                components[current->node] = current->class;
                elist_append(component, (void *)current->node);
                /* Add neighbours to stack */
                if (elist_length(E + current->node) > 0) {
                    tmp = (NodeClass *)
                          xmalloc(elist_length(E + current->node) * sizeof(NodeClass));
                    for (i = 0; i < elist_length(E + current->node); i++) {
                        j = (int)elist_get(E + current->node, i);
                        /* Check whether a decision has been already been mode for
                         * this neighbour to avoid sequences of superfluous pushes and
                         * pops - without this check the checks would still be carried
                         * out in the recursive calls but in some situations multiple
                         * times.
                         */
                        if ((components[j] <= 0) && (components[j] != current->class)) {
                            tmp[n].node = (int)elist_get(E + current->node, i);
                            tmp[n].class = current->class;
                            Push(stack, tmp + n);
                            n++;
                        }
                    }
                }
                _coalesce_cande_recursion(stack, component, g, E, C, components, f);
                /* Clean up */
                if (elist_length(E + current->node) > 0) {
                    for (i = 0; i < n; i++)
                        Pop(stack);
                    free(tmp);
                }
                if (current->class - 1 == current->node)
                    /* This node initiated a new component */
                    elist_destroy(component);
                else
                    /* Did not! */
                    elist_removelast(component);
            }

            /* Leave current node out of current component if component is
             * not node's own.
             */
            if (current->node != current->class - 1) {
                components[current->node] = -current->class;
                _coalesce_cande_recursion(stack, component, g, E, C, components, f);
            }

            /* Restore old components value for current node */
            components[current->node] = old;
        }
        else
            /* Decision already made for current, continue with rest of
             * nodes on stack.
             */
            _coalesce_cande_recursion(stack, component, g, E, C, components, f);

        /* Restore current to stack */
        Push(stack, current);
    }
    
}


/* Enumerate all partitions of sequences in g into sets that
 * constitute connected components in the graph with an edge between
 * two sequences if they have entangled ancestral material and
 * constitute cliques in the graph with an edge between two sequences
 * if they are compatible. For each partition, coalesce sequences that
 * are in the same set and apply f to the resulting
 * HistoryFragment. It is the responsibility of the calling function
 * to free memory used for the HistoryFragment.
 */
static void _coalesce_compatibleandentangled_map
(Genes *g, void (*f)(Genes *))
{
    int i, j;
    EList *component = elist_make();
    EList *E = (EList *)xmalloc(g->n * sizeof(EList));
    int **C = (int **)xmalloc(g->n * sizeof(int *));
    LList *stack = MakeLList();
    int *components = (int *)xcalloc(g->n, sizeof(int));
    NodeClass *tmp = xmalloc(g->n * sizeof(NodeClass));

    /* Start by constructing edge list for each node of compatible and
     * entangled graphs.
     */
    for (i = 0; i < g->n; i++) {
        elist_init(E + i);
        C[i] = (int *)xcalloc(g->n, sizeof(int));
        /* Also use loop to insert each node in stack with itself as class */
        tmp[i].node = i;
        tmp[i].class = i + 1;
        Enqueue(stack, tmp + i);
    }
    for (i = 0; i < g->n; i++)
        for (j = i + 1; j < g->n; j++)
            if (compatible(g, i, j)) {
                C[i][j] = C[j][i] = 1;
                /* There is no point in pursuing an entanglement edge if the
                 * two concerned sequences are not compatible; hence we put
                 * the edge addition to E inside the check of compatibility.
                 */
                if (entangled(g, i, j)) {
                    elist_append(E + i, (void *)j);
                    elist_append(E + j, (void *)i);
                }
            }

    /* Initiate enumeration of all splits into connected components */
    _coalesce_cande_recursion(stack, component, g, E, C, components, f);

    /* Clean up */
    for (i = 0; i < g->n; i++) {
        elist_free(E + i);
        free(C[i]);
    }
    free(E);
    free(C);
    elist_destroy(component);
    free(components);
    free(tmp);
    DestroyLList(stack);
}

/* Enumerate all partitions of sequences in g into sets that
 * constitute connected components in the graph with an edge between
 * two sequences if they have entangled ancestral material and
 * constitute cliques in the graph with an edge between two sequences
 * if they are compatible. For each partition, coalesce sequences that
 * are in the same set and insert all resulting HistoryFragments in a
 * list that is returned.
 */
static EList *_coalesce_compatibleandentangled_states;
static void _coalesce_compatibleandentangled_f(Genes *g)
{
    HistoryFragment *f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));

    f->g = g;
    f->event = eventlist;
    elist_append(_coalesce_compatibleandentangled_states, f);
}
static EList *_coalesce_compatibleandentangled(Genes *g)
{
    _coalesce_compatibleandentangled_states = elist_make();
    _coalesce_compatibleandentangled_map(g, _coalesce_compatibleandentangled_f);
    return _coalesce_compatibleandentangled_states;
}

/* Recursion actually implementing the branch&bound procedure of
 * beagle. The hash table t is used for the dynamic programming and
 * branches known to lead to more recombinations than target are cut
 * off. If try_coalesces is true, all sensible coalesces are tried;
 * otherwise, only events involving at least one split are tried. The
 * rationale behind this, is that calls from an invocation where we
 * try all sensible coalesces should not themselves try all sensible
 * coalesces. If reusable is true, the values with which states are
 * inserted in the hash table are such that the hash table can be
 * reused. Return value just states whether target can be met (or
 * bettered). The number of ancestral states visited is reduced by
 * only attempting splits that will allow coalesces that are in some
 * sense maximal.
 */
static int beagle_recursion(Genes *g, HashTable *t, int target,
                            int try_coalesces)
{
    int i, j, n;
    Index *start, *end;
    EList *coalesced, *prefix, *postfix, *infix, *overlap;
    BeagleSplitInformation *splits;
#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif
    PackedGenes *p;

    if (try_coalesces) {
        /* We have just imploded genes, but we still need to pursue paths
         * coalescing compatible sequences where neither is subsumed in the
         * other but where the ancestral material is still entangled.
         */
        coalesced = _coalesce_compatibleandentangled(g);
        n = elist_length(coalesced);
        if (check_for_bottom(coalesced)) {
#ifdef ENABLE_VERBOSE
            if (v) {
                printf("Reached by coalescing some compatible sequences from:\n");
                output_genes_indexed(g, NULL);
            }
            set_verbose(v);
#endif
            return 1;
        }
    }
    else
        n = 0;

    i = 0;
    if (target > 0) {
        /* There is room for at least one split */
        /* Determine interesting recombination ranges */
        start = maximumsubsumedprefixs(g);
        end = maximumsubsumedpostfixs(g);

        /* Try all sensible events with one split */
#ifdef ENABLE_VERBOSE
        set_verbose(v - 1);
#endif
        prefix = maximal_prefix_coalesces(g, start, end);
        if (!exact_randomise && check_for_bottom(prefix)) {
            free(start);
            free(end);
            if (try_coalesces)
                free_elist_elements(coalesced);
#ifdef ENABLE_VERBOSE
            if (v) {
                printf("Reached by coalescing a compatible prefix from:\n");
                output_genes_indexed(g, NULL);
            }
            set_verbose(v);
#endif
            return 1;
        }

        postfix = maximal_postfix_coalesces(g, start, end);
        /* If a random solution is required, we need to randomise the splits */
        if (exact_randomise) {
            elist_extend(postfix, prefix);
            permute_elist(postfix);
            n += elist_length(postfix);
        }
        else
            n += elist_length(prefix) + elist_length(postfix);

        if (check_for_bottom(postfix)) {
            free(start);
            free(end);
            if (try_coalesces)
                free_elist_elements(coalesced);
            if (!exact_randomise)
                free_elist_elements(prefix);
#ifdef ENABLE_VERBOSE
            if (v) {
                if (exact_randomise)
                    printf("Reached by coalescing a compatible prefix or postfix from:\n");
                else
                    printf("Reached by coalescing a compatible postfix from:\n");
                output_genes_indexed(g, NULL);
            }
            set_verbose(v);
#endif
            return 1;
        }

        if (target > 1) {
            /* Try all sensible events with two splits */
            infix = maximal_infix_coalesces(g, start, end);
            if (!exact_randomise && check_for_bottom(infix)) {
                free(start);
                free(end);
                if (try_coalesces)
                    free_elist_elements(coalesced);
                free_elist_elements(prefix);
                free_elist_elements(postfix);
#ifdef ENABLE_VERBOSE
                if (v) {
                    printf("Reached by coalescing a compatible infix from:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
            overlap = maximal_overlap_coalesces(g, start, end);
            /* If a random solution is required, we need to randomise the splits */
            if (exact_randomise) {
                elist_extend(overlap, infix);
                permute_elist(overlap);
                n += elist_length(overlap);
            }
            else
                n += elist_length(infix) + elist_length(overlap);

            if (check_for_bottom(overlap)) {
                free(start);
                free(end);
                if (try_coalesces)
                    free_elist_elements(coalesced);
                free_elist_elements(postfix);
                if (!exact_randomise) {
                    free_elist_elements(prefix);
                    free_elist_elements(infix);
                }
#ifdef ENABLE_VERBOSE
                if (v) {
                    if (exact_randomise)
                        printf("Reached by coalescing compatible overlaps or a compatible infix from:\n");
                    else
                        printf("Reached by coalescing compatible overlaps from:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }

            if (n > 0) {
                splits = (BeagleSplitInformation *)
                         xmalloc(n * sizeof(BeagleSplitInformation));
                if (!exact_randomise
                        && transfer2splitinformation(splits, infix, 2, &i, t, target)) {
                    free(start);
                    free(end);
                    if (try_coalesces)
                        free_elist_elements(coalesced);
                    free_elist_elements(prefix);
                    free_elist_elements(postfix);
                    free_elist_elements(overlap);
#ifdef ENABLE_VERBOSE
                    if (v) {
                        printf("Found resolved ancestral state in hash table, reached from coalescing\ncompatible overlaps in:\n");
                        output_genes_indexed(g, NULL);
                    }
                    set_verbose(v);
#endif
                    return 1;
                }
                if (transfer2splitinformation(splits, overlap, 2, &i, t, target)) {
                    free(start);
                    free(end);
                    if (try_coalesces)
                        free_elist_elements(coalesced);
                    if (!exact_randomise)
                        free_elist_elements(prefix);
                    free_elist_elements(postfix);
#ifdef ENABLE_VERBOSE
                    if (v) {
                        if (exact_randomise)
                            printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible infix or compatible overlaps in:\n");
                        else
                            printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible infix in:\n");
                        output_genes_indexed(g, NULL);
                    }
                    set_verbose(v);
#endif
                    return 1;
                }
            }
            else {
                if (!exact_randomise)
                    elist_destroy(infix);
                elist_destroy(overlap);
            }
        }
        else if (n > 0) {
            splits = (BeagleSplitInformation *)
                     xmalloc(n * sizeof(BeagleSplitInformation));
        }
        if (n > 0) {
            if (!exact_randomise
                    && transfer2splitinformation(splits, prefix, 1, &i, t, target)) {
                free(start);
                free(end);
                if (try_coalesces)
                    free_elist_elements(coalesced);
                free_elist_elements(postfix);
#ifdef ENABLE_VERBOSE
                if (v) {
                    printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible prefix in:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
            if (transfer2splitinformation(splits, postfix, 1, &i, t, target)) {
                free(start);
                free(end);
                if (try_coalesces)
                    free_elist_elements(coalesced);
#ifdef ENABLE_VERBOSE
                if (v) {
                    if (exact_randomise)
                        printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible prefix or postfix in:\n");
                    else
                        printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible postfix in:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
        }
        else {
            if (!exact_randomise)
                elist_destroy(prefix);
            elist_destroy(postfix);
        }

        /* Clean up */
        free(start);
        free(end);
    }
    else if (try_coalesces && (elist_length(coalesced) > 0))
        splits = (BeagleSplitInformation *)
                 xmalloc(elist_length(coalesced) * sizeof(BeagleSplitInformation));

    if (try_coalesces) {
        if (elist_length(coalesced) > 0) {
            if (transfer2splitinformation(splits, coalesced, 0, &i, t, target)) {
#ifdef ENABLE_VERBOSE
                if (v) {
                    printf("Found resolved ancestral state in hash table, reached by coalescing some\ncompatible sequences in:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
        }
        else
            elist_destroy(coalesced);
    }

    /* We now have all the ancestral states that could sensibly lead to
     * the current ancestral state g stored in splits. Sort them and
     * pursue each branch in sorted order.
     */
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif

    fflush(stdout);
    if (i > 0) {
        /* ...and there are actually some */
        if (exact_randomise)
            permute(splits, i);
        else
            merge_sort(splits, i, sizeof(BeagleSplitInformation),
                       (int (*)(void *, void *))compareboundandancestralmaterial);
#ifdef ENABLE_VERBOSE
        if (v > 1)
            for (j = 0; j < i; j++) {
                output_genes_indexed(splits[j].g->g, NULL);
                printf("is #%d; splits %d, and ancestral material %d\n",
                       j, splits[j].splits, splits[j].am);
            }
#endif
        for (j = 0; j < i; j++) {
            g = splits[j].g->g;
            if (beagle_recursion(g, t, target - splits[j].splits,
                                 (splits[j].splits > 0))) {
                /* We found a branch leading to at most target recombinations */
#ifdef ENABLE_VERBOSE
                if (v) {
                    printf("Obtained from %d splits in:\n", splits[j].splits);
                    output_genes_indexed(g, NULL);
                }
#endif
#ifdef DEBUG
                if (ancestral_state_trace != NULL) {
                    /* Insert current ancestral state as one visited in the
                     * minimum history we are in the process of returning back
                     * from.
                     */
                    p = pack_genes(g);
                    hashtable_insert(p, NULL, ancestral_state_trace);
                }
#endif
                if (reusable) {
                    p = pack_genes(g);
                    /* We know this set of sequences is present in the hash
                     * table, as it was at the latest inserted in
                     * transfer2splitinformation.
                     */
                    hashtable_update((void *)p,
                                     (void *)(splits[j].splits - target - 1), t, NULL);
                    free_packedgenes(p);
                }
                if (eventlist != NULL)
                    Prepend(splits[j].g->event, eventlist);
                free_genes(splits[j].g->g);
                free(splits[j].g);
                j++;
                for (; j < i; j++) {
                    free_genes(splits[j].g->g);
                    if (eventlist != NULL) {
                        while(Length(splits[j].g->event) > 0)
                            free(Pop(splits[j].g->event));
                        DestroyLList(splits[j].g->event);
                    }
                    free(splits[j].g);
                }
                free(splits);
                return 1;
            }
            else if (reusable) {
                p = pack_genes(g);
                hashtable_update((void *)p,
                                 (void *)(2 * (target - splits[j].splits)
                                          + try_coalesces + 1), t, NULL);
                free_packedgenes(p);
            }

            free_genes(splits[j].g->g);
            if (eventlist != NULL) {
                while(Length(splits[j].g->event) > 0)
                    free(Pop(splits[j].g->event));
                DestroyLList(splits[j].g->event);
            }
            free(splits[j].g);
        }
    }
    if (n > 0)
        free(splits);

    /* We cannot explain g with at most target recombinations. Return
     * this informative fact.
     */
    return 0;
}

#ifdef ENABLE_VERBOSE
static void print_gene(Genes *g, va_list args)
{
    output_genes_indexed(g, NULL);
}
#endif

#ifdef HAPLOTYPE_BLOCKS
/* Return maximum of a and b */
static int _intmax(int a, int b)
{
    return (a > b ? a : b);
}
#endif

/* Use branch&bound plus dynamic programming techniques to reconstruct
 * a history for g requiring a minimum number of recombinations. If
 * eventlist is not NULL, a list of the events leading to this number
 * of recombinations is compiled in eventlist. If haploblocks is not
 * NULL, local minimum number of recombinations are stored in
 * haploblocks; haploblocks is assumed to be a table initialised to 0s
 * - entry i, j is set to the minimum number of recombinations
 * required in the region from site i to site i + j + 1. It is assumed
 * that lower is a valid lower bound on the number of recombinations,
 * and the search is terminated if it is established that no history
 * with at most upper recombinations exists (and a number larger than
 * upper is returned). If global variable reusable is true, t should
 * not be NULL and it is used as initial hash table, otherwise a local
 * hash table is allocated and later deallocated for keeping track of
 * ancestral states encountered.
 */
static int beagle_core(Genes *g, FILE *print_progress, int lower, int upper,
                       HashTable *t)
{
    Genes *h;
    int table_size, bound = 0, r1 = 0, r2 = 0;
    void *lookup;
    LList *implode, *tmp = eventlist;
    Event *e;
#ifdef HAPLOTYPE_BLOCKS
    int i, j, n, localbound, oldreusable = reusable;
    LList *l, *tmprep = representativeness;
    LListCounter *tmprepcount = representativeness_counter;
    SuperColumn *c;
#endif
    PackedGenes *p;
#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif

    /* Check whether we can rule out the necessity for recombinations */
    if (no_recombinations_required(g))
        return 0;

    /* We cannot handle upper bounds larger than half of INT_MAX */
    if (upper >= INT_MAX / 2)
        upper = INT_MAX / 2 - 1;

    /* Check whether sequence set is already present in hash table */
    if (reusable) {
        p = pack_genes(g);
        if (hashtable_lookup((void *)p, t, &lookup)) {
            bound = (int)lookup;
            free_packedgenes(p);
            if (bound < 0) {
                /* We know the true minimum for this sequence set */
                if (!skip_lookup) {
                    /* And we are going to use it */
                    if ((eventlist != NULL) && (-bound <= upper)) {
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = LOOKUP;
                        e->event.lookup = -bound - 1;
                    }
                    return -bound - 1;
                }
            }
            else
                /* We have a lower bound for this sequence set */
                if (2 * upper + 1 < bound)
                    return upper + 1;
        }
        else
            free_packedgenes(p);
    }

    /* Compute a good lower bound on the number of recombinations */
    if (eventlist != NULL)
        eventlist = MakeLList();
#ifdef HAPLOTYPE_BLOCKS
    if (haploblocks != NULL) {
        n = g->length;
        representativeness = MakeLList();
        for (i = 0; i < g->length; i++) {
            c = (SuperColumn *)xmalloc(sizeof(SuperColumn));
            c->left = c->right = i;
            Enqueue(representativeness, (void *)c);
        }
        representativeness_counter = MakeCounter(representativeness, FIRST);
    }
#endif
    /* Create working copy of g */
    g = copy_genes(g);
    implode_genes(g);

#ifdef DEBUG
    /* Insert initial ancestral state as one that was visited */
    if ((ancestral_state_trace != NULL) && (g->n > 0)) {
        p = pack_genes(g);
        hashtable_insert(p, NULL, ancestral_state_trace);
    }
#endif

    implode = eventlist;
    eventlist = tmp;
#ifdef HAPLOTYPE_BLOCKS
    if (haploblocks != NULL) {
        l = representativeness;
        representativeness = NULL;
        free(representativeness_counter);
    }
#endif
    if (!no_recombinations_required(g)) {
        bound = LOWERBOUND(g);

        if (!reusable) {
            /* Determine size of and allocate hash table */
            table_size = msb((g->n - 3) * g->length) + bound;
            t = new_packedgeneshashtable(table_size);
        }

#ifdef HAPLOTYPE_BLOCKS
        if (haploblocks != NULL)
            reusable = 1;
#endif
#ifdef ENABLE_VERBOSE
        set_verbose(v);
#endif
        /* Determine minimum number of recombinations required for g */
        if (bound < lower)
            bound = lower;
        if (upper < 0)
            upper = INT_MAX;
        if (print_progress != NULL)
            fprintf(print_progress, "At least %d recombination%s required\n", bound,
                    (bound != 1 ? "s" : ""));
        p = pack_genes(g);
        if (hashtable_update((void *)p, (void *)(bound << reusable),
                             t, NULL) < 0) {
            r1 = 1;
        }
        for (; (bound <= upper) && !beagle_recursion(g, t, bound, 1);)
        {
            bound++;
#ifdef ENABLE_VERBOSE
            if (v)
                printf("%d states examined with bound %d\n", hashtable_size(t),
                       bound - 1);
            else if (print_progress != NULL)
                fprintf(print_progress, "At least %d recombinations required\n",
                        bound);
#else
            if (print_progress != NULL)
                fprintf(print_progress, "At least %d recombinations required\n",
                        bound);
#endif
            hashtable_update((void *)p, (void *)(bound << reusable),
                             t, NULL);
        }

        if (!r1)
            free_packedgenes(p);
        if ((eventlist != NULL) && (bound > upper))
            /* We failed to find a valid history */
            while (Length(implode) > 0)
                free(Pop(implode));
#ifdef ENABLE_VERBOSE
        if (v) {
            printf("%d states examined with true minimum (%d)\n", hashtable_size(t),
                   bound);
            if (v > 1) {
                printf("Hash modulo was %ld and number of collisions %d\n", t->size,
                       hashtable_collisions(t));
                hashtable_printlargestbucket(t, (void (*)(void *, va_list))print_gene);
                printf("One largest bucket was %d\n", hashtable_largestbucket(t));
            }
        }
#endif

#ifdef HAPLOTYPE_BLOCKS
        if (haploblocks != NULL) {
            /* We already computed minimum number of recombinations for full
             * sequence set.
             */
            haploblocks[0][g->length - 2] = bound;
            /* Do not remember events from local bound computations */
            eventlist = NULL;
            /* Run through regions in order of increasing length */
            for (i = 2; i < g->length; i++)
                for (j = 0; j <= g->length - i; j++) {
                    h = copy_region(g, j, j + i);
                    implode_genes(h);
                    if (!no_recombinations_required(h)) {
                        p = pack_genes(h);
                        r2 = 0;
                        lookup = (void *)0;
                        if (!hashtable_lookup((void *)p, t, &lookup)
                                || ((int)lookup >= 0)) {
                            localbound = (int)lookup / 2;
                            if (i > 2)
                                localbound = _intmax(localbound,
                                                     _intmax(haploblocks[j][i - 3],
                                                             haploblocks[j + 1][i - 3]));
                            if (hashtable_update((void *)p, (void *)(localbound << reusable),
                                                 t, NULL) < 0)
                                r2 = 1;
                            for (; (localbound <= upper)
                                    && !beagle_recursion(h, t, localbound, 1);)
                                hashtable_update((void *)p, (void *)(++localbound << reusable),
                                                 t, NULL);
                            haploblocks[j][i - 2] = localbound;
#ifdef ENABLE_VERBOSE
                            if (v)
                                printf("%d recombinations required between site %d and site %d\n",
                                       localbound, j, j + i - 1);
#endif
                        }
                        else
                            haploblocks[j][i - 2] = -(int)lookup - 1;
                        if (!r2)
                            free_packedgenes(p);
                    }
                    free_genes(h);
                }

            /* So far we have worked on imploded sequence set - reconstruct
             * local bounds for original sequence set.
             */
            explode_local(haploblocks, l, n);
            /* Clean up */
            while (Length(l) > 0)
                free(Pop(l));
            DestroyLList(l);
            /* Restore eventlist */
            eventlist = tmp;
            reusable = oldreusable;
        }
#endif
        /* Clean up */
        if (!reusable)
            hashtable_destroy(t, (void (*)(void *))free_packedgenes,
                              NULL, (void (*)(void *))free);
    }
#ifdef HAPLOTYPE_BLOCKS
    else if (haploblocks != NULL) {
        while (Length(l) > 0)
            DestroyLList(Pop(l));
        DestroyLList(l);
    }
    representativeness = tmprep;
    representativeness_counter = tmprepcount;
#endif

    if (eventlist != NULL) {
        if (bound <= upper)
            Prepend(implode, eventlist);
        else
            while (Length(eventlist) > 0)
                free(Pop(eventlist));
    }

    free_genes(g);

    return bound;
}

/* Compute minimum number of recombinations required by any history for g */
int beagle(Genes *g, FILE *print_progress)
{
    return beagle_core(g, print_progress, 0, INT_MAX, NULL);
}

/* Compute minimum number of recombinations required by any history
 * for g, assuming lower is a lower bound on this number and
 * terminating the search once it has been established that upper does
 * not suffice.
 */
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper)
{
    return beagle_core(g, print_progress, lower, upper, NULL);
}

/* Compute minimum number of recombinations required by any history
 * for g, using existing hash table t for checking and storing
 * ancestral states encountered.
 */
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t)
{
    int n;

    reusable = 1;
    n = beagle_core(g, print_progress, 0, INT_MAX, t);
    reusable = 0;

    return n;
}

/* Compute minimum number of recombinations required by any history
 * for g, assuming lower is a lower bound on this number and
 * terminating the search once it has been established that upper does
 * not suffice; t is used as hash table for checking and storing
 * ancestral states encountered.
 */
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
                            int upper, HashTable *t)
{
    int n;

    reusable = 1;
    n = beagle_core(g, print_progress, lower, upper, t);
    reusable = 0;

    return n;
}

/* Find a random evolutionary history with at most r recombinations.
 * If r is less than the minimum number of recombinations required
 * for g, NULL is returned. If t is not NULL it is assumed to be a
 * reused and reusable hash table of ancestral configurations for
 * computing minimum number of recombinations required for this set of
 * genes. The evolutionary history inferred is returned as an llist of
 * Events.
 */

LList *beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t)
{
    int reuse = 1, i, j, *permutation;
    Genes *h;
    LList *tmp = eventlist, *history;
    EList *configurations, *events;
    Event *e;
    eventlist = NULL;

    /* Allocate hash table for ancestral configurations encountered if
     * none such was provided.
     */
    if (t == NULL) {
        t = beagle_allocate_hashtable(g, -1);
        reuse = 0;
    }

    /* Determine number of recombinations required for data set */
    if (beagle_reusable_bounded(g, print_progress, 0, r, t) > r) {
        eventlist = tmp;
        return NULL;
    }

    /* Find a history with at most bound recombinations and store events
     * in history.
     */
    history = MakeLList();
    events = elist_make();
    h = copy_genes(g);
    while (h->n > 1) {
        /* Non-segregating sites confuse the backtrack - remove */
        /* Removal of nonsegregating sites will eliminate all sequences in
         * a deterministic manner if all sites are eliminated. So remember
         * the number of sites so that we can coalesce them in a random
         * manner if all sites are removed.
         */
        j = h->n;
        remove_nonsegregating(h);
        /* If this results in all sites being eliminated, we are done (but,
         * possibly, for a few coalescences).
         */
        if (h->length == 0) {
            for (i = j; i > 1; i--) {
                e = (Event *)xmalloc(sizeof(Event));
                e->type = COALESCENCE;
                e->event.c.s1 = unbiased_random(i);
                j = unbiased_random(i - 1);
                if (j >= e->event.c.s1)
                    e->event.c.s2 = j + 1;
                else {
                    e->event.c.s2 = e->event.c.s1;
                    e->event.c.s1 = j;
                }
                Enqueue(history, e);
            }
            break;
        }
        /* Determine configurations reachable by mutation */
        configurations = force_mutation(h, events);
        /* Determine configurations reachable by coalescence */
        elist_extend(configurations, force_coalesce(h, events));
        /* Determine configurations reachable by recombinations */
        for (i = 0; i < h->n; i++)
            elist_extend(configurations, force_split(h, i, events));
        /* Go through configurations in randomly permuted order and choose
         * first that does not exceed recombination allowance.
         */
        permutation = (int *)xmalloc(elist_length(configurations) * sizeof(int));
        for (i = 0; i < elist_length(configurations); i++)
            permutation[i] = i;
        for (i = 0; i < elist_length(configurations); i++) {
            j = unbiased_random(elist_length(configurations) - i);
            free_genes(h);
            h = (Genes *)elist_get(configurations, permutation[i + j]);
            e = elist_get(events, permutation[i + j]);
            if(beagle_reusable_bounded(h, NULL, 0, r - (e->type == RECOMBINATION), t)
                    <= r - (e->type == RECOMBINATION)) {
                /* Found next configuration */
                Enqueue(history, e);
                if (e->type == RECOMBINATION)
                    /* Used one more recombination from our total allowance */
                    r--;
                /* Free all configurations and events not yet inspected */
                permutation[i + j] = permutation[i];
                for (i++; i < elist_length(configurations); i++) {
                    free_genes(elist_get(configurations, permutation[i]));
                    free(elist_get(events, permutation[i]));
                }
                /* Clean up */
                elist_destroy(configurations);
                events->count = 0;
                free(permutation);
                break;
            }
            else {
                /* This direction would lead to excess recombinations */
                free(elist_get(events, permutation[i + j]));
                permutation[i + j] = permutation[i];
                if (i == configurations->count - 1) {
                    fprintf(stderr,
                            "Error in randomised search. Please email error report\n");
                    exit(1);
                }
            }
        }
    }

    /* Clean up */
    elist_destroy(events);
    free_genes(h);
    eventlist = tmp;

    return history;
}

/* Allocate a hash table for storing ancestral states. The table_size
 * should be the logarithm of the number of buckets in the hash table
 * - if no valid number is provided (i.e. table_size is non-positive)
 * and a set of sequences is provided, a reasonable table size is
 * estimated.
 */
HashTable *beagle_allocate_hashtable(Genes *g, int table_size)
{
    HashTable *t;

    /* Determine size of and allocate hash table */
    if ((g != NULL) && (table_size <= 0)) {
        if (no_recombinations_required(g))
            table_size = 3;
        else
            table_size = msb((g->n - 3) * g->length) + LOWERBOUND(g);
    }
    else if (table_size <= 0)
        table_size = 10;

    t = new_packedgeneshashtable(table_size);

    return t;
}

/* Deallocate a hash table used by beagle and all the elements stored in it */
void beagle_deallocate_hashtable(HashTable *t)
{
    hashtable_destroy(t, (void (*)(void *))free_packedgenes, NULL,
                      (void (*)(void *))free);
}

/* Functions for interfacing with lower bound computations */
static int _greedy_rmin, _greedy_hk;
static Genes *_greedy_currentstate;
/* Initialise parameters for hashing function call information into m bins */
static int *_greedy_initparam(unsigned long m)
{
    int i, *p = (int *)xmalloc(5 * sizeof(int));

    p[0] = m;
    for (i = 1; i < 5; i++)
        p[i] = xrandom() % m;

    return p;
}

/* Compute minimum number of recombinations needed for current state */
static void *_noexp_rmin()
{
    /* Set up hash table for common use for all beagle invocations */
    if (_greedy_beaglereusable == NULL) {
        _greedy_beaglereusable = beagle_allocate_hashtable(_greedy_currentstate, -1);
    }
    if (_greedy_rmin < 0) {
        /* We haven't computed r_min for this configuration yet */
        _greedy_rmin = beagle_reusable(_greedy_currentstate, NULL,
                                       _greedy_beaglereusable);
    }
    
}


/* Compute haplotype lower bound with the heuristic parameters
 * specified by p.
 */
static int _hb(Genes *g)
{
    int i;
    void **a = (void **)xmalloc(4 * sizeof(void *));

    /* Create array specifying this function call */
    for (i = 0; i < 3; i++) {
        a[i + 1] = (void *)INT_MAX;
    }
    a[0] = (void *)_hb;

    /* Check whether this function call has previously been invoked; if
     * not, compute and store value.
     */
    if (!hashtable_lookup(a, _greedy_functioncalls, (void **)&i)) {
        i = haplotype_bound_genes(g);
        hashtable_insert(a, (void *)i, _greedy_functioncalls);
    }

    /* Store lower bound in expression */

    return i;
}

/* Compute lower bound from local exact minimum number of
 * recombinations combined using the composite method.
 */
static int _eagl(Genes *g)
{
    void **a = (void **)xmalloc(2 * sizeof(void *));
    int b, **B;
    Sites *s;

    /* Create array for specifying this function call */
    a[0] = (void *)_eagl;
    a[1] = (void *)(int)(10);

    /* Check whether this function call has previously been invoked; if
     * not, compute and store value.
     */
    if (!hashtable_lookup(a, _greedy_functioncalls, (void **)&b)) {
        s = genes2sites(g);
        B = hudson_kaplan_local(s);
        free_sites(s);
        b = eagl(g, 10, B, NULL, NULL);
        hashtable_insert(a, (void *)b, _greedy_functioncalls);
    }

    return b;
}



/* Hash the function call information stored in elm using the
 * parameters in p.
 */
static unsigned long _greedy_hash(void **elm, int *p)
{
    int i;
    unsigned long v = 0;

    if (*elm == _hb)
        for (i = 1; i <= 3; i++)
            v += p[i + 1] * (int)elm[i];
    else
        v = p[1] + p[2] * (int)elm[1];

    return v % p[0];
}

/* Compare the two function calls a and b */
static int _greedy_compare(void **a, void **b)
{
    int i;

    if (*a != *b)
        return 0;

    if (*a == _hb) {
        for (i = 1; i <= 3; i++)
            if ((int)a[i] != (int)b[i])
                return 0;
    }
    else
        return (int)a[1] == (int)b[1];

    return 1;
}

/* Set current ancestral state to g, update am, seq and len to reflect
 * g, and remove information from previous ancestral state from cache.
 */
static double _am, _seq, _len;
static void _reset_builtins(Genes *g)
{
    _greedy_currentstate = g;

    _am = ancestral_material(g);
    _seq = g->n;
    _len = g->length;

    _greedy_rmin = _greedy_hk = -1;
    if (_greedy_functioncalls != NULL) {
        hashtable_cleanout(_greedy_functioncalls, free, NULL);
    }
    else
        _greedy_functioncalls = hashtable_new
                                (6, (unsigned long (*)(void *, void *))_greedy_hash,
                                 (int (*)(void *, void *))_greedy_compare,
                                 (void *(*)(unsigned long))_greedy_initparam);
}

/* Update global quantities with contribution from g and free any
 * events that may be stored in eventlist.
 */
static int _choice_fixed;
HistoryFragment *_greedy_choice;
static double _minam, _maxam, _minseq, _maxseq, _minlen, _maxlen;
static Action ac;
static void __update(Genes *g)
{
    int am = ancestral_material(g);

    if (am < _minam)
        _minam = am;
    else if (am > _maxam)
        _maxam = am;
    if (g->n < _minseq)
        _minseq = g->n;
    else if (g->n > _maxseq)
        _maxseq = g->n;
    if (g->length < _minlen)
        _minlen = g->length;
    else if (g->length > _maxlen)
        _maxlen = g->length;
}

static void _update(Genes *g)
{
    /* Check whether we have found a path to the MRCA */
    if (!_choice_fixed && no_recombinations_required(g)) {
        /* Found a path to the MRCA - choose it */
        _choice_fixed = 1;
        _greedy_choice = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));
        _greedy_choice->g = g;
        _greedy_choice->event = eventlist;
        _greedy_choice->recombinations = _recombinations;
        _greedy_choice->elements = elements;
        _greedy_choice->sites = sites;
        _greedy_choice->action = ac;
        
        return;
    }

    /* Update global quantities */
    if (!_choice_fixed)
        __update(g);

    /* Free memory used for g and event leading to it */
    if (eventlist != NULL) {
        while (Length(eventlist) != 0)
            free(Pop(eventlist));
        DestroyLList(eventlist);
    }
    free_genes(g);
}

/* Score computation for each state in the neighbourhood */
static double sc_min = DBL_MAX, sc_max = 0;
static double prev_lb = 0, current_lb = 0, _lb;
double scoring_function(Genes *g) {
    
    double sc;
    double lb;
    int sign;
    
    // If we have already reached the end, we still cycle through all the possible
    // choices of last step and select the cheapest.
    // We set the score to -(cost of step) if it resolves the last incompatibility,
    // otherwise set the score to -(very big number). The random_select function will
    // pick the move with the least negative score in this case, as needed.
    if(_choice_fixed) {
        sign = (Temp < 0) - (Temp > 0) - (Temp == 0);
        if(no_recombinations_required(g)) {
            sc = sign * _recombinations;
        }
        else {
            sc = sign * DBL_MAX;
        }
    }
    // If we have not reached the end, score the move as usual. 
    else {
        if(_maxam < 75) {
            _noexp_rmin();
            lb = _greedy_rmin;
        }
//         else if(_am < 150) {
//             lb = _eagl(g);
//         }
        else if(_maxam < 200){
            lb = _hb(g);
        }
        else {
            lb = hudson_kaplan_genes(g);
        }
        
        _lb = lb;
        
        sc = (_recombinations + lb) * _maxam + _am;
        if(sc < sc_min) {
            sc_min = sc;
        }
        if(sc > sc_max) {
            sc_max = sc;
        }
    }
    
    return sc;
}

/* Once scores have been computed, renormalise and apply annealing */
double score_renormalise(Genes *g, double sc) {
    
    int sign;
    
    if(_choice_fixed) {
        sign = (Temp < 0) - (Temp > 0) - (Temp == 0);
        if(no_recombinations_required(g)) {
            sc = sign * _recombinations;
        }
        else {
            sc = sign * DBL_MAX;
        }
    }
    else {
        if(sc_max != sc_min) {
            if(Temp != -1) {
                sc = exp(Temp * (1 - (sc - sc_min)/(sc_max - sc_min)));
            } 
        }
        else {
            sc = 1;
        }
    }
    
    return sc;
}

/* Store HistoryFragments of possible predecessors in predecessors */
static EList *_predecessors = NULL;
static void _store(Genes *g)
{
    HistoryFragment *f;
    
    /* Wrap configuration and events leading to it in a HistoryFragment */
    f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));
    f->event = eventlist;
    f->g = g;
    f->recombinations = _recombinations;
    f->elements = elements;
    f->sites = sites;
    f->action = ac;
    if(elements != NULL && g->n != 0 && g->n != elist_length(elements)) {
        fprintf(stderr, "Error: number of sequence labels in elements [%d] not equal to current size of dataset [%d]. Event type: %.1f", elist_length(elements), g->n, _recombinations);
        exit(0);
    }
    if(elements != NULL && g->length > 0 && g->length != elist_length(sites)) {
        fprintf(stderr, "Error: number of site labels in sites not equal to current size of dataset.");
        exit(0);
    }
    if (!_choice_fixed && no_recombinations_required(g)) {
        /* Found a path to the MRCA - choose it */
        _choice_fixed = 1;
//         _greedy_choice = f;
    }

    elist_append(_predecessors, f);
    __update(g);

}

static int (*_choice_function)(double);


/* Update the lookup list of SE/RM and recombination numbers
 * This is of length rec_max, and keeps track of the maximum number of RM events seen
 * for each given number of recombinations already proposed. For example, if we have
 * seen a solution with 5 recombinations and 10 RMs, and the current solution reaches
 * 5 recombinations and 10 RMs but has not yet resolved all incompatibilities, then
 * this solution will be sub-optimal and can be abandoned.
 */
void update_lookup(EList *lku, int index, int bd) {
    int i, j, k;
    // Let S = number of SE + RM in the solution
    // Let R = number of recombinations in the solution
    // Then lookup[R] = S, lookup[R + 1 : R + 2*S] <= S, lookup[R + 2*S : end] = 0
    k = (lku->count - 1 > index + 2*bd ? index + 2*bd : lku->count - 1);
    elist_change(lku, index, (void *)bd);
    for(i = index + 1; i <= k; i++) {
        j = (int)elist_get(lku, i);
        if(j > bd) {
            elist_change(lku, i, (void *)bd);
        }
    }
    for(i = k+1; i < lku->count; i++) {
        elist_change(lku, i, (void *)0);
    }
}

/* Main function of kwarg implementing neighbourhood search.
 */
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double), void (*reset)(void), int ontheflyselection)
{
    int global, i, nbdsize = 0, total_nbdsize = 0, seflips = 0, rmflips = 0, recombs = 0, preds, bad_soln = 0;
    double r = 0;
    Index *start, *end;
    LList *tmp = eventlist;
    double printscore = 0;
    HistoryFragment *f;
    void (*action)(Genes *);
    const char *names[5];
    names[0] = "Coalescence";
    names[1] = "Sequencing error";
    names[2] = "Recurrent mutation";
    names[3] = "Single recombination";
    names[4] = "Double recombination";
    double *score_array;
    
    #ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
    #endif
    
    
    /* Create working copy of g */
    g = copy_genes(g);
    
    if(howverbose > 0) {
        fprintf(print_progress, "Input data:\n");
        if(howverbose == 2) {
            output_genes(g, print_progress, NULL);
        }
        fprintf(print_progress, "%d sequences with %d sites\n", g->n, g->length);
    }

    // Reduce the dataset
    implode_genes(g);
    if(howverbose > 0) {
        printf("%d sequences with %d sites after reducing\n", g->n, g->length);
    }
    if(lookup != NULL) {
        if((int)elist_get(lookup, 0) == INT_MAX)
            update_lookup(lookup, 0, g->n * g->length);
    }

    global = 1;
    
    /* Repeatedly choose an event back in time, until data set has been
     * explained.
     */
    _choice_function = select;
    if (!ontheflyselection && global)
        _predecessors = elist_make();
    if ((_choice_fixed = no_recombinations_required(g)) != 0)
        /* Data set can be explained without recombinations */
        free_genes(g);
    
     while (!_choice_fixed) {
        /* Reset statistics of reachable configurations */
        _minam = _minseq = _minlen = INT_MAX;
        _maxam = _maxseq = _maxlen = 0;
        _greedy_choice = NULL;
        nbdsize = 0;
        preds = 0;

        /* Determine interesting recombination ranges */
        start = maximumsubsumedprefixs(g);
        end = maximumsubsumedpostfixs(g);
     

        action = _store;
        
        /* We have just imploded genes, but we still need to pursue paths
            * coalescing compatible sequences where neither is subsumed in the
            * other but where the ancestral material is still entangled.
            */
        if(howverbose > 0) {
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
            fprintf(print_progress, "Searching possible predecessors:\n");
        }
        no_events = 0;
        _recombinations = 0;
        ac = COAL;
        preds = 0;
        nbdsize = 0;
        
        _coalesce_compatibleandentangled_map(g, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Coalescing entangled: ", preds);
            }
        
        if(se_cost != -1) {
            no_events = 1;
            _recombinations = se_cost;
            ac = SE;
            
            seqerror_flips(g, action);
                preds = elist_length(_predecessors) - nbdsize;
                nbdsize = elist_length(_predecessors);
                if(howverbose > 0) {
                    fprintf(print_progress, "%-40s %3d\n", "Sequencing errors: ", preds);
                }
        }
        
        if(rm_cost != -1) {
            no_events = 1;
            _recombinations = rm_cost;
            ac = RM;
            
            recmut_flips(g, action);

            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Recurrent mutations: ", preds);
            }
        }
        
        /* Try all sensible events with one split */
        if(r_cost != -1) {
            no_events = 1;
            _recombinations = r_cost;
            ac = RECOMB1;
            
            maximal_prefix_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Prefix recombinations: ", preds);
            }
            
            maximal_postfix_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Postfix recombinations: ", preds);
            }
        }
            
        /* Try all sensible events with two splits */
        if(rr_cost != -1) {
            no_events = 2;
            _recombinations = rr_cost;
            ac = RECOMB2;
            
            maximal_infix_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (infix): ", preds);
            }
            
            maximal_overlap_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (overlap): ", preds);
            }
        }
        
        if(howverbose > 0) {
            fprintf(print_progress, "%-40s %3d\n", "Finished constructing predecessors.", elist_length(_predecessors));
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
        }
        
        
        /* Finalise choice and prepare for next iteration */
        free_genes(g);
            /* Still looking for path to MRCA */
            if (!ontheflyselection && global) {
                /* So far we have only enumerated putative predecessors -
                 * score these and choose one.
                 */
                
                // Set the tracking lists to NULL for the score computation, and destroy the old elements/sites
                eventlist = NULL;
                elist_destroy(elements);
                elements = NULL;
                elist_destroy(sites);
                sites = NULL;
                reset();
                
                nbdsize = elist_length(_predecessors); // number of predecessors we score
                if(nbdsize == 0) {
                    fprintf(stderr, "No neighbours left to search but MRCA not reached.");
                }
                total_nbdsize = total_nbdsize + nbdsize;
                
                // Calculate all the scores and store in an array
                // Update sc_min and sc_max for renormalising the score later
                score_array = malloc(elist_length(_predecessors) * sizeof(double));
                if(!_choice_fixed) {
                    sc_min = DBL_MAX, sc_max = 0;
                    for (i = 0; i < elist_length(_predecessors); i++) {
                        f = (HistoryFragment *)elist_get(_predecessors, i);
                        _reset_builtins(f->g); // set f to be _greedy_currentstate
                        _recombinations = f->recombinations;
                        // Calculate all the scores and update the min and max
                        score_array[i] = scoring_function(f->g);
                    }
                }
                
                // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
                for (i = 0; i < elist_length(_predecessors); i++) {
                    f = (HistoryFragment *)elist_get(_predecessors, i);
                    _reset_builtins(f->g); // set _greedy_currentstate to be f->g
                    // Bug fix: need to update _recombinations otherwise this will always be 2
                    _recombinations = f->recombinations; 
                    printscore = score_renormalise(f->g, score_array[i]);
                    if (print_progress != NULL && howverbose == 2) {
                        fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i+1, f->recombinations);
                        output_genes(f->g, print_progress, NULL);
                        print_elist(f->elements, "Sequences: ");
                        print_elist(f->sites, "Sites: ");
                        fprintf(print_progress, "Predecessor score: %.0f \n\n", 
                                (printscore == -DBL_MAX ? -INFINITY : (printscore == DBL_MAX ? INFINITY : printscore)));
                        fflush(print_progress);
                    }
                    if (select(printscore)) {
                        // compute score and check if better than that of _greedy_choice
                        /* If so, discard old choice */
                        if (_greedy_choice != NULL) {
                            free_genes(_greedy_choice->g);
                            if (_greedy_choice->event != NULL) {
                                while (Length(_greedy_choice->event) != 0)
                                    free(Pop(_greedy_choice->event));
                                DestroyLList(_greedy_choice->event);
                            }
                            if(_greedy_choice->elements != NULL)
                                elist_destroy(_greedy_choice->elements);
                            if(_greedy_choice->sites != NULL)
                                elist_destroy(_greedy_choice->sites);
                            free(_greedy_choice);
                        }
                        /* Set f to be new choice */
                        _greedy_choice = f;
                    }
                    else {
                            /* Discard f */
                            free_genes(f->g);
                            if (f->event != NULL) {
                                while (Length(f->event) != 0)
                                    free(Pop(f->event));
                                DestroyLList(f->event);
                            }
                            if(f->elements != NULL) {
                                elist_destroy(f->elements);
                            }
                            if(f->sites != NULL) {
                                elist_destroy(f->sites);
                            }
                            free(f);
                    }
                }
                
                free(score_array);
                
                eventlist = tmp;
                elist_empty(_predecessors, NULL); // this should now be empty
            }
            
            g = _greedy_choice->g;
            elements = _greedy_choice->elements;
            sites = _greedy_choice->sites;
            
            switch(_greedy_choice->action) {
                case COAL:
                    break;
                case SE:
                    seflips = seflips + _greedy_choice->recombinations/se_cost;
                    break;
                case RM:
                    rmflips = rmflips + _greedy_choice->recombinations/rm_cost;
                    break;
                case RECOMB1:
                    recombs++;
                    break;
                case RECOMB2:
                    recombs += 2;
                    break;
            }
            
            if (print_progress != NULL && howverbose == 2) {
                fprintf(print_progress, "%s completed at cost of %.3f.\n", names[_greedy_choice->action], _greedy_choice->recombinations);
                fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
                fprintf(print_progress, "Current data:\n");
                output_genes(_greedy_choice->g, print_progress, NULL);
                fflush(print_progress);
            }
            if (print_progress != NULL && howverbose == 1) {
                fprintf(print_progress, "%s at cost %.3f \n", names[_greedy_choice->action], _greedy_choice->recombinations);
                fflush(print_progress);
            }
            /* Predecessor and events leading to it are stored in _greedy_choice */
//         }
        
        
        if (eventlist != NULL) {
            Append(eventlist, _greedy_choice->event);
        }
        
        r += _greedy_choice->recombinations;
        
        /* Clean up */
        free(_greedy_choice);
        free(start);
        free(end);
        if(_choice_fixed) {
            free_genes(g);
        }
        
        // Can abandon the run if the number of recombinations already exceeds rec_max
        if(recombs > rec_max) {
            bad_soln = 1;
            break;
        }
        
        // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
        // seen in earlier solutions.
        if(rec_max != INT_MAX && lookup != NULL) {
            if(seflips + rmflips > (int)elist_get(lookup, recombs)) {
                bad_soln = 1;
                break;
            }
        }
        
    }
    
    // If we exited the loop because of a sub-optimal solution, record this
    if(bad_soln) {
        if(reference > 0){
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", reference, r_seed, Temp, se_cost, rm_cost, r_cost, rr_cost, total_nbdsize);
        } else {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", r_seed, Temp, se_cost, rm_cost, r_cost, rr_cost, total_nbdsize);
        }
    }
    else {
    // Otherwise, record the result
        if (print_progress != NULL && howverbose > 0) {
            fprintf(print_progress, "\nTotal number of states considered: %d\n", total_nbdsize);
            fprintf(print_progress, "Total event cost: %.1f\n", r);
            if(reference > 0) {
                fprintf(print_progress, "%10s %13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Ref", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
                    "SE", "RM", "R", "N_states", "Time");
            } else {
                fprintf(print_progress, "%13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost", "SE", "RM", "R", "N_states", "Time");
            }
        }
        if(reference > 0) {
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", reference, r_seed, Temp, se_cost, rm_cost, r_cost, rr_cost, seflips, rmflips, recombs, total_nbdsize);
        } else {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", r_seed, Temp, se_cost, rm_cost, r_cost, rr_cost, seflips, rmflips, recombs, total_nbdsize);
        }
        if(lookup != NULL) {
            if(seflips + rmflips < (int)elist_get(lookup, recombs)){
                // If found a better bound r < rec_max for Rmin, update.
                if(seflips + rmflips == 0) {
                    rec_max = recombs;
                }
                update_lookup(lookup, recombs, seflips + rmflips);
            }
        }
    }
    
    elist_destroy(_predecessors);
    

    
    return r;

}

