/************************************************************************
 * 
 * enumerate.c: Implementation of functions to enumerate ancestral
 * states in histories with at most a given number of recombinations
 * under the infinite sites assumption for an SNP data set.
 *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "enumerate.h"
#include "elist.h"
#include "common.h"
#include "bounds.h"
#include "exact.h"
#include "bitfunctions.h"
#include "mergesort.h"
#include "hashtable.h"

/* Select method for computing a lower bound on the number of recombinations */
#ifdef GREVEN_HAPLOTYPE
#define LOWERBOUND(g) haplotype_bound_genes(g)
#elif defined(GREVEN_HAPLOTYPEHEURISTIC)
#define LOWERBOUND(g) haplotype_heuristic_genes(g, INT_MAX, INT_MAX, 1)
#else
#define LOWERBOUND(g) hudson_kaplan_genes(g)
#endif

/* Declaration of enumerate_recursion as it is invoked from
 * handle_ancestralstate.
 */
static void enumerate_recursion(Genes *g, int n, int lower, HashTable *t,
				HashTable *states);
/* Check whether g is a novel ancestral state that can be explained
 * with at most n recombinations. If this is the case, insert it in
 * states and pursue branches leading back in time from g.
 */
static void handle_ancestralstate(Genes *g, int n, int lower, HashTable *t,
				 HashTable *states)
{
  void *lookup;
#ifdef ENABLE_VERBOSE
  int v = verbose();
  set_verbose(v - 1);
#endif
  PackedGenes *p = pack_genes(g);

  if (lower < 0)
    lower = 0;

  if (!hashtable_lookup((void *)p, states, &lookup)){
    /* We reached a novel ancestral state - determine whether we
     * can extend the history without using more than n further
     * recombinations.
     */
    if (beagle_reusable_bounded(g, NULL, lower, n, t) > n){
      free_genes(g);
      free_packedgenes(p);
      return;
    }

    /* The ancestral state g is visitable */
    hashtable_insert((void *)p, (void *)n, states);
#ifdef ENABLE_VERBOSE
    if (v){
      printf("New visitable ancestral state:\n");
      output_genes_indexed(g, NULL);
    }
    set_verbose(v);
#endif
  }
  else if (n <= (int)lookup){
    /* Been here, done that */
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif
    free_genes(g);
    free_packedgenes(p);
    return;
  }
  else{
    /* Been here before, but using more recombinations */
    hashtable_update((void *)p, (void *)n, states, NULL);
    free_packedgenes(p);
  }

#ifdef ENABLE_VERBOSE
  set_verbose(v);
#endif
  enumerate_recursion(g, n, lower, t, states);

  free_genes(g);
}

/* The recursion actually enumerating ancestral states in histories
 * ending with g and going back in time using at most n
 * recombinations; lower is a guaranteed lower bound on the number of
 * recombinations required to explain g, t is an auxiliary hash table
 * used for determining the minimum number of recombinations required
 * for a given ancestral state, and states is the hash table the
 * enumerated ancestral states is stored in.
 */
static void enumerate_recursion(Genes *g, int n, int lower, HashTable *t,
				HashTable *states)
{
  int i, j;
  Genes *h;
  EList *forced;

  /* Try all possible coalesces */
  for (i = 0; i < g->n - 1; i++)
    for (j = i + 1; j < g->n; j++)
      if (compatible(g, i, j)){
	h = copy_genes(g);
	coalesce(h, i, j);
#ifdef ENABLE_VERBOSE
	if (verbose() > 1){
	  printf("Coalescing sequences %d and %d", i, j);
	  output_genes_indexed(h, NULL);
	}
#endif
	handle_ancestralstate(h, n, lower, t, states);
      }

  /* Try all possible mutations */
#ifdef ENABLE_VERBOSE
  set_verbose(verbose() - 1);
#endif
  forced = force_mutation(g, NULL);
#ifdef ENABLE_VERBOSE
  set_verbose(verbose() + 1);
#endif
  for (i = 0; i < elist_length(forced); i++){
    h = elist_get(forced, i);
    handle_ancestralstate(h, n, lower, t, states);
  }
  elist_destroy(forced);

  if (n > 0)
    /* Try all possible recombinations */
    for (i = 0; i < g->n; i++){
#ifdef ENABLE_VERBOSE
      set_verbose(verbose() - 1);
#endif
      forced = force_split(g, i, NULL);
#ifdef ENABLE_VERBOSE
      set_verbose(verbose() + 1);
#endif
      for (j = 0; j < elist_length(forced); j++){
	h = elist_get(forced, j);
	handle_ancestralstate(h, n - 1, lower - 1, t, states);
      }
      elist_destroy(forced);
    }
}

/* Enumerate the different ancestral states visitable in histories
 * with at most n recombinations.
 */
HashTable *enumerate_absolute(Genes *g, int n)
{
  HashTable *t, *states = NULL;
  Genes *h;
  int table_size, bound;
#ifdef ENABLE_VERBOSE
  int v = verbose();
  set_verbose(v - 1);
#endif

  /* Compute a good lower bound on the number of recombinations */
  h = copy_genes(g);
  implode_genes(h);
  if (no_recombinations_required(h))
    bound = 0;
  else{
    reallocate_genes(h);
    bound = LOWERBOUND(h);
  }
  if (bound > n){
    /* More recombinations required than we have available */
    free_genes(h);
    return NULL;
  }

  /* Determine sizes of and allocate hash tables */
  table_size = msb(h->n * h->length) + n;
  t = new_packedgeneshashtable(table_size);
  /* Compute minimum number of recombinations required */
  if (!no_recombinations_required(h))
    bound = beagle_reusable_bounded(h, NULL, bound, n, t);
  free_genes(h);

  if (bound <= n){
    table_size = msb(g->n * g->length) + n;
    states = new_packedgeneshashtable(table_size);

    /* Enumerate all ancestral states visitable with at most bound
     * recombinations.
     */
    hashtable_insert((void *)pack_genes(g), (void *)bound, states);
    enumerate_recursion(g, n, bound, t, states);
  }
  /* Clean up */
  hashtable_destroy(t, (void (*)(void *))free_packedgenes, NULL,
		    (void (*)(void *))free);

  return states;
}

/* Compute the number of different ancestral states visitable in
 * histories with at most n recombinations.
 */
int count_absolute(Genes *g, int n)
{
  HashTable *states = enumerate_absolute(g,n);
  int i = 0;

  if (states != NULL){
    i = hashtable_size(states);
    /* Clean up */
    hashtable_destroy(states, (void (*)(void *))free_packedgenes, NULL,
		      (void (*)(void *))free);
  }

  return i;
}

/* Enumerate the different ancestral states visitable in histories
 * with at most n recombinations more than the minimum number of
 * recombinations required to explain the sequences in g.
 */
HashTable *enumerate_relative(Genes *g, int n)
{
  HashTable *t, *states;
  Genes *h;
  int table_size, bound;
#ifdef ENABLE_VERBOSE
  int v = verbose();
  set_verbose(v - 1);
#endif

  if (n < 0)
    /* Cannot explain g with a negative number of recombinations */
    return NULL;

  /* Compute a good lower bound on the number of recombinations */
  h = copy_genes(g);
  implode_genes(h);
  if (no_recombinations_required(h))
    bound = 0;
  else{
    reallocate_genes(h);
    bound = LOWERBOUND(h);
  }

  /* Determine sizes of and allocate hash tables */
  table_size = msb(h->n * h->length) + n;
  t = new_packedgeneshashtable(table_size);
  /* Compute minimum number of recombinations required */
  if (!no_recombinations_required(h))
    bound = beagle_reusable_bounded(h, NULL, bound, INT_MAX, t);
  free_genes(h);


  table_size = msb(g->n * g->length) + n;
  states = new_packedgeneshashtable(table_size);

  /* Enumerate all ancestral states visitable with at most bound
   * recombinations.
   */
  hashtable_insert((void *)pack_genes(g), (void *)(bound + n), states);

  enumerate_recursion(g, n + bound, bound, t, states);
  /* Clean up */
  hashtable_destroy(t, (void (*)(void *))free_packedgenes, NULL,
		    (void (*)(void *))free);

  return states;
}

/* Compute the number of different ancestral states visitable in
 * histories with at most n recombinations more than the minimum
 * number of recombinations required to explain the sequences in g.
 */
int count_relative(Genes *g, int n)
{
  HashTable *states = enumerate_relative(g, n);
  int i = 0;

  if (states != NULL){
    i = hashtable_size(states);
    hashtable_destroy(states, (void (*)(void *))free_packedgenes, NULL,
		      (void (*)(void *))free);
  }

  return i;
}
