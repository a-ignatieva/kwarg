/* common.h: header file for common.c; definition of common data types.
 */
#ifndef _COMMON_H
#define _COMMON_H

#include "llist.h"
#include "gene.h"
#include "arg.h"
/* Function prototypes; a brief explanation of each function should be
 * provided prior to its implementation in common.c.
 */
#ifdef ENABLE_VERBOSE
int verbose();
void set_verbose(int v);
#endif

#ifdef HAPLOTYPE_BLOCKS
typedef struct _SuperColumn {
  int left;
  int right;
} SuperColumn;

extern LList *representativeness;
extern LListCounter *representativeness_counter;
extern int **haploblocks;
void explode_local(int **local, LList *r, int n);
#endif
extern LList *eventlist;
extern EList *elements;
extern EList *sites;
extern EList *lookup;
extern int seq_numbering;
extern double se_cost;
extern double rm_cost;
extern double r_cost;
extern double rr_cost;
extern int howverbose;
extern double _recombinations;
extern int no_events;
extern int gc_enabled;
extern double Temp;
extern double r_seed;
extern int rec_max, rm_max;
extern long int x2seed;
extern long int xseed;
extern int counter;
extern int reference;
extern HashTable *_greedy_functioncalls, *_greedy_beaglereusable;
#ifdef DEBUG
extern HashTable *ancestral_state_trace;
#endif

void *xmalloc(int n);
void *xcalloc(int m, int n);
void *xrealloc(void *oldadr, int n);
#define XRAND_MAX RAND_MAX
void initialise_xrandom();
void initialise_x2random(double seed);
long int xrandom();
long int x2random();
char *i2a(int n);
void pretty_print(FILE *fp, char *s, int l, int i);
void print_option(FILE *fp, char *option, char *description, int l, int i);
void remove_element(int *array, int index, int array_length);
void delete_i(int *array, int i, int array_length);
void delete_by_value(int *array, int v, int array_length);
int max_value(int *array, int array_length);
void print_elist(EList *e, char *comment);
void set_array(double *a1, double *a2, int a2_length, int b);

#endif
