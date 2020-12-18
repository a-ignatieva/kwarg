/***************************************************************************
*
*    common.c: Implementation of functions that should be commonly available
*    to all source files.
*
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <unistd.h>

#include "common.h"
#include "gene.h"
#include "llist.h"
#include "arg.h"

#ifdef ENABLE_VERBOSE
static int verbosity = 0;
int verbose()
{
  return (verbosity > 0 ? verbosity : 0);
}

void set_verbose(int v)
{
  verbosity = v;
}
#endif
#ifdef HAPLOTYPE_BLOCKS
LList *representativeness = NULL;
LListCounter *representativeness_counter;
int **haploblocks = NULL;
#endif
double se_cost;
double rm_cost;
double r_cost;
double rr_cost;
LList *eventlist;
EList *elements;
EList *sites;
EList *lookup;
int seq_numbering;
int howverbose = 0;
int no_events = 0;
double _recombinations;
int gc_enabled = 0;
double Temp = 1;
double r_seed;
int rec_max;
int counter = 0;
int reference = 0;
HashTable *_greedy_functioncalls = NULL, *_greedy_beaglereusable = NULL;
#ifdef DEBUG
/* Define structure for storing trace of ancestral states as we return
 * from a successful path.
 */
HashTable *ancestral_state_trace = NULL;
#endif

/* xmalloc(n): Allocate n bytes of memory, checking for successful allocation.
 */
void *xmalloc(int n)
{
  void *adr;

  if (n <= 0){
    fprintf(stderr,
	    "Erroneous memory allocation - please email error report\n");
    return NULL;
  }

  if ((adr = malloc(n)) == NULL){
    fprintf(stderr, "Unable to allocate sufficient amount of memory\n");
    exit(2);
  }

  return adr;
}

/* xcalloc(m, n): Allocate n words of size m of memory intialised to
 * zero, checking for successful allocation.
 */
void *xcalloc(int m, int n)
{
  void *adr;

  if ((n <= 0) || (m <= 0)){
    fprintf(stderr,
	    "Erroneous memory allocation - please email error report\n");
    return NULL;
  }

  if ((adr = calloc(m, n)) == NULL){
    fprintf(stderr, "Unable to allocate sufficient amount of memory\n");
    exit(2);
  }

  return adr;
}

/* xrealloc(n): Resize memory allocated at oldadr to n bytes of memory,
 * checking for successful allocation.
 */
void *xrealloc(void *oldadr, int n)
{
  void *adr;

  if (n <= 0){
    fprintf(stderr,
	    "Erroneous memory allocation - please email error report\n");
    return NULL;
  }

  if ((adr = realloc(oldadr, n)) == NULL){
    fprintf(stderr, "Unable to allocate sufficient amount of memory\n");
    exit(2);
  }

  return adr;
}

/* Interface to random number generator. This interface is used
 * throughout, so only this needs to be replaced to change random
 * number generator.
 */
/* initialise_xrandom(): Initialise random number generator, using the
 * current time.
 */
long int x2seed;
void initialise_x2random(double seed)
{
#ifndef DEBUG
    if(seed == 0) {
        r_seed = (double)time(NULL) + (double)xrandom();
    }
    else {
        r_seed = seed;
    }

  srandom(r_seed);
#else
  /* Make sure random sequence is the same for every run */
  srandom(1);
#endif
    
}

long int xseed;
void initialise_xrandom()
{
    xseed = (double)time(NULL) + counter;
}

/* xrandom(): Return (pseudo-)random number between 0 and XRAND_MAX 
 */
long int x2random()
{
    long int r = random();
    counter++;
//     printf("%li ", r);
//     fflush(stdout);
    return r;
//     x2seed = (x2seed * 1664519 + 1013904229) % XRAND_MAX;
//     return x2seed;
    
}

/* Simple LCG, only used for initialising hash tables.
 */
long int xrandom()
{
    xseed = (xseed * 1664525 + 1013904223) % XRAND_MAX;
    return xseed;
}

/* Convert n to string */
char *i2a(int n)
{
  char *s, *t;

  if (n > 0)
    s = t = xmalloc(((int)log10(n) + 2) * sizeof(char));
  else if (n == 0)
    s = t = xmalloc(2 * sizeof(char));
  else{
    s = xmalloc(((int)log10(n) + 2) * sizeof(char));
    s[0] = '-';
    t = s + 1;
  }

  sprintf(t, "%d", n);

  return s;
}

/* Output s to fp with line length l and indentation i */
void pretty_print(FILE *fp, char *s, int l, int i)
{
  int j;
  char *last = s + strlen(s) - 1;

  if (fp == NULL)
    fp = stdout;

  /* Remove initial stretches of white space */
  while (isspace(*s))
    s++;
  if (s > last)
    /* s contains only white space */
    return;
  while (isspace(*last))
    last--;

  while (s <= last){
    /* Output indentation */
    for (j = 0; j < i; j++)
      fputc(' ', fp);

    /* Is there a new line within reach? */
    for (j = 0; (j <= l - i - 2) && (s + j <= last); j++)
      if (s[j] == '\n'){
	fwrite(s, sizeof(char), j + 1, fp);
	s += j + 1;
	break;
      }
    if ((j <= l - i - 2) && (s + j <= last))
      /* New line encountered - continue with next line */
      continue;

    if (s + l - i <= last){
      /* Find good place for next line break */
      for (j = l - i - 2; j > 0; j--){
	if ((s[j] == '-') && (s[j - 1] != ' ')){
	  /* Break after hyphen */
	  j++;
	  break;
	}
	if ((s[j] == ' ') && (s[j + 1] != '-')){
	  /* Break at space, unless it borders a dash */
	  while ((j > 0) && (s[j] == ' '))
	    j--;
	  if ((j > 0) && (s[j - 1] != '-'))
	    break;
	}
      }
      if (j == 0){
	/* No good line break found - look for acceptable line break */
	for (j = l - i - 2; (j > 0) && (s[j] != ' '); j--);
	while (s[j] == ' ')
	  j--;
	if (j == 0)
	  j = l - i - 2;
      }
    }
    else
      /* Rest of text fits on one line */
      j = last - s;
    fwrite(s, sizeof(char), j + 1, fp);
    fputc('\n', fp);
    s += j + 1;
    /* Remove initial stretch of white space */
    while (isspace(*s))
      s++;
  }
}

/* Print an option description to fp with line length l and subsequent
 * line indentation i (length of option if i negative).
 */
void print_option(FILE *fp, char *option, char *description, int l, int i)
{
  int n, m, j;

  if (fp == NULL)
    fp = stdout;

  /* Remove initial stretches of white space */
  while (isspace(*option))
    option++;
  n = strlen(option);
  if (n > 0)
    while (isspace(option[n - 1]))
      n--;
  while (isspace(*description))
    description++;
  m = strlen(description);
  if (m > 0)
    while (isspace(description[m - 1]))
      m--;

  /* Output option */
  fputc(' ', fp);
  fwrite(option, sizeof(char), n, fp);
  fputc(' ', fp);

  if (n + 2 + m <= l){
    /* Whole description fits on first line - output it */
    fwrite(description, sizeof(char), m, fp);
    fputc('\n', fp);
  }
  else{
    /* Check whether there is a new line within reach */
    for (j = 0; j <= l - n - 2; j++)
      if (description[j] == '\n'){
	fwrite(description, sizeof(char), j + 1, fp);
	pretty_print(fp, description + j + 1, l, (i < 0 ? n + 2 : i));
	return;
      }

    /* Output first line of description */
    /* Find good place for next line break */
    for (j = l - n - 3; j > 0; j--){
      if ((description[j] == '-') && (description[j - 1] != ' ')){
	/* Break after hyphen */
	j++;
	break;
      }
      if ((description[j] == ' ') && (description[j + 1] != '-')){
	/* Break at space, unless it borders a dash */
	while ((j > 0) && (description[j] == ' '))
	  j--;
	if ((j > 0) && (description[j - 1] != '-'))
	  break;
      }
    }
    if (j == 0){
      /* No good line break found - look for acceptable line break */
      for (j = l - n - 4; (j > 0) && (description[j] != ' '); j--);
      while (description[j] == ' ')
	j--;
      if (j == 0)
	j = l - n - 4;
    }
    fwrite(description, sizeof(char), j + 1, fp);
    fputc('\n', fp);
    description += j + 1;
    /* Output remainder of description */
    pretty_print(fp, description, l, (i < 0 ? n + 2 : i));
  }
}

#ifdef HAPLOTYPE_BLOCKS
/* From a matrix local of local bounds on recombinations in an
 * imploded data set and a list r describing what sites in the
 * original data set a site in the imploded data set corresponds to,
 * construct, in local, the matrix of local bounds on recombinations
 * in the original data set. It is assumed that local has been
 * allocated large enough to hold the expanded matrix.
 */
void explode_local(int **local, LList *r, int n)
{
  int i, j, a, b;
  SuperColumn *f, *l, *f2;
  LListCounter *first, *last;

  first = MakeCounter(r, LAST);
  last = MakeCounter(r, LAST);
  f = (SuperColumn *)Prev(first);
  /* Run through all pairs of sites in imploded sequence set */
  while (f != NULL){
    f2 = (SuperColumn *)Prev(first);
    a = (f2 != NULL ? f2->right : -1);
    b = n;
    l = (SuperColumn *)Prev(last);
    while (l != f){
      /* Run through all pairs of sites in original sequence set the
       * pair of sites in the imploded sequence set corresponds to.
       */
      for (i = f->right; i > a; i--)
	for (j = b - 1; j >= l->left; j--)
	  local[i][j - i - 1] = local[GetPosition(first) + 1]
	    [GetPosition(last) - GetPosition(first) - 2];
      b = l->left;
      l = Prev(last);
    }
    /* Make sure to reset entries corresponding to bounds between
     * sites in the same collapsed position to 0.
     */
    for (i = b - 2; i > a; i--)
      for (j = b - 1; j > i; j--)
	local[i][j - i - 1] = 0;
    f = f2;
    InitCounter(last, r, LAST);
  }
  DestroyCounter(first);
  DestroyCounter(last);
}
#endif

void remove_element(int *array, int index, int array_length)
{
    int i;
    for(i = index; i < array_length - 1; i++) 
        array[i] = array[i + 1];
    array[array_length - 1] = 0;
}

void delete_i(int *array, int i, int array_length)
{
    int j;
    for(j = 0; j < array_length; j++) 
        if(array[j] == i) {
            remove_element(array, j, array_length);
            break;
        }
}

void delete_by_value(int *array, int v, int array_length)
{
    int j = 0;
    while(j < array_length){
        if(array[j] == v) {
            remove_element(array, j, array_length);
        }
        else {
            j++;
        }
    }
}

void print_elist(EList *e, char *comment) {
    int i;
    int p;
    if(comment != NULL) {
        printf("%s", comment);
    }
    for(i=0; i < e->count; i++) {
        p = (int)(elist_get(e, i));
        if(p > 0) {
            printf("%d ", p);
        }
        else {
            printf("X ");
        }
    }
    printf("\n");
}

void set_array(double *a1, double *a2, int a2_length, int b) {
    int i;
    
    for(i = 0; i < a2_length; i++) {
        a1[b + i] = a2[i];
    }
}

