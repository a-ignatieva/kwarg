/************************************************************************ 
 * 
 * topological_sort.c: Implementation of topological sort of an
 * equation system given by matrix Ap, Ai, Ax, where the matrix
 * representation is as described for umfpack. The matrix columns are
 * reordered according to the sort and an array with the sizes of the
 * strongly connected components is returned.
 * 
 ************************************************************************/

#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "topological_sort.h"
#include "common.h"
#include "mergesort.h"
#include "bitfunctions.h"

static int min(int a, int b)
{
  return (a < b ? a : b);
}

/* Visit nodes in depth first order */
static int dfs(int v, int *n, int *Ap, int *Ai, int *dfsnum, int **stack,
	       int *newindex, int *oldindex, int *next, int **size)
{
  int i, low, *tmp;

  if (dfsnum[v] != 0)
    return dfsnum[v];

  *((*stack)++) = v;
  low = dfsnum[v] = ++(*n);
  /* Run through edges leaving v */
  for (i = Ap[v]; i < Ap[v + 1]; i++)
    low = min(low, dfs(Ai[i], n, Ap, Ai, dfsnum, stack, newindex, oldindex,
		       next, size));

  /* Check for component */
  if (low == dfsnum[v]){
    /* Remember old top of stack */
    tmp = *stack;
    /* Pop elements from stack until we reach v */
    while ((i = *(--(*stack))) != v){
      /* Insert the popped node in the reindexing arrays */
      oldindex[*next] = i;
      newindex[i] = (*next)++;
      dfsnum[i] = INT_MAX;
    }
    /* Insert last popped node (v) in the reindexing arrays */
    oldindex[*next] = i;
    newindex[i] = (*next)++;
    /* Number of elements popped is old minus new top of stack */
    *((*size)++) = tmp - *stack;
    dfsnum[v] = INT_MAX;
  }

  return low;
}

/* Determine the strongly connected components in the equation system
 * given by n, Ap, Ai, where the representation is as described for
 * umfpack. The variables component and size are set to point to
 * arrays containing the nodes grouped according to component (where
 * the components are inversely topologically sorted, i.e. there are
 * no outgoing edges from the first component and no incoming edges
 * into the last component) and the size of each component. The return
 * value is the number of connected components.
 */
int strongly_connected_components(int n, int *Ap, int *Ai, int **newindex,
				  int **oldindex, int **accsize)
{
  int *stack = (int *)xmalloc(n * sizeof(int)),
    *tmpstack = stack, *tmpsize;
  int *dfsnum = (int *)xcalloc(n, sizeof(int));
  int i, j = 0, index = 0;

  /* Allocate component and size arrays */
  *newindex = (int *)xmalloc(n * sizeof(int));
  *oldindex = (int *)xmalloc(n * sizeof(int));
  tmpsize = *accsize = (int *)xmalloc(n * sizeof(int));

  /* Find connected components, starting from each node in turn */
  for (i = 0; i < n; i++)
    dfs(i, &j, Ap, Ai, dfsnum, &tmpstack, *newindex, *oldindex, &index,
	&tmpsize);

  /* Number of components is the number of elements inserted into the
   * size array.
   */
  i = tmpsize - *accsize;
  /* Compact size array and compute accumulated sizes */
  *accsize = (int *)xrealloc(*accsize, (i + 1) * sizeof(int));
  for (j = 1; j < i; j++)
    (*accsize)[j] += (*accsize)[j - 1];
  memmove(*accsize + 1, *accsize, i * sizeof(int));
  **accsize = 0;

  /* Clean up */
  free(stack);
  free(dfsnum);

  return i;
}

/* We need to be able to access array holding new indeces for sorting
 * variables.
 */
int *tmpnew, *tmpAi;
/* Compare the new index of two old variables */
int compare_indeces(int *a, int *b)
{
  return (tmpnew[tmpAi[*a]] < tmpnew[tmpAi[*b]]);
}

/* Reorder the equation system given by n, Ap, Ai, Ax, bx, where the
 * representation is as described for umfpack, according to a
 * topological sort of the strongly connected components. The return
 * value is an array holding the size of each of the strongly
 * connected components.
 */
int topological_sort(int n, int *Ap, int *Ai, double *Ax, double *bx,
		      int **newindex, int **oldindex, int **accsize)
{
  int i, j, m,
    *Rp = (int *)xmalloc((n + 1) * sizeof(int)),
    *Ri = (int *)xmalloc(Ap[n] * sizeof(int)),
    *sort = (int *)xmalloc(n * sizeof(int)),
    *bucket = NULL;
  double *Rx = (double *)xmalloc(Ap[n] * sizeof(double)),
    *rx = (double *)xmalloc(Ap[n] * sizeof(double));

  /* Start by finding strongly connected components */
  m = strongly_connected_components(n, Ap, Ai, newindex, oldindex, accsize);

  /* Create permuted equation system in R */
  tmpnew = *newindex;
  Rp[0] = 0;
  for (i = 0; i < n; i++){
    /* Construct column for variable i in the new order */
    m = Ap[(*oldindex)[i] + 1] - Ap[(*oldindex)[i]];
    Rp[i + 1] = Rp[i] + m;
    /* Sort variables in column according to their new index */
    /* Check whether an O(m log(m)) mergesort or O(n) bucket sort is faster */
    if (m * msb(m) >= n){
      /* Do bucket sort */
      if (bucket == NULL)
	bucket = (int *)xcalloc(n, sizeof(int));
      else
	for (j = 0; j < n; j++)
	  bucket[j] = 0;
      for (j = 0; j < m; j++)
	bucket[(*newindex)[Ai[Ap[(*oldindex)[i]] + j]]] = 1;
      for (j = 1; j < n; j++)
	bucket[j] += bucket[j - 1];
      for (j = 0; j < m; j++)
	sort[bucket[(*newindex)[Ai[Ap[(*oldindex)[i]] + j]]] - 1] = j;
    }
    else{
      /* Do merge sort */
      for (j = 0; j < m; j++)
	sort[j] = j;
      tmpAi = Ai + Ap[(*oldindex)[i]];
      merge_sort(sort, m, sizeof(int),
		 (int (*)(void *, void *))compare_indeces);
    }
    /* Set the entries of Ri and Rx relating to this column */
    for (j = 0; j < m; j++){
      Ri[Rp[i] + j] = (*newindex)[Ai[Ap[(*oldindex)[i]] + sort[j]]];
      Rx[Rp[i] + j] = Ax[Ap[(*oldindex)[i]] + sort[j]];
      rx[i] = bx[(*oldindex)[i]];
    }
  }

  /* Move everything back to A */
  memcpy(Ap, Rp, (n + 1) * sizeof(int));
  memcpy(Ai, Ri, Ap[n] * sizeof(int));
  memcpy(Ax, Rx, Ap[n] * sizeof(double));
  memcpy(bx, rx, n * sizeof(double));

  /* Clean up */
  if (bucket != NULL)
    free(bucket);
  free(Rp);
  free(Ri);
  free(Rx);
  free(rx);
  free(sort);

  return m;
}
