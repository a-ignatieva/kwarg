/*********************************************************************** 
 * 
 * mergesort.c: Implementation of merge sort using a temporary array
 * of the same size as the data to be sorted. The data is assumed to
 * be an array of void pointers.
 *
 ***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mergesort.h"

/* xmalloc(n): Tries to allocate n bytes of memory; if unsuccessful an
 * error is reported.
 */
static void *xmalloc(size_t n)
{
  void *adr = malloc(n);

  if (adr == NULL){
    fprintf(stderr, "Alzheimer's error\n");
    exit(1);
  }

  return adr;
}

/* merge(s1, e1, s2, e2, dest, size, less_than): Merge the two arrays
 * of elements, each element requiring size bytes of storage, starting
 * at locations s1 and s2 into an array starting at dest. The elements
 * are compared using less_than and e1 (e2) specifies the address
 * immediately after the end of array 1 (array 2).
 */
static void merge(void *s1, void *e1, void *s2, void *e2, void *dest,
		  size_t size, int (*less_than)(void *, void *))
{
  for (;;){
    /* By testing whether the next element in array 2 is smaller than
     * the next element in array 1, rather than the other way round,
     * we achieve in place sorting (assuming that array 1 is always to
     * the left of array 2 in the full data set).
     */
    if (less_than(s2, s1)){
      /* Move element from array 2 to destination */
      memcpy(dest, s2, size);
      dest += size;
      s2 += size;
      if (s2 == e2){
	/* Finished with array 2, copy rest of array 1 to dest */
	memcpy(dest, s1, e1 - s1);
	return;
      }
    }
    else{
      /* Move element from array 1 to destination */
      memcpy(dest, s1, size);
      dest += size;
      s1 += size;
      if (s1 == e1){
	/* Finished with array 1, copy rest of array 2 to dest */
	memcpy(dest, s2, e2 - s2);
	return;
      }
    }
  }
}

/* merge_sort(data, n, size, less_than): Sort the n elements, each of
 * requiring size bytes of storage, in the array starting at data
 * according to the less_than comparison function; the data array is
 * overwritten with the sorted array.
 */
void merge_sort(void *data, size_t n, size_t size,
		int (*less_than)(void *, void *))
{
  void *tmp = (void *)xmalloc(n * size), *tmp2;
  int i, istmp = 0, l = size;

  /* Intervals and array lengths are measured in multiples of size to
   * avoid constantly having to multiply by size.
   */
  n = n * size;
  while (l < n){
    /* Merge intervals of length l */
    for (i = 0; i + 2 * l < n; i += 2 * l){
      /* Merge intervals [i, i + l) and [i + l, i + 2l) */
      merge(data + i, data + i + l, data + i + l, data + i + 2 * l, tmp + i,
	    size, less_than);
    }
    if (i + l >= n)
      /* Only one interval left */
      memcpy(tmp + i, data + i, n - i);
    else
      /* Two intervals, of which the latter is no longer than l, left */
      merge(data + i, data + i + l, data + i + l, data + n, tmp + i, 
	    size, less_than);
    /* Double interval length and swap data and tmp */
    l *= 2;
    tmp2 = data;
    data = tmp;
    tmp = tmp2;
    /* Update whether current data is original or temporary array */
    istmp = !istmp;
  }

  if (istmp){
    /* Final merge was into temporary array - copy sorted array to data */
    /* The variable data actually points to the temporary array, while
     * tmp points to the original data array.
     */
    memcpy(tmp, data, n);
    /* Free temporary array */
    free(data);
  }
  else
    /* Free temporary array */
    free(tmp);
}
