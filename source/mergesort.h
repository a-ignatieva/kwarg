/* mergesort.h: Prototype for mergesort function */

#ifndef _MERGESORT_H
#define _MERGESORT_H

#include <stdlib.h>

void merge_sort(void *data, size_t n, size_t size,
		int (*less_than)(void *, void *));

#endif
