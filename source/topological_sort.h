/*******************************************************************

    topological_sort.h

    Description of methods for sorting an equation system according
    to its strongly connected components

********************************************************************/

#ifndef _TOPOLOGICAL_SORT_H
#define _TOPOLOGICAL_SORT_H
int strongly_connected_components(int n, int *Ap, int *Ai, int **newindex,
				  int **oldindex, int **accsize);
int topological_sort(int n, int *Ap, int *Ai, double *Ax, double *bx,
		     int **newindex, int **oldindex, int **accsize);
#endif
