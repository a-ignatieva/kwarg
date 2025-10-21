#ifndef _BOUNDS_H
#define _BOUNDS_H

#include "gene.h"

int yun(Genes *g, int best, KwargContext *ctx);
int **hudson_kaplan_local(Sites *s);
int hudson_kaplan(Sites *s);
int hudson_kaplan_genes(Genes *g, KwargContext *ctx);
int **haplotype_heuristic_local(Sites *s, int maxsetsize,
				int maxintervallength,
				int subsetincreasethreshold);
int haplotype_heuristic(Sites *s, int maxsetsize, int maxintervallength,
			int subsetincreasethreshold, KwargContext *ctx);
int haplotype_heuristic_genes(Genes *g, int maxsetsize, int maxintervallength,
			      int subsetincreasethreshold, KwargContext *ctx);
int **haplotype_bound_local(Sites *s);
int haplotype_bound(Sites *s, KwargContext *ctx);
int haplotype_bound_genes(Genes *g, KwargContext *ctx);
int **eagl_local(Genes *g, int max_length, int **B,
		 int (*per_region_action)(void),
		 int (*per_increment_action)(int, int, int **), KwargContext *ctx);
int eagl(Genes *g, int max_length, int **B,
	 int (*per_region_action)(void),
	 int (*per_increment_action)(int, int, int **), KwargContext *ctx);
int local2global(int n, int **B, KwargContext *ctx);
void local2locals(int n, int **B);

#endif
