#ifndef _EXACT_H
#define _EXACT_H

#include <stdio.h>
#include "gene.h"
#include "common.h"
#include "hashtable.h"

extern int exact_randomise;
int beagle(Genes *g, FILE *print_progress, KwargContext *ctx);
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper, KwargContext *ctx);
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t, KwargContext *ctx);
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
			    int upper, HashTable *t, KwargContext *ctx);
LList *beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t, KwargContext *ctx);
HashTable *beagle_allocate_hashtable(Genes *g, int table_size, KwargContext *ctx);
void beagle_deallocate_hashtable(HashTable *t);
double scoring_function(Genes *g, KwargContext *ctx);
double score_renormalise(Genes *g, double sc, KwargContext *ctx);
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double, KwargContext*),
               void (*reset)(KwargContext*), int ontheflyselection, int reference, KwargContext *ctx);
#endif
