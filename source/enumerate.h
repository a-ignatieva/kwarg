#ifndef _ENUMERATE_H
#define _ENUMERATE_H

#include "gene.h"
#include "hashtable.h"

HashTable *enumerate_absolute(Genes *g, int n);
int count_absolute(Genes *g, int n);
HashTable *enumerate_relative(Genes *g, int n);
int count_relative(Genes *g, int n);

#endif
