#ifndef _HASHTABLE_H
#define _HASHTABLE_H

#include <stdio.h>
#include <stdarg.h>
#include "elist.h"

typedef struct _HashTable {
  unsigned long size;                    /* Size of hash table */
  unsigned long (*hash)(void *, void *); /* Hash function */
  int (*compare)(void *, void *);        /* Function to compare two elements */
  void *parameters;                      /* Parameters to hash function */
  EList *table;                          /* Actual table */
} HashTable;

/* Prototypes */
void hashtable_init(int bits, HashTable *table,
		    unsigned long (*hash)(void *, void *),
		    int (*compare)(void *, void *),
		    void *(*initialise_parameters)(unsigned long));
HashTable *hashtable_new(int bits, unsigned long (*hash)(void *, void *),
			 int (*compare)(void *, void *),
			 void *(*initialise_parameters)(unsigned long));
void hashtable_free(HashTable *t, void (*free_key)(void *),
		    void (*free_value)(void *),
		    void (*free_parameters)(void *));
void hashtable_destroy(HashTable *t, void (*free_key)(void *),
		    void (*free_value)(void *),
		    void (*free_parameters)(void *));
void hashtable_cleanout(HashTable *t, void (*free_key)(void *),
			void (*free_value)(void *));
int hashtable_lookup(void *elm, HashTable *t, void **value);
void *hashtable_lookuprepresentative(void *elm, HashTable *t, void **value);
void hashtable_insert(void *elm, void *value, HashTable *t);
int hashtable_update(void *elm, void *value, HashTable *t,
		     int (*update)(void *, void *));
int hashtable_size(HashTable *t);
void hashtable_map(HashTable *t, void (*f)(void *, void *, va_list), ...);
int hashtable_collisions(HashTable *t);
int hashtable_largestbucket(HashTable *t);
void hashtable_printlargestbucket(HashTable *t,
				  void (*print_elm)(void *, va_list));

#endif
