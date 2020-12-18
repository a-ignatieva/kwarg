/************************************************************************
 *
 * hashtable.c: Generic functions for maintaining a set as a hash
 * table. The hashing is assumed to be modulo some prime number, so a
 * function for finding a suitable prime is supplied. But the hash
 * function itself necessarily has to be supplied by the invoking
 * function.
 *
 ************************************************************************/

#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "hashtable.h"
#include "elist.h"
#include "common.h"

/* Find the largest integer m with m * m <= n */
static unsigned long intsqrt(unsigned long n)
{
  unsigned long i = (unsigned long)sqrt(n);

  while (i * i < n)
    i++;

  while (i * i > n)
    i--;

  return i;
}

/* Find a prime number close to n. No intelligent primality testing
 * is used, but as we are only looking for primes with at most
 * sizeof(long) bits, this shouldn't be a concern.
 */
static unsigned long prime(unsigned long n)
{
  unsigned long i;

  /* We only search numbers that we know are not divisible with 2 or
   * 3; to do this we start from a number we know has remainder 1
   * modulo 6.
   */
  n = n - (n % 6) + 1;

  for (;;n += 6){
        /* Check whether n is prime */
        for (i = 5; i <= intsqrt(n); i += 6){
            if (n % i == 0)
                /* n is not a prime */
                break;
            if (n % (i + 4) == 0)
                /* n is not a prime */
                break;
        }
        if (i > intsqrt(n + 1))
        /* n is a prime */
        break;
        /* Check whether n + 4 is prime */
        for (i = 5; i <= intsqrt(n + 4); i += 6){
            if ((n + 4) % i == 0)
                /* n + 4 is not a prime */
                break;
            if ((n + 4) % (i + 4) == 0)
                /* n + 4 is not a prime */
                break;
        }
        if (i > intsqrt(n + 4)){
        /* n + 4 is a prime */
        n += 4;
        break;
        }
  }

  return n;
}

/* Return a prime number roughly between 2^bits and 2^(bits + 1) */
static unsigned long modulo(int bits)
{
  if (bits >= sizeof(unsigned long) * CHAR_BIT)
    /* Requested modulo too large for return type */
    bits = sizeof(unsigned long) * CHAR_BIT - 1;

  return prime(((unsigned long)1 << bits) + (xrandom() & (((unsigned long)1 << bits) - 1)));
}

/* Initialise a new hash table mapping to a universe of size between
 * 2^bits and 2^(bits + 1).
 */
void hashtable_init(int bits, HashTable *table,
		    unsigned long (*hash)(void *, void *),
		    int (*compare)(void *, void *),
		    void *(*initialise_parameters)(unsigned long))
{
  unsigned long i;

  /* Determine hash table size */
  table->size = modulo(bits);

  /* Initialise hash table */
  table->parameters = initialise_parameters(table->size);
  table->hash = hash;
  table->compare = compare;
  table->table = (EList *)xmalloc(table->size * sizeof(EList));
  for (i = 0; i < table->size; i++)
    elist_init(table->table + i);
}

/* Create a new hash table mapping to a universe of size between
 * 2^bits and 2^(bits + 1).
 */
HashTable *hashtable_new(int bits, unsigned long (*hash)(void *, void *),
			 int (*compare)(void *, void *),
			 void *(*initialise_parameters)(unsigned long))
{
  unsigned long i;
  HashTable *new = (HashTable *)xmalloc(sizeof(HashTable));

  /* Determine hash table size */
  new->size = modulo(bits);

  /* Initialise hash table */
  new->parameters = initialise_parameters(new->size);
  new->hash = hash;
  new->compare = compare;
  new->table = (EList *)xmalloc(new->size * sizeof(EList));
  for (i = 0; i < new->size; i++)
    elist_init(new->table + i);

  return new;
}

typedef struct _HashTablePair {
  void *key;
  void *value;
} HashTablePair;

typedef void (*_HASHTABLE_V_Func_VP)(void *);

static void free_hashtablepair(HashTablePair *entry, va_list args)
{
  void (*free_key)(void *) = va_arg(args, _HASHTABLE_V_Func_VP);
  void (*free_value)(void *) = va_arg(args, _HASHTABLE_V_Func_VP);
  

  if (free_key != NULL)
    free_key(entry->key);
  if (free_value != NULL)
    free_value(entry->value);
  free(entry);
}

/* Deallocate memory used for hash table t, except for the HashTable
 * structure t itself; if free_key or free_value are non-NULL they are
 * applied to the key, respectively value, of each key/value pair
 * stored in the hash table. If free_parameters is non-NULL it is
 * applied to the parameters stored for the hash function.
 */
void hashtable_free(HashTable *t, void (*free_key)(void *),
		    void (*free_value)(void *),
		    void (*free_parameters)(void *))
{
  int i;

  for (i = 0; i < t->size; i++)
    if (t->table[i].list != NULL){
      if ((free_key != NULL) || (free_value != NULL))
	elist_map(t->table + i,
		  (void (*)(void *, va_list))free_hashtablepair,
		  free_key, free_value);
      free(t->table[i].list);
    }
  free(t->table);
  if (free_parameters != NULL)
    free_parameters(t->parameters);
}

/* Deallocate memory used for hash table t; if free_key or free_value
 * are non-NULL they are applied to the key, respectively value, of
 * each key/value pair stored in the hash table. If free_parameters is
 * non-NULL it is applied to the parameters stored for the hash
 * function.
 */
void hashtable_destroy(HashTable *t, void (*free_key)(void *),
		       void (*free_value)(void *),
		       void (*free_parameters)(void *))
{
  hashtable_free(t, free_key, free_value, free_parameters);
  free(t);
}

/* Remove all the elements stored in t; if free_genes or free_keys are
 * non-NULL they are applied to the key, respectively value, of each
 * key/value pair stored in the hash table.
 */
void hashtable_cleanout(HashTable *t, void (*free_key)(void *),
			void (*free_value)(void *))
{
  int i;

  for (i = 0; i < t->size; i++)
    if (t->table[i].count != 0){
        if ((free_key != NULL) || (free_value != NULL))
            elist_map(t->table + i,
            (void (*)(void *, va_list))free_hashtablepair,
            free_key, free_value);
        t->table[i].count = 0;
    }
}

/* Determine elm's index in elist, if present. Otherwise return -1 */
static int find_index(void *elm, EList *elist, int (*compare)(void *, void *))
{
  int i;
  HashTablePair *entry;

  for (i = 0; i < elist_length(elist); i++){
    entry = (HashTablePair *)elist_get(elist, i);
    if (compare(elm, entry->key))
      return i;
  }

  return -1;
}

/* Check whether elm is already present in t. If value is not NULL, set
 * it equal to elm's value if elm is present in t.
 */
int hashtable_lookup(void *elm, HashTable *t, void **value)
{
  unsigned long h = t->hash(elm, t->parameters);
  int i = find_index(elm, t->table + h, t->compare);
  HashTablePair *entry;

  if (i < 0)
    /* g is not present in t */
    return 0;
  else{
    if (value != NULL){
      entry = elist_get(t->table + h, i);
      *value = entry->value;
    }
    return 1;
  }
}

/* Check whether elm is already present in t. If elm is present,
 * return the representative of elm stored in t; otherwise return
 * NULL. If value is not NULL, set it equal to elm's value if elm is
 * present in t.
 */
void *hashtable_lookuprepresentative(void *elm, HashTable *t, void **value)
{
  unsigned long h = t->hash(elm, t->parameters);
  int i = find_index(elm, t->table + h, t->compare);
  HashTablePair *entry;

  if (i < 0)
    /* g is not present in t */
    return NULL;
  else{
    if (value != NULL){
      entry = elist_get(t->table + h, i);
      *value = entry->value;
    }
    return entry->key;
  }
}

/* Insert elm into t with value associated. It is assumed that elm is not
 * already present in t.
 */
void hashtable_insert(void *elm, void *value, HashTable *t)
{
  HashTablePair *entry = (HashTablePair *)xmalloc(sizeof(HashTablePair));
  unsigned long h = t->hash(elm, t->parameters);

  entry->key = elm;
  entry->value = value;
  elist_append(t->table + h, entry);
}

/* If elm is already present in t, change its associated value to
 * value if update is NULL or returns true when invoked with the old
 * and the new value. If elm is not present, insert g with value
 * associated.  Returns 1 if elm was present and the associated value
 * was changed, 0 if elm was present but left unmodified, and -1 if
 * elm was inserted.
 */
int hashtable_update(void *elm, void *value, HashTable *t,
		     int (*update)(void *, void *))
{
  unsigned long h = t->hash(elm, t->parameters);
  int i = find_index(elm, t->table + h, t->compare);
  HashTablePair *entry;

  if (i < 0){
    /* g is not already present in t */
    hashtable_insert(elm, value, t);
    return -1;
  }
  else{
    /* g is present in t */
    entry = elist_get(t->table + h, i);
    if ((update == NULL) || update(entry->value, value)){
      /* Update g's associated value */
      entry->value = value;
      return 1;
    }
    else
      return 0;
  }
}

/* Determine the number of elements stored in t */
int hashtable_size(HashTable *t)
{
  int i, size = 0;

  for (i = 0; i < t->size; i++)
    size += t->table[i].count;

  return size;
}

/* Apply f to all key, value pairs in t */
void hashtable_map(HashTable *t, void (*f)(void *, void *, va_list), ...)
{
  int i, j;
  HashTablePair *entry;
  va_list args;

  for (i = 0; i < t->size; i++)
    for (j = 0; j < elist_length(t->table + i); j++){
      entry = elist_get(t->table + i, j);
      va_start(args, *f);
      f(entry->key, entry->value, args);
      va_end(args);
    }
}

/* Determine number of collisions in t */
int hashtable_collisions(HashTable *t)
{
  int i, collisions = 0;

  for (i = 0; i < t->size; i++)
    if (t->table[i].count > 1)
      collisions += t->table[i].count * (t->table[i].count - 1) / 2;

  return collisions;
}

/* Determine most number of elements sharing hash value */
int hashtable_largestbucket(HashTable *t)
{
  int i, largest = t->table[0].count;

  for (i = 1; i < t->size; i++)
    if (t->table[i].count > largest)
      largest = t->table[i].count;

  return largest;
}

/* Print the genes in (one of) the largest bucket(s) */
void hashtable_printlargestbucket(HashTable *t,
				  void (*print_elm)(void *, va_list))
{
  int i, l = hashtable_largestbucket(t);

  for (i = 0; t->table[i].count != l; i++);

  elist_map(t->table + i, print_elm);
}
