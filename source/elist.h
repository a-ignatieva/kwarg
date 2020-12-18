/*******************************************************************

    llist.h
  
    Description of an expandable list

********************************************************************/

#ifndef _ELIST_H
#define _ELIST_H

#include <stdarg.h>

/* Datatypes */
typedef struct _EList {
  unsigned int count;  /* Number of elements currently in list */
  unsigned int size;   /* Maximum number of elements list can hold */
  void **list;         /* Actual list of elements */
} EList;

/* Prototypes */
EList *elist_make();
void elist_init(EList *elist);
void elist_empty(EList *elist, void (*release)(void *));
void elist_free(EList *elist);
void elist_destroy(EList *elist);
void elist_append(EList *elist, void *elm);
void elist_extend(EList *l1, EList *l2);
void elist_safeextend(EList *l1, EList *l2);
void *elist_removelast(EList *elist);
void *elist_remove(EList *elist, unsigned int index);
void *elist_change(EList *elist, unsigned int index, void *elm);
void elist_swap(EList *elist, int i, int j);
void *elist_get(EList *elist, unsigned int index);
void elist_map(EList *elist, void (*f)(void *, va_list), ...);
int elist_length(EList *elist);

#endif
