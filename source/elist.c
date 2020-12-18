/*******************************************************************
 *   
 *   elist.c: Implementation of an expandable list with amortised constant time
 *   operations.
 *   
 ********************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "elist.h"

/* xmalloc(n): Allocate n bytes of memory, checking for successful allocation.
 */
static void *xmalloc(int n)
{
    void *adr;
    if ((adr = malloc(n)) == NULL){
        fprintf(stderr, "Unable to allocate sufficient amount of memory\n");
        exit(2);
    }
    
    return adr;
}

/* xrealloc(n): Resize memory allocated at oldadr to n bytes of memory,
 * checking for successful allocation.
 */
static void *xrealloc(void *oldadr, int n)
{
    void *adr;
    if ((adr = realloc(oldadr, n)) == NULL){
        fprintf(stderr, "Unable to allocate sufficient amount of memory\n");
        exit(2);
    }
    
    return adr;
}

/* Create an empty EList */
EList *elist_make()
{
    EList *elist = (EList *)xmalloc(sizeof(EList));
    elist_init(elist);
    
    return elist;
}

/* Initialise an EList structure */
void elist_init(EList *elist)
{
    elist->count = 0;
    elist->size = 0;
    elist->list = NULL;
}

/* Empty elist, applying release (if not NULL) to each element stored
 * in elist.
 */
void elist_empty(EList *elist, void (*release)(void *))
{
    int i;
    
    if (release != NULL)
        for (i = 0; i < elist->count; i++)
            release(elist->list[i]);
        
        elist->count = 0;
}

/* Free memory used by an EList structure, except for the EList
 * structure of elist itself.
 */
void elist_free(EList *elist)
{
    if (elist->list != NULL)
        free(elist->list);
    elist_init(elist);
}

/* Free memory used by an EList structure */
void elist_destroy(EList *elist)
{
    elist_free(elist);
    free(elist);
}

/* Insert a new element elm at the end of EList */
void elist_append(EList *elist, void *elm)
{
    if (elist->count == elist->size){
        /* We need to expand the list */
        if (elist->size)
            elist->size *= 2;
        else
            elist->size = 1;
        elist->list = (void **)xrealloc(elist->list, elist->size * sizeof(void *));
    }
    
    elist->list[elist->count++] = elm;
}

/* Append l2 to l1, destroying l2 in the process */
void elist_extend(EList *l1, EList *l2)
{
    int i;
    
    if (l2->count > 0){
        /* l2 does contain elements that needs to be moved */
        if (l1->count + l2->count > l1->size){
            /* And we will need more space */
            l1->size = 2 * (l1->count + l2->count);
            l1->list = (void **)xrealloc(l1->list, l1->size * sizeof(void *));
        }
        /* Now move elements from l2 to l1 */
        for (i = 0; i < l2->count; i++)
            l1->list[l1->count + i] = l2->list[i];
        l1->count += l2->count;
    }
    
    /* Destroy l2 */
    elist_destroy(l2);
}

/* Append l2 to l1, NOT destroying l2 in the process */
void elist_safeextend(EList *l1, EList *l2)
{
    int i;
    
    if (l2->count > 0){
        /* l2 does contain elements that needs to be moved */
        if (l1->count + l2->count > l1->size){
            /* And we will need more space */
            l1->size = 2 * (l1->count + l2->count);
            l1->list = (void **)xrealloc(l1->list, l1->size * sizeof(void *));
        }
        /* Now move elements from l2 to l1 */
        for (i = 0; i < l2->count; i++)
            l1->list[l1->count + i] = l2->list[i];
        l1->count += l2->count;
    }
    
}

/* Remove last element from elist and return it */
void *elist_removelast(EList *elist)
{
    if (elist->count)
        return elist->list[--elist->count];
    else
        return NULL;
}

void *elist_deletelast(EList *elist) {   
    void *elm;
    
    if (elist->count) {
        elm = elist->list[elist->count - 1];
        elist_remove(elist, elist->count - 1);
    }
    else
        elm = NULL;
    
    return elm;
}

/* Remove element at index and return it. This function takes time
 * proportional to the number of elements after the one removed.
 */
// AI bug fix: was not working properly (memmove using incorrect pointers).
void *elist_remove(EList *elist, unsigned int index)
{
    void *elm;
    
    if (index < elist->count){
        elm = elist->list[index];
        memmove(elist->list + index, elist->list + index + 1, (elist->count - index - 1)*sizeof(void *));
        elist->count--;
    }
    else
        elm = NULL;
    
    return elm;
}

/* Change the element at index to elm, returning the old element */
void *elist_change(EList *elist, unsigned int index, void *elm)
{
    void *oldelm;
    
    if (index < elist->count){
        oldelm = elist->list[index];
        elist->list[index] = elm;
    }
    else
        oldelm = NULL;
    
    return oldelm;
}

/* Swap elements in positions i and j - it is assumed that i and j are
 * valid indeces in elist.
 */
void elist_swap(EList *elist, int i, int j)
{
    void *tmp;
    
    if (i != j){
        tmp = elist->list[i];
        elist->list[i] = elist->list[j];
        elist->list[j] = tmp;
    }
}

/* Get element at index i - it is assumed that i is a valid index in elist */
void *elist_get(EList *elist, unsigned int i)
{
    return elist->list[i];
}

/* Apply f to all elements in elist, with extra arguments passed as an
 * va_list as second argument to f.
 */
void elist_map(EList *elist, void (*f)(void *, va_list), ...)
{
    unsigned int i;
    va_list args;
    
    for (i = 0; i < elist->count; i++){
        va_start(args, *f);
        (*f)(elist->list[i], args);
        va_end(args);
    }
}

int elist_length(EList *elist)
{
    return elist->count;
}
