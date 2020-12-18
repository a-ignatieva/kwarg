/*******************************************************************
*
*    llist.c
*  
*    Implementation of a general purpose linked list with a stack and
*    queue interface.
*		
*    Christian Storm (cstorm@daimi.aau.dk), May 1995
*
*    Revision History:
*
*    January 2001. Modified LListMap to take a variable list of
*                  arguments, instead of a void pointer; minor
*		           change to InitCounter to allow the use of FIRST
*		           and LAST.
*		           Rune Lyngsø (rlyngsoe@brics.dk)
*
*    March 2005. Added functions Append and Prepend
*                Rune Lyngsø (lyngsoe@stats.ox.ac.uk)
*
*    April 2007. Added MergeSort
*                Rune Lyngsø (lyngsoe@stats.ox.ac.uk)
*
********************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <stdarg.h>
#include "llist.h"

/* Wrapper for malloc; checks whether memory was available. */
static void *xmalloc(int n)
{
  void *adr;
  if ((adr = malloc(n)) == NULL){
    fprintf(stderr, "Virtual memory exhausted.\n");
    exit(1);
  }

  return adr;
}

/* Create an empty list. */
LList *MakeLList()
{
  /* Allocate space and initialize a 'LList' */
  LList *llist = (LList *)xmalloc(sizeof(LList));
  llist->first = (LListNode *)xmalloc(sizeof(LListNode));
  llist->last = (LListNode *)xmalloc(sizeof(LListNode));
  llist->count = 0;

  /* Initialize the `stopblocks' */
  (llist->first)->prev = NULL;
  (llist->first)->next = llist->last;
  (llist->first)->elm = NULL;
  (llist->last)->prev = llist->first;
  (llist->last)->next = NULL;
  (llist->last)->elm = NULL;
  
  return llist;
}

/* Initialize an empty list. */
void InitLList(LList *llist)
{
  /* Allocate space and initialize a 'LList' */
  llist->first = (LListNode *)xmalloc(sizeof(LListNode));
  llist->last = (LListNode *)xmalloc(sizeof(LListNode));
  llist->count = 0;

  /* Initialize the `stopblocks' */
  (llist->first)->prev = NULL;
  (llist->first)->next = llist->last;
  (llist->first)->elm = NULL;
  (llist->last)->prev = llist->first;
  (llist->last)->next = NULL;
  (llist->last)->elm = NULL;
}

/* Destroys the LList-structure. That is, all except the elements
 * possibly pointed to by the elm-field in the LListNodes. The
 * function makes it very easy to lose all contact with these
 * elements ... use with care.  
 */
void DestroyLList(LList *llist)
{
  LListNode *temp;

  /* free all the LListNodes */
  while ((temp = llist->first->next) != NULL)
    {
      llist->first->next = temp->next;
      free(temp);
    }
  free(llist->first);

  /* free the LList */
  free(llist);
}

/* Inserts a new node to the right of `lnode', thus it would be a BIG
 * mistake to call with `lnode' equal to `llist->last'.  
 */
void Insert(LList *llist, LListNode *lnode, void *elm)
{
  LListNode *newnode = (LListNode *)xmalloc(sizeof(LListNode));

  /* Initialize the new node */
  newnode->elm = elm;
  newnode->prev = lnode;
  newnode->next = lnode->next;

  /* Update the list */
  (newnode->prev)->next = newnode;
  (newnode->next)->prev = newnode;
  (llist->count)++;
}

/* Removes the node `lnode'. It shall only be used the remove nodes
 * inserted using `Insert'.  
 */
void *Remove(LList *llist, LListNode *lnode)
{
  void *elm = lnode->elm;
  
  (lnode->prev)->next = lnode->next;
  (lnode->next)->prev = lnode->prev;
  (llist->count)--;
  free(lnode);

  return elm;
}

/* Functions implementing a stack- and queue-interface to the linked
 * list using the generic `Insert' and `Remove' functions.
 */
void Push(LList *llist, void *elm)
{
  Insert(llist, llist->first, elm);
}

void *Pop(LList *llist)
{
  void *elm = NULL;

  if ((llist->first)->next != llist->last)
    /* The list is not empty */
    elm = Remove(llist, (llist->first)->next);

  return elm;
}

void *Top(LList *llist)
{
  /* If list is empty we return the elm member of the terminal
   * `stopblock', which is NULL.
   */
  return llist->first->next->elm;
}

void Enqueue(LList *llist, void *elm)
{
  Insert(llist, (llist->last)->prev, elm);
}

void *Dequeue(LList *llist)
{
  void *elm = NULL;

  if ((llist->last)->prev != llist->first)
    /* The list is not empty */
    elm = Remove(llist, (llist->last)->prev);

  return elm;
}

/* Append llist2 to llist1, destroying llist2 in the process. */
void Append(LList *llist1, LList *llist2)
{
  llist1->count += llist2->count;
  llist1->last->prev->next = llist2->first->next;
  llist2->first->next->prev = llist1->last->prev;
  free(llist1->last);
  free(llist2->first);
  llist1->last = llist2->last;
  free(llist2);
}

/* Prepend llist1 to llist2, destroying llist1 in the process. */
void Prepend(LList *llist1, LList *llist2)
{
  llist2->count += llist1->count;
  llist1->last->prev->next = llist2->first->next;
  llist2->first->next->prev = llist1->last->prev;
  free(llist1->last);
  free(llist2->first);
  llist2->first = llist1->first;
  free(llist1);
}

/* Return last element of llist */
void *Bottom(LList *llist)
{
  /* If list is empty we return the elm member of the initial
   * `stopblock', which is NULL.
   */
  return llist->last->prev->elm;
}

/* Return element at position index (with first position numbered 0)
 * in llist, or NULL if llist contains at most index elements.
 */
void *GetByIndex(LList *llist, unsigned int index)
{
  int i;
  LListNode *current;

  if (index >= llist->count)
    return NULL;

  if (index < llist->count / 2){
    /* Position closest to start of llist */
    current = llist->first;
    for (i = 0; i <= index; i++)
      current = current->next;
  }
  else{
    /* Position closest to end of llist */
    current = llist->last;
    for (i = llist->count; i > index; i--)
      current = current->prev;
  }

  return current->elm;
}

/* Functions used to iterate through a linked list. */
void LListMap(LList *llist, void (*f)(void *, va_list), ...)
{
  LListNode *current = llist->first;
  va_list args;

  while ((current = current->next) != llist->last){
    va_start(args, *f);
    (*f)(current->elm, args);
    va_end(args);
  }
}

/* Allocate and return a new LListCounter initialised to LList llist
 * and position pos.
 */
LListCounter *MakeCounter(LList *llist, int pos)
{
  LListCounter *lcounter = (LListCounter *)xmalloc(sizeof(LListCounter));

  InitCounter(lcounter, llist, pos);
  
  return lcounter;
}

/* Init LListCounter lcounter to LList llist and position pos */
void InitCounter(LListCounter *lcounter, LList *llist, int pos)
{
  if (pos == LAST)
    pos = llist->count;

  if (pos <= llist->count - 1 - pos)
    {
      /* The node 'pos' is not closest to the last node */
      lcounter->pos = -1;
      lcounter->current = llist->first;
    }
  else
    {
      /* The node 'pos' is closest to the last node */
      lcounter->pos = llist->count;
      lcounter->current = llist->last;
    }
  SetCounter(lcounter, pos);
}

/* Free memory used by an LListCounter */
void DestroyCounter(LListCounter *lcounter)
{
  free(lcounter);
}

/* Positions are from -1 to the number of elements n in the llist; -1
 * sets the counter to the initial `stopblock', n sets the counter to
 * the terminal `stopblock', and a position i in the range from 0 to n
 * - 1 sets the counter to the i'th member of the llist (counting from
 * 0). FIRST (LAST) will set the counter to the initial (terminal)
 * `stopblock'.
 */
void *SetCounter(LListCounter *lcounter, int pos)
{
  if (pos == LAST)
    {
      /* Move to the end of the llist */
      while ((lcounter->current)->next != NULL)
	{
	  (lcounter->pos)++;
	  lcounter->current = (lcounter->current)->next;
	}
    }
  else if (pos > lcounter->pos)
    {
      /* Move 'pos - lcounter->pos' nodes to the right */
      while (lcounter->pos < pos && (lcounter->current)->next != NULL)
	{
	  (lcounter->pos)++;
	  lcounter->current = (lcounter->current)->next;
	}
    }
  else if (pos < lcounter->pos)
    {
      /* Move 'lcounter->pos - pos' nodes to the left */
      while (lcounter->pos > pos && (lcounter->current)->prev != NULL)
	{
	  (lcounter->pos)--;
	  lcounter->current = (lcounter->current)->prev;
	}
    }

  return (lcounter->current)->elm;
}

/* Return element at current position */
void *Current(LListCounter *lcounter)
{
  return lcounter->current->elm;
}

/* Change element in current position, returning the old element */
void *ChangeCurrent(LListCounter *lcounter, void *elm)
{
  void *oldelm = lcounter->current->elm;

  lcounter->current->elm = elm;

  return oldelm;
}

/* Return the current position of the counter */
int GetPosition(LListCounter *lcounter)
{
  if (lcounter->current->next == NULL)
    return LAST;
  else
    return lcounter->pos;
}

/* Update counter to point to next element and return this */
void *Next(LListCounter *lcounter)
{
  lcounter->pos += 1;
  lcounter->current = lcounter->current->next;
  return lcounter->current->elm;
}

/* Update counter to point to previous element and return this */
void *Prev(LListCounter *lcounter)
{
  lcounter->pos -= 1;
  lcounter->current = lcounter->current->prev;
  return lcounter->current->elm;
}

/* Removes the current node of a counter. It shall only be used when
 * the counter is not at FIRST or LAST. The counter is set to the
 * preceeding element.
 */
void *RemoveMoveLeft(LList *llist, LListCounter *lcounter)
{
  LListNode *tmp = lcounter->current;
  void *elm = tmp->elm;

  Prev(lcounter);
  Remove(llist, tmp);
  return elm;
}

/* Removes the current node of a counter. It shall only be used when
 * the counter is not at FIRST or LAST. The counter is set to the
 * following element.
 */
void *RemoveMoveRight(LList *llist, LListCounter *lcounter)
{
  LListNode *tmp = lcounter->current;
  void *elm = tmp->elm;

  lcounter->current = lcounter->current->next;
  Remove(llist, tmp);

  return elm;
}

/* returns number of items in llist */
int Length(LList *llist)
{
  return llist->count;
}

static LListNode *_merge(LListNode *a, LListNode *b,
		   int (*less_than)(void *, void *))
{
  LListNode *new, *tmp;

  /* Find first element in merged list */
  if (less_than(a->elm, b->elm)){
    new = tmp = a;
    a = a->next;
  }
  else{
    new = tmp = b;
    b = b->next;
  }

  /* Continue merging as long as both lists are non-empty */
  while ((a != NULL) && (b != NULL))
    if (less_than(a->elm, b->elm)){
      tmp->next = a;
      a = a->next;
    }
    else{
      tmp->next = b;
      b = b->next;
    }

  /* Append remaining non-empty list to merged list */
  if (a != NULL)
    tmp->next = a;
  else
    tmp->next = b;

  return new;
}

/* Sort llist in ascending order according to less_than, reusing the nodes */
void MergeSort(LList *llist, int (*less_than)(void *, void *))
{
  LListNode *tmp, **sublist;
  int i, n = Length(llist);

  /* Empty list (or for that matter list containing one element) is
   * already sorted.
   */
  if (n <= 1)
    return;

  /* Set up structure for holding sorted sublists */
  sublist = (LListNode **)xmalloc(n * sizeof(LListNode *));
  sublist[0] = llist->first->next;
  for (i = 1; i < n; i++){
    sublist[i] = (sublist[i - 1])->next;
    sublist[i - 1]->next = NULL;
  }
  sublist[n - 1]->next = NULL;

  /* Do the rounds of merging */
  while (n > 1){
    for (i = 0; i < n / 2; i++)
      sublist[i] = _merge(sublist[2 * i], sublist[2 * i + 1], less_than);
    /* With an odd number of lists the last one is not merged with other
     * lists - we still need to move it to its proper place though.
     */
    if (n & 1)
      sublist[n / 2] = sublist[n - 1];
  }

  /* First element in tmp should now be the beginning of the sorted
   * list; reattach the flanking stop blocks and update the pointers
   * to the previous element of the list.
   */
  llist->first->next = sublist[0];
  tmp = llist->first;
  while (tmp->next != NULL){
    tmp->next->prev = tmp;
    tmp = tmp->next;
  }
  llist->last->prev = tmp;

  /* Clean up */
  free(sublist);
}
