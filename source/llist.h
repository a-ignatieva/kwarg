/*******************************************************************

    llist.h
  
    Description of a general purpose linked list with a stack and
    queue interface etc.
		
    Christian Storm (cstorm@daimi.aau.dk), May 1995

********************************************************************/

#ifndef _LLIST_H
#define _LLIST_H

#include <stdarg.h>

#define FIRST -1
#define LAST  -2

/* Datatypes */
typedef struct tagLListNode
{
  void *elm;			
  struct tagLListNode *prev;
  struct tagLListNode *next;
} LListNode;

typedef struct tagLList
{
  int count;
  struct tagLListNode *first;
  struct tagLListNode *last;
} LList;

typedef struct tagLListCounter
{
  int pos;
  struct tagLListNode *current;
} LListCounter;

/* Prototypes */
LList *MakeLList();
void InitLList(LList *llist);
void DestroyLList(LList *llist);
void Insert(LList *llist, LListNode *lnode, void *elm);
void *Remove(LList *llist, LListNode *lnode);
void Push(LList *llist, void *elm);
void *Pop(LList *llist);
void *Top(LList *llist);
void Enqueue(LList *llist, void *elm);
void *Dequeue(LList *llist);
void Append(LList *llist1, LList *llist2);
void Prepend(LList *llist1, LList *llist2);
void *Bottom(LList *llist);
void *GetByIndex(LList *llist, unsigned int index);
void LListMap(LList *llist, void (*f)(void *, va_list), ...);
LListCounter *MakeCounter(LList *llist, int pos);
void InitCounter(LListCounter *lcounter, LList *llist, int pos);
void DestroyCounter(LListCounter *lcounter);
void *SetCounter(LListCounter *lcounter, int pos);
int GetPosition(LListCounter *lcounter);
void *Current(LListCounter *lcounter);
void *ChangeCurrent(LListCounter *lcounter, void *elm);
void *Next(LListCounter *lcounter);
void printLList(LList *llist);
void *Prev(LListCounter *lcounter);
void *RemoveMoveLeft(LList *llist, LListCounter *lcounter);
void *RemoveMoveRight(LList *llist, LListCounter *lcounter);
int Length(LList *llist);
void MergeSort(LList *llist, int (*less_than)(void *, void *));

#endif
