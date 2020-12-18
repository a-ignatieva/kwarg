/*******************************************************************
*
*    mystring.c
*  
*    Implementation of functions for a concatenable string type.
*
********************************************************************/

#include <stdlib.h>
#include <string.h>
#include "mystring.h"
#include "llist.h"
#include "common.h"

/* mystring_init(): Return an empty string */
MyString *mystring_init()
{
  return (MyString *)MakeLList();
}

/* mystring_str2mystr(s): Return a MyString representation of s */
MyString *mystring_str2mystr(char *s)
{
  int i = 0;
  MyString *t = mystring_init();

  /* Transfer characters to linked list one at a time */
  while (s[i] != '\0')
    Enqueue((LList *)t, (void *)s[i++]);

  return t;
}

/* mystring_len(s): Return length of s */
int mystring_len(MyString *s)
{
  return Length((LList *)s);
}

/* mystring_mystr2str(s): Convert s to a string */
char *mystring_mystr2str(MyString *s)
{
  LListCounter *lcounter = MakeCounter((LList *)s, FIRST);
  void *c;
  char *t = (char *)xmalloc((mystring_len(s) + 1) * sizeof(char));
  int i = 0;

  /* Copy characters to string one at a time */
  while ((c = Next(lcounter)) != NULL)
    t[i++] = (char)c;
  /* Null terminate */
  t[i] = '\0';

  /* Clean up */
  free(lcounter);

  return t;
}

/* mystring_concat(s, t): Concatenate s and t into s, deleting t in
 * the process.
 */
MyString *mystring_concat(MyString *s, MyString *t)
{
  Append((LList *)s, (LList *)t);

  return s;
}

/* mystring_prepend(s, t): Prepend string t to s */
MyString *mystring_prepend(MyString *s, char *t)
{
  int i = strlen(t);

  for (; i > 0; Push((LList *)s, (void *)t[--i]));

  return s;
}

/* mysstring_append(s, t): Append string t to s */
MyString *mystring_append(MyString *s, char *t)
{
  int i;

  for (i = 0; t[i] != '\0'; i++)
    Enqueue((LList *)s, (void *)t[i]);

  return s;
}

/* mystring_addfront(s, c): Add c to front of s */
MyString *mystring_addfront(MyString *s, char c)
{
  Push((LList *)s, (void *)c);

  return s;
}

/* mystring_addend(s, c): Add c to end of s */
MyString *mystring_addend(MyString *s, char c)
{
  Enqueue((LList *)s, (void *)c);

  return s;
}

/* mystring_copy(s): Create a copy of s */
MyString *mystring_copy(MyString *s)
{
  LListCounter *lcounter = MakeCounter((LList *)s, FIRST);
  MyString *t = mystring_init();
  void *c;

  /* Copy characters of s one at a time */
  while ((c = Next(lcounter)) != NULL)
    Enqueue((LList *)t, c);

  /* Clean up */
  free(lcounter);

  return t;
}

void mystring_free(MyString *s)
{
  DestroyLList((LList *)s);
}
