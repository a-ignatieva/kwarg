/*******************************************************************

    mystring.h
  
    Description of functions for a concatenable string type

********************************************************************/

#ifndef _MYSTRING_H
#define _MYSTRING_H

#include "llist.h"

typedef LList MyString;

MyString *mystring_init();
int mystring_len(MyString *s);
MyString *mystring_str2mystr(char *s);
char *mystring_mystr2str(MyString *s);
MyString *mystring_concat(MyString *s, MyString *t);
MyString *mystring_prepend(MyString *s, char *t);
MyString *mystring_append(MyString *s, char *t);
MyString *mystring_addfront(MyString *s, char c);
MyString *mystring_addend(MyString *s, char c);
MyString *mystring_copy(MyString *s);
void mystring_free(MyString *s);

#endif
