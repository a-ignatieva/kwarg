/************************************************************************
 * 
 * getwordsize.c: A simple program that just prints the number of bits
 * in an unsigned long and in an unsigned int.
 * 
 ************************************************************************/

#include <stdio.h>
#include <limits.h>

int main(int argc, char **argv)
{
  printf("%d\n", CHAR_BIT * (int)sizeof(unsigned long));
  printf("%d\n", CHAR_BIT * (int)sizeof(unsigned int));

  return 0;
}
