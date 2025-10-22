#ifndef _BITFUNCTIONS_H
#define _BITFUNCTIONS_H

#define BLOCKSIZE 64
#define TERNARY_BLOCKSIZE 20

int msb(unsigned long a);
int lsb(unsigned long a);
unsigned int weight(unsigned long a);
unsigned int mulblocksize(unsigned int a);
unsigned int divblocksize(unsigned int a);
unsigned int modblocksize(unsigned int a);
#endif
