#include "bitfunctions.h"

/* msb(a): finds index, going from 0 to the word size - 1, of most significant
 * bit that is 1 in a; if all bits are 0, -1 is returned.
 */
int msb(unsigned long a)
{
  if (a){
    /* Most significant bit is in 0 - 63 */
    if (a & 0xffffffff00000000)
      /* Most significant bit is in 32 - 63 */
      if (a & 0xffff000000000000)
        /* Most significant bit is in 48 - 63 */
        if (a & 0xff00000000000000)
          /* Most significant bit is in 56 - 63 */
          if (a & 0xf000000000000000)
            /* Most significant bit is in 60 - 63 */
            if (a & 0xc000000000000000)
              /* Most significant bit is in 62 - 63 */
              if (a & 0x8000000000000000)
                return 63;
              else
                return 62;
            else
              /* Most significant bit is in 60 - 61 */
              if (a & 0x2000000000000000)
                return 61;
              else
                return 60;
          else
            /* Most significant bit is in 56 - 59 */
            if (a & 0xc00000000000000)
              /* Most significant bit is in 58 - 59 */
              if (a & 0x800000000000000)
                return 59;
              else
                return 58;
            else
              /* Most significant bit is in 56 - 57 */
              if (a & 0x200000000000000)
                return 57;
              else
                return 56;
        else
          /* Most significant bit is in 48 - 55 */
          if (a & 0xf0000000000000)
            /* Most significant bit is in 52 - 55 */
            if (a & 0xc0000000000000)
              /* Most significant bit is in 54 - 55 */
              if (a & 0x80000000000000)
                return 55;
              else
                return 54;
            else
              /* Most significant bit is in 52 - 53 */
              if (a & 0x20000000000000)
                return 53;
              else
                return 52;
          else
            /* Most significant bit is in 48 - 51 */
            if (a & 0xc000000000000)
              /* Most significant bit is in 50 - 51 */
              if (a & 0x8000000000000)
                return 51;
              else
                return 50;
            else
              /* Most significant bit is in 48 - 49 */
              if (a & 0x2000000000000)
                return 49;
              else
                return 48;
      else
        /* Most significant bit is in 32 - 47 */
        if (a & 0xff0000000000)
          /* Most significant bit is in 40 - 47 */
          if (a & 0xf00000000000)
            /* Most significant bit is in 44 - 47 */
            if (a & 0xc00000000000)
              /* Most significant bit is in 46 - 47 */
              if (a & 0x800000000000)
                return 47;
              else
                return 46;
            else
              /* Most significant bit is in 44 - 45 */
              if (a & 0x200000000000)
                return 45;
              else
                return 44;
          else
            /* Most significant bit is in 40 - 43 */
            if (a & 0xc0000000000)
              /* Most significant bit is in 42 - 43 */
              if (a & 0x80000000000)
                return 43;
              else
                return 42;
            else
              /* Most significant bit is in 40 - 41 */
              if (a & 0x20000000000)
                return 41;
              else
                return 40;
        else
          /* Most significant bit is in 32 - 39 */
          if (a & 0xf000000000)
            /* Most significant bit is in 36 - 39 */
            if (a & 0xc000000000)
              /* Most significant bit is in 38 - 39 */
              if (a & 0x8000000000)
                return 39;
              else
                return 38;
            else
              /* Most significant bit is in 36 - 37 */
              if (a & 0x2000000000)
                return 37;
              else
                return 36;
          else
            /* Most significant bit is in 32 - 35 */
            if (a & 0xc00000000)
              /* Most significant bit is in 34 - 35 */
              if (a & 0x800000000)
                return 35;
              else
                return 34;
            else
              /* Most significant bit is in 32 - 33 */
              if (a & 0x200000000)
                return 33;
              else
                return 32;
    else
      /* Most significant bit is in 0 - 31 */
      if (a & 0xffff0000)
        /* Most significant bit is in 16 - 31 */
        if (a & 0xff000000)
          /* Most significant bit is in 24 - 31 */
          if (a & 0xf0000000)
            /* Most significant bit is in 28 - 31 */
            if (a & 0xc0000000)
              /* Most significant bit is in 30 - 31 */
              if (a & 0x80000000)
                return 31;
              else
                return 30;
            else
              /* Most significant bit is in 28 - 29 */
              if (a & 0x20000000)
                return 29;
              else
                return 28;
          else
            /* Most significant bit is in 24 - 27 */
            if (a & 0xc000000)
              /* Most significant bit is in 26 - 27 */
              if (a & 0x8000000)
                return 27;
              else
                return 26;
            else
              /* Most significant bit is in 24 - 25 */
              if (a & 0x2000000)
                return 25;
              else
                return 24;
        else
          /* Most significant bit is in 16 - 23 */
          if (a & 0xf00000)
            /* Most significant bit is in 20 - 23 */
            if (a & 0xc00000)
              /* Most significant bit is in 22 - 23 */
              if (a & 0x800000)
                return 23;
              else
                return 22;
            else
              /* Most significant bit is in 20 - 21 */
              if (a & 0x200000)
                return 21;
              else
                return 20;
          else
            /* Most significant bit is in 16 - 19 */
            if (a & 0xc0000)
              /* Most significant bit is in 18 - 19 */
              if (a & 0x80000)
                return 19;
              else
                return 18;
            else
              /* Most significant bit is in 16 - 17 */
              if (a & 0x20000)
                return 17;
              else
                return 16;
      else
        /* Most significant bit is in 0 - 15 */
        if (a & 0xff00)
          /* Most significant bit is in 8 - 15 */
          if (a & 0xf000)
            /* Most significant bit is in 12 - 15 */
            if (a & 0xc000)
              /* Most significant bit is in 14 - 15 */
              if (a & 0x8000)
                return 15;
              else
                return 14;
            else
              /* Most significant bit is in 12 - 13 */
              if (a & 0x2000)
                return 13;
              else
                return 12;
          else
            /* Most significant bit is in 8 - 11 */
            if (a & 0xc00)
              /* Most significant bit is in 10 - 11 */
              if (a & 0x800)
                return 11;
              else
                return 10;
            else
              /* Most significant bit is in 8 - 9 */
              if (a & 0x200)
                return 9;
              else
                return 8;
        else
          /* Most significant bit is in 0 - 7 */
          if (a & 0xf0)
            /* Most significant bit is in 4 - 7 */
            if (a & 0xc0)
              /* Most significant bit is in 6 - 7 */
              if (a & 0x80)
                return 7;
              else
                return 6;
            else
              /* Most significant bit is in 4 - 5 */
              if (a & 0x20)
                return 5;
              else
                return 4;
          else
            /* Most significant bit is in 0 - 3 */
            if (a & 0xc)
              /* Most significant bit is in 2 - 3 */
              if (a & 0x8)
                return 3;
              else
                return 2;
            else
              /* Most significant bit is in 0 - 1 */
              if (a & 0x2)
                return 1;
              else
                return 0;

  }
  /* If all bits are 0, return -1 */
  return -1;
}

/* lsb(a): finds index, going from 0 to the word size - 1, of least significant
 * bit that is 1 in a; if all bits are 0, -1 is returned.
 */
int lsb(unsigned long a)
{
  if (a){
    /* Least significant bit is in 0 - 63 */
    if (a & 0xffffffff)
      /* Least significant bit is in 0 - 31 */
      if (a & 0xffff)
        /* Least significant bit is in 0 - 15 */
        if (a & 0xff)
          /* Least significant bit is in 0 - 7 */
          if (a & 0xf)
            /* Least significant bit is in 0 - 3 */
            if (a & 0x3)
              /* Least significant bit is in 0 - 1 */
              if (a & 0x1)
                return 0;
              else
                return 1;
            else
              /* Least significant bit is in 2 - 3 */
              if (a & 0x4)
                return 2;
              else
                return 3;
          else
            /* Least significant bit is in 4 - 7 */
            if (a & 0x30)
              /* Least significant bit is in 4 - 5 */
              if (a & 0x10)
                return 4;
              else
                return 5;
            else
              /* Least significant bit is in 6 - 7 */
              if (a & 0x40)
                return 6;
              else
                return 7;
        else
          /* Least significant bit is in 8 - 15 */
          if (a & 0xf00)
            /* Least significant bit is in 8 - 11 */
            if (a & 0x300)
              /* Least significant bit is in 8 - 9 */
              if (a & 0x100)
                return 8;
              else
                return 9;
            else
              /* Least significant bit is in 10 - 11 */
              if (a & 0x400)
                return 10;
              else
                return 11;
          else
            /* Least significant bit is in 12 - 15 */
            if (a & 0x3000)
              /* Least significant bit is in 12 - 13 */
              if (a & 0x1000)
                return 12;
              else
                return 13;
            else
              /* Least significant bit is in 14 - 15 */
              if (a & 0x4000)
                return 14;
              else
                return 15;
      else
        /* Least significant bit is in 16 - 31 */
        if (a & 0xff0000)
          /* Least significant bit is in 16 - 23 */
          if (a & 0xf0000)
            /* Least significant bit is in 16 - 19 */
            if (a & 0x30000)
              /* Least significant bit is in 16 - 17 */
              if (a & 0x10000)
                return 16;
              else
                return 17;
            else
              /* Least significant bit is in 18 - 19 */
              if (a & 0x40000)
                return 18;
              else
                return 19;
          else
            /* Least significant bit is in 20 - 23 */
            if (a & 0x300000)
              /* Least significant bit is in 20 - 21 */
              if (a & 0x100000)
                return 20;
              else
                return 21;
            else
              /* Least significant bit is in 22 - 23 */
              if (a & 0x400000)
                return 22;
              else
                return 23;
        else
          /* Least significant bit is in 24 - 31 */
          if (a & 0xf000000)
            /* Least significant bit is in 24 - 27 */
            if (a & 0x3000000)
              /* Least significant bit is in 24 - 25 */
              if (a & 0x1000000)
                return 24;
              else
                return 25;
            else
              /* Least significant bit is in 26 - 27 */
              if (a & 0x4000000)
                return 26;
              else
                return 27;
          else
            /* Least significant bit is in 28 - 31 */
            if (a & 0x30000000)
              /* Least significant bit is in 28 - 29 */
              if (a & 0x10000000)
                return 28;
              else
                return 29;
            else
              /* Least significant bit is in 30 - 31 */
              if (a & 0x40000000)
                return 30;
              else
                return 31;
    else
      /* Least significant bit is in 32 - 63 */
      if (a & 0xffff00000000)
        /* Least significant bit is in 32 - 47 */
        if (a & 0xff00000000)
          /* Least significant bit is in 32 - 39 */
          if (a & 0xf00000000)
            /* Least significant bit is in 32 - 35 */
            if (a & 0x300000000)
              /* Least significant bit is in 32 - 33 */
              if (a & 0x100000000)
                return 32;
              else
                return 33;
            else
              /* Least significant bit is in 34 - 35 */
              if (a & 0x400000000)
                return 34;
              else
                return 35;
          else
            /* Least significant bit is in 36 - 39 */
            if (a & 0x3000000000)
              /* Least significant bit is in 36 - 37 */
              if (a & 0x1000000000)
                return 36;
              else
                return 37;
            else
              /* Least significant bit is in 38 - 39 */
              if (a & 0x4000000000)
                return 38;
              else
                return 39;
        else
          /* Least significant bit is in 40 - 47 */
          if (a & 0xf0000000000)
            /* Least significant bit is in 40 - 43 */
            if (a & 0x30000000000)
              /* Least significant bit is in 40 - 41 */
              if (a & 0x10000000000)
                return 40;
              else
                return 41;
            else
              /* Least significant bit is in 42 - 43 */
              if (a & 0x40000000000)
                return 42;
              else
                return 43;
          else
            /* Least significant bit is in 44 - 47 */
            if (a & 0x300000000000)
              /* Least significant bit is in 44 - 45 */
              if (a & 0x100000000000)
                return 44;
              else
                return 45;
            else
              /* Least significant bit is in 46 - 47 */
              if (a & 0x400000000000)
                return 46;
              else
                return 47;
      else
        /* Least significant bit is in 48 - 63 */
        if (a & 0xff000000000000)
          /* Least significant bit is in 48 - 55 */
          if (a & 0xf000000000000)
            /* Least significant bit is in 48 - 51 */
            if (a & 0x3000000000000)
              /* Least significant bit is in 48 - 49 */
              if (a & 0x1000000000000)
                return 48;
              else
                return 49;
            else
              /* Least significant bit is in 50 - 51 */
              if (a & 0x4000000000000)
                return 50;
              else
                return 51;
          else
            /* Least significant bit is in 52 - 55 */
            if (a & 0x30000000000000)
              /* Least significant bit is in 52 - 53 */
              if (a & 0x10000000000000)
                return 52;
              else
                return 53;
            else
              /* Least significant bit is in 54 - 55 */
              if (a & 0x40000000000000)
                return 54;
              else
                return 55;
        else
          /* Least significant bit is in 56 - 63 */
          if (a & 0xf00000000000000)
            /* Least significant bit is in 56 - 59 */
            if (a & 0x300000000000000)
              /* Least significant bit is in 56 - 57 */
              if (a & 0x100000000000000)
                return 56;
              else
                return 57;
            else
              /* Least significant bit is in 58 - 59 */
              if (a & 0x400000000000000)
                return 58;
              else
                return 59;
          else
            /* Least significant bit is in 60 - 63 */
            if (a & 0x3000000000000000)
              /* Least significant bit is in 60 - 61 */
              if (a & 0x1000000000000000)
                return 60;
              else
                return 61;
            else
              /* Least significant bit is in 62 - 63 */
              if (a & 0x4000000000000000)
                return 62;
              else
                return 63;

  }
  /* If all bits are 0, return -1 */
  return -1;
}

/* weight(a): counts the number of bits in a that are 1 */
unsigned int weight(unsigned long a)
{
  /* Adding blocks of size 1 */
  a = (a & 0x5555555555555555) + (a >> 1 & 0x5555555555555555);
  /* Adding blocks of size 2 */
  a = (a & 0x3333333333333333) + (a >> 2 & 0x3333333333333333);
  /* Adding blocks of size 4 */
  a = (a & 0x0f0f0f0f0f0f0f0f) + (a >> 4 & 0x0f0f0f0f0f0f0f0f);
  /* Adding blocks of size 8 */
  a = (a & 0x00ff00ff00ff00ff) + (a >> 8 & 0x00ff00ff00ff00ff);
  /* Adding blocks of size 16 */
  a = (a & 0x0000ffff0000ffff) + (a >> 16 & 0x0000ffff0000ffff);

  /* Adding blocks of size at most 32 */
  return (a & 0x00000000ffffffff) + (a >> 32);
}

/* Computes a divided by BLOCKSIZE */
unsigned int mulblocksize(unsigned int a)
{
  return a << 6;
}

/* Computes a divided by BLOCKSIZE */
unsigned int divblocksize(unsigned int a)
{
  return a >> 6;
}

/* Computes a modulo BLOCKSIZE */
unsigned int modblocksize(unsigned int a)
{
  return a & 63;
}
