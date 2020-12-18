#!/usr/bin/python
from os import popen
from sys import stdout
from math import floor, log

# Convert a string with four bits to the corresponding hexadecimal character
def fourbit(s):
  a = 0
  for i in range(0, min(len(s), 4)):
    if s[-1 - i] == '1':
      a = a + (1 << i)
  return "%x" % a
# Convert a string with a binary representation of a pattern to a string with
# a hexadecimal representation of the same pattern
def bin2hex(s):
  t = ""
  while s != "":
    t = fourbit(s[-4:])[-1] + t
    s = s[:-4]
  return "0x" + t

# Generate code for finding most significant bit in time logarithmic in
# wordsize
def generate_msbcode(low, high, prefix = "", f = stdout):
  if (low < high):
    # We still need to check more ranges
    f.write(prefix + "/* Most significant bit is in " + repr(low) + " - " + repr(high) +  " */\n")
    # Generate test pattern
    mean = int((low + high) / 2 + 1)
    f.write(prefix + "if (a & " + bin2hex((high - mean + 1) * "1" + mean * "0") + ")\n")
    generate_msbcode(mean, high, prefix + "  ", f)
    f.write(prefix + "else\n")
    generate_msbcode(low, mean - 1, prefix + "  ", f)
  else:
    f.write(prefix + "return " + repr(low) + ";\n")

def generate_lsbcode(low, high, prefix = "", f = stdout):
  if (low < high):
    # We still need to check more ranges
    f.write(prefix + "/* Least significant bit is in " + repr(low) + " - " + repr(high) + " */\n")
    # Generate test pattern
    mean = int((low + high) / 2 + 1)
    f.write(prefix + "if (a & " + bin2hex((mean - low) * "1" + low * "0") + ")\n")
    generate_lsbcode(low, mean - 1, prefix + "  ", f)
    f.write(prefix + "else\n")
    generate_lsbcode(mean, high, prefix + "  ", f)
  else:
    f.write(prefix + "return " + repr(low) + ";\n")

def generate_weightcode(wordsize, prefix = "", f = stdout):
  blocksize = 1
  while (2 * blocksize < wordsize):
    # Build pattern
    s = blocksize * "0" + blocksize * "1"
    while len(s) < wordsize:
      s = s + s
    s = bin2hex(s[:wordsize])
    # Output code
    f.write(prefix + "/* Adding blocks of size " + repr(blocksize) + " */\n")
    f.write(prefix + "a = (a & " + s + ") + (a >> " + repr(blocksize) + " & " + s + ");\n")
    blocksize = blocksize * 2
  # One final sum to compute
  # Build pattern
  s = ""
  while len(s) < wordsize:
    s = blocksize * "0" + blocksize * "1" + s
  s = bin2hex(s[:wordsize])
  # Output code
  f.write("\n" + prefix + "/* Adding blocks of size at most " + repr(blocksize) + " */\n")
  f.write(prefix + "return (a & " + s + ") + (a >> " + repr(blocksize) + ");\n")

def generate_mulcode(wordsize, prefix = "", f = stdout):
  i = 0
  tmp = wordsize
  while tmp > 1:
    if (tmp & 1) == 1:
      # Wordsize is not a power of two so we use standard multiplication
      f.write(prefix + 'return a * ' + repr(wordsize) + ';\n')
      return
    tmp = tmp >> 1
    i = i + 1
  # Wordsize is a power of two
  f.write(prefix + 'return a << ' + repr(i) + ';\n')

def generate_divcode(wordsize, prefix = "", f = stdout):
  i = 0
  tmp = wordsize
  while tmp > 1:
    if (tmp & 1) == 1:
      # Wordsize is not a power of two so we use standard division
      f.write(prefix + 'return a / ' + repr(wordsize) + ';\n')
      return
    tmp = tmp >> 1
    i = i + 1
  # Wordsize is a power of two
  f.write(prefix + 'return a >> ' + repr(i) + ';\n')

def generate_modcode(wordsize, prefix = "", f = stdout):
  i = 0
  tmp = wordsize
  while tmp > 1:
    if (tmp & 1) == 1:
      # Wordsize is not a power of two so we use standard modulo
      f.write(prefix + 'return a % ' + repr(wordsize) + ';\n')
      return
    tmp = tmp >> 1
    i = i + 1
  # Wordsize is a power of two
  f.write(prefix + 'return a & ' + repr(wordsize - 1) + ';\n')

file = popen("./getwordsize")
longsize = int(file.readline())
intsize = int(file.readline())
file.close()
file = open("bitfunctions.h", "w")
file.write("#ifndef _BITFUNCTIONS_H\n#define _BITFUNCTIONS_H\n\n#define BLOCKSIZE %d\n#define TERNARY_BLOCKSIZE %d\n\nint msb(unsigned long a);\nint lsb(unsigned long a);\nunsigned int weight(unsigned long a);\nunsigned int mulblocksize(unsigned int a);\nunsigned int divblocksize(unsigned int a);\nunsigned int modblocksize(unsigned int a);\n#endif\n" % (longsize, int(floor(intsize * log(2) / log(3)))))
file.close()
file = open("bitfunctions.c", "w")
file.write('#include "bitfunctions.h"\n\n/* msb(a): finds index, going from 0 to the word size - 1, of most significant\n * bit that is 1 in a; if all bits are 0, -1 is returned.\n */\nint msb(unsigned long a)\n{\n  if (a){\n')
generate_msbcode(0, longsize - 1, "    ", file)
file.write('\n  }\n  /* If all bits are 0, return -1 */\n  return -1;\n}\n')
file.write('\n/* lsb(a): finds index, going from 0 to the word size - 1, of least significant\n * bit that is 1 in a; if all bits are 0, -1 is returned.\n */\nint lsb(unsigned long a)\n{\n  if (a){\n')
generate_lsbcode(0, longsize - 1, "    ", file)
file.write('\n  }\n  /* If all bits are 0, return -1 */\n  return -1;\n}\n')
file.write('\n/* weight(a): counts the number of bits in a that are 1 */\nunsigned int weight(unsigned long a)\n{\n')
generate_weightcode(longsize, "  ", file)
file.write('}\n\n/* Computes a divided by BLOCKSIZE */\nunsigned int mulblocksize(unsigned int a)\n{\n')
generate_mulcode(longsize, "  ", file)
file.write('}\n\n/* Computes a divided by BLOCKSIZE */\nunsigned int divblocksize(unsigned int a)\n{\n')
generate_divcode(longsize, "  ", file)
file.write('}\n\n/* Computes a modulo BLOCKSIZE */\nunsigned int modblocksize(unsigned int a)\n{\n')
generate_modcode(longsize, "  ", file)
file.write('}\n')
file.close()
