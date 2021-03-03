# Tools to use
CC = gcc

# Flags to the C compiler
# If DEBUG is 'yes', the executable is not optimised and compiled with
# information used by the gnu debugging and profiling tools
ifeq (${DEBUG}, yes)
CFLAGS = -Wall -g -pg -DDEBUG
else
CFLAGS = -O3 -w
endif
# If ENABLE_VERBOSE is defined, code for outputting information about the
# progress of the computation at several levels of detail is included
ifdef ENABLE_VERBOSE
CFLAGS := ${CFLAGS} -DENABLE_VERBOSE
endif
# If BEAGLE_NOHAPLOTYPEBLOCKS is not defined, code for computing all local
# minimum number of recombinations is included
ifndef BEAGLE_NOHAPLOTYPEBLOCKS
CFLAGS := ${CFLAGS} -DHAPLOTYPE_BLOCKS
endif
# If BEAGLE_HAPLOTYPE is defined, the full composite haplotype bound will be
# used in the branch & bound computation of the minimum number of
# recombinations required.
ifdef BEAGLE_HAPLOTYPE
CFLAGS := ${CFLAGS} -DBEAGLE_HAPLOTYPE
else
# If BEAGLE_HUDSONKAPLAN is not defined, the composite haplotype bound
# computed using the third heuristic (see HAPLOTYPE_HEURISTICS) with a
# required increase in the number of subsets of at least two is used to
# obtain lower bounds in the branch an bound computation of the minimum
# number of recombinations required. If BEAGLE_HUDSONKAPLAN is defined,
# Hudson & Kaplan's lower bound is used.
ifndef BEAGLE_HUDSONKAPLAN
CFLAGS := ${CFLAGS} -DBEAGLE_HAPLOTYPEHEURISTIC
endif
endif
# The branch & bound computation of the minimum number of recombinations
# required uses a hash table to store ancestral states already encountered to
# avoid tracing back from the same ancestral state more than once. I.e. a sort
# of dynamic programming to ensure that we only go through the same
# configuration at most once. If BEAGLE_STORELEAVES is not defined, the
# ancestral states that are eliminated directly from the lower bound are not
# stored in the hash table. This leads to a reduction in the amount of
# memory used to store ancestral states already visited, but increases the
# computational time as it is usually faster to look up an ancestral state in
# the hash table than to recompute the lower bound.
ifndef BEAGLE_STORELEAVES
CFLAGS := ${CFLAGS} -DBEAGLE_DONOTSTORELEAVES
endif
# A minimum number of recombinations computation similar to the branch and
# bound computation of the minimum number of recombinations is part of the
# code to enumerate all ancestral states visitable in minimum, or near minimum,
# histories. This always uses the reduced state space approach, but whether
# the ancestral states that are eliminated directly from the lower bound are
# stored in the hash table is determined by ENUMERATE_STORELEAVES.
ifndef ENUMERATE_STORELEAVES
CFLAGS := ${CFLAGS} -DENUMERATE_DONOTSTORELEAVES
endif
# The lower bound method used for computing minimum number of recombinations
# required in the visitable ancestral states enumeration is chosen in the
# same manner as for the branch and bound computation of the minimum number
# of recombinations required
ifdef ENUMERATE_HAPLOTYPE
CFLAGS := ${CFLAGS} -DENUMERATE_HAPLOTYPE
else
ifndef ENUMERATE_HUDSONKAPLAN
CFLAGS := ${CFLAGS} -DENUMERATE_HAPLOTYPEHEURISTIC
endif
endif
# Flags to the linker
ifeq (${PROFILE}, yes)
LDFLAGS := ${LDFLAGS} -pg
LDLIBS := -lm
ifneq (${DEBUG}, yes)
CFLAGS := ${CFLAGS} -pg
endif
else
ifeq (${DEBUG}, yes)
LDFLAGS := ${LDFLAGS}
LDLIBS := -lm
else
LDFLAGS := ${LDFLAGS} -O3
LDLIBS := -lm
endif
endif

SRCFILES = beagle.c kwarg.c simplify.c flip.c gene.c llist.c elist.c hashtable.c getwordsize.c common.c bounds.c exact.c mergesort.c enumerate.c arg.c backtrack.c topological_sort.c mystring.c
HDRFILES = gene.h llist.h elist.h hashtable.h common.h bounds.h exact.h mergesort.h enumerate.h arg.h backtrack.h topological_sort.h mystring.h
EXECUTABLES = beagle kwarg simplify flip getwordsize
DATA = kreitman_snp.txt
SCRIPTS = generatebitfunctions.py
INSTRUCTIONS = 
AUTOSRCFILES = bitfunctions.c
AUTOHDRFILES = bitfunctions.h

all: beagle kwarg simplify flip

beagle: beagle.o gene.o bounds.o exact.o llist.o elist.o hashtable.o common.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

backtrack.o: backtrack.c backtrack.h arg.h common.h

arg.o: arg.c arg.h llist.h gene.h common.h bitfunctions.h

mystring.o: mystring.c mystring.h llist.h

kwarg: kwarg.o gene.o bounds.o common.o exact.o llist.o elist.o hashtable.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

simplify: simplify.o gene.o bounds.o exact.o llist.o elist.o hashtable.o common.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

flip: flip.o gene.o bounds.o exact.o llist.o elist.o hashtable.o common.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

beagle.o: beagle.c gene.h bounds.h common.h exact.h Makefile

kwarg.o: kwarg.c gene.h bounds.h common.h exact.h backtrack.h Makefile

simplify.o: simplify.c gene.h bounds.h exact.h common.h backtrack.h Makefile

flip.o: flip.c gene.h bounds.h exact.h common.h backtrack.h Makefile

gene.o: gene.c gene.h common.h bitfunctions.h llist.h elist.h mergesort.h Makefile

bounds.o: bounds.c bounds.h exact.h gene.h common.h mergesort.h Makefile

exact.o: exact.c exact.h gene.h common.h mergesort.h bounds.h llist.h hashtable.h bitfunctions.h

enumerate.o: enumerate.c enumerate.h gene.h elist.h common.h bounds.h exact.h bitfunctions.h mergesort.h

greedy.o: greedy.c gene.h common.h bounds.h

llist.o: llist.c llist.h Makefile

elist.o: elist.c elist.h Makefile

hashtable.o: hashtable.c hashtable.h elist.h common.h Makefile

common.o: common.c common.h llist.h Makefile

bitfunctions.o: bitfunctions.c bitfunctions.h Makefile

bitfunctions.c bitfunctions.h: generatebitfunctions.py getwordsize Makefile
	python3 generatebitfunctions.py; rm -f getwordsize

getwordsize: getwordsize.c Makefile
	${CC} ${CFLAGS} getwordsize.c -o getwordsize

mergesort.o: mergesort.c mergesort.h Makefile

topological_sort.o: topological_sort.c topological_sort.h common.h Makefile

archive: ${SRCFILES} ${HDRFILES} ${DATA} ${SCRIPTS} ${INSTRUCTIONS} Makefile
	tar cvf kwarg.`date +%d%m%y`.tar ${SRCFILES} ${HDRFILES} ${DATA} ${SCRIPTS} ${INSTRUCTIONS} Makefile

# make clean removes all builds and all intermediates
.PHONY: clean
clean:
	rm -f *.o *~ core ${EXECUTABLES} ${AUTOSRCFILES} ${AUTOHDRFILES}

# make semiclean removes all intermediates, like object files, but leaves
# the executables
.PHONY: semiclean
semiclean:
	rm -f *.o ${AUTOSRCFILES} ${AUTOHDRFILES}

