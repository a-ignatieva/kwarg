/*******************************************************************

    arg.h
  
    Description of functions for building and outputting an ancestral
    recombination graph

********************************************************************/

#ifndef _ARG_H
#define _ARG_H

#include "llist.h"
#include "mystring.h"

typedef enum {ARGDOT,    /* Output ARG as ARG in DOT format */
	      ARGGDL,    /* Output ARG as ARG in GDL format */
	      ARGGML,    /* Output ARG as ARG in GML format */
	      TREEDOT,   /* Output ARG as list of trees in DOT format */
	      TREEGDL,   /* Output ARG as list of trees in GDL format */
	      TREEGML,   /* Output ARG as list of trees in GML format */
	      TREENEWICK /* Output ARG as list of trees in Newick format */
} ARGFormat;

typedef enum {ARGSAMPLE,        /* Node represents a sampled sequence */
	      ARGCOALESCENCE,   /* Node represents a coalescence */
	      ARGRECOMBINATION, /* Node represents a recombination */
	      ARGANCESTOR       /* Node represents an ancestral sequence */
} ARGNodeType;

typedef struct _ARGEdge {
  int target;
  LList *mutations;
} ARGEdge;

typedef struct _ARGNode {
  ARGNodeType type;
  char *label;
  char *sequence;
  union {
    ARGEdge one;
    struct {
      ARGEdge prefix;
      ARGEdge suffix;
      int position;
    } two;
  } predecessor;
} ARGNode;

typedef struct _ARG {
  int n;
  ARGNode *nodes;
} ARG;

typedef enum { ARGNONE,          /* Do not label nodes */
	       ARGBOTH,          /* Label nodes with both label and sequence */
	       ARGSEQUENCE,      /* Label nodes only with sequence */
	       ARGLABEL,         /* Label nodes only with label */
	       ARGSEQUENCEFIRST, /* Label nodes with sequence, or with label
				  * if sequence is NULL. */
	       ARGLABELFIRST     /* Label nodes with label, or with sequence
				  * if label is NULL. */
} ARGLabels;

/* Postpone this inclusion until data types have been defined */
#include "gene.h"

/* Prototypes */
ARG *arg_new();
void arg_destroy(ARG *arg);
int arg_addnode(ARG *arg, ARGNodeType type, char *label, char *s, ...);
ARGNode *arg_getnode(ARG *arg, int i);
void arg_finalise(ARG *arg);
void arg_output(ARG *arg, AnnotatedGenes *a, FILE *fp, ARGFormat format,
		ARGLabels nodelabels, int annotate_edges, int generate_id,
		...);
#endif
