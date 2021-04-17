/*******************************************************************
*
*    arg.c :Implementation of functions for building and outputting an
*    ancestral recombination graph
*
********************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include "bitfunctions.h"
#include "llist.h"
#include "arg.h"
#include "common.h"

/* arg_new(): create an empty ARG */
ARG *arg_new()
{
  ARG *arg = (ARG *)xmalloc(sizeof(ARG));

  arg->nodes = (ARGNode *)xmalloc(sizeof(ARGNode));
  arg->n = 0;

  return arg;
}

/* Destroy arg and free all memory used by it */
void arg_destroy(ARG *arg)
{
  int i;

  for (i = 0; i < arg->n; i++){
    switch(arg->nodes[i].type){
    case ARGSAMPLE:
    case ARGCOALESCENCE:
      /* This node has one outgoing edge */
      while (Length(arg->nodes[i].predecessor.one.mutations) > 0)
	Pop(arg->nodes[i].predecessor.one.mutations);
      DestroyLList(arg->nodes[i].predecessor.one.mutations);
      break;
    case ARGRECOMBINATION:
      /* This node has two outgoing edges */
      while (Length(arg->nodes[i].predecessor.two.prefix.mutations) > 0)
	Pop(arg->nodes[i].predecessor.two.prefix.mutations);
      DestroyLList(arg->nodes[i].predecessor.two.prefix.mutations);
      while (Length(arg->nodes[i].predecessor.two.suffix.mutations) > 0)
	Pop(arg->nodes[i].predecessor.two.suffix.mutations);
      DestroyLList(arg->nodes[i].predecessor.two.suffix.mutations);
      break;
    case ARGANCESTOR:
      /* This node has no outgoing edges */
      break;
    }
    if (arg->nodes[i].label != NULL)
      free(arg->nodes[i].label);
    if (arg->nodes[i].sequence != NULL)
      free(arg->nodes[i].sequence);
  }
  free(arg->nodes);
  free(arg);
}

/* Add a type node to arg, with label and s as label and sequence. If
 * type is ARGRECOMBINATION, an extra argument specifying the first
 * site of the suffix contribution to the recombinant should be
 * provided. Return value is index of new node. The strings supplied
 * for label and sequence will be deallocated by arg_destroy.
 */
int arg_addnode(ARG *arg, ARGNodeType type, char *label, char *s,...)
{
  va_list args;

  if (weight(arg->n) == 1)
    /* Ran out of space for nodes - double size of array */
    arg->nodes = (ARGNode *)xrealloc(arg->nodes, 2 * arg->n * sizeof(ARGNode));

  arg->nodes[arg->n].type = type;
  arg->nodes[arg->n].label = label;
  arg->nodes[arg->n].sequence = s;

  switch(type){
  case ARGSAMPLE:
  case ARGCOALESCENCE:
    /* These types of nodes have one predecessor */
    arg->nodes[arg->n].predecessor.one.mutations = MakeLList();
    arg->nodes[arg->n].predecessor.one.target = -1;
    break;
  case ARGRECOMBINATION:
    /* These types of nodes have two predecessors */
    arg->nodes[arg->n].predecessor.two.prefix.mutations = MakeLList();
    arg->nodes[arg->n].predecessor.two.prefix.target = -1;
    arg->nodes[arg->n].predecessor.two.suffix.target = -1;
    va_start(args, s);
    arg->nodes[arg->n].predecessor.two.position = (int)va_arg(args, int);
    va_end(args);
    arg->nodes[arg->n].predecessor.two.suffix.mutations = MakeLList();
    break;
  case ARGANCESTOR:
    /* These types of nodes have no predecessors */
    break;
  }

  return (arg->n)++;
}

/* Return node with index i in arg */
ARGNode *arg_getnode(ARG *arg, int i)
{
  return arg->nodes + i;
}

/* Finalise arg, changing nodes with non-terminated outgoing edges to
 * ANCESTOR type.
 */
void arg_finalise(ARG *arg)
{
  int i;

  for (i = 0; i < arg->n; i++){
    switch(arg->nodes[i].type){
    case ARGSAMPLE:
    case ARGCOALESCENCE:
      if (arg->nodes[i].predecessor.one.target == -1){
	/* Sanity check */
	if (Length(arg->nodes[i].predecessor.one.mutations) > 0){
	  fprintf(stderr, "Something is wrong in ARG finalisation - please email error report\n");
	  exit(1);
	}
	DestroyLList(arg->nodes[i].predecessor.one.mutations);
	arg->nodes[i].type = ARGANCESTOR;
      }
      break;
    case ARGRECOMBINATION:
      if (arg->nodes[i].predecessor.two.prefix.target == -1)
	arg->nodes[i].predecessor.two.prefix.target
	  = arg_addnode(arg, ARGANCESTOR, NULL, NULL);
      if (arg->nodes[i].predecessor.two.suffix.target == -1)
	arg->nodes[i].predecessor.two.suffix.target
	  = arg_addnode(arg, ARGANCESTOR, NULL, NULL);
      break;
    case ARGANCESTOR:
      break;
    }
  }
}

/* Output labels in llist to fp, separated by semicolons */
static void output_edgelabels(LList *llist, AnnotatedGenes *a, int maiden,
			      FILE *fp)
{
  LListCounter *lcounter = MakeCounter(llist, FIRST);
  int i;
  void *p;

  for (i = 0; i < Length(llist); i++){
    p = Next(lcounter);
    if (!maiden)
      fprintf(fp, ";");
    if((int)p == -INT_MAX) {
        if (a->positions != NULL)
            fprintf(fp, "*%s", (char *)GetByIndex(a->positions, 0));
        else
            fprintf(fp, "*%d", 1);
    }
    else {
        if((int)p < 0) {
            if (a->positions != NULL)
                fprintf(fp, "*%s", (char *)GetByIndex(a->positions, -(int)p));
            else
                fprintf(fp, "*%d", -(int)p + 1);
        }
        else {
            if (a->positions != NULL)
            fprintf(fp, "%s", (char *)GetByIndex(a->positions, (int)p));
            else
            fprintf(fp, "%d", (int)p + 1);
        }
    }
    maiden = 0;
  }

  DestroyCounter(lcounter);
}

static void update_descendant(int (*a)[2], int d)
{
  if ((*a)[0] == -1)
    (*a)[0] = d;
  else
    (*a)[1] = d;
}

static void transfer_mutations(LList *from, LList *to, int startsite,
			       int endsite)
{
  LListCounter *lcounter = MakeCounter(from, FIRST);
  int i, site;

  for (i = 0; i < Length(from); i++){
    site = (int)Next(lcounter);
    if ((site >= startsite) && (site < endsite))
      Enqueue(to, (void *)site);
  }
  free(lcounter);
}

static void visit(ARG *arg, ARG *tree, int current, int (*descendant)[2],
		  int *tree_node, int startsite, int endsite)
{
  switch(arg->nodes[current].type){
  case ARGSAMPLE:
    update_descendant(descendant + arg->nodes[current].predecessor.one.target,
		      current);
    transfer_mutations
      (arg->nodes[current].predecessor.one.mutations,
       tree->nodes[tree_node[current]].predecessor.one.mutations,
       startsite, endsite);
    visit(arg, tree, arg->nodes[current].predecessor.one.target, descendant,
	  tree_node, startsite, endsite);
    break;
  case ARGCOALESCENCE:
  case ARGANCESTOR:
    if (descendant[current][1] == -1)
      /* We still need to visit the second descendant of this node */
      return;
    if (tree_node[descendant[current][0]] != -1){
      if (tree_node[descendant[current][1]] != -1){
	/* Both descendants are part of tree */
	memcpy(tree->nodes + tree->n, arg->nodes + current, sizeof(ARGNode));
	tree->nodes[tree_node[descendant[current][0]]].predecessor.one.target
	  = tree->n;
	tree->nodes[tree_node[descendant[current][1]]].predecessor.one.target
	  = tree->n;
	if (arg->nodes[current].type != ARGANCESTOR){
	  tree->nodes[tree->n].predecessor.one.mutations = MakeLList();
	  transfer_mutations(arg->nodes[current].predecessor.one.mutations,
			     tree->nodes[tree->n].predecessor.one.mutations,
			     startsite, endsite);
	  update_descendant
	    (descendant + arg->nodes[current].predecessor.one.target, current);
	}
	tree_node[current] = (tree->n)++;
      }
      else{
	/* Only first descendant is part of tree */
	if (arg->nodes[current].type == ARGANCESTOR){
	  /* First descendant is common ancestor of all samples */
	  tree->nodes[tree_node[descendant[current][0]]].type = ARGANCESTOR;
	  DestroyLList
	    (tree->nodes[tree_node[descendant[current][0]]].predecessor.one.mutations);
	}
	else{
	  /* Descendant of ancestor to current node in tree is first
	   * descendant to current node.
	   */
	  update_descendant
	    (descendant + arg->nodes[current].predecessor.one.target,
	     descendant[current][0]);
	  /* Copy relevant mutations from edge to predecessor in ARG
	   * to edge from first descendant in marginal tree.
	   */
	  transfer_mutations(arg->nodes[current].predecessor.one.mutations,
			     tree->nodes[tree_node[descendant[current][0]]].predecessor.one.mutations,
			     startsite, endsite);
	}
      }
    }
    else{
      if (tree_node[descendant[current][1]] != -1){
	/* Only second descendant is part of tree */
	if (arg->nodes[current].type == ARGANCESTOR){
	  /* Second descendant is common ancestor of all samples */
	  tree->nodes[tree_node[descendant[current][1]]].type = ARGANCESTOR;
	  DestroyLList
	    (tree->nodes[tree_node[descendant[current][1]]].predecessor.one.mutations);
	}
	else{
	  /* Descendant of ancestor to current node in tree is second
	   * descendant to current node.
	   */
	  update_descendant
	    (descendant + arg->nodes[current].predecessor.one.target,
	     descendant[current][1]);
	  /* Copy relevant mutations from edge to predecessor in ARG
	   * to edge from first descendant in marginal tree.
	   */
	  transfer_mutations(arg->nodes[current].predecessor.one.mutations,
			     tree->nodes[tree_node[descendant[current][1]]].predecessor.one.mutations,
			     startsite, endsite);
	}
      }
      else
	/* None of the descendants are part of the tree */
	update_descendant
	  (descendant + arg->nodes[current].predecessor.one.target, current);
    }
    if (arg->nodes[current].type != ARGANCESTOR)
      visit(arg, tree, arg->nodes[current].predecessor.one.target, descendant,
	    tree_node, startsite, endsite);
    break;
  case ARGRECOMBINATION:
    if (startsite < arg->nodes[current].predecessor.two.position){
      /* Predecessor providing prefix is predecessor at this site */
      update_descendant
	(descendant + arg->nodes[current].predecessor.two.prefix.target,
	 descendant[current][0]);
      /* Copy relevant mutations from edge to predecessor in ARG
       * to edge from first descendant in marginal tree.
       */
      if (tree_node[descendant[current][0]] != -1)
	transfer_mutations
	  (arg->nodes[current].predecessor.two.prefix.mutations,
	   tree->nodes[tree_node[descendant[current][0]]].predecessor.one.mutations,
	   startsite, endsite);
      visit(arg, tree, arg->nodes[current].predecessor.two.prefix.target,
	    descendant, tree_node, startsite, endsite);
      /* Visit non-predecessor */
      update_descendant
	(descendant + arg->nodes[current].predecessor.two.suffix.target,
	 current);
      visit(arg, tree, arg->nodes[current].predecessor.two.suffix.target,
	    descendant, tree_node, startsite, endsite);
    }
    else{
      /* Predecessor providing suffix is predecessor at this site */
      update_descendant
	(descendant + arg->nodes[current].predecessor.two.suffix.target,
	 descendant[current][0]);
      /* Copy relevant mutations from edge to predecessor in ARG
       * to edge from first descendant in marginal tree.
       */
      if (tree_node[descendant[current][0]] != -1)
	transfer_mutations
	  (arg->nodes[current].predecessor.two.suffix.mutations,
	   tree->nodes[tree_node[descendant[current][0]]].predecessor.one.mutations,
	   startsite, endsite);
      visit(arg, tree, arg->nodes[current].predecessor.two.suffix.target,
	    descendant, tree_node, startsite, endsite);
      /* Visit non-predecessor */
      update_descendant
	(descendant + arg->nodes[current].predecessor.two.prefix.target,
	 current);
      visit(arg, tree, arg->nodes[current].predecessor.two.prefix.target,
	    descendant, tree_node, startsite, endsite);
    }
    break;
  }
}

static void check_positions(void *site, va_list args)
{
  int i, *quote = va_arg(args, int *);
  LList *labels = va_arg(args, LList *);
  char *label = GetByIndex(labels, (int)site);

  /* If we have already seen a single quote there is no point in
   * looking for more.
   */
  if (!*quote)
    for (i = 0; i < strlen(label); i++)
      if (label[i] == '\'')
	*quote = 1;
}

static void newick_visit(ARG *tree, int node, AnnotatedGenes *a,
			 ARGLabels nodelabels, int annotate_edges,
			 int generate_id, MyString **parts, MyString *subtree)
{
  int i, j, single_quotes;
  char *label = tree->nodes[node].label;
  LListCounter *lcounter;

  if ((tree->nodes[node].type == ARGSAMPLE) || (parts[node] != NULL)){
    /* Last visit to this node */
    if (tree->nodes[node].type != ARGSAMPLE){
      /* Finish construction of bifucating node */
      mystring_addend(parts[node], ',');
      if (annotate_edges)
	/* Add space between nodes to improve readability */
	mystring_addend(parts[node], ' ');
      mystring_concat(parts[node], subtree);
      mystring_addend(parts[node], ')');
    }

    /* Add label to node */
    if (((label != NULL)
	 || ((tree->nodes[node].type == ARGSAMPLE) && (generate_id)))
	&& (nodelabels != ARGNONE) && (nodelabels != ARGSEQUENCE)
	&& ((nodelabels != ARGSEQUENCEFIRST)
	    || (tree->nodes[node].sequence == NULL))){
      /* The label of the node is used for actual labelling */
      /* If label contains single quotes we need to escape them */
      single_quotes = 0;
      if (label != NULL){
	/* Genuine label */
	for (i = 0; label[i] != '\0'; i++)
	  if (label[i] == '\'')
	    single_quotes++;
      }
      else
	/* Label is just sample index */
	label = i2a(node + 1);

      if (single_quotes > 0){
	/* Label did contain single quotes */
	j = single_quotes;
	label = (char *)xmalloc((i + j + 1) * sizeof(char));
	for (; i + j >= 0; i--){
	  label[i + j] = label[i];
	  if (label[i] == '\'')
	    label[i + --j] = '\'';
	}
      }

      /* Add label annotation to Newick representation */
      if (((nodelabels == ARGBOTH) && (tree->nodes[node].sequence != NULL))
	  || (single_quotes > 0))
	/* Node label needs to be quoted */
	mystring_addend(parts[node], '\'');
      mystring_append(parts[node], tree->nodes[node].label);
      if ((nodelabels == ARGBOTH) && (tree->nodes[node].sequence != NULL)){
	/* Both label and sequence should be used as annotation */
	mystring_addend(parts[node], ';');
	mystring_append(parts[node], tree->nodes[node].sequence);
	mystring_addend(parts[node], '\'');
      }
      else if (single_quotes > 0)
	/* End label quote */
	mystring_addend(parts[node], '\'');

      /* Clean up */
      if (tree->nodes[node].label == NULL)
	free(label);
    }
    else if ((tree->nodes[node].sequence != NULL) &&
	     (nodelabels != ARGNONE) && (nodelabels != ARGLABEL)){
      /* The sequence of the node is used for labelling */
      mystring_append(parts[node], tree->nodes[node].sequence);
    }

    if (tree->nodes[node].type != ARGANCESTOR){
      if (annotate_edges){
	/* Add mutations on edge from node to predecessor */
	if (Length(tree->nodes[node].predecessor.one.mutations) > 0){
	  mystring_addend(parts[node], ':');
	  lcounter = MakeCounter(tree->nodes[node].predecessor.one.mutations,
				 FIRST);
	  /* Check whether list of mutations needs to be quoted */
	  single_quotes = 0;
	  if (a->positions != NULL)
	    LListMap(tree->nodes[node].predecessor.one.mutations,
		     (void (*)(void *, va_list))check_positions,
		     &single_quotes, a->positions);
	  if (single_quotes)
	    /* It does */
	    mystring_addend(parts[node], '\'');
	  /* Add each mutation to 'length' of this node's ancestral edge */
	  for (i = 0; i < Length(tree->nodes[node].predecessor.one.mutations);
	       i++){
	    /* Separate mutations by semicolon */
	    if (i > 0)
	      mystring_addend(parts[node], ';');
	    /* Get next mutation and add it to length label */
	    j = (int)Next(lcounter);
	    if (a->positions != NULL){
	      if (single_quotes){
		/* Quote the single quotes in each site label */
		label = (char *)GetByIndex(a->positions, j);
		for (j = 0; j < strlen(label); j++){
		  mystring_addend(parts[node], label[j]);
		  if (label[j] == '\'')
		    mystring_addend(parts[node], '\'');
		}
	      }
	      else
		mystring_append
		  (parts[node], (char *)GetByIndex(a->positions, j));
	    }
	    else{
	      label = i2a(j + 1);
	      mystring_append(parts[node], label);
	      free(label);
	    }
	  }
	  free(lcounter);
	}
      }

      /* Recurse on parent node */
      if(tree->nodes[node].predecessor.one.target < tree->n) {
        newick_visit(tree, tree->nodes[node].predecessor.one.target, a, nodelabels, annotate_edges, generate_id, parts, parts[node]);
      }
    }
  }
  else{
    /* First visit to bifucating node - initialise construction */
    parts[node] = mystring_str2mystr("(");
    mystring_concat(parts[node], subtree);
  }
}

static void newick_output(ARG *tree, AnnotatedGenes *a, FILE *fp,
			  ARGLabels nodelabels, int annotate_edges,
			  int generate_id, char *label)
{
  int i;
  MyString **parts = xmalloc(tree->n * sizeof(MyString));
  char *s;

  /* Initialise parts */
  for (i = 0; i < tree->n; i++)
    if (tree->nodes[i].type == ARGSAMPLE)
      parts[i] = mystring_init();
    else
      parts[i] = NULL;

  /* Generate newick representation */
  for (i = 0; i < tree->n; i++)
    if (tree->nodes[i].type == ARGSAMPLE)
      newick_visit(tree, i, a, nodelabels, annotate_edges, generate_id, parts,
		   NULL);

  /* Output Newick representation */
  /* Ancestor should always be last node in ARG structure */
  s = mystring_mystr2str(parts[tree->n - 1]);
  fputs(s, fp);
  if (label != NULL)
    fputs(label, fp);
  putc(';', fp);
  putc('\n', fp);

  /* Clean up */
  /* All MyStrings in parts should have been destroyed by
   * concatenations to other MyStrings, except for the one at the
   * ancestor.
   */
  mystring_free(parts[tree->n - 1]);
  free(parts);
  free(s);
}

static void marginal_trees(ARG *arg, AnnotatedGenes *a, FILE *fp,
			   ARGFormat format, ARGLabels nodelabels,
			   int annotate_edges, int generate_id, int intervals,
			   char *label)
{
  int i, j, next, *starts, *tree_node, (*descendant)[2];
  ARG *tree;
  char *t;
  MyString *s;

  if (intervals){
    /* Start by determining the recombination free intervals */
    starts = (int *)xcalloc(a->g->length, sizeof(int));
    starts[0] = 1;
    for (i = 0; i < arg->n; i++)
      if (arg->nodes[i].type == ARGRECOMBINATION)
	starts[arg->nodes[i].predecessor.two.position] = 1;
    /* Convert Booleans to index of next interval start */
    next = a->g->length;
    for (i = a->g->length - 1; i >= 0; i--)
      if (starts[i]){
	starts[i] = next;
	next = i;
      }
  }
  else{
    /* Each position should be treated as its own interval */
    starts = (int *)xmalloc(a->g->length * sizeof(int));
    for (i = 0; i < a->g->length; i++)
      starts[i] = i + 1;
  }

  /* Generate and output evolutionary tree(s) for each recombination
   * free interval.
   */
  tree = (ARG *)xmalloc(sizeof(ARG));
  tree->nodes = xmalloc((2 * a->g->n - 1) * sizeof(ARGNode));
  tree->n = 0;
  i = 0;
  tree_node = (int *)xmalloc(arg->n * sizeof(int));
  descendant = (int (*)[2])xmalloc(arg->n * sizeof(int[2]));
  /* We need to copy all sample nodes to the marginal tree first to
   * have generate_id have the desired effect.
   */
  for (j = 0; j < arg->n; j++)
    if (arg->nodes[j].type == ARGSAMPLE){
      memcpy(tree->nodes + tree->n, arg->nodes + j, sizeof(ARGNode));
      tree->nodes[tree->n].predecessor.one.mutations = MakeLList();
      tree_node[j] = (tree->n)++;
    }

  /* Continue until we reach the end of the sequences */
  while (i < a->g->length){
    /* Wipe all non-sample nodes from tree */
    tree->n = a->g->n;
    /* Initialise arrays of descendants of a node and for mapping from
     * ARG nodes to marginal tree nodes.
     */
    for (j = 0; j < arg->n; j++)
      if (arg->nodes[j].type != ARGSAMPLE){
	descendant[j][0] = descendant[j][1] = -1;
	tree_node[j] = -1;
      }

    /* Start traversal back to root from each sample node */
    for (j = 0; j < arg->n; j++)
      if (arg->nodes[j].type == ARGSAMPLE)
	visit(arg, tree, j, descendant, tree_node, i, starts[i]);

    /* Copy sequence to tree nodes */
    for (j = 0; j < arg->n; j++)
      if ((tree_node[j] != -1) && (arg->nodes[j].sequence != NULL)){
	/* Copy relevant part of the sequence from ARG node to tree node */
	tree->nodes[tree_node[j]].sequence
	  = (char *)xmalloc((starts[i] - i + 1) * sizeof(char));
	strncpy(tree->nodes[tree_node[j]].sequence, arg->nodes[j].sequence,
		starts[i] - i);
	tree->nodes[tree_node[j]].sequence[starts[i] - i] = '\0';
      }
    /* Output current tree */
    /* GDL: graph can vaere et element i en anden graph; title: text angiver
     *      titel paa graf
     * GML: */
    /* Generate label specifying site or region tree is valid for */
    if (a->positions != NULL){
      s = mystring_str2mystr(GetByIndex(a->positions, i));
      if (i < starts[i] - 1){
	/* Tree valid for interval rather than single site */
	mystring_addend(s, '-');
	mystring_append(s, GetByIndex(a->positions, starts[i] - 1));
      }
    }
    else{
      t = i2a(i + 1);
      s = mystring_str2mystr(t);
      free(t);
      if (i < starts[i] - 1){
	/* Tree valid for interval rather than single site */
	mystring_addend(s, '-');
	t = i2a(starts[i]);
	mystring_append(s, t);
	free(t);
      }
    }
    t = mystring_mystr2str(s);
    mystring_free(s);
    if (format == TREENEWICK)
      newick_output(tree, a, fp, nodelabels, annotate_edges, generate_id, t);
    else if (format == TREEDOT)
      arg_output(tree, a, fp, ARGDOT, nodelabels, annotate_edges, generate_id);
    else if (format == TREEGDL)
      arg_output(tree, a, fp, ARGGDL, nodelabels, annotate_edges, generate_id);
    else if (format == TREEGML)
      arg_output(tree, a, fp, ARGGML, nodelabels, annotate_edges, generate_id);
    /* Clean up */
    free(t);
    /* Clear everything but non-site specific information of sample
     * nodes from marginal tree data structure.
     */
    for (j = 0; j < tree->n; j++){
      free(tree->nodes[j].sequence);
      if (tree->nodes[j].type != ARGANCESTOR){
	while (Length(tree->nodes[j].predecessor.one.mutations) > 0)
	  Pop(tree->nodes[j].predecessor.one.mutations);
	if (tree->nodes[j].type == ARGCOALESCENCE)
	  DestroyLList(tree->nodes[j].predecessor.one.mutations);
      }
    }
    i = starts[i];
  }

  /* Clean up */
  for (j = 0; j < tree->n; j++)
    if (tree->nodes[j].type == ARGSAMPLE)
      DestroyLList(tree->nodes[j].predecessor.one.mutations);
  free(tree->nodes);
  free(tree);
  free(starts);
  free(tree_node);
  free(descendant);
}

/* Output arg to file fp (stdout if fp is NULL). Sequence and site
 * information is taken from a, graph description language (DOT, GDL
 * or GML) is specified by format, node labelling (see definition of
 * ARGLabels in arg.h) is specified by nodelabels; if annotate_edges
 * is true edges are labelled with mutations, and if generate_id is
 * true sample nodes with no label in a are assigned their index in a
 * as label. If format is one of the formats prefixed with TREE an
 * extra argument should be provided, specifying whether one tree
 * for every site (if argument is false) or one tree for every
 * recombination free interval (if argument is true) should be output.
 */
void arg_output(ARG *arg, AnnotatedGenes *a, FILE *fp, ARGFormat format,
		ARGLabels nodelabels, int annotate_edges, int generate_id, ...)
{
  int i, maiden;
  va_list args;
  int intervals;

  /* Check whether marginal trees rather than the arg should be output */
  switch(format){
  case TREEDOT:
  case TREEGDL:
  case TREEGML:
  case TREENEWICK:
    /* Optional argument should specify whether a tree should be
     * output for every site, or only for every recombination free
     * interval.
     */
    va_start(args, generate_id);
    intervals = (int)va_arg(args, int);
    va_end(args);
    marginal_trees(arg, a, fp, format, nodelabels, annotate_edges, generate_id,
		   intervals, NULL);
    return;
  }

  if (fp == NULL)
    fp = stdout;

  /* Output prelude */
  switch(format){
  case ARGDOT:
    fprintf(fp, "digraph ARG {\n  { rank = same;");
    for (i = 0; i < a->g->n; i++)
      fprintf(fp, " %d;", i);
    fprintf(fp, " }\n");
    break;
  case ARGGDL:
    fprintf(fp, "graph: {\n");
    break;
  case ARGGML:
    fprintf(fp, "graph [\n  directed 1\n");
    break;
  }
  /* Output nodes and their edges */
  for (i = 0; i < arg->n; i++){
    /* Output node i */
    maiden = 1;
    switch(format){
    case ARGDOT:
      fprintf(fp, "  %d [label=\"", i);
      break;
    case ARGGDL:
      fprintf(fp, "  node: { title: \"%d\" label: \"", i);
      break;
    case ARGGML:
      fprintf(fp, "  node [\n    id %d\n    label \"", i);
      break;
    }
    /* Generate node label */
    if (((arg->nodes[i].label != NULL)
	 || ((arg->nodes[i].type == ARGSAMPLE) && generate_id))
	&& ((nodelabels == ARGLABEL) || (nodelabels == ARGBOTH)
	    || (nodelabels == ARGLABELFIRST)
	    || ((nodelabels == ARGSEQUENCEFIRST)
		&& (arg->nodes[i].sequence == NULL)))){
      if (arg->nodes[i].label != NULL)
	fprintf(fp, "%s", arg->nodes[i].label);
      else
	fprintf(fp, "%d", i + 1);
      maiden = 0;
    }
    if ((arg->nodes[i].sequence != NULL)
	&& ((nodelabels == ARGSEQUENCE) || (nodelabels == ARGBOTH)
	    || (nodelabels == ARGSEQUENCEFIRST)
	    || ((nodelabels == ARGLABELFIRST) && maiden))){
      if (!maiden)
	fprintf(fp, "; %s", arg->nodes[i].sequence);
      else
	fprintf(fp, "%s", arg->nodes[i].sequence);
      maiden = 0;
    }
    if (maiden){
      /* No node label printed */
      if ((arg->nodes[i].type == ARGRECOMBINATION)
	  && (arg->nodes[i].label != NULL))
	/* If a recombination node is still unlabeled, label it (label
	 * should be recombination point) if possible.
	 */
	fprintf(fp, "%s\"", arg->nodes[i].label);
      else{
	switch(format){
	case ARGDOT:
	  fprintf(fp, "\",shape=point");
	  break;
	case ARGGDL:
	  fprintf(fp, "\" scaling: 0.0");
	  break;
	case ARGGML:
	  fprintf(fp, " \"");
	  break;
	}
      }
    }
    else
      fprintf(fp, "\"");
    switch(format){
    case ARGDOT:
      fprintf(fp, ",color=");
      break;
    case ARGGDL:
      fprintf(fp, " shape: circle bordercolor: ");
      break;
    case ARGGML:
      fprintf(fp, "\n    graphics [\n      outline \"");
      break;
    }
    switch(arg->nodes[i].type){
    case ARGSAMPLE:
      fprintf(fp, "red");
      break;
    case ARGCOALESCENCE:
      fprintf(fp, "green");
      break;
    case ARGRECOMBINATION:
      fprintf(fp, "blue");
      break;
    case ARGANCESTOR:
      fprintf(fp, "yellow");
      break;
    }
    /* End node description */
    switch(format){
    case ARGDOT:
      fprintf(fp, "];\n");
      break;
    case ARGGDL:
      if (i < a->g->n)
	fprintf(fp, " vertical_order: maxlevel");
      fprintf(fp, " }\n");
      break;
    case ARGGML:
      fprintf(fp,
	      "\"\n    ]\n    vgj [\n      type \"Oval\"\n      labelPosition \"in\"\n    ]\n  ]\n");
      break;
    }
    /* Output edges out of this node */
    switch(arg->nodes[i].type){
    case ARGSAMPLE:
    case ARGCOALESCENCE:
      /* One outgoing edge */
      switch(format){
      case ARGDOT:
	fprintf(fp, "  %d -> %d", arg->nodes[i].predecessor.one.target, i);
	break;
      case ARGGDL:
	fprintf(fp, "  edge: { sourcename: \"%d\" targetname: \"%d\"",
		arg->nodes[i].predecessor.one.target, i);
	break;
      case ARGGML:
	fprintf(fp,
		"  edge [\n    linestyle \"solid\"\n    source %d\n    target %d",
		arg->nodes[i].predecessor.one.target, i);
	break;
      }
      /* Output edge label */
      if (annotate_edges && (arg->nodes[i].predecessor.one.mutations != NULL)
	  && (Length(arg->nodes[i].predecessor.one.mutations) > 0)){
	switch(format){
	case ARGDOT:
	  fprintf(fp, " [label=\"");
	  break;
	case ARGGDL:
	  fprintf(fp, " label: \"");
	  break;
	case ARGGML:
	  fprintf(fp, "\n    label \"");
	  break;
	}
	output_edgelabels(arg->nodes[i].predecessor.one.mutations, a, 1, fp);
	fprintf(fp, "\"");
	if (format == ARGDOT)
	  fprintf(fp, "]");
      }
      switch(format){
      case ARGDOT:
	fprintf(fp, ";\n");
	break;
      case ARGGDL:
	fprintf(fp, " }\n");
	break;
      case ARGGML:
	fprintf(fp, "\n  ]\n");
	break;
      }
      break;
    case ARGRECOMBINATION:
      /* Two outgoing edges */
      /* Prefix edge */
      switch(format){
      case ARGDOT:
	fprintf(fp, "  %d -> %d [label=\"P",
		arg->nodes[i].predecessor.two.prefix.target, i);
	break;
      case ARGGDL:
	fprintf(fp,
		"  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"P",
		arg->nodes[i].predecessor.two.prefix.target, i);
	break;
      case ARGGML:
	fprintf(fp,
		"  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"P",
		arg->nodes[i].predecessor.two.prefix.target, i);
	break;
      }
      if (annotate_edges
	  && (arg->nodes[i].predecessor.two.prefix.mutations != NULL)
	  && (Length(arg->nodes[i].predecessor.two.prefix.mutations) > 0))
	output_edgelabels(arg->nodes[i].predecessor.two.prefix.mutations, a, 0,
			  fp);
      switch(format){
      case ARGDOT:
	fprintf(fp, "\"]\n");
	break;
      case ARGGDL:
	fprintf(fp, "\" }\n");
	break;
      case ARGGML:
	fprintf(fp, "\"\n  ]\n");
	break;
      }
      /* Suffix edge */
      switch(format){
      case ARGDOT:
	fprintf(fp, "  %d -> %d [label=\"S",
		arg->nodes[i].predecessor.two.suffix.target, i);
	break;
      case ARGGDL:
	fprintf(fp,
		"  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"S",
		arg->nodes[i].predecessor.two.suffix.target, i);
	break;
      case ARGGML:
	fprintf(fp,
		"  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"S",
		arg->nodes[i].predecessor.two.suffix.target, i);
	break;
      }
      if (annotate_edges
	  && (arg->nodes[i].predecessor.two.suffix.mutations != NULL)
	  && (Length(arg->nodes[i].predecessor.two.suffix.mutations) > 0))
	output_edgelabels(arg->nodes[i].predecessor.two.suffix.mutations, a, 0,
			  fp);
      switch(format){
      case ARGDOT:
	fprintf(fp, "\"]\n");
	break;
	case ARGGDL:
	fprintf(fp, "\" }\n");
	break;
      case ARGGML:
	fprintf(fp, "\"\n  ]\n");
	break;
      }
      break;
    case ARGANCESTOR:
      /* No outgoing edges */
      break;
    }
  }
  /* Output postlude */
  switch(format){
  case ARGDOT:
  case ARGGDL:
    fprintf(fp, "}\n");
    break;
  case ARGGML:
    fprintf(fp, "]\n");
    break;
  }
}
