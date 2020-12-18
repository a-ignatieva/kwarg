#ifndef _GENE_H
#define _GENE_H

#include <stdio.h>

#include "elist.h"
#include "llist.h"
#include "hashtable.h"

/* Datatypes */
typedef struct _Gene {
  unsigned long *type;      /* Type bit vector */
  unsigned long *ancestral; /* Ancestral site bit vector */
} Gene;

typedef struct _Genes {
  int n;      /* Number of sequences */
  int length; /* Length of sequences */
  Gene *data; /* Sequences */
} Genes;

typedef struct _AnnotatedGenes {
  Genes *g;         /* Data */
  LList *positions; /* Site labels */
  LList *sequences; /* Sequence labels */
} AnnotatedGenes;

typedef enum { GENE_ANY, GENE_BEAGLE, GENE_FASTA } Gene_Format;
typedef enum { GENE_BINARY, GENE_NUCLEIC, GENE_AMINO } Gene_SeqType;

typedef struct _PackedGenes {
  int n;               /* Number of sequences */
  int length;          /* Length of sequences */
  unsigned int *data; /* Sequences */
} PackedGenes;

typedef struct _Site {
  unsigned long *type;      /* Type bit vector */
  unsigned long *ancestral; /* Ancestral sequence bit vector */
} Site;

typedef struct _Sites {
  int n;      /* Number of sequences */
  int length; /* Length of sequences */
  Site *data; /* Sites */
} Sites;

typedef struct _Index {
  int index;
  int block;
} Index;

/* Postpone this inclusion until data types have been defined */
#include "backtrack.h"

/* Global variable for specifying whether the common ancestral
 * sequence is known or not.
 */
extern int gene_knownancestor;
/* Prototypes */
void free_genes(Genes *g);
Genes *make_genes();
void free_annotatedgenes(AnnotatedGenes *g);
void free_sites(Sites *s);
Sites *make_sites();
void free_gene(Gene *g);
void free_site(Site *s);
Genes *new_genes(int n, int length, char **s, Gene_SeqType t);
AnnotatedGenes *read_genes(char *fname, Gene_Format f, Gene_SeqType t);
void output_genes(Genes *g, FILE *fp, char *comment);
void output_labelled_genes(Genes *g, FILE *fp, LList *labels);
void output_genes_indexed(Genes *s, FILE *fp);
void output_annotatedgenes(AnnotatedGenes *a, FILE *fp, char *comment);
void add_gene(Genes *g, Gene *new, ...);
void add_site(Sites *s, Site *new, ...);
Gene *get_gene(Genes *g, int i);
Site *get_site(Sites *s, int i);
Genes *select_genes(Genes *g, EList *l);
Sites *select_sites(Sites *s, EList *l);
char get_genes_character(Genes *g, int seq, int site);
char get_sites_character(Sites *s, int seq, int site);
void set_genes_character(Genes *g, int seq, int site, char c);
void set_sites_character(Sites *s, int seq, int site, char c);
void swap_genes(Genes *g, int a, int b);
void add_ancestral(Genes *g);
void add_ancestral_sites(Sites *s);
Genes *copy_genes(Genes *g);
Sites *copy_sites(Sites *s);
Genes *copy_allbutone(Genes *g, int a);
Genes *copy_region(Genes *g, int a, int b);
void remove_gene(Genes *g, int a);
void remove_annotatedgene(AnnotatedGenes *g, int a);
void reallocate_genes(Genes *g);
Genes *random_genes(int n, int m);
Sites *genes2sites(Genes *g);
Genes *sites2genes(Sites *s);
char **genes2strings(Genes *g);
Genes *strings2genes(char **s);
Gene *string2gene(char *s);
Site *string2site(char *s);
char *gene2string(Gene *g, int length);
char *site2string(Site *s, int n);
char **genes2string(Genes *g);
void output_sites(Sites *s, FILE *fp, char *comment);
void output_sites_indexed(Sites *s, FILE *fp);
int remove_siamesetwins(Genes *g);
int remove_uninformative(Genes *g);
int remove_nonsegregating(Genes *g);
int coalesce_subsumed(Genes *g);
int reduce_coalesce(Genes *g, int *elements);
int implode_genes(Genes *g);
int no_recombinations_required(Genes *g);
void force_safeevents(Genes *g);
int force_mutations(Genes *g);
int mutate(Genes *g, int pos, int mutant);
EList *force_mutation(Genes *g, EList *event);
int segregating_site(Genes *g, int i);
int conflicting_sites(Genes *g);
int next_ancestral(Genes *g, int a, int i);
int identical(Genes *g, int a, int b);
int compatible(Genes *g, int a, int b);
EList *incompatible_sites(Genes *g, int a, int b);
int subsumed_sequence(Genes *g, int a, int b);
int subsumed_site(Sites *s, int a, int b);
EList *ancestral_sites(Genes *g, int a);
EList *zero_sequences(Sites *s, int i);
EList *one_sequences(Sites *s, int i);
int find_safe_coalescence(Genes *g, int a);
int entangled(Genes *g, int a, int b);
void coalesce(Genes *g, int a, int b);
EList *force_coalesce(Genes *g, EList *event);
int coalescence_amleft(Genes *g, int a, int b);
void split(Genes *g, int a, int i);
EList *force_split(Genes *g, int a, EList *event);
void split_coalesceprefix(Genes *g, int a, int index, int block, int b);
void split_coalescepostfix(Genes *g, int a, int index, int block, int b);
void splitafter_coalescepostfix(Genes *g, int a, int index, int block, int b);
void split_removeprefix(Genes *g, int a, int index, int block);
void split_removepostfix(Genes *g, int a, int index, int block);
int individual_ancestral_material(Genes *g, int i);
int ancestral_material(Genes *g);
int individual_all_ancestral(Genes *g, int i);
int all_ancestral(Genes *g);
int informative_ancestral_material(Genes *g);
int ancestral_material_overlap(Genes *g);
int pairwise_am_overlap(Genes *g, int a, int b);
int minimum_compatiblechops(Genes *g, int a);
int minimum_subsumedchops(Genes *g, int a);
Index *maximumsubsumedprefixs(Genes *g);
Index *maximumsubsumedpostfixs(Genes *g);
Index *maximumsubsumedprefix(Genes *g, int s);
Index *maximumsubsumedpostfix(Genes *g, int s);
EList *maximal_prefix_coalesces(Genes *g, Index *a, Index *b);
void maximal_prefix_coalesces_map(Genes *g, Index *a, Index *b,
				    void (*f)(Genes *));
EList *maximal_postfix_coalesces(Genes *g, Index *a, Index *b);
void maximal_postfix_coalesces_map(Genes *g, Index *a, Index *b,
				     void (*f)(Genes *));
void seqerror_flips(Genes* g, void (*f)(Genes *));
void recmut_flips(Genes* g, void (*f)(Genes *));
EList *maximal_infix_coalesces(Genes *g, Index *a, Index *b);
void maximal_infix_coalesces_map(Genes *g, Index *a, Index *b,
				  void (*f)(Genes *));
EList *maximal_overlap_coalesces(Genes *g, Index *a, Index *b);
void maximal_overlap_coalesces_map(Genes *g, Index *a, Index *b,
				     void (*f)(Genes *));
int compare_sequences(Genes *g, int a, int b);
int compare_sites(Sites *s, int a, int b);
int compare_genes(Genes *g, Genes *h);
PackedGenes *pack_genes(Genes *g);
void free_packedgenes(PackedGenes *p);
Genes *unpack_genes(PackedGenes *p);
int compare_packedgenes(PackedGenes *g, PackedGenes *h);
HashTable *new_geneshashtable(int bits);
void init_geneshashtable(HashTable *table, int bits);
HashTable *new_packedgeneshashtable(int bits);
void init_packedgeneshashtable(HashTable *table, int bits);
#endif
