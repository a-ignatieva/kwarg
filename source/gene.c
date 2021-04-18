/***************************************************************************
*
*    gene.c: Implementation of functions to read, print, and mutate an
*    SNP data set.
*
*    Implementation of functions for building and outputting an
*    ancestral recombination graph
*
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>

#include "gene.h"
#include "llist.h"
#include "common.h"
#include "bitfunctions.h"
#include "hashtable.h"
#include "mergesort.h"
#include "backtrack.h"

/* Private type declarations */
typedef enum { _GENE_SAME = 0, _GENE_OPPOSITE = 1, _GENE_MATCH = 2, _GENE_LEFT = 4, _GENE_RIGHT = 8 } _Gene_TwinTypeComponents;
typedef enum { _GENE_INCOMPATIBLE = 0, _GENE_IDENTICAL = _GENE_SAME | _GENE_MATCH, _GENE_CONJUGATE = _GENE_OPPOSITE | _GENE_MATCH, _GENE_LEFTIDENTICAL = _GENE_SAME | _GENE_LEFT, _GENE_LEFTCONJUGATE = _GENE_OPPOSITE | _GENE_LEFT, _GENE_RIGHTIDENTICAL = _GENE_SAME | _GENE_RIGHT, _GENE_RIGHTCONJUGATE = _GENE_OPPOSITE | _GENE_RIGHT } _Gene_TwinType;
typedef struct __SiameseBlock {
    int start;  /* Start of block of Siamese sites */
    int end;    /* End of block of Siamese sites */
    int master; /* Master site of block */
} _SiameseBlock;

/* Global variable determining whether the most recent common ancestor
 * is assumed known.
 */
int gene_knownancestor = 0;

/* Free memory used by a genes data structure */
void free_genes(Genes *g)
{
    int i;

    for (i = 0; i < g->n; i++) {
        free(g->data[i].type);
        free(g->data[i].ancestral);
    }
    free(g->data);
    free(g);
}

/* Create a Genes data structure containing no sequences. This should
 * be used with care as some (most) functions assume non-empty data
 * set. One exception is add_gene, so used in combination with this
 * function it can be used to build a data set.
 */
Genes *make_genes()
{
    Genes *g = (Genes *)xmalloc(sizeof(Genes));
    g->data = (Gene *)xmalloc(sizeof(Gene));
    g->data[0].type = (unsigned long *)xcalloc(1, sizeof(unsigned long));
    g->data[0].ancestral = (unsigned long *)xcalloc(1, sizeof(unsigned long));
    g->n = 0;
    g->length = 0;

    return g;
}

/* Free memory used by an annotatedgenes data structure */
void free_annotatedgenes(AnnotatedGenes *g)
{
    if (g->positions != NULL) {
        while (Length(g->positions) > 0)
            free(Pop(g->positions));
        DestroyLList(g->positions);
    }
    if (g->sequences != NULL) {
        while (Length(g->sequences) > 0)
            free(Pop(g->sequences));
        DestroyLList(g->sequences);
    }
    free_genes(g->g);
    free(g);
}

/* Free memory used by a sites data structure */
void free_sites(Sites *s)
{
    int i;

    for (i = 0; i < s->length; i++) {
        free(s->data[i].type);
        free(s->data[i].ancestral);
    }
    free(s->data);
    free(s);
}

/* Create a Sites data structure containing no sites. This should be
 * used with care as some (most) functions assume non-empty data
 * set. One exception is add_site, so used in combination with this
 * function it can be used to build a data set.
 */
Sites *make_sites()
{
    Sites *s = (Sites *)xmalloc(sizeof(Sites));
    s->data = (Site *)xmalloc(sizeof(Site));
    s->data[0].type = (unsigned long *)xcalloc(1, sizeof(unsigned long));
    s->data[0].ancestral = (unsigned long *)xcalloc(1, sizeof(unsigned long));
    s->n = 0;
    s->length = 0;

    return s;
}

/* Free memory used by a gene data structure */
void free_gene(Gene *g)
{
    free(g->type);
    free(g->ancestral);
}

/* Free memory used by a site data structure */
void free_site(Site *s)
{
    free(s->type);
    free(s->ancestral);
}

/* Initialise an already allocated gene structure g requiring n sites,
 * where the type is given by s.
 */
static void init_gene(Gene *g, int n, char *s, char *anc, Gene_SeqType t)
{
    int i, j;
    unsigned long index;
    char c, a;

    g->type = (unsigned long *)xcalloc(divblocksize(n - 1) + 1,
                                       sizeof(unsigned long));
    g->ancestral = (unsigned long*)xcalloc(divblocksize(n - 1) + 1,
                                           sizeof(unsigned long));

    index = 1;
    j = 0;
    /* Initialise the sites */
    for (i = 0; i < n; i++) {
        switch(t) {
        case GENE_BINARY:
            if (s[i] == '1')
                /* Type in this position is 1 */
                g->type[j] |= index;
            else if (s[i] != '0')
                /* This site is not resolved so shouldn't be considered ancestral */
                g->ancestral[j] |= index;
            break;
        case GENE_NUCLEIC:
            c = tolower(s[i]);
            a = tolower(anc[i]);
            switch(a) {
            case 'a':
            case 'c':
            case 'g':
            case 't':
            case 'u':
                /* Ancestral sequence is resolved for this site */
                switch(c) {
                case 'a':
                case 'c':
                case 'g':
                case 't':
                case 'u':
                    if ((c != a) && ((c != 't') || (a != 'u')) && ((c != 'u') || (a != 't')))
                        g->type[j] |= index;
                    break;
                default:
                    /* If symbol is not a fully resolved nucleic acid, leave it
                     * unresolved.
                     */
                    g->ancestral[j] |= index;
                }
                break;
            default:
                /* Ancestral sequence as yet undetermined in this site */
                switch(c) {
                case 'a':
                case 'c':
                case 'g':
                case 't':
                case 'u':
                    /* Take ancestral base from this sequence */
                    anc[i] = c;
                    break;
                default:
                    /* Sequence unresolved in this site */
                    g->ancestral[j] |= index;
                }
            }
            break;
        case GENE_AMINO:
            c = tolower(s[i]);
            a = tolower(anc[i]);
            switch(a) {
            case 'a':
            case 'c':
            case 'd':
            case 'e':
            case 'f':
            case 'g':
            case 'h':
            case 'i':
            case 'k':
            case 'l':
            case 'm':
            case 'n':
            case 'p':
            case 'q':
            case 'r':
            case 's':
            case 't':
            case 'v':
            case 'w':
            case 'y':
                /* Ancestral sequence is resolved for this site */
                switch(c) {
                case 'a':
                case 'c':
                case 'd':
                case 'e':
                case 'f':
                case 'g':
                case 'h':
                case 'i':
                case 'k':
                case 'l':
                case 'm':
                case 'n':
                case 'p':
                case 'q':
                case 'r':
                case 's':
                case 't':
                case 'v':
                case 'w':
                case 'y':
                    if (c != a)
                        g->type[j] |= index;
                    break;
                default:
                    /* Symbol is not a resolved amino acid */
                    g->ancestral[j] |= index;
                }
                break;
            default:
                switch(c) {
                case 'a':
                case 'c':
                case 'd':
                case 'e':
                case 'f':
                case 'g':
                case 'h':
                case 'i':
                case 'k':
                case 'l':
                case 'm':
                case 'n':
                case 'p':
                case 'q':
                case 'r':
                case 's':
                case 't':
                case 'v':
                case 'w':
                case 'y':
                    /* Take ancestral amino acid from this sequence */
                    anc[i] = c;
                    break;
                default:
                    /* Sequence is not resolved in this site */
                    g->ancestral[j] |= index;
                }
            }
        }
        index <<= 1;
        if (index == 0) {
            /* The ancestral information created is complement of the actual
             * ancestral inforamtion.
             */
            g->ancestral[j] = ~g->ancestral[j];
            /* Reached end of block - increment block number and reset */
            index = 1;
            j++;
        }
    }
    if (index != 1)
        /* Convert ancestral information in last block to its complement */
        g->ancestral[j] = ~g->ancestral[j] & (index - 1);
}

/* Allocate and initialise new set of genes structure with length
 * sites and n genes of types given by s.
 */
Genes *new_genes(int n, int length, char **s, Gene_SeqType t)
{
    int i;
    Genes *g = (Genes *)xmalloc(sizeof(Genes));

    /* Set up basic structure */
    g->n = n;
    g->length = length;
    g->data = (Gene *)xmalloc(n * sizeof(Gene));

    /* Copy genes to structure */
    for (i = 0; i < n; i++)
        init_gene(g->data + i, g->length, s[i], s[0], t);

    return g;
}

/* Read data file into string */
static char *read_file(FILE *fp)
{
    char *data = (char *)xmalloc(sizeof(char));
    long int i = 0, n = 1;

    /* Read file, character by character */
    while ((data[i] = fgetc(fp)) != EOF)
        if (++i == n) {
            /* Running out of space - double it */
            data = xrealloc(data, 2 * n * sizeof(char));
            n *= 2;
        }

    /* Terminate string */
    data[i] = '\0';
    data = xrealloc(data, (i + 1) * sizeof(char));

    return data;
}

/* Return next token from s. The `\' character escapes the following symbol */
typedef enum { CHAR, ESCCHAR } Gene_Token_Type;
typedef struct _Gene_Token {
    Gene_Token_Type type;
    char value;
} Gene_Token;

static Gene_Token *get_token(char **s)
{
    Gene_Token *token;

    if (**s == '\0')
        /* Reached end of the string */
        return NULL;

    if (**s == '\\') {
        /* Escaped character */
        (*s)++;
        if (**s == '\0')
            /* Reached end of the string (in a rather weird way) */
            return NULL;
        else if (**s == '\n') {
            /* Escape character glues two lines together */
            (*s)++;
            return get_token(s);
        }
        token = (Gene_Token *)xmalloc(sizeof(Gene_Token));
        token->type = ESCCHAR;
        switch(**s) {
        case 'n':
            token->value = '\n';
            break;
        case 'b':
            token->value = '\b';
            break;
        case 't':
            token->value = '\t';
            break;
        default:
            token->value = **s;
        }
    }
    else {
        /* Normal character */
        token = (Gene_Token *)xmalloc(sizeof(Gene_Token));
        token->type = CHAR;
        token->value = **s;
    }

    (*s)++;
    return token;
}

/* Return next line from s as an llist of tokens */
static LList *get_line(char **s)
{
    Gene_Token *t;
    LList *line = MakeLList();

    for (;;) {
        t = get_token(s);
        if (t == NULL) {
            /* Reached end of file */
            if (Length(line) == 0) {
                DestroyLList(line);
                return NULL;
            }
            else
                return line;
        }
        else if ((t->type == CHAR) && (t->value == '\n')) {
            /* Reached end of line */
            free(t);
            if (Length(line) > 0)
                return line;
        }
        else
            Enqueue(line, t);
    }
}

/* Convert an llist of tokens to a string, removing white space and
 * comment. If the string is empty, NULL is returned.
 */
static char *get_sequence(LList *s)
{
    char *t = NULL;
    int i, n;
    LListCounter *c = MakeCounter(s, FIRST);
    Gene_Token *token = (Gene_Token *)Next(c);

    while (token != NULL) {
        if (token->type == CHAR) {
            /* Next character was not escaped */
            if (token->value == '#') {
                /* Remainder of line is comment */
                n = Length(s);
                for (i = GetPosition(c); i < n; i++)
                    free(Dequeue(s));
                break;
            }
            else if (isspace(token->value)) {
                /* Ignore white space */
                free(RemoveMoveRight(s, c));
                token = (Gene_Token *)Current(c);
                continue;
            }
        }
        token = (Gene_Token *)Next(c);
    }

    /* Transfer characters to string */
    n = Length(s);
    if (n > 0) {
        t = (char *)xmalloc((n + 1) * sizeof(char));
        for (i = 0; i < n; i++) {
            token = (Gene_Token *)Pop(s);
            t[i] = token->value;
            free(token);
        }
        t[n] = '\0';
    }

    /* Clean up */
    DestroyLList(s);
    DestroyCounter(c);

    return t;
}

/* Interpret s as a white space separated list of labels, and return
 * an llist containing these labels as strings.
 */
static LList *get_labels(LList *s)
{
    char *t;
    int i, n;
    LListCounter *c = MakeCounter(s, FIRST);
    Gene_Token *token = (Gene_Token *)Next(c);
    LList *labels = MakeLList();

    while (token != NULL) {
        /* Remove stretch of white spaces */
        do {
            if ((!isspace(token->value)) || (token->type == ESCCHAR))
                /* Next character is not white space */
                break;
            free(RemoveMoveRight(s, c));
            token = (Gene_Token *)Current(c);
        } while (token != NULL);
        if (token == NULL)
            /* Reached end of line */
            break;

        /* Find stretch of non-white space characters */
        do {
            if ((isspace(token->value)) && (token->type == CHAR))
                /* Next token is white space */
                break;
            token = (Gene_Token *)Next(c);
        } while (token != NULL);
        /* Transfer this stretch to a string */
        n = GetPosition(c);
        if (n < 0)
            n = Length(s);
        t = (char *)xmalloc((n + 1) * sizeof(char));
        for (i = 0; i < n; i++) {
            token = (Gene_Token *)Pop(s);
            t[i] = token->value;
            free(token);
        }
        t[n] = '\0';
        Enqueue(labels, t);
        InitCounter(c, s, FIRST);
        token = (Gene_Token *)Next(c);
    }

    /* Clean up */
    DestroyLList(s);
    DestroyCounter(c);

    return labels;
}

/* Free the llists of g, c, and the n strings in sequences */
static void abort_parse(AnnotatedGenes *g, LListCounter *c, char **sequences,
                        int n, char *t)
{
    int i;

    while (Length(g->positions) > 0)
        free(Pop(g->positions));
    DestroyLList(g->positions);
    while (Length(g->sequences)> 0)
        free(Pop(g->sequences));
    DestroyLList(g->sequences);
    if (c != NULL)
        DestroyCounter(c);
    if (n > 0) {
        for (i = 0; i < n; i++)
            free(sequences[i]);
        free(sequences);
    }
    if (t != NULL)
        free(t);
    free(g);
}

/* Try to parse string s as a data set in original beagle format. If
 * this is not possible, NULL is returned.
 */
static AnnotatedGenes *parse_genes_beagle(char *s, Gene_SeqType seqtype)
{
    LList *line = get_line(&s);
    LListCounter *c;
    AnnotatedGenes *g = (AnnotatedGenes *)xmalloc(sizeof(AnnotatedGenes));
    Gene_Token *token;
    int tagcount = 2, i, tag, l, m, n = 0;
    static char tags[][11] = { ">", "positions:" };
    char *t = NULL, *sequence, **sequences;

    if (line == NULL) {
        free(g);
        return NULL;
    }
    c = MakeCounter(line, FIRST);
    g->positions = MakeLList();
    g->sequences = MakeLList();
    while (line != NULL) {
        token = (Gene_Token *)Top(line);
        if ((token->value == '#') && (token->type == CHAR)) {
            /* Entire line is comment */
            free(Pop(line));
            /* Remove initial white spaces */
            while (Length(line) > 0) {
                token = (Gene_Token *)Top(line);
                if (!isspace(token->value) || (token->type == ESCCHAR))
                    break;
                free(Pop(line));
            }
            /* Check whether it contains a sequence label or a list of site labels */
            for (tag = 0; tag < tagcount; tag++) {
                InitCounter(c, line, FIRST);
                for (i = 0; i < strlen(tags[tag]); i++) {
                    token = (Gene_Token *)Next(c);
                    if ((tolower(token->value) != tags[tag][i])
                            || (token->type == ESCCHAR))
                        break;
                }
                if (i >= strlen(tags[tag]))
                    /* Matched tag */
                    break;
            }
            switch(tag) {
            case 0:
                /* Matched sequence label tag */
                for (i = 0; i < strlen(tags[0]); i++)
                    free(Pop(line));
                /* Remove initial and trailing white space stretches */
                while (Length(line) > 0) {
                    token = (Gene_Token *)Top(line);
                    if (!isspace(token->value) || (token->type == ESCCHAR))
                        break;
                    free(Pop(line));
                }
                while (Length(line) > 0) {
                    token = (Gene_Token *)Bottom(line);
                    if (!isspace(token->value) || (token->type == ESCCHAR))
                        break;
                    free(Dequeue(line));
                }
                /* Transfer remaining symbols to string */
                if (t != NULL)
                    /* Free old label */
                    free(t);
                t = (char *)xmalloc((Length(line) + 1) * sizeof(char));
                for (i = 0; Length(line) > 0; i++) {
                    token = Pop(line);
                    t[i] = token->value;
                    free(token);
                }
                t[i] = '\0';
                DestroyLList(line);
                break;
            case 1:
                /* Matched site labels tag */
                for (i = 0; i < strlen(tags[1]); i++)
                    free(Pop(line));
                /* Get site labels */
                Append(g->positions, get_labels(line));
                break;
            default:
                /* No tags - ignore entire line */
                while (Length(line) > 0)
                    free(Pop(line));
                DestroyLList(line);
            }
        }
        else {
            /* Line is not all comment - must be a sequence */
            sequence = get_sequence(line);
            if (sequence != NULL) {
                /* And it was non-empty */
                if (n == 0) {
                    /* First sequence encountered */
                    sequences = (char **)xmalloc(sizeof(char *));
                    sequences[0] = sequence;
                    m = n = 1;
                    l = strlen(sequence);
                }
                else {
                    if (l != strlen(sequence)) {
                        /* Sequences are not of identical length - abort */
                        abort_parse(g, c, sequences, n, t);
                        free(sequence);
                        return NULL;
                    }
                    if (m == n) {
                        /* We are running out of space - double it */
                        m *= 2;
                        sequences = (char **)xrealloc(sequences, m * sizeof(char *));
                    }
                    sequences[n++] = sequence;
                }
                if (t != NULL) {
                    /* This sequence had a label */
                    for (i = Length(g->sequences) + 1; i < n; i++)
                        /* Make sure all previous sequences have an associated
                         * label by assigning unlabeled sequences the NULL label.
                         */
                        Enqueue(g->sequences, NULL);
                    if (strlen(t) == 0) {
                        /* Empty sequence label */
                        free(t);
                        Enqueue(g->sequences, NULL);
                    }
                    else
                        Enqueue(g->sequences, t);
                    t = NULL;
                }
            }
        }
        line = get_line(&s);
    }

    /* Check for consistency, etc. */
    if ((n == 0) || ((Length(g->positions) > 0) && (Length(g->positions) != l))) {
        /* either no sequences, or number of site labels does not match
         * number of sites.
         */
        abort_parse(g, c, sequences, n, t);
        return NULL;
    }
    if (Length(g->sequences) > 0)
        /* Make sure all sequences have an associated label by assigning
         * unlabeled sequences the NULL label.
         */
        for (i = Length(g->sequences); i < n; i++)
            Enqueue(g->sequences, NULL);
    else {
        DestroyLList(g->sequences);
        g->sequences = NULL;
    }
    if (Length(g->positions) == 0) {
        DestroyLList(g->positions);
        g->positions = NULL;
    }

    g->g = new_genes(n, l, sequences, seqtype);

    /* Clean up */
    DestroyCounter(c);
    for (i = 0; i < n; i++)
        free(sequences[i]);
    free(sequences);
    if (t != NULL)
        free(t);

    return g;
}

/* Add the next l non-white space characters from s to the array of
 * strings in sequences.
 */
static void add_sequence(char *s, int *l, int *n, int *m, char ***sequences,
                         AnnotatedGenes *g, char *sequence)
{
    int i = 0;

    if (*n == 0) {
        /* First sequence to be added */
        *sequences = (char **)xmalloc(sizeof(char *));
        *m = 1;
    }
    else if (*n == *m) {
        /* We are running out of space - double it */
        (*m) *= 2;
        *sequences = (char **)xrealloc(*sequences, *m * sizeof(char *));
    }

    /* Copy sequence to sequences */
    (*sequences)[*n] = (char *)xmalloc((*l + 1) * sizeof(char));
    while (i < *l) {
        if (!isspace(*s))
            (*sequences)[*n][i++] = *s;
        s++;
    }
    (*sequences)[*n][i] = '\0';
    (*n)++;

    /* Insert sequence label in llist of sequence labels */
    Enqueue(g->sequences, sequence);

    /* Reset character count */
    *l = 0;
}

/* Try to parse string s as a data set in fasta format. If this is not
 * possible, NULL is returned.
 */
static AnnotatedGenes *parse_genes_fasta(char *s, Gene_SeqType seqtype)
{
    int tagcount = 1, i, j, tag, k = 0, l = 0, m, n = 0;
    AnnotatedGenes *g = (AnnotatedGenes *)xmalloc(sizeof(AnnotatedGenes));
    static char tags[][11] = { "positions:" };
    char *sequence = NULL, *label, **sequences, *t;

    /* First line should be first sequence label */
    if (*s != '>') {
        free(g);
        return NULL;
    }

    g->sequences = MakeLList();
    g->positions = MakeLList();
    while (*s != '\0') {
        switch (*s) {
        case '>':
            /* Beginning of new sequence */
            if (sequence != NULL) {
                /* Store previous sequence */
                if (l == 0)
                    l = k;
                if ((k == 0) || (k != l)) {
                    /* No sequence for last sequence label or sequences not of
                     * equal length.
                     */
                    abort_parse(g, NULL, sequences, n, sequence);
                    return NULL;
                }
                add_sequence(t, &k, &n, &m, &sequences, g, sequence);
            }
            /* Start by scanning through leading > and white space symbols */
            while (*s == '>')
                s++;
            while ((isspace(*s)) && (*s != '\n'))
                s++;
            /* Determine sequence label */
            if (*s == '\0') {
                /* End of file - last sequence is empty */
                abort_parse(g, NULL, sequences, n, sequence);
                return NULL;
            }
            if (*s == '\n')
                /* No sequence name */
                sequence = NULL;
            else {
                for (i = 0; (s[i] != '\n') && (s[i] != '\0'); i++);
                /* Removing trailing stretch of white space */
                for (j = i - 1; isspace(s[j]); j--);
                sequence = (char *)xcalloc(j + 2, sizeof(char));
                strncpy(sequence, s, j + 1);
                s += i;
            }
            /* Record (possible) start of sequence data */
            t = s + 1;
            break;
        case ';':
            /* Comment line */
            if (k > 0) {
                /* Comments should precede sequences */
                abort_parse(g, NULL, sequences, n, sequence);
                return NULL;
            }
            /* Remove initial stretch of ; characters */
            while (*s == ';')
                s++;
            /* Remove white spaces */
            while (isspace(*s) && (*s != '\n'))
                s++;
            if (*s == '\n') {
                /* Empty comment line */
                /* Record (possible) start of sequence data */
                t = s + 1;
                break;
            }
            /* Check whether comment line contains special tag */
            for (tag = 0; tag < tagcount; tag++) {
                for (i = 0; i < strlen(tags[tag]); i++)
                    if (tolower(s[i]) != tags[tag][i])
                        break;
                if (i >= strlen(tags[tag])) {
                    /* Matched tag */
                    s += i;
                    break;
                }
            }
            switch(tag) {
            case 0:
                /* Matched site labels tag */
                while ((*s != '\n') && (*s != '\0')) {
                    /* Skip stretches of white space */
                    while (isspace(*s) && (*s != '\n'))
                        s++;
                    if ((*s == '\n') || (*s == '\0'))
                        break;
                    /* Determine next label */
                    for (i = 0; !isspace(s[i]) && (s[i] != '\0'); i++);
                    label = (char *)xcalloc(i + 1, sizeof(char));
                    strncpy(label, s, i);
                    Enqueue(g->positions, label);
                    s += i;
                }
                break;
            default:
                /* Ordinary comment line - skip */
                while ((*s != '\n') && (*s != '\0'))
                    s++;
            }
            /* Record (possible) start of sequence data */
            t = s + 1;
            break;
        default:
            /* Line containts (part of) the sequence - count number of symbols */
            while ((*s != '\n') && (*s != '\0')) {
                s++;
                k++;
            }
        }
        /* End of line - skip stretches of new lines */
        while (*s == '\n')
            s++;
    }

    /* Store last sequence */
    if ((k == 0) || ((n > 0) && (l != k))) {
        /* No sequence for last sequence label or sequences not of
         * equal length.
         */
        abort_parse(g, NULL, sequences, n, sequence);
        return NULL;
    }
    add_sequence(t, &k, &n, &m, &sequences, g, sequence);

    /* Check for consistency */
    if ((Length(g->positions) > 0) && (Length(g->positions) != l)) {
        /* Number of site labels does not match number of sites */
        abort_parse(g, NULL, sequences, n, NULL);
        return NULL;
    }
    if (Length(g->positions) == 0) {
        DestroyLList(g->positions);
        g->positions = NULL;
    }

    /* Transfer sequences to a Genes structure */
    g->g = new_genes(n, l, sequences, seqtype);

    /* Clean up */
    for (i = 0; i < n; i++)
        free(sequences[i]);
    free(sequences);

    return g;
}

/* Read a set of genes from file given by fname (stdin if fname is NULL) */
AnnotatedGenes *read_genes(char *fname, Gene_Format f, Gene_SeqType t)
{
    char *s;
    FILE *fp;
    AnnotatedGenes *g = NULL;

    /* Set up file for input */
    if (fname == NULL)
        fp = stdin;
    else {
        if ((fp = fopen(fname, "r")) == NULL) {
            fprintf(stderr, "Cannot open file %s\n", fname);
            exit(1);
        }
    }

    /* Read and store file */
    s = read_file(fp);
    if ((f == GENE_ANY) || (f == GENE_BEAGLE))
        g = parse_genes_beagle(s, t);
    if ((g == NULL) && ((f == GENE_ANY) || (f == GENE_FASTA)))
        g = parse_genes_fasta(s, t);

    /* Clean up */
    if (fname != NULL)
        fclose(fp);
    free(s);

    /* Did we succeed? */
    if (g == NULL) {
        /* No! */
        fprintf(stderr, "Could not read data from %s", (fname != NULL ? fname :
                "standard input"));
        switch(f) {
        case GENE_ANY:
            fprintf(stderr, ".\n");
            break;
        case GENE_BEAGLE:
            fprintf(stderr, " in original beagle format.\n");
            break;
        case GENE_FASTA:
            fprintf(stderr, " in fasta format.\n");
            break;
        }
        exit(1);
    }

    return g;
}

/* Output gene data in g to file fp, prefixing each sequence with its
 * label from labels if labels is not NULL. Return value is the
 * indentation of sequence data due to labels.
 */
static int _output_genes(Genes *g, FILE *fp, LList *labels)
{
    int i, j, l, block;
    unsigned long index;
    LListCounter *c;
    char *label;

    if (labels != NULL) {
        /* Find length of longest sequence label */
        l = 0;
        c = MakeCounter(labels, FIRST);
        while ((label = (char *)Next(c)) != NULL) {
            i = strlen(label);
            l = (i > l ? i : l);
        }
        InitCounter(c, labels, FIRST);
    }

    /* Output genes */
    for (i = 0; i < g->n; i++) {
        if (labels != NULL) {
            /* First output label */
            label = Next(c);
            if (fprintf(fp, "%s", label) < 0) {
                fprintf(stderr, "Error outputting data\n");
                exit(1);
            }
            /* Output alignment spaces */
            for (j = strlen(label); j <= l; j++)
                if (fputc((int)' ', fp) == EOF) {
                    fprintf(stderr, "Error outputting data\n");
                    exit(1);
                }
        }
        index = 1;
        block = 0;
        for (j = 0; j < g->length; j++) {
            if (g->data[i].ancestral[block] & index) {
                if (g->data[i].type[block] & index) {
                    /* Gene i has type 1 at site j */
                    if (fputc((int)'1', fp) == EOF) {
                        fprintf(stderr, "Error outputting data\n");
                        exit(1);
                    }
                }
                else
                    /* Gene i has type 0 at site j */
                    if (fputc((int)'0', fp) == EOF) {
                        fprintf(stderr, "Error outputting data\n");
                        exit(1);
                    }
            }
            else {
                /* Gene i does not have ancestral data at site j */
                if (fputc((int)'x', fp) == EOF) {
                    fprintf(stderr, "Error outputting data\n");
                    exit(1);
                }
            }
            index <<= 1;
            if (index == 0) {
                index = 1;
                block++;
            }
        }
        /* Finished outputting gene i */
        if (fputc((int)'\n', fp) == EOF) {
            fprintf(stderr, "Error outputting data\n");
            exit(1);
        }
    }

    if (labels != NULL) {
        /* Clean up */
        DestroyCounter(c);
        return l + 1;
    }

    return 0;
}

/* Output gene data in g to file fp (stdout if fp is NULL), prefixed
 * with comment if not NULL.
 */
void output_genes(Genes *g, FILE *fp, char *comment)
{
    /* Set up file for output */
    if (fp == NULL)
        fp = stdout;

    /* Output genes */
    if (comment != NULL)
        fprintf(fp, "%s", comment);
    _output_genes(g, fp, NULL);
}

/* Output gene data in g to file fp (stdout if fp is NULL), prefixed
 * with comment if not NULL.
 */
void output_labelled_genes(Genes *g, FILE *fp, LList *labels)
{
    /* Set up file for output */
    if (fp == NULL)
        fp = stdout;
    
    _output_genes(g, fp, labels);
}

/* Output gene data in g to file fp (stdout if fp is NULL), followed
 * by a line of site indeces.
 */
void output_genes_indexed(Genes *g, FILE *fp)
{
    int i;

    /* Set up file for output */
    if (fp == NULL)
        fp = stdout;

    output_genes(g, fp, NULL);
    for (i = 0; i < g->length; i++)
        fprintf(fp, "%d", i % 10);
    fprintf(fp, "\n");
}

/* Output gene data in a to file fp (stdout if fp is NULL), prefixed
 * with comment if not NULL.
 */
void output_annotatedgenes(AnnotatedGenes *a, FILE *fp, char *comment)
{
    int i, j, l, m;
    char *label;
    LListCounter *c;

    /* Set up file for output */
    if (fp == NULL)
        fp = stdout;

    /* Output genes */
    if (comment != NULL)
        fprintf(fp, "%s", comment);
    l = _output_genes(a->g, fp, a->sequences);

    /* Output labels */
    if (a->positions != NULL) {
        c = MakeCounter(a->positions, FIRST);
        while ((label = (char *)Next(c)) != NULL)
            if (strlen(label) > 0)
                break;
        if (label != NULL) {
            /* At least one label was not the empty string */
            i = 0;
            m = 1;
            while (m) {
                m = 0;
                /* Start by emitting as many spaces as the sequence labels took up */
                for (j = 0; j < l; j++)
                    if (fputc((int)' ', fp) == EOF) {
                        fprintf(stderr, "Error outputting data\n");
                        exit(1);
                    }
                /* Now emit next symbol of the sequence labels */
                InitCounter(c, a->positions, FIRST);
                while ((label = (char *)Next(c)) != NULL) {
                    j = strlen(label);
                    if (i < j) {
                        if (i < j - 1)
                            /* This label will require at least one more iteration */
                            m = 1;
                        j = label[i];
                    }
                    else
                        /* This label does not contain any more symbols */
                        j = ' ';
                    if (fputc(j, fp) == EOF) {
                        fprintf(stderr, "Error outputting data\n");
                        exit(1);
                    }
                }
                if (fputc('\n', fp) == EOF) {
                    fprintf(stderr, "Error outputting data\n");
                    exit(1);
                }
                i++;
            }
        }
        DestroyCounter(c);
    }
}

/* Add the gene specified by the gene new to g. The new gene is
 * assumed to have the same length as the genes already present in
 * g. If g does not contain any genes, an extra argument specifying
 * the number of sites in new is required.
 */
void add_gene(Genes *g, Gene *new, ...)
{
    int blocks;
    va_list args;

    /* Determine length of sequences and corresponding number of blocks */
    if (g->n == 0) {
        va_start(args, new);
        g->length = (int)va_arg(args, int);
        va_end(args);
    }
    blocks = divblocksize(g->length - 1) + 1;

    /* Modify basic structure */
    g->n++;
    g->data = (Gene *)xrealloc(g->data, g->n * sizeof(Gene));

    /* Copy new gene to structure */
    g->data[g->n - 1].type
        = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(g->data[g->n - 1].type, new->type,
           blocks * sizeof(unsigned long));
    g->data[g->n - 1].ancestral
        = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(g->data[g->n - 1].ancestral, new->ancestral,
           blocks * sizeof(unsigned long));
}

/* Add the site specified by the site new to s. The new site is
 * assumed to have the same length as the sites already present in
 * s. If s does not contain any sites, an extra argument specifying
 * the number of genes in new is required.
 */
void add_site(Sites *s, Site *new, ...)
{
    int blocks;

    va_list args;
    /* Determine number of sequences and corresponding number of blocks */
    if (s->length == 0) {
        va_start(args, new);
        s->n = (int)va_arg(args, int);
        va_end(args);
    }
    blocks = divblocksize(s->n - 1) + 1;

    /* Modify basic structure */
    s->length++;
    s->data = (Site *)xrealloc(s->data, s->length * sizeof(Site));

    /* Copy new site to structure */
    s->data[s->length - 1].type
        = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(s->data[s->length - 1].type, new->type,
           blocks * sizeof(unsigned long));
    s->data[s->length - 1].ancestral
        = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(s->data[s->length - 1].ancestral, new->ancestral,
           blocks * sizeof(unsigned long));
}

/* Return the i'th sequence of g as a Gene data structure */
Gene *get_gene(Genes *g, int i)
{
    int blocks = divblocksize(g->length - 1) + 1;
    Gene *new = (Gene *)xmalloc(sizeof(Gene));

    new->type = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(new->type, g->data[i].type, blocks * sizeof(unsigned long));
    new->ancestral = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(new->ancestral, g->data[i].ancestral, blocks * sizeof(unsigned long));

    return new;
}

/* Return the i'th site of s as a Site data structure */
Site *get_site(Sites *s, int i)
{
    int blocks = divblocksize(s->n - 1) + 1;
    Site *new = (Site *)xmalloc(sizeof(Site));

    new->type = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(new->type, s->data[i].type, blocks * sizeof(unsigned long));
    new->ancestral = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    memcpy(new->ancestral, s->data[i].ancestral, blocks * sizeof(unsigned long));

    return new;
}

/* Construct a Genes data set that contains the subset of sequences in
 * g specified by l. All indeces in l are assumed to be less than
 * g->n.
 */
Genes *select_genes(Genes *g, EList *l)
{
    int i, blocks = divblocksize(g->length - 1) + 1;
    Genes *h = xmalloc(sizeof(Genes));
    h->n = elist_length(l);

    if (h->n == 0) {
        /* Empty list of sequences */
        h->length = 0;
        h->data = (Gene *)xmalloc(sizeof(Gene));
        h->data[0].type = (unsigned long *)xcalloc(1, sizeof(unsigned long));
        h->data[0].ancestral = (unsigned long *)xcalloc(1, sizeof(unsigned long));
    }
    else {
        /* Add specified sequences to h */
        h->length = g->length;
        h->data = (Gene *)xmalloc(h->n * sizeof(Gene));
        for (i = 0; i < h->n; i++) {
            h->data[i].type
                = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
            memcpy(h->data[i].type, g->data[(int)elist_get(l, i)].type,
                   blocks * sizeof(unsigned long));
            h->data[i].ancestral
                = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
            memcpy(h->data[i].ancestral, g->data[(int)elist_get(l, i)].ancestral,
                   blocks * sizeof(unsigned long));
        }
    }

    return h;
}

/* Construct a Sites data set that contains the subset of sites in s
 * specified by l. All indeces in l are assumed to be less than
 * s->length.
 */
Sites *select_sites(Sites *s, EList *l)
{
    int i, blocks = divblocksize(s->n - 1) + 1;
    Sites *t = xmalloc(sizeof(Sites));
    t->length = elist_length(l);

    if (t->length == 0) {
        /* Empty list of sequences */
        t->n = 0;
        t->data = (Site *)xmalloc(sizeof(Site));
        t->data[0].type = (unsigned long *)xcalloc(1, sizeof(unsigned long));
        t->data[0].ancestral = (unsigned long *)xcalloc(1, sizeof(unsigned long));
    }
    else {
        /* Add specified sites to t */
        t->n = s->n;
        t->data = (Site *)xmalloc(elist_length(l) * sizeof(Gene));
        for (i = 0; i < t->n; i++) {
            t->data[i].type
                = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
            memcpy(t->data[i].type, s->data[(int)elist_get(l, i)].type,
                   blocks * sizeof(unsigned long));
            t->data[i].ancestral
                = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
            memcpy(t->data[i].ancestral, s->data[(int)elist_get(l, i)].ancestral,
                   blocks * sizeof(unsigned long));
        }
    }

    return t;
}

/* Return character in site site of sequence seq of g - if character
 * is ancestral either 0 or 1 is returned depending on type, if
 * character is non-ancestral 2 is returned. It is assumed that seq is
 * less than g->n and site is less than g->length.
 */
char get_genes_character(Genes *g, int seq, int site)
{
    int block = divblocksize(site);
    int index = modblocksize(site);

    if ((g->data[seq].ancestral[block] & ((unsigned long)1 << index)) != 0)
        /* Character is ancestral */
        return ((g->data[seq].type[block] >> index) & 1);
    else
        /* Character is non-ancestral */
        return 2;
}

/* Return character in site site of sequence seq of g - if character
 * is ancestral either 0 or 1 is returned depending on type, if
 * character is non-ancestral 2 is returned. It is assumed that seq is
 * less than g->n and site is less than g->length.
 */
char get_sites_character(Sites *s, int seq, int site)
{
    int block = divblocksize(seq);
    int index = modblocksize(seq);

    if ((s->data[site].ancestral[block] & ((unsigned long)1 << index)) != 0)
        /* Character is ancestral */
        return ((s->data[site].type[block] >> index) & 1);
    else
        /* Character is non-ancestral */
        return 2;
}

/* Set character in site site of sequence seq of g to c, where 0 and 1
 * denotes ancestral types and any other value of c corresponds to a
 * non-ancestral type. It is assumed that seq is less than g->n and
 * site is less than g->length.
 */
void set_genes_character(Genes *g, int seq, int site, char c)
{
    int block = divblocksize(site);
    int index = modblocksize(site);

    /* Set type for character */
    if (c == 1)
        /* Ancestral type 1 */
        g->data[seq].type[block] |= (unsigned long)1 << index;
    else
        /* Non-ancestral type or ancestral type 0 */
        g->data[seq].type[block] &= ~((unsigned long)1 << index);

    /* Set ancestrality for character */
    if (c <= 1)
        /* Ancestral type */
        g->data[seq].ancestral[block] |= (unsigned long)1 << index;
    else
        /* Non-ancestral type */
        g->data[seq].ancestral[block] &= ~((unsigned long)1 << index);
}

/* Set character in site site of sequence seq of s to c, where 0 and 1
 * denotes ancestral types and any other value of c corresponds to a
 * non-ancestral type. It is assumed that seq is less than s->n and
 * site is less than s->length.
 */
void set_sites_character(Sites *s, int seq, int site, char c)
{
    int block = divblocksize(seq);
    int index = modblocksize(seq);

    /* Set type for character */
    if (c == 1)
        /* Ancestral type 1 */
        s->data[site].type[block] |= (unsigned long)1 << index;
    else
        /* Non-ancestral type or ancestral type 0 */
        s->data[site].type[block] &= ~((unsigned long)1 << index);

    /* Set ancestrality for character */
    if (c <= 1)
        /* Ancestral type */
        s->data[site].ancestral[block] |= (unsigned long)1 << index;
    else
        /* Non-ancestral type */
        s->data[site].ancestral[block] &= ~((unsigned long)1 << index);
}

/* Swap sequences a and b */
void swap_genes(Genes *g, int a, int b)
{
    unsigned long *tmp = g->data[a].type;

    g->data[a].type = g->data[b].type;
    g->data[b].type = tmp;
    tmp = g->data[a].ancestral;
    g->data[a].ancestral = g->data[b].ancestral;
    g->data[b].ancestral = tmp;
}

/* Add the all 0 sequence as the last sequence to g */
void add_ancestral_genes(Genes *g)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    /* Modify basic structure */
    g->n++;
    g->data = (Gene *)xrealloc(g->data, g->n * sizeof(Gene));

    /* Insert the all 0 sequence */
    g->data[g->n - 1].type
        = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    g->data[g->n - 1].ancestral
        = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    for (i = 0; i < blocks - 1; i++) {
        g->data[g->n - 1].type[i] = 0;
        g->data[g->n - 1].ancestral[i] = ~0;
    }
    g->data[g->n - 1].type[i] = 0;
    g->data[g->n - 1].ancestral[i] = ((unsigned long)2 << modblocksize(g->length - 1)) - 1;
}

/* Add the all 0 sequence as the last sequence to s */
void add_ancestral_sites(Sites *s)
{
    int i, blocks = divblocksize(s->n) + 1, index = modblocksize(s->n);

    /* Modify basic structure */
    s->n++;

    /* Insert the all 0 sequence */
    if (index == 0)
        /* We need to extend the bit arrays for each size */
        for (i = 0; i < s->length; i++) {
            s->data[i].type
                = (unsigned long *)xrealloc(s->data[i].type,
                                            blocks * sizeof(unsigned long));
            s->data[i].type[blocks - 1] = 0;
            s->data[i].ancestral
                = (unsigned long *)xrealloc(s->data[i].ancestral,
                                            blocks * sizeof(unsigned long));
            s->data[i].ancestral[blocks - 1] = 1;
        }
    else
        /* Room for new sequence in existing bit arrays */
        for (i = 0; i < s->length; i++)
            s->data[i].ancestral[blocks - 1] |= (unsigned long)1 << index;
}

/* Make a copy of a structure containing a set of genes */
Genes *copy_genes(Genes *g)
{
    Genes *new = (Genes *)xmalloc(sizeof(Genes));
    int i, blocks = divblocksize(g->length - 1) + 1;

    new->length = g->length;
    new->n = g->n;
    new->data = (Gene *)xmalloc(new->n * sizeof(Gene));

    for (i = 0; i < new->n; i++) {
        new->data[i].type
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].type, g->data[i].type, blocks * sizeof(unsigned long));
        new->data[i].ancestral
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].ancestral, g->data[i].ancestral,
               blocks * sizeof(unsigned long));
    }

    return new;
}

/* Make a copy of a structure containing a set of sites */
Sites *copy_sites(Sites *s)
{
    Sites *new = (Sites *)xmalloc(sizeof(Sites));
    int i, blocks = divblocksize(s->n - 1) + 1;

    new->length = s->length;
    new->n = s->n;
    new->data = (Site *)xmalloc(new->length * sizeof(Site));

    for (i = 0; i < new->length; i++) {
        new->data[i].type
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].type, s->data[i].type, blocks * sizeof(unsigned long));
        new->data[i].ancestral
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].ancestral, s->data[i].ancestral,
               blocks * sizeof(unsigned long));
    }

    return new;
}

/* Make a copy of a structure containing a set of genes, except for
 * sequence a.
 */
Genes *copy_allbutone(Genes *g, int a)
{
    Genes *new = (Genes *)xmalloc(sizeof(Genes));
    int i, blocks = divblocksize(g->length - 1) + 1;

    new->length = g->length;
    new->n = g->n - 1;
    new->data = (Gene *)xmalloc(new->n * sizeof(Gene));

    for (i = 0; i < a; i++) {
        new->data[i].type
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].type, g->data[i].type, blocks * sizeof(unsigned long));
        new->data[i].ancestral
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].ancestral, g->data[i].ancestral,
               blocks * sizeof(unsigned long));
    }
    for (; i < new->n; i++) {
        new->data[i].type
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].type, g->data[i + 1].type,
               blocks * sizeof(unsigned long));
        new->data[i].ancestral
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        memcpy(new->data[i].ancestral, g->data[i + 1].ancestral,
               blocks * sizeof(unsigned long));
    }

    return new;
}

/* Make a copy of a structure containing a set of genes in the region
 * starting at site a and ending just before site b.
 */
Genes *copy_region(Genes *g, int a, int b)
{
    Genes *new = (Genes *)xmalloc(sizeof(Genes));
    int i, j, blocks = divblocksize(b - a - 1) + 1, start = divblocksize(a),
              shift = modblocksize(a), gblocks = divblocksize(g->n - 1) + 1;
    unsigned long pattern = ((unsigned long)1 << modblocksize(b - a)) - 1;

    if ((new->length = b - a) < 0) {
        /* Cannot copy a region of negative size */
        free (new);
        return NULL;
    }
    new->n = g->n;
    new->data = (Gene *)xmalloc(new->n * sizeof(Gene));

    for (i = 0; i < new->n; i++) {
        /* Copy sequence i region */
        new->data[i].type
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        new->data[i].ancestral
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        for (j = 0; j < blocks - 1; j++) {
            new->data[i].type[j] = g->data[i].type[j + start] >> shift;
            new->data[i].ancestral[j] = g->data[i].ancestral[j + start] >> shift;
            if (shift != 0) {
                new->data[i].type[j]
                |= g->data[i].type[j + start + 1] << (BLOCKSIZE - shift);
                new->data[i].ancestral[j]
                |= g->data[i].ancestral[j + start + 1] << (BLOCKSIZE - shift);
            }
        }
        /* Last block - handle with care */
        new->data[i].type[j] = g->data[i].type[j + start] >> shift;
        new->data[i].ancestral[j] = g->data[i].ancestral[j + start] >> shift;
        if ((shift != 0) && (j + start < gblocks - 1)) {
            new->data[i].type[j]
            |= g->data[i].type[j + start + 1] << (BLOCKSIZE - shift);
            new->data[i].ancestral[j]
            |= g->data[i].ancestral[j + start + 1] << (BLOCKSIZE - shift);
        }
        if (pattern != 0) {
            new->data[i].type[j] &= pattern;
            new->data[i].ancestral[j] &= pattern;
        }
    }

    return new;
}

/* Remove sequence a from g */
void remove_gene(Genes *g, int a)
{
    if ((a >= 0) && (a < g->n)) {
        g->n--;
        free(g->data[a].type);
        free(g->data[a].ancestral);
        memmove(g->data + a, g->data + a + 1, (g->n - a) * sizeof(Gene));
    }
}

/* Remove sequence a from g */
void remove_annotatedgene(AnnotatedGenes *g, int a)
{
    LListCounter *lcounter;

    if ((a >= 0) && (a < g->g->n)) {
        remove_gene(g->g, a);
        if (g->positions != NULL) {
            lcounter = MakeCounter(g->positions, a);
            RemoveMoveRight(g->positions, lcounter);
            free(lcounter);
        }
        if (g->sequences != NULL) {
            lcounter = MakeCounter(g->sequences, a);
            RemoveMoveRight(g->sequences, lcounter);
            free(lcounter);
        }
    }
}

/* Reallocate the memory used by a gene - some operations may reduce
 * the size of our data set, so by reallocating we can reduce the
 * amount of memory used to store a particular set of genes.
 */
void reallocate_genes(Genes *g)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    g->data = (Gene *)xrealloc(g->data, g->n * sizeof(Gene));
    for (i = 0; i < g->n; i++) {
        g->data[i].type
            = (unsigned long *)xrealloc(g->data[i].type,
                                        blocks * sizeof(unsigned long));
        g->data[i].ancestral
            = (unsigned long *)xrealloc(g->data[i].ancestral,
                                        blocks * sizeof(unsigned long));
    }
}

/* Create a random set of n genes of length m */
Genes *random_genes(int n, int m)
{
    int i, j, bits = 1, blocks = divblocksize(m - 1) + 1;
    unsigned long filter = 1, data;
    Genes *g = (Genes *)xmalloc(sizeof(Genes));

    /* Set up basic structure */
    g->n = n;
    g->length = m;
    g->data = (Gene *)xmalloc(n * sizeof(Gene));

    /* Determine how many unbiased bits we obtain from rand, and create
     * a filter for obtaining them.
     */
    while (filter & XRAND_MAX) {
        bits++;
        filter <<= 1;
    }
    if (bits > BLOCKSIZE) {
        bits = BLOCKSIZE;
        filter = ~0;
    }
    else
        filter--;

    /* Create new genes */
    for (i = 0; i < n; i++) {
        /* Allocate memory for gene i */
        g->data[i].type = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
        g->data[i].ancestral
            = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        /* Set type to random value */
        for (j = 0; j < m; j += bits) {
            data = xrandom();
            g->data[i].type[divblocksize(j)] |= data << modblocksize(j);
            if ((modblocksize(j) + bits > BLOCKSIZE)
                    && (divblocksize(j) < blocks - 1))
                /* This block of random bits extends into the next block in g */
                g->data[i].type[divblocksize(j) + 1]
                    = (j == 0 ? 0 : data >> (BLOCKSIZE - modblocksize(j)));
        }
        /* Set gene to have ancestral material in all positions */
        for (j = 0; j < blocks - 1; j++)
            g->data[i].ancestral[j] = ~0;
        if (modblocksize(m) == 0)
            g->data[i].ancestral[j] = ~0;
        else {
            g->data[i].ancestral[j] = ((unsigned long)1 << modblocksize(m)) - 1;
            g->data[i].type[j] &= g->data[i].ancestral[j];
        }
    }

    return g;
}

/* Convert gene to a site oriented structure where the data is
 * comprised of columns instead of sequences.
 */
Sites *genes2sites(Genes *g)
{
    int i, j, block, blocks = divblocksize(g->n - 1) + 1;
    unsigned long index;
    Sites *s = (Sites *)xmalloc(sizeof(Sites));

    s->n = g->n;
    s->length = g->length;
    s->data = (Site *)xmalloc(g->length * sizeof(Site));
    for (i = 0; i < s->length; i++) {
        /* Construct site i */
        s->data[i].type
            = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
        s->data[i].ancestral
            = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
        /* Collect data from all sequences */
        index = 1;
        block = 0;
        for (j = 0; j < s->n; j++) {
            if (g->data[j].type[divblocksize(i)] & (unsigned long)1 << modblocksize(i))
                s->data[i].type[block] |= index;
            if (g->data[j].ancestral[divblocksize(i)] & (unsigned long)1 << modblocksize(i))
                s->data[i].ancestral[block] |= index;
            index <<= 1;
            if (index == 0) {
                index = 1;
                block += 1;
            }
        }
    }

    return s;
}

/* Convert a site oriented structure where the data is represented as
 * columns to a sequence oriented structure.
 */
Genes *sites2genes(Sites *s)
{
    int i, j, block, blocks = divblocksize(s->length - 1) + 1;
    unsigned long index;
    Genes *g = (Genes *)xmalloc(sizeof(Genes));

    g->n = s->n;
    g->length = s->length;
    g->data = (Gene *)xmalloc(g->n * sizeof(Gene));
    for (i = 0; i < g->n; i++) {
        /* Construct sequence i */
        g->data[i].type
            = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
        g->data[i].ancestral
            = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
        /* Collect data from all sites */
        index = 1;
        block = 0;
        for (j = 0; j < g->length; j++) {
            /* Handle site j */
            if (s->data[j].type[divblocksize(i)] & (unsigned long)1 << modblocksize(i))
                g->data[i].type[block] |= index;
            if (s->data[j].ancestral[divblocksize(i)] & (unsigned long)1 << modblocksize(i))
                g->data[i].ancestral[block] |= index;
            index <<= 1;
            if (index == 0) {
                index = 1;
                block += 1;
            }
        }
    }

    return g;
}

/* Return the gene set as an array of strings */
char **genes2string(Genes *g)
{
    int i, j;
    char **s = NULL;

    if (g->n > 0) {
        /* Allocate array for the sequences */
        s = (char **)xmalloc(g->n * sizeof(char *));

        /* Now convert the sequences in g to strings */
        for (i = 0; i < g->n; i++) {
            s[i] = (char *)xcalloc(g->length + 1, sizeof(char));
            for (j = 0; j < g->length; j++)
                if (g->data[i].type[divblocksize(j)] & (unsigned long)1 << modblocksize(j))
                    s[i][j] = '1';
                else if (g->data[i].ancestral[divblocksize(j)] & (unsigned long)1 << modblocksize(j))
                    s[i][j] = '0';
                else
                    s[i][j] = 'x';
        }
    }

    return s;
}


/* Output site data in s to file fp (stdout if fp is NULL), prefixed
 * with comment if not NULL.
 */
void output_sites(Sites *s, FILE *fp, char *comment)
{
    int i, j, block;
    unsigned long index;

    /* Set up file for output */
    if (fp == NULL)
        fp = stdout;

    /* Output genes */
    if (comment != NULL)
        fprintf(fp, "%s", comment);
    index = 1;
    block = 0;
    for (i = 0; i < s->n; i++) {
        for (j = 0; j < s->length; j++) {
            if (s->data[j].ancestral[block] & index) {
                if (s->data[j].type[block] & index) {
                    /* Gene i has type 1 at site j */
                    if (fputc((int)'1', fp) == EOF) {
                        fprintf(stderr, "Error outputting data\n");
                        exit(1);
                    }
                }
                else
                    /* Gene i has type 0 at site j */
                    if (fputc((int)'0', fp) == EOF) {
                        fprintf(stderr, "Error outputting data\n");
                        exit(1);
                    }
            }
            else {
                /* Gene i does not have ancestral data at site j */
                if (fputc((int)'x', fp) == EOF) {
                    fprintf(stderr, "Error outputting data\n");
                    exit(1);
                }
            }
        }
        /* Finished outputting gene i */
        if (fputc((int)'\n', fp) == EOF) {
            fprintf(stderr, "Error outputting data\n");
            exit(1);
        }
        index <<= 1;
        if (index == 0) {
            index = 1;
            block++;
        }
    }
}

/* Output site data in s to file fp (stdout if fp is NULL), followed
 * by a line of site indeces.
 */
void output_sites_indexed(Sites *s, FILE *fp)
{
    int i;

    /* Set up file for output */
    if (fp == NULL)
        fp = stdout;

    output_sites(s, NULL, NULL);
    for (i = 0; i < s->length; i++)
        fprintf(fp, "%d", i % 10);
    fprintf(fp, "\n");
}

/* Add data from g's block from index low to high to new starting from
 * index. If new becomes full, copy back to nblock in g and start a
 * new block in new. index and nblock are updated to reflect the changes.
 */
static void transfer_data(Genes *g, int block, int low, int high,
                          unsigned long **new, int *index, int *nblock)
{
    int i;
    /* Filter picks out the region that should be copied */
    unsigned long tmp, filter;

    filter = ((unsigned long)2 << high) - ((unsigned long)1 << low);

    if (high - low + 1 + *index >= BLOCKSIZE) {
        /* We have enough data to fill a block, and data does not need to
         * be shifted right as the above inequality implies that index >= low.
         */
        for (i = 0; i < g->n; i++) {
            /* Copy type data */
            tmp = g->data[i].type[block] & filter;
            if (low == *index) {
                /* If low = index, the standard shift to construct new initial
                 * value of new would be equal to BLOCKSIZE.
                 */
                g->data[i].type[*nblock] = new[0][i] | tmp;
                new[0][i] = 0;
            }
            else {
                g->data[i].type[*nblock] = new[0][i] | tmp << (*index - low);
                new[0][i] = tmp >> (low + BLOCKSIZE - *index);
            }
            /* Copy ancestral data */
            tmp = g->data[i].ancestral[block] & filter;
            if (low == *index) {
                /* If low = index, the standard shift to construct new initial
                 * value of new would be equal to BLOCKSIZE.
                 */
                g->data[i].ancestral[*nblock] = new[1][i] | tmp;
                new[1][i] = 0;
            }
            else {
                g->data[i].ancestral[*nblock] = new[1][i] | tmp << (*index - low);
                new[1][i] = tmp >> (low + BLOCKSIZE - *index);
            }
        }
        /* Update indeces - we have added high - low + 1 new sites and
         * copied an entire block to g.
         */
        *index += high - low + 1 - BLOCKSIZE;
        *nblock += 1;
    }
    else {
        /* Not enough data to fill block */
        if (*index > low)
            /* Data needs to be shifted left */
            for (i = 0; i < g->n; i++) {
                /* Copy type data */
                new[0][i] |= (g->data[i].type[block] & filter) << (*index - low);
                /* Copy ancestral data */
                new[1][i] |= (g->data[i].ancestral[block] & filter) << (*index - low);
            }
        else
            /* Data needs to be shifted right */
            for (i = 0; i < g->n; i++) {
                /* Copy type data */
                new[0][i] |= (g->data[i].type[block] & filter) >> (low - *index);
                /* Copy ancestral data */
                new[1][i] |= (g->data[i].ancestral[block] & filter) >> (low - *index);
            }
        *index += high - low + 1;
    }
}

/* Copy site given by from_index and from_block to site given by
 * to_index and to_block, assuming that the destination site is
 * subsumed in the source site. If conjugate is true, the copying is
 * of the conjugate of the source site.
 */
static void upgrade_subsumed(Genes *g, int to_index, int to_block,
                             int from_index, int from_block, int conjugate)
{
    int i;
    unsigned long to_pattern = (unsigned long)1 << to_index,
                  from_pattern = (unsigned long)1 << from_index;

    for (i = 0; i < g->n; i++)
        if (((g->data[i].ancestral[to_block] & to_pattern) == 0)
                && ((g->data[i].ancestral[from_block] & from_pattern) != 0)) {
            /* No ancestral material in destination site but ancestral
             * material in source site for sequence i.
             */
            g->data[i].ancestral[to_block] |= to_pattern;
            if (conjugate ^ ((g->data[i].type[from_block] & from_pattern) == 0))
                /* New type is 0 */
                g->data[i].type[to_block] &= ~to_pattern;
            else
                /* New type is 1 */
                g->data[i].type[to_block] |= to_pattern;
        }
}

/* Determine whether the site indicated by block and index and its
 * successor form a siamese pair. Return type indicates type of
 * twinning (identical or conjugate) and strict subsumation, if
 * present.
 */
static _Gene_TwinType siamese_site(int index, int block, Genes *g)
{
    int i;
    unsigned long pattern, ancestral,
             subsumation = 3; /* Initially both sites could be subsumed in each other */

    if (index != BLOCKSIZE - 1) {
        /* Both sites are in the current block */
        pattern = (unsigned long)1 << index;
        for (i = 0; i < g->n; i++) {
            ancestral = (g->data[i].ancestral[block] >> index) & 3;
            if (ancestral != 3) {
                /* At most one site carries ancestral material in this sequence */
                if (ancestral != 0) {
                    if ((ancestral & subsumation) == 0) {
                        /* Neither site is subsumed in the other */
                        index += 1;
                        return _GENE_INCOMPATIBLE;
                    }
                    /* Record which site subsumed the other in this sequences */
                    subsumation = ancestral;
                }
            }
            else {
                /* Both sites carry ancestral material in this sequence */
                if (((g->data[i].type[block] >> 1) ^ g->data[i].type[block])
                        & pattern) {
                    /* One site carries a 0, the other a 1 */
                    if (gene_knownancestor)
                        /* The two sites are different */
                        return _GENE_INCOMPATIBLE;
                    /* Look for conjugate columns */
                    for (i++; i < g->n; i++) {
                        ancestral = (g->data[i].ancestral[block] >> index) & 3;
                        if (ancestral != 3) {
                            /* At most one site carries ancestral material in this
                             * sequence.
                             */
                            if (ancestral != 0) {
                                /* Ancestral material in sites for sequence i */
                                if ((ancestral & subsumation) == 0)
                                    /* Neither site is subsumed in the other */
                                    return _GENE_INCOMPATIBLE;
                                /* Record which site subsumed the other in this sequences */
                                subsumation = ancestral;
                            }
                        }
                        else /* Both sites carry ancestral material */
                            if ((((g->data[i].type[block] >> 1) ^ g->data[i].type[block])
                                    & pattern) == 0)
                                /* Sites are not conjugate for sequence i */
                                return _GENE_INCOMPATIBLE;
                    }
                    /* Found no counterexample to conjugate columns */
                    switch(subsumation) {
                    case 1:
                        /* Right column is subsumed in left column */
                        return _GENE_LEFTCONJUGATE;
                    case 2:
                        /* Left column is subsumed in right column */
                        return _GENE_RIGHTCONJUGATE;
                    case 3:
                        /* Columns always share ancestral material */
                        return _GENE_CONJUGATE;
                    }
                }
                else {
                    /* Sites carry identical type - look for identical columns */
                    for (i++; i < g->n; i++) {
                        ancestral = (g->data[i].ancestral[block] >> index) & 3;
                        if (ancestral != 3) {
                            /* At most one site carries ancestral material in this
                             * sequence.
                             */
                            if (ancestral != 0) {
                                /* Ancestral material in sites for sequence i */
                                if ((ancestral & subsumation) == 0)
                                    /* Neither site is subsumed in the other */
                                    return _GENE_INCOMPATIBLE;
                                /* Record which site subsumed the other in this sequences */
                                subsumation = ancestral;
                            }
                        }
                        else if ((((g->data[i].type[block] >> 1) ^ g->data[i].type[block])
                                  & pattern) != 0)
                            /* Sites are not identical for sequence i */
                            return _GENE_INCOMPATIBLE;
                    }
                    /* Found no counterexample to identical columns */
                    switch(subsumation) {
                    case 1:
                        /* Right column is subsumed in left column */
                        return _GENE_LEFTIDENTICAL;
                    case 2:
                        /* Left column is subsumed in right column */
                        return _GENE_RIGHTIDENTICAL;
                    case 3:
                        /* Columns always share ancestral material */
                        return _GENE_IDENTICAL;
                    }
                }
            }
        }
        /* One column carries only non-ancestral material, or there are no
         * sequeces.
         */
        switch(subsumation) {
        case 1:
            /* Right column is subsumed in left column */
            return _GENE_LEFTIDENTICAL;
        case 2:
            /* Left column is subsumed in right column */
            return _GENE_RIGHTIDENTICAL;
        case 3:
            /* Sanity check - if there is at least one sequence in g we
             * should never reach this point.
             */
            if (g->n > 0) {
                fprintf(stderr, "ERROR: Something wrong in siamese_site function, please email error report\n");
                exit(1);
            }
            return _GENE_IDENTICAL;
        }
    }
    else {
        /* Sites span a block boundary */
        for (i = 0; i < g->n; i++) {
            ancestral = g->data[i].ancestral[block] >> (BLOCKSIZE - 1);
            if ((ancestral & g->data[i].ancestral[block + 1] & 1) == 0) {
                /* Not ancestral material in both sites for sequence i */
                if ((ancestral | (g->data[i].ancestral[block + 1] & 1)) != 0) {
                    /* Ancestral material in one of the sites */
                    if ((subsumation != 3) && (subsumation != ancestral))
                        /* Neither site is subsumed in the other */
                        return _GENE_INCOMPATIBLE;
                    else
                        /* Record subsumation */
                        subsumation = ancestral;
                }
            }
            else {
                /* Ancestral material in both sites */
                if ((g->data[i].type[block] >> (BLOCKSIZE - 1))
                        ^ (g->data[i].type[block + 1] & 1)) {
                    if (gene_knownancestor)
                        /* The two sites are not identical */
                        return _GENE_INCOMPATIBLE;
                    /* Look for conjugate columns */
                    for (i++; i < g->n; i++) {
                        ancestral = g->data[i].ancestral[block] >> (BLOCKSIZE - 1);
                        if ((ancestral & g->data[i].ancestral[block + 1] & 1) == 0) {
                            /* Not ancestral material in both sites for sequence i */
                            if ((ancestral | (g->data[i].ancestral[block + 1] & 1)) != 0) {
                                /* Ancestral material in one of the sites */
                                if ((subsumation != 3) && (subsumation != ancestral))
                                    /* Neither site is subsumed in the other */
                                    return _GENE_INCOMPATIBLE;
                                else
                                    /* Record subsumation */
                                    subsumation = ancestral;
                            }
                        }
                        else {
                            /* Ancestral material in both sites for sequence i */
                            if (!((g->data[i].type[block] >> (BLOCKSIZE - 1))
                                    ^ (g->data[i].type[block + 1] & 1)))
                                /* The two sites are not conjugate */
                                return _GENE_INCOMPATIBLE;
                        }
                    }
                    /* Found no counterexample to conjugate columns */
                    switch(subsumation) {
                    case 1:
                        /* Right column is subsumed in left column */
                        return _GENE_LEFTCONJUGATE;
                    case 0:
                        /* Left column is subsumed in right column */
                        return _GENE_RIGHTCONJUGATE;
                    case 3:
                        /* Columns always share ancestral material */
                        return _GENE_CONJUGATE;
                    }
                }
                else {
                    /* Look for identical columns */
                    for (i++; i < g->n; i++) {
                        ancestral = g->data[i].ancestral[block] >> (BLOCKSIZE - 1);
                        if ((ancestral & g->data[i].ancestral[block + 1] & 1) == 0) {
                            /* Not ancestral material in both sites for sequence i */
                            if ((ancestral | (g->data[i].ancestral[block + 1] & 1)) != 0) {
                                /* Ancestral material in one of the sites */
                                if ((subsumation != 3) && (subsumation != ancestral))
                                    /* Neither site is subsumed in the other */
                                    return _GENE_INCOMPATIBLE;
                                else
                                    /* Record subsumation */
                                    subsumation = ancestral;
                            }
                        }
                        else {
                            if ((g->data[i].type[block] >> (BLOCKSIZE - 1))
                                    ^ (g->data[i].type[block + 1] & 1))
                                /* The two sites are not identical */
                                return _GENE_INCOMPATIBLE;
                        }
                    }
                    /* Found no counterexample to identical columns */
                    switch(subsumation) {
                    case 1:
                        /* Right column is subsumed in left column */
                        return _GENE_LEFTIDENTICAL;
                    case 0:
                        /* Left column is subsumed in right column */
                        return _GENE_RIGHTIDENTICAL;
                    case 3:
                        /* Columns always share ancestral material */
                        return _GENE_IDENTICAL;
                    }
                }
            }
        }
        /* One column carries only non-ancestral material, or there are no
         * sequeces.
         */
        switch(subsumation) {
        case 1:
            /* Right column is subsumed in left column */
            return _GENE_LEFTIDENTICAL;
        case 0:
            /* Left column is subsumed in right column */
            return _GENE_RIGHTIDENTICAL;
        case 3:
            /* Sanity check - if there is at least one sequence in g we
             * should never reach this point.
             */
            if (g->n > 0) {
                fprintf(stderr, "ERROR: Something wrong in siamese_site function, please email error report\n");
                exit(1);
            }
            return _GENE_IDENTICAL;
        }
    }

    /* We should never reach this point */
    fprintf(stderr, "ERROR: Something wrong in siamese_site function, please email error report\n");
    exit(1);
    return _GENE_IDENTICAL;
}

/* Determine whether the sites indicated by nblock, nindex and block,
 * index form a siamese pair, and, if it is present, subsumation.
 */
static _Gene_TwinType siamese_sites(int nindex, int nblock, int index, int block,
                                    Genes *g)
{
    int i;
    unsigned long npattern = (unsigned long)1 << nindex,
                  pattern = (unsigned long)1 << index,
                  ancestral,
                  subsumation = 3; /* Record subsumation with 3 for both ways, npattern
		      * for nindex, nblock subsuming, and 0 for index, block
		      * subsuming. */

    for (i = 0; i < g->n; i++)
        if (((ancestral = (g->data[i].ancestral[nblock] & npattern)) != 0)
                || ((g->data[i].ancestral[block] & pattern) != 0)) {
            /* Sequence i has ancestral material in sites */
            if ((ancestral == 0)
                    || ((g->data[i].ancestral[block] & pattern) == 0)) {
                /* But only in one of the sites */
                if ((subsumation != 3) && (subsumation != ancestral))
                    /* Neither site is subsumed in the other */
                    return _GENE_INCOMPATIBLE;
                /* Record subsumation */
                subsumation = ancestral;
            }
            else
                /* Sequence i has ancestral material in both sites */
                if (((g->data[i].type[nblock] >> nindex) & 1)
                        != ((g->data[i].type[block] >> index) & 1)) {
                    if (gene_knownancestor)
                        /* Sites are not identical */
                        return _GENE_INCOMPATIBLE;
                    /* Look for conjugate sites */
                    for (i++; i < g->n; i++)
                        if (((ancestral =(g->data[i].ancestral[nblock] & npattern)) != 0)
                                || ((g->data[i].ancestral[block] & pattern) != 0)) {
                            /* Sequence i has ancestral material in sites */
                            if ((ancestral == 0)
                                    || ((g->data[i].ancestral[block] & pattern) == 0)) {
                                /* But only in one of the sites */
                                if ((subsumation != 3) && (subsumation != ancestral))
                                    /* Neither site is subsumed in the other */
                                    return _GENE_INCOMPATIBLE;
                                /* Record subsumation */
                                subsumation = ancestral;
                            }
                            else if (!(((g->data[i].type[nblock] >> nindex) & 1)
                                       ^ ((g->data[i].type[block] >> index) & 1)))
                                /* Sites are not conjugate */
                                return _GENE_INCOMPATIBLE;
                        }
                    /* No counterexample to conjugate columns found */
                    switch(subsumation) {
                    case 0:
                        /* Left column is subsumed in right column */
                        return _GENE_RIGHTCONJUGATE;
                    case 3:
                        /* Columns always share ancestral material */
                        return _GENE_CONJUGATE;
                    default:
                        /* Right column is subsumed in left column */
                        return _GENE_LEFTCONJUGATE;
                    }
                }
                else {
                    /* Look for identical sites */
                    for (i++; i < g->n; i++)
                        if (((ancestral =(g->data[i].ancestral[nblock] & npattern)) != 0)
                                || ((g->data[i].ancestral[block] & pattern) != 0)) {
                            /* Sequence i has ancestral material in sites */
                            if ((ancestral == 0)
                                    || ((g->data[i].ancestral[block] & pattern) == 0)) {
                                /* But only in one of the sites */
                                if ((subsumation != 3) && (subsumation != ancestral))
                                    /* Neither site is subsumed in the other */
                                    return _GENE_INCOMPATIBLE;
                                /* Record subsumation */
                                subsumation = ancestral;
                            }
                            else if (((g->data[i].type[nblock] >> nindex) & 1)
                                     ^ ((g->data[i].type[block] >> index) & 1))
                                /* Sites are not identical */
                                return _GENE_INCOMPATIBLE;
                        }
                    /* No counterexample to identical columns found */
                    switch(subsumation) {
                    case 0:
                        /* Left column is subsumed in right column */
                        return _GENE_RIGHTIDENTICAL;
                    case 3:
                        /* Columns always share ancestral material */
                        return _GENE_IDENTICAL;
                    default:
                        /* Right column is subsumed in left column */
                        return _GENE_LEFTIDENTICAL;
                    }
                }
        }
    /* One column didn't contain ancestral material, or there are no
     * sequences in g.
     */
    switch(subsumation) {
    case 0:
        /* Left column is subsumed in right column */
        return _GENE_RIGHTIDENTICAL;
    case 3:
        /* Sanity check - if there is at least one sequence in g we
         * should never reach this point.
         */
        if (g->n > 0) {
            fprintf(stderr, "ERROR: Something wrong in siamese_site function, please email error report\n");
            exit(1);
        }
        return _GENE_IDENTICAL;
    default:
        /* Right column is subsumed in left column */
        return _GENE_RIGHTIDENTICAL;
    }
}

/* Remove Siamese twins: if two neighbouring sites are identical (or
 * conjugate if the ancestral type is unknown) we can safely merge
 * them. This also applies if one site is subsumed in the other, but
 * identity or conjugation holds for all sequences where both carry
 * ancestral material. Return value tells whether any sites were
 * merged.
 */
int remove_siamesetwins(Genes *g)
{
    int i, j, n = 0, k, kk, p, s, index = 0, mindex, mblock, pindex, pblock, cindex, cblock;
    _SiameseBlock *state;
    unsigned long *tmp[2];
    _Gene_TwinType types = _GENE_INCOMPATIBLE, oldtypes;
    Event *e;
#ifdef HAPLOTYPE_BLOCKS
    SuperColumn *c, *d;
#endif

    /* Skip initial stretch of non-siamese sites */
    while ((index < g->length - 1) &&
            !(types = siamese_site(modblocksize(index), divblocksize(index), g)))
        index++;

    if (types) {
        /* There are Siamese twins in g */
        state = (_SiameseBlock *)xmalloc((1 + g->length / 2)
                                         * sizeof(_SiameseBlock));
        for (;;) {
            state[n].start = index;
            if (types & _GENE_RIGHT)
                /* Left site strictly subsumed in right site */
                state[n].master = index + 1;
            else
                /* Right site subsumed in left site */
                state[n].master = index;

            /* Explore extension of block of Siamese sites to the right */
            mindex = modblocksize(state[n].master);
            mblock = divblocksize(state[n].master);
            index += 2;
            while ((index < g->length) &&
                    (types = siamese_sites(mindex, mblock, modblocksize(index),
                                           divblocksize(index), g))) {
                oldtypes = types;
                if (types & _GENE_RIGHT) {
                    /* Previous master site strictly subsumed in this site */
                    state[n].master = index;
                    mindex = modblocksize(index);
                    mblock = divblocksize(index);
                }
                index++;
            }
            /* Found rightmost site in Siamese block */
            state[n].end = index - 1;

            /* If leftmost site is not master site, then there may be more
             * Siamese sites to the left.
             */
            if (state[n].master !=  state[n].start) {
                /* We prefer master site at first or last site of block, since
                 * it cannot be first let's start by checking if it can be last.
                 */
                if ((index < g->length) && (index > state[n].start + 2))
                    /* Determine type of last Siamese call */
                    types = oldtypes;

                if ((types & _GENE_LEFT) == 0)
                    /* Master site can be last site */
                    state[n].master = state[n].end;

                /* Now check for extra Siamese sites to the left */
                index = state[n].start - 1;
                while ((index >= 0) &&
                        siamese_sites(mindex, mblock, modblocksize(index),
                                      divblocksize(index), g)) {
                    /* We do not need to check for subsumation direction, as any
                     * site to the left of where we first identified Siameseness
                     * cannot be a master site.
                     */
                    if ((n > 0) && (index == state[n - 1].end)) {
                        /* This block includes previous block - if master site of
                         * current block is a Siamese twin with master site of
                         * previous block we can incorporate the previous block
                         * into the current block.
                         */
                        if (siamese_sites(mindex, mblock, modblocksize(state[n - 1].master),
                                          divblocksize(state[n - 1].master), g)) {
                            index = state[n - 1].start - 1;
                            state[n - 1].end = state[n].end;
                            /* As before, we do not need to check subsumation direction */
                            state[n - 1].master = state[n].master;
                            n--;
                        }
                        else
                            /* Cannot merge the two blocks, so this block cannot be
                             * extended further to the left.
                             */
                            break;
                    }
                    else
                        /* Site not part of other block, so just continue with
                         * previous site.
                         */
                        index--;
                }

                if (index < 0)
                    /* Loop terminated by reaching beginning of sequences */
                    state[n].start = 0;
                else
                    /* Reached site not included in Siamese block */
                    state[n].start = index + 1;
                index = state[n].end + 1;
            }
            n++;

            /* Continue scanning for next Siamese block */
            while ((index < g->length - 1) &&
                    !(types =
                          siamese_site(modblocksize(index), divblocksize(index), g)))
                index++;
            if (index >= g->length - 1)
                /* Reached end of sequences */
                break;
        }

        /* Replace Siamese blocks with master site */
        /* Allocate structures for compacting the information in g */
        tmp[0] = (unsigned long *)xcalloc(sizeof(unsigned long), g->n);
        tmp[1] = (unsigned long *)xcalloc(sizeof(unsigned long), g->n);

        /* Insert dummy Siamese block at end */
        state[n].start = g->length;
        state[n].end = state[n].master = g->length + 1;
        n++;

        /* Determine first point where sites are being collapsed (i.e. removed) */
        mblock = divblocksize(state[0].start + (state[0].start == state[0].master));
        index = modblocksize(state[0].start + (state[0].start == state[0].master));
        mindex = 0;
        if (index > 0)
            /* Transfer part of block prior to Siamese sites to tmp */
            transfer_data(g, mblock, 0, index - 1, tmp, &mindex, &mblock);
        if ((state[0].start != state[0].master) &&
                (state[0].end != state[0].master))
            /* Also transfer master site to tmp */
            transfer_data(g, divblocksize(state[0].master),
                          modblocksize(state[0].master),
                          modblocksize(state[0].master), tmp, &mindex, &mblock);

        /* Now transfer all the sites between blocks of Siamese sites */
        for (i = 1; i < n; i++) {
            /* Block and index of first site to transfer */
            pblock = divblocksize(state[i - 1].end
                                  + (state[i - 1].master != state[i - 1].end));
            pindex = modblocksize(state[i - 1].end
                                  + (state[i - 1].master != state[i - 1].end));
            /* Block and index of last site to transfer */
            cblock = divblocksize(state[i].start
                                  - (state[i].master != state[i].start));
            cindex = modblocksize(state[i].start
                                  - (state[i].master != state[i].start));
            if ((pblock == cblock) && (pindex <= cindex))
                /* All sites to be transferred in same block */
                transfer_data(g, cblock, pindex, cindex, tmp, &mindex, &mblock);
            else if (state[i].start > state[i - 1].end + 1) {
                /* Sites to be transferred span multiple blocks */
                /* First transfer initial block */
                transfer_data(g, pblock, pindex, BLOCKSIZE - 1, tmp, &mindex, &mblock);
                /* Then transfer full blocks */
                for (j = pblock + 1; j < cblock; j++)
                    transfer_data(g, j, 0, BLOCKSIZE - 1, tmp, &mindex, &mblock);
                /* Finally transfer last block */
                transfer_data(g, cblock, 0, cindex, tmp, &mindex, &mblock);
            }
		// WARNING bug fix (AI) 
		// start
		else if (state[i-1].end + 1 == state[i].master) {
		// Master sites are straddling two sides of a block split
		// Need to transfer both
        // First transfer initial block 
                transfer_data(g, pblock, pindex, BLOCKSIZE - 1, tmp, &mindex, &mblock);
                /* Then transfer full blocks */
                for (j = pblock + 1; j < cblock; j++)
                    transfer_data(g, j, 0, BLOCKSIZE - 1, tmp, &mindex, &mblock);
                /* Finally transfer last block */
                transfer_data(g, cblock, 0, cindex, tmp, &mindex, &mblock);
		}
		// end
            if ((state[i].master != state[i].start)
                    && (state[i].master != state[i].end))
                /* We need to transfer master site from this Siamese block */
                transfer_data(g, divblocksize(state[i].master),
                              modblocksize(state[i].master),
                              modblocksize(state[i].master), tmp, &mindex, &mblock);
        }
        /* Remove dummy */
        n--;

        /* Update g to reflect changes */
        if (mindex > 0) {
            /* Copy data in tmp back into g */
            for (i = 0; i < g->n; i++) {
                g->data[i].type[mblock] = tmp[0][i];
                g->data[i].ancestral[mblock] = tmp[1][i];
            }
            g->length = mulblocksize(mblock) + mindex;
        }
        else
            g->length = mulblocksize(mblock);

        if (eventlist != NULL)
            /* Insert collapse of Siamese twins into list of events */
            for (i = n - 1; i >= 0; i--) {
                /* Collapse into master */
                for (j = state[i].master; j < state[i].end; j++) {
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = COLLAPSE;
                    e->event.collapse = state[i].master;
                    Enqueue(eventlist, e);
                }
                for (j = state[i].master - 1; j >= state[i].start; j--) {
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = COLLAPSE;
                    e->event.collapse = j;
                    Enqueue(eventlist, e);
                }
            }
        
        // Updated list of sites, labelling the master sites with -(number of columns collapsed).
        if(sites != NULL) {    
            k = 0;
            for (i = 0; i < n ; i++) {
                s = -(int)elist_get(sites, state[i].master);
                if(s <= 0) {
                    s = 1;
                }
                for (j = state[i].start; j < state[i].master; j++) {
//                     printf("Merging Siamese twin columns %d -> %d\n", j, state[i].master);
                    p = (int)elist_get(sites, j - k);
                    if(p < 0) {
                        s = s - p;
                    } else {
                        s++;
                    }
                    elist_remove(sites, j - k);
                    k++;
                }
                kk = k;
                for (j = state[i].master + 1; j <= state[i].end; j++) {
//                     printf("Merging Siamese twin columns %d <- %d\n", state[i].master, j);
                    p = (int)elist_get(sites, j - k);
                    if(p < 0) {
                        s = s - p;
                    } else {
                        s++;
                    }
                    elist_remove(sites, j - k);
                    k++;
                }
                // Master site is now at position state[i].start
                // Change this to be -(number of columns collapsed)
//                 printf("%d %d\n", state[i].master-kk, -s);
                elist_change(sites, state[i].master-kk, (void *)(-s));
            }
        }

#ifdef ENABLE_VERBOSE
        if (verbose())
            /* Print information about Siamese columns collapsed */
            for (i = 0; i < n ; i++) {
                for (j = state[i].start; j < state[i].master; j++)
                    printf("Merging Siamese twin columns %d and %d\n", j,
                           state[i].master);
                for (j = state[i].master + 1; j <= state[i].end; j++)
                    printf("Merging Siamese twin columns %d and %d\n", state[i].master,
                           j);
            }
#endif
#ifdef HAPLOTYPE_BLOCKS
        if (representativeness != NULL)
            /* Update list of blocks of collapsed sites to reflect the
             * Siamese blocks just collapsed.
             */
            for (i = n - 1; i >= 0; i--) {
                /* Find start position of this Siamese block in list of block
                 * of sites that have already been collapsed.
                 */
                c = (SuperColumn *)SetCounter(representativeness_counter,
                                              state[i].start);
                /* Merge with block for other sites in this Siamese block */
                for (j = state[i].start; j < state[i].end; j++) {
                    d = (SuperColumn *)Next(representativeness_counter);
                    c->right = d->right;
                    RemoveMoveLeft(representativeness, representativeness_counter);
                }
            }
#endif



        /* Clean up */
        free(tmp[0]);
        free(tmp[1]);
        free(state);

        /* Some Siamese twins were collapsed */
        return 1;
    }

    /* No Siamese twins were found */
    return 0;
}

/* Remove all uninformative sites, i.e. sites where at most one of the
 * sequences with ancestral material differs. Return value tells
 * whether any columns were uninformative.
 */
int remove_uninformative(Genes *g)
{
    int i, j, first, nindex = -1, nblock,
                     blocks = divblocksize(g->length - 1) + 1;
    unsigned long *tmp[2], one0, two0, one1, two1, current, index;
    Event *e;
    int n = 0;
#ifdef HAPLOTYPE_BLOCKS
    int m = 0;
#endif
    
    for (i = 0; i < blocks; i++) {
        /* Compile type information from all genes for block i */
        one1 = g->data[0].type[i];
        one0 = ~g->data[0].type[i] & g->data[0].ancestral[i];
        two0 = two1 = 0;
        for (j = 1; j < g->n; j++) {
            two1 |= one1 & g->data[j].type[i];
            one1 |= g->data[j].type[i];
            current = ~g->data[j].type[i] & g->data[j].ancestral[i];
            two0 |= one0 & current;
            one0 |= current;
        }
        /* two1 (two0) now has a 1 in all positions where there are at
         * least two 1s (0s) in the data.
         */
        if (gene_knownancestor)
            /* When wild type is known we only care about mutant count */
            current = two1;
        else
            /* Site is informative if there is at least two occurrences of
             * each type.
             */
            current = two0 & two1;

        /* Remove uninformative columns in block i */
        if (((i < blocks - 1) && ~current)
                || ((i >= blocks - 1)
                    && (current != ((unsigned long)2 << modblocksize(g->length - 1)) - 1))) {
            /* Of which there are some */
            if (nindex < 0) {
                /* But this is the first block containing uninformative columns */
                nblock = i;
                nindex = 0;
                /* Set up structures for building new gene data */
                tmp[0] = (unsigned long *)xcalloc(g->n, sizeof(unsigned long));
                tmp[1] = (unsigned long *)xcalloc(g->n, sizeof(unsigned long));
            }
            j = 0;
            index = 1;
            while (j < BLOCKSIZE) {
                /* Skip uninformative columns */
                while (~current & index) {
#ifdef ENABLE_VERBOSE
                    if (verbose()) {
                        if (mulblocksize(i) + j < g->length)
                            printf("Removing uninformative column %d\n",
                                   mulblocksize(i) + j);
                        else {
                            j = BLOCKSIZE;
                            index = 0;
                        }
                    }
#endif

                    if ((eventlist != NULL)
                            && ((gene_knownancestor ? one1 : one0 & one1) & index)) {
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = SUBSTITUTION;
                        e->event.s.seq = -1;
                        e->event.s.site = mulblocksize(i) + j - n;
                        Enqueue(eventlist, e);
                    }
                    if(sites != NULL) {
                        if (mulblocksize(i) + j < g->length) {
//                             printf("Removing site labelled %d by deleting element number %d\n", mulblocksize(i) + j, mulblocksize(i) + j - n);
                            elist_remove(sites, mulblocksize(i) + j - n);
                        }
                        else {
                                j = BLOCKSIZE;
                                index = 0;
                        }
                    }
                    n++;
#ifdef HAPLOTYPE_BLOCKS
                    if ((representativeness != NULL)
                            && (mulblocksize(i) + j < g->length)) {
                        /* Remove uninformative sites */
                        SetCounter(representativeness_counter, mulblocksize(i) + j - m);
                        free(RemoveMoveLeft(representativeness,
                                            representativeness_counter));
                        m++;
                    }
#endif
                    j++;
                    index <<= 1;
                }

                if (index) {
                    /* We haven't reached end of block */
                    first = j;

                    /* Find block of informative columns */
                    while (current & index) {
                        j++;
                        index <<= 1;
                    }

                    /* Copy current block of informative sites */
                    transfer_data(g, i, first, j - 1, tmp, &nindex, &nblock);
                }
            }
        }
        else if (nindex >= 0) {
            /* No uninformative columns in block i; but we have started
             * copying - continue this process.
             */
            if (i < blocks - 1)
                transfer_data(g, i, 0, BLOCKSIZE - 1, tmp, &nindex, &nblock);
            else
                transfer_data(g, i, 0, modblocksize(g->length - 1), tmp, &nindex,
                              &nblock);
        }
    }

    /* Update gene */
    if (nindex >= 0) {
        /* The data contained at least one uninformative column */
        if (nindex > 0) {
            /* Copy data in tmp back into g */
            for (i = 0; i < g->n; i++) {
                g->data[i].type[nblock] = tmp[0][i];
                g->data[i].ancestral[nblock] = tmp[1][i];
            }
            g->length = mulblocksize(nblock) + nindex;
        }
        else
            g->length = mulblocksize(nblock);

        if (g->length == 0) {
            /* Everything was removed */
            if (eventlist != NULL)
                for (i = 1; i < g->n; i++) {
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = COALESCENCE;
                    e->event.c.s1 = 0;
                    e->event.c.s2 = 1;
                    Enqueue(eventlist, e);
                }
            for (i = 0; i < g->n; i++) {
                free(g->data[i].type);
                free(g->data[i].ancestral);
            }
            g->n = 0;
        }
        /* Clean up */
        free(tmp[0]);
        free(tmp[1]);

        return 1;
    }
    else
        return 0;
}


/* Remove all non-segregating sites, i.e. sites where none of the
 * sequences with ancestral material differs. Return value tells
 * whether any columns were non-segregating.
 */
int remove_nonsegregating(Genes *g)
{
    int i, j, first, nindex = -1, nblock,
                     blocks = divblocksize(g->length - 1) + 1;
    unsigned long *tmp[2], one0, one1, current, index;
    Event *e;
#ifdef HAPLOTYPE_BLOCKS
    int m = 0;
#endif

    for (i = 0; i < blocks; i++) {
        /* Compile type information from all genes for block i */
        one1 = g->data[0].type[i];
        one0 = ~g->data[0].type[i] & g->data[0].ancestral[i];
        for (j = 1; j < g->n; j++) {
            one1 |= g->data[j].type[i];
            one0 |= ~g->data[j].type[i] & g->data[j].ancestral[i];
        }
        /* one1 (one0) now has a 1 in all positions where there are at
         * least one 1 (0) in the data.
         */
        if (gene_knownancestor)
            /* When wild type is known we only care about mutant count */
            current = one1;
        else
            /* Site is segregating if there is at least one occurrence of
             * each type.
             */
            current = one0 & one1;

        /* Remove non-segregating columns in block i */
        if (((i < blocks - 1) && ~current)
                || ((i >= blocks - 1)
                    && (current != ((unsigned long)2 << modblocksize(g->length - 1)) - 1))) {
            /* Of which there are some */
            if (nindex < 0) {
                /* But this is the first block containing uninformative columns */
                nblock = i;
                nindex = 0;
                /* Set up structures for building new gene data */
                tmp[0] = (unsigned long *)xcalloc(g->n, sizeof(unsigned long));
                tmp[1] = (unsigned long *)xcalloc(g->n, sizeof(unsigned long));
            }
            j = 0;
            index = 1;
            while (j < BLOCKSIZE) {
                /* Skip non-segregating columns */
                while (~current & index) {
#ifdef ENABLE_VERBOSE
                    if (verbose()) {
                        if (mulblocksize(i) + j < g->length)
                            printf("Removing non-segregating column %d\n",
                                   mulblocksize(i) + j);
                        else {
                            j = BLOCKSIZE;
                            index = 0;
                        }
                    }
#endif
#ifdef HAPLOTYPE_BLOCKS
                    if ((representativeness != NULL)
                            && (mulblocksize(i) + j < g->length)) {
                        /* Remove non-segregating sites */
                        SetCounter(representativeness_counter, mulblocksize(i) + j - m);
                        free(RemoveMoveLeft(representativeness,
                                            representativeness_counter));
                        m++;
                    }
#endif
                    j++;
                    index <<= 1;
                }

                if (index) {
                    /* We haven't reached end of block */
                    first = j;

                    /* Find block of segregating columns */
                    while (current & index) {
                        j++;
                        index <<= 1;
                    }

                    /* Copy current block of informative sites */
                    transfer_data(g, i, first, j - 1, tmp, &nindex, &nblock);
                }
            }
        }
        else if (nindex >= 0) {
            /* We have started copying - continue this process */
            if (i < blocks - 1)
                transfer_data(g, i, 0, BLOCKSIZE - 1, tmp, &nindex, &nblock);
            else
                transfer_data(g, i, 0, modblocksize(g->length - 1), tmp, &nindex,
                              &nblock);
        }
    }

    /* Update gene */
    if (nindex >= 0) {
        /* The data contained at least one non-segregating column */
        if (nindex > 0) {
            /* Copy data in tmp back into g */
            for (i = 0; i < g->n; i++) {
                g->data[i].type[nblock] = tmp[0][i];
                g->data[i].ancestral[nblock] = tmp[1][i];
            }
            g->length = mulblocksize(nblock) + nindex;
        }
        else
            g->length = mulblocksize(nblock);

        if (g->length == 0) {
            /* Everything was removed */
            if (eventlist != NULL)
                for (i = 1; i < g->n; i++) {
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = COALESCENCE;
                    e->event.c.s1 = 0;
                    e->event.c.s2 = 1;
                    Enqueue(eventlist, e);
                }
            for (i = 0; i < g->n; i++) {
                free(g->data[i].type);
                free(g->data[i].ancestral);
            }
            g->n = 0;
        }
        /* Clean up */
        free(tmp[0]);
        free(tmp[1]);

        return 1;
    }
    else
        return 0;
}

/* Coalesce all compatible sequences where the ancestral material of
 * one is a subset of the ancestral material of the other. Return
 * value is number of coalesces performed.
 */
int coalesce_subsumed(Genes *g)
{
    int i, j, k, changes = 0, conflicts, blocks = divblocksize(g->length - 1) + 1;
    Event *e;

//     printf("=====================================\n");
//     printf("coalesce_subsumed on data:\n");
//     output_genes(g, NULL, NULL);
//     printf("=====================================\n");

    /* Look for coalescences */
    for (i = 0; i < g->n; i++)
        if (g->data[i].type != NULL) {
            for (j = i + 1; j < g->n; j++)
                if (g->data[j].type != NULL) {
                    /* Check whether sequence j is subsumed in sequence i, or
                     * sequence i is subsumed in sequence j.
                     */
                    conflicts = 0;
                    /* Check types */
                    for (k = 0; k < blocks; k++)
                        if ((g->data[i].type[k] ^ g->data[j].type[k])
                                & (g->data[i].ancestral[k] & g->data[j].ancestral[k])) {
                            /* We found a conflict */
                            break;
                        }
                    if (!(k < blocks)) {
                        /* Check ancestral material */
                        for (k = 0; k < blocks; k++) {
                            if (~g->data[i].ancestral[k] & g->data[j].ancestral[k]) {
                                /* Sequence j has ancestral material in a place where
                                 * sequence i doesn't
                                 */
                                if (conflicts == 1)
                                    /* Already encountered a situtation where sequence i
                                     * had ancestral material in a place where sequence
                                     * j didn't.
                                     */
                                    break;
                                else
                                    /* Record event */
                                    conflicts = -1;
                            }
                            if (~g->data[j].ancestral[k] & g->data[i].ancestral[k]) {
                                /* Sequence i has ancestral material in a place where
                                 * sequence j doesn't
                                 */
                                if (conflicts == -1)
                                    /* Already encountered a situtation where sequence j
                                     * had ancestral material in a place where sequence
                                     * i didn't.
                                     */
                                    break;
                                else
                                    /* Record event */
                                    conflicts = 1;
                            }
                        }
                        if (!(k < blocks)) {
                            /* We can assimilate one sequence into the other */
                            changes++;
                            #ifdef ENABLE_VERBOSE
                            if (verbose())
                                printf("Coalescing sequences %d and %d\n", i, j);
                            #endif
//                             printf("Coalescing sequences %d and %d\n", i, j);
                            if(elements != NULL) {
                                elist_change(elements, i, (void *)(-1));
                                elist_change(elements, j, (void *)(-1));
                            }
                            
                            if (conflicts == -1) {
                                /* Sequence j cannot be assimilated into sequence i,
                                 * so it must be the other way round.
                                 */
                                free(g->data[i].type);
                                free(g->data[i].ancestral);
                                g->data[i].type = NULL;
                                /* We just coalesced sequence i into sequence j, so we
                                 * should look no further for sequences that can
                                 * coalesce with sequence i.
                                 */
                                break;
                            }
                            else {
                                /* Sequence j can be assimilated into sequence i */
                                free(g->data[j].type);
                                free(g->data[j].ancestral);
                                g->data[j].type = NULL;
                            }
                        }
                    }
                }
        }

    if (changes > 0) {
        
        /* We did coalesce some sequences - compact list of sequences */
        for (i = 0; g->data[i].type != NULL; i++);
        if(elements != NULL) {
            elist_remove(elements, i);
        }
        if (eventlist != NULL) {
            e = (Event *)xmalloc(sizeof(Event));
            e->type = REMOVE;
            e->event.remove = i;
            Enqueue(eventlist, e);
        }
        for (j = i + 1; j < g->n; j++)
            if (g->data[j].type != NULL) {
                /* This sequence was not assimilated - move it to its new position */
                memcpy(&(g->data[i]), &(g->data[j]), sizeof(Gene));
                i++;
            }
            else  {
                if(elements != NULL) {
                    elist_remove(elements, i);
                }
                if (eventlist != NULL) {
                e = (Event *)xmalloc(sizeof(Event));
                e->type = REMOVE;
                e->event.remove = i;
                Enqueue(eventlist, e);
                }
            }
        /* Update g */
        g->n = i;
    }
    
    return changes;
}

/* Perform all Siamese twin collapses, removal of uninformative sites
 * and coalesces that are guaranteed not to increase the number of
 * recombinations required. The data structure g is modified to
 * reflect events, and the number of events is returned.
 */
int implode_genes(Genes *g)
{
    int n = g->n, m = g->length, change = 1, tmp;

#ifdef ENABLE_VERBOSE
//     if (verbe())
        printf("Imploding genes:\n");
#endif

    while (change) {
#ifdef ENABLE_VERBOSE
        if (verbose())
            output_genes_indexed(g, NULL);
#endif
        /* Look for uninformative sites */
        change = remove_uninformative(g);
        if (g->n == 0) break;
        /* Look for removable Siamese twins */
        /* Collapsing Siamese twins should not allow new mutations, and
         * coalesces are handled presently.
         */
        remove_siamesetwins(g); // collapses identical neighbouring cols
        /* Look for safe coalesces */
        tmp = coalesce_subsumed(g);
        change |= tmp;
    }

    return (n - g->n) + (m - g->length);
}

/* Simple, fast checks for whether recombinations may be required to
 * explain g. This is not an exhaustive check, so even if it returns
 * False, g may be explained without recombinations.
 */
int no_recombinations_required(Genes *g)
{
    int i;
    Event *e;
    if (g->n < (gene_knownancestor ? 3 : 4)) {
        if (eventlist != NULL) {
            if (g->length > 0)
                remove_uninformative(g);
            for (i = 1; i < g->n; i++) {
                e = (Event *)xmalloc(sizeof(Event));
                e->type = COALESCENCE;
                e->event.c.s1 = 0;
                e->event.c.s2 = 1;
                Enqueue(eventlist, e);
            }
        }
        /* Too few sequences to form segregating sites */
        return 1;
    }

    if (!ancestral_material_overlap(g)) {
        if (eventlist != NULL)
            for (i = 1; i < g->n; i++) {
                e = (Event *)xmalloc(sizeof(Event));
                e->type = COALESCENCE;
                e->event.c.s1 = 0;
                e->event.c.s2 = 1;
                Enqueue(eventlist, e);
            }
        /* Only one sequence carrying ancestral material in each site */
        return 1;
    }

    return 0;
}

/* Perform all mutations and coalesces that are guaranteed not to
 * increase the number of recombinations required. The data structure
 * g may be modified to reflect events.
 */
void force_safeevents(Genes *g)
{
    int change = 1, tmp;

#ifdef ENABLE_VERBOSE
    if (verbose()) printf("Forcing safe events:\n");
#endif
    while (change) {
        /* Look for mutations */
        change = force_mutations(g);
        /* Look for safe coalesces */
        tmp = coalesce_subsumed(g);
        change |= tmp;
    }
}

/* Force mutations everywhere possible, but do not remove columns as
 * is done in remove_uninformative. No events are registered in eventlist.
 */
int force_mutations(Genes *g)
{
    int i, j, n = 0, blocks = divblocksize(g->length - 1) + 1;
    unsigned long one0, two0, one1, two1, current;

    for (i = 0; i <= blocks - 1; i++) {
        /* Compile type information from all genes for block i */
        one1 = g->data[0].type[i];
        one0 = ~g->data[0].type[i] & g->data[0].ancestral[i];
        two0 = two1 = 0;
        for (j = 1; j < g->n; j++) {
            two1 |= one1 & g->data[j].type[i];
            one1 |= g->data[j].type[i];
            current = ~g->data[j].type[i] & g->data[j].ancestral[i];
            two0 |= one0 & current;
            one0 |= current;
        }
        /* two1 (two0) now has a 1 in all positions where there are at
         * least two 1s (0s) in the data.
         */
        /* Determine sites where there is exactly one 1, and where there
         * is exactly one 0 and more than one 1.
         */
        one1 &= ~two1;
        one0 &= ~two0 & two1;

        /* Force mutations */
        if (one1) {
            for (j = 0; j < g->n; j++)
                g->data[j].type[i] &= ~one1;
            n += weight(one1);
#ifdef ENABLE_VERBOSE
            if (verbose())
                printf("Mutated %d 1s to 0s in block %d\n", weight(one1), i);
#endif
        }
        if (one0 && gene_knownancestor) {
            for (j = 0; j < g->n; j++)
                g->data[j].type[i]
                    = (g->data[j].type[i] | one0) & g->data[j].ancestral[i];
            n += weight(one0);
#ifdef ENABLE_VERBOSE
            if (verbose())
                printf("Mutated %d 0s to 1s in block %d\n", weight(one0), i);
#endif
        }
    }

    return n;
}

/* For each mutation possible in g, force it and insert the resulting
 * set of sequences in an EList that is returned. A character can be
 * mutated if it is the only occurence of that type in its site and
 * there is at least one occurrence of the other type. If event is not
 * NULL the recombinations are appended to this EList in the same
 * order as the resulting configurations are stored in the EList
 * returned.
 */
EList *force_mutation(Genes *g, EList *event)
{
    int i, j, k, w0, w1, index, block, blocks = divblocksize(g->n - 1) + 1;
    Sites *s = genes2sites(g);
    Genes *h;
    EList *forced = elist_make();
    Event *e;

    for (i = 0; i < g->length; i++) {
        /* Check whether we can force a mutation in site i */
        w1 = weight(s->data[i].type[0]);
        w0 = weight(~s->data[i].type[0] & s->data[i].ancestral[0]);
        for (j = 1; j < blocks; j++) {
            w1 += weight(s->data[i].type[j]);
            w0 += weight(~s->data[i].type[j] & s->data[i].ancestral[j]);
        }
        if ((w1 == 1) && (gene_knownancestor || (w0 > 0))) {
            /* Exactly one occurrence of 1 in this site - force its mutation */
            index = modblocksize(i);
            block = divblocksize(i);
            h = copy_genes(g);
            /* Find sequence containing the 1 */
            for (j = 0; (k = lsb(s->data[i].type[j])) < 0; j++);
            /* Now force the mutation */
#ifdef ENABLE_VERBOSE
            if (verbose()) {
                printf("Mutating 1 in sequence %d, site %d to a 0\n",
                       k + mulblocksize(j), index + mulblocksize(block));
                output_genes_indexed(g, NULL);
            }
#endif
            h->data[k + mulblocksize(j)].type[block] &= ~((unsigned long)1 << index);
            elist_append(forced, h);
            if (event != NULL) {
                /* Insert corresponding event in list of events */
                e = (Event *)xmalloc(sizeof(Event));
                e->type = SUBSTITUTION;
                e->event.s.seq = k + mulblocksize(j);
                e->event.s.site = i;
                elist_append(event, e);
            }
        }
        if (!gene_knownancestor && (w0 == 1) && (w1 > 0)) {
            /* Exactly one occurrence of 1 in this site - force its mutation */
            index = modblocksize(i);
            block = divblocksize(i);
            h = copy_genes(g);
            /* Find sequence containing the 0 */
            for (j = 0; (k = lsb(~s->data[i].type[j] & s->data[i].ancestral[j])) < 0;
                    j++);
            /* Now force the mutation */
#ifdef ENABLE_VERBOSE
            if (verbose()) {
                printf("Mutating 0 in sequence %d, site %d to a 1\n",
                       k + mulblocksize(j), index + mulblocksize(block));
                output_genes_indexed(g, NULL);
            }
#endif
            h->data[k + mulblocksize(j)].type[block] |= ((unsigned long)1 << index);
            elist_append(forced, h);
            if (event != NULL) {
                /* Insert corresponding event in list of events */
                e = (Event *)xmalloc(sizeof(Event));
                e->type = SUBSTITUTION;
                e->event.s.seq = k + mulblocksize(j);
                e->event.s.site = i;
                elist_append(event, e);
            }
        }
    }

    /* Clean up */
    free_sites(s);

    return forced;
}

/* Find a character in position pos that can undergo a substitution,
 * substitute it, and return the sequence undergoing the
 * subtitution. If no substitution is possible in position pos -1 is
 * returned. If mutant is 0 or 1, it is assumed to be the mutant type
 * (i.e. the mutation will be of a site carrying type mutant); for any
 * other value of mutant, either type can be the mutant type (with a
 * preference for 1 being the mutant type). If the mutant type is known,
 * mutant has no effect.
 */
int mutate(Genes *g, int pos, int mutant)
{
    int i, block = divblocksize(pos), zero = -1, one = -1;
    unsigned long index = (unsigned long)1 << modblocksize(pos);

    if (gene_knownancestor)
        /* We know the mutant type */
        mutant = 1;

    /* Find substitution */
    for (i = 0; i < g->n; i++)
        if ((g->data[i].ancestral[block] & index) != 0) {
            /* Sequence i carries ancestral material in this position */
            if ((g->data[i].type[block] & index) != 0) {
                /* Sequence i has type 1 in this position */
                if (mutant != 0) {
                    /* And 1s can be mutated */
                    if (one != -1)
                        one = -2;
                    else
                        one = i;
                }
            }
            else if (mutant != 1) {
                /* Sequence i has type 0 in this position, and 0s can be mutated */
                if (zero != -1)
                    zero = -2;
                else
                    zero = i;
            }
        }

    /* Perform substitution */
    if (one >= 0) {
        g->data[one].type[block] &= ~index;
        return one;
    }
    else if (zero >= 0) {
        g->data[zero].type[block] |= index;
        return zero;
    }

    return -1;
}

/* Determines whether site i in sequence set g is segregating */
int segregating_site(Genes *g, int i)
{
    int j, block = divblocksize(i);
    unsigned long index = (unsigned long)1 << modblocksize(i);
    int k;

    if (!gene_knownancestor) {
        for (j = 0; j < g->n; j++)
            /* Find first occurrence of sequence carrying ancestral material
             * in this site.
             */
            if ((g->data[j].ancestral[block] & index) != 0) {
                /* Check whether all other sequences carrying ancestral material
                 * in this site agrees.
                 */
                for (k = j + 1; k < g->n; k++)
                    if (((g->data[j].type[block] ^ g->data[k].type[block])
                            & g->data[k].ancestral[block] & index) != 0)
                        /* Sequence k doesn't */
                        return 1;
            }
    }
    else
        for (j = 0; j < g->n; j++)
            if ((g->data[j].type[block] & index) != 0)
                /* Mutant type encountered */
                return 1;

    /* No evidence of segregation */
    return 0;
}

/* Determines whether there exists a pair of conflicting sites in g */
int conflicting_sites(Genes *g)
{
    Sites *s = genes2sites(g);
    unsigned long type00, type01, type10, type11, filter;
    int i, j, k;

    for (i = 0; i < s->length - 1; i++) {
        /* Find conflicting sites with site i as first site */
        for (j = i + 1; j < s->length; j++) {
            /* Compare sites i and j */
            type01 = type10 = type11 = 0;
            type00 = gene_knownancestor;
            for (k = 0; k < divblocksize(s->n - 1) + 1; k++) {
                filter = s->data[i].ancestral[k] & s->data[j].ancestral[k];
                type00 |= ~s->data[i].type[k] & ~s->data[j].type[k] & filter;
                type01 |= ~s->data[i].type[k] & s->data[j].type[k] & filter;
                type10 |= s->data[i].type[k] & ~s->data[j].type[k] & filter;
                type11 |= s->data[i].type[k] & s->data[j].type[k] & filter;
            }
            if (type00 && type01 && type10 && type11) {
                /* Sites are in conflict */
                free_sites(s);
                return 1;
            }
        }
    }

    free_sites(s);
    return 0;
}

/* Return index of first site after site i carrying ancestral material
 * in sequence a. If no such site exists, -1 is returned. Negative
 * values of i causes the search to start from the beginning of the
 * sequence.
 */
int next_ancestral(Genes *g, int a, int i)
{
    int block, blocks = divblocksize(g->length - 1) + 1;
    unsigned long pattern;

    if (i >= g->length - 1)
        return -1;
    else if (i < 0) {
        pattern = g->data[a].ancestral[0];
        block = 0;
    }
    else {
        block = divblocksize(i + 1);
        i = modblocksize(i + 1);
        if (i == 0)
            pattern = ~0;
        else
            pattern = ~(((unsigned long)1 << i) - 1);
        pattern &= g->data[a].ancestral[block];
    }

    while (((i = lsb(pattern)) == -1) && (++block < blocks))
        pattern = g->data[a].ancestral[block];

    if (i >= 0)
        return mulblocksize(block) + i;
    else
        return -1;
}

/* Determines whether sequences a and b are identical */
int identical(Genes *g, int a, int b)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    /* Check for identity */
    for (i = 0; i < blocks; i++)
        if (((g->data[a].ancestral[i] ^ g->data[b].ancestral[i]) != 0)
                || ((g->data[a].type[i] ^ g->data[b].type[i]) != 0))
            /* Conflict */
            return 0;

    /* No conflicts found */
    return 1;
}

/* Determines whether sequences a and b are compatible for coalescence */
int compatible(Genes *g, int a, int b)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    /* Check for compatibility */
    for (i = 0; i < blocks; i++)
        if ((g->data[a].ancestral[i] & g->data[b].ancestral[i]
                & (g->data[a].type[i] ^ g->data[b].type[i])) != 0)
            /* Conflict in site where both sequences have ancestral material */
            return 0;

    /* No conflicts found */
    return 1;
}

/* Returns list of sites where sequences a and b are incompatible */
EList *incompatible_sites(Genes *g, int a, int b)
{
    int i, j, k, blocks = divblocksize(g->length - 1) + 1;
    EList *l = elist_make();
    unsigned long incomp;

    for (i = 0; i < blocks; i++)
        /* Compute bit vector representation of incompatibilities in block i */
        if ((incomp = g->data[a].ancestral[i] & g->data[b].ancestral[i]
                      & (g->data[a].type[i] ^ g->data[b].type[i])) != 0) {
            /* Convert these positions to indeces and insert in list */
            /* First position is found in logarithmic time, then we do a linear
             * scan.
             */
            k = mulblocksize(i);
            j = lsb(incomp);
            elist_append(l, (void *)(k + j));
            incomp ^= (unsigned long)1 << j;
            for (j++; incomp != 0; j++)
                if ((((unsigned long)1 << j) & incomp) != 0) {
                    elist_append(l, (void *)(k + j));
                    incomp ^= (unsigned long)1 << j;
                }
        }

    return l;
}

/* Determines whether sequence a is subsumed in sequence b */
int subsumed_sequence(Genes *g, int a, int b)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    /* Check for subsumation */
    for (i = 0; i < blocks; i++)
        if (((~g->data[b].ancestral[i] | (g->data[a].type[i] ^ g->data[b].type[i]))
                & g->data[a].ancestral[i]) != 0)
            /* Either b does not carry ancestral material or there is a
             * conflict in a site where a carries ancestral material.
             */
            return 0;

    /* No conflicts found */
    return 1;
}

/* Determines whether site a is subsumed in site b */
int subsumed_site(Sites *s, int a, int b)
{
    int i, blocks = divblocksize(s->length - 1) + 1;

    /* Check for subsumation */
    for (i = 0; i < blocks; i++)
        if (((~s->data[b].ancestral[i] | (s->data[a].type[i] ^ s->data[b].type[i]))
                & s->data[a].ancestral[i]) != 0)
            /* Either b does not carry ancestral material or there is a
             * conflict in a sequence where a carries ancestral material.
             */
            return 0;

    /* No conflicts found */
    return 1;
}

/* Return list of indeces of ancestral sites in sequence a */
EList *ancestral_sites(Genes *g, int a)
{
    int i, block = 0;
    unsigned long index = 1;
    EList *sites = elist_make();

    for (i = 0; i < g->length; i++) {
        if (g->data[a].ancestral[block] & index)
            elist_append(sites, (void *)i);
        if ((index <<= 1) == 0) {
            /* End of block - continue with next block */
            block++;
            index = 1;
        }
    }

    return sites;
}

/* Return list of indeces of sequences carrying a 0 in site i */
EList *zero_sequences(Sites *s, int i)
{
    int h, j, k, blocks = divblocksize(s->n - 1) + 1;
    EList *l = elist_make();
    unsigned long zeros;

    for (h = 0; h < blocks; h++)
        /* Compute bit vector representation of sequences carrying 0 in block i */
        if ((zeros = ~s->data[i].type[h] & s->data[i].ancestral[h]) != 0) {
            /* Convert these positions to indeces and insert in list */
            /* First position is found in logarithmic time, then we do a linear
             * scan.
             */
            k = mulblocksize(h);
            j = lsb(zeros);
            elist_append(l, (void *)(k + j));
            zeros ^= (unsigned long)1 << j;
            for (j++; zeros != 0; j++)
                if ((((unsigned long)1 << j) & zeros) != 0) {
                    elist_append(l, (void *)(k + j));
                    zeros ^= (unsigned long)1 << j;
                }
        }

    return l;
}

/* Return list of indeces of sequences carrying a 1 in site i */
EList *one_sequences(Sites *s, int i)
{
    int h, j, k, blocks = divblocksize(s->n - 1) + 1;
    EList *l = elist_make();
    unsigned long ones;

    for (h = 0; h < blocks; h++)
        /* Copy bit vector representation of sequences carrying 1 in block i */
        if ((ones = s->data[i].type[h]) != 0) {
            /* Convert these positions to indeces and insert in list */
            /* First position is found in logarithmic time, then we do a linear
             * scan.
             */
            k = mulblocksize(h);
            j = lsb(ones);
            elist_append(l, (void *)(k + j));
            ones ^= (unsigned long)1 << j;
            for (j++; ones != 0; j++)
                if ((((unsigned long)1 << j) & ones) != 0) {
                    elist_append(l, (void *)(k + j));
                    ones ^= (unsigned long)1 << j;
                }
        }

    return l;
}

/* Determine whether b subsumes a in all sites given by segregating;
 * the data should consist of n blocks.
 */
static int subsumed_segregating(unsigned long *atype,
                                unsigned long *aancestral,
                                unsigned long *btype,
                                unsigned long *bancestral,
                                unsigned long *segregating, int n)
{
    int i;

    for (i = 0; i < n; i++)
        if ((segregating[i] & aancestral[i]
                & (~bancestral[i] | (atype[i] ^ btype[i]))) != 0)
            /* Non-subsumation detected */
            return 0;

    return 1;
}

/* Find a sequence that subsumes sequence a in all segregating sites;
 * if no such sequence exists, -1 is returned.
 */
int find_safe_coalescence(Genes *g, int a)
{
    int i, j, blocks = divblocksize(g->length - 1) + 1;
    unsigned long type, ancestral, *segregating
        = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));

    /* Determine set of segregating sites */
    for (i = 0; i < blocks; i++) {
        if(!gene_knownancestor) {
            ancestral = g->data[0].ancestral[i];
            type = g->data[0].type[i];
            for (j = 1; j < g->n; j++) {
                segregating[i] |= ancestral & g->data[j].ancestral[i] & (type ^ g->data[j].type[i]);
                type |= g->data[j].type[i];
                ancestral |= g->data[j].ancestral[i];
            }
        }
        else {
            // If know root state then just looking for sites that have at least one 1
            for(j = 0; j < g->n; j++) {
                segregating[i] |= g->data[j].ancestral[i] & g->data[j].type[i];
            }
        }
    }

    /* Find sequence subsuming a in all segregating sites */
    for (j = 0; j < g->n; j++)
        if ((j != a) &&
                (subsumed_segregating(g->data[a].type, g->data[a].ancestral,
                                      g->data[j].type, g->data[j].ancestral,
                                      segregating, blocks))) {
            free(segregating);
            return j;
        }

    free(segregating);
    return -1;
}

/* Determines whether a coalescence will entangle material from
 * sequences a and b. It is assumed that the sequences have positive
 * length.
 */
int entangled(Genes *g, int a, int b)
{
    int i, aflank, bflank, blocks = divblocksize(g->length - 1) + 1;

    /* Determine leftmost occurrence of ancestral material in
     * sequences a and b.
     */
    for (i = 0; i < blocks; i++) {
        aflank = lsb(g->data[a].ancestral[i]);
        if (aflank >= 0) {
            aflank += mulblocksize(i);
            break;
        }
    }
    for (i = 0; i < blocks; i++) {
        bflank = lsb(g->data[b].ancestral[i]);
        if (bflank >= 0) {
            bflank += mulblocksize(i);
            break;
        }
    }
    if (aflank == bflank)
        /* Entanglement */
        return 1;
    else if (aflank < bflank) {
        /* Leftmost occurrence of ancestral material in a is to the left
         * of leftmost occurrence of ancestral material in b; determine
         * rightmost occurrence of ancestral material in a.
         */
        for (i = blocks - 1; i >= 0; i--) {
            aflank = msb(g->data[a].ancestral[i]);
            if (aflank >= 0) {
                aflank += mulblocksize(i);
                break;
            }
        }
        if (aflank < bflank)
            /* No entanglement possible */
            return 0;
        else
            /* Entanglement */
            return 1;
    }
    else {
        /* Leftmost occurrence of ancestral material in b is to the left
         * of leftmost occurrence of ancestral material in a; determine
         * rightmost occurrence of ancestral material in b.
         */
        for (i = blocks - 1; i >= 0; i--) {
            bflank = msb(g->data[b].ancestral[i]);
            if (bflank >= 0) {
                bflank += mulblocksize(i);
                break;
            }
        }
        if (bflank < aflank)
            /* No entanglement possible */
            return 0;
        else
            /* Entanglement */
            return 1;
    }
}

/* Coalesce the two sequences a and b. The resulting sequence will be
 * left in a, while b will be replaced by the last sequence in the
 * data set. It is assumed that a and b are compatible.
 */
void coalesce(Genes *g, int a, int b)
{
    int i, blocks = divblocksize(g->length - 1) + 1, j;
    
    if(elements != NULL) {
        elist_swap(elements, b, elements->count - 1);
        elist_removelast(elements);
    }

    for (i = 0; i < blocks; i++) {
        /* Copy all type information from b to a */
        g->data[a].type[i] |= g->data[b].type[i];
        /* Coalescent will have ancestral material anywhere either a or b
         * has ancestral material.
         */
        g->data[a].ancestral[i] |= g->data[b].ancestral[i];
    }

    /* Clean up */
    free(g->data[b].type);
    free(g->data[b].ancestral);

    /* Update g */
    g->n -= 1;
    if (b != g->n)
        /* Copy last sequence to b */
        memcpy(&(g->data[b]), &(g->data[g->n]), sizeof(Gene));

}



/* For each pair of compatible sequences of g, coalesce them and
 * insert the resulting set of sequences in an EList that is returned.
 * If event is not NULL the coalesces are appended to this EList
 * in the same order as the resulting configurations are stored in the
 * EList returned.
 */
EList *force_coalesce(Genes *g, EList *event)
{
    int i, j;
    Genes *h;
    EList *forced = elist_make();
    Event *e;

    for (i = 0; i < g->n - 1; i++)
        for (j = i + 1; j < g->n; j++)
            if (compatible(g, i, j)) {
#ifdef ENABLE_VERBOSE
                if (verbose()) {
                    printf("Coalescing sequences %d and %d\n", i, j);
                    output_genes_indexed(g, NULL);
                }
#endif
                h = copy_genes(g);
                coalesce(h, i, j);
                elist_append(forced, h);
                if (event != NULL) {
                    /* Insert corresponding event in list of events */
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = COALESCENCE;
                    e->event.c.s1 = i;
                    e->event.c.s2 = j;
                    elist_append(event, e);
                }
            }

    return forced;
}

/* Computes the amount of ancestral material left after coalescing the
 * two compatible sequences a and b and performing all safe
 * evolutionary events. It runs implode on the data set obtained by
 * the coalescence and reports the amount of ancestral material
 * left.
 */
int coalescence_amleft(Genes *g, int a, int b)
{
    Genes *h = copy_genes(g);
    int n;
    LList *tmp = eventlist;
    eventlist = NULL;

    coalesce(h, a, b);
    implode_genes(h);
    n = ancestral_material(h);

    /* Clean up */
    free_genes(h);
    eventlist = tmp;

    return n;
}

/* Split sequence a into two sequences before site given by index and
 * block. The postfix is inserted as last sequence.
 */
static void _split(Genes *g, int a, int index, int block)
{
    int j, blocks = divblocksize(g->length - 1) + 1;
    unsigned long filter;

    g->data = (Gene *)xrealloc(g->data, (g->n + 1) * sizeof(Gene));
    g->data[g->n].type
        = (unsigned long*)xcalloc(blocks, sizeof(unsigned long));
    g->data[g->n].ancestral
        = (unsigned long*)xcalloc(blocks, sizeof(unsigned long));

    if (index != 0) {
        /* Split happens inside block */
        filter = ((unsigned long)1 << index) - 1;
        g->data[g->n].type[block] = g->data[a].type[block] & ~filter;
        g->data[g->n].ancestral[block] = g->data[a].ancestral[block] & ~filter;
        g->data[a].type[block] = g->data[a].type[block] & filter;
        g->data[a].ancestral[block] = g->data[a].ancestral[block] & filter;
        /* Now this block has been correctly initialised */
        j = block + 1;
    }
    else
        /* We need to copy block to the new sequence */
        j = block;

    for (; j < blocks; j++) {
        /* Split happened before this block */
        g->data[g->n].ancestral[j] = g->data[a].ancestral[j];
        g->data[g->n].type[j] = g->data[a].type[j];
        g->data[a].type[j] = g->data[a].ancestral[j] = 0;
    }

    /* Update gene to contain one more sequence */
    g->n += 1;
    
    if(elements != NULL) {
        if((int)(elist_get(elements, a)) != -1) {
            elist_append(elements, (void *)seq_numbering);
            seq_numbering++;
        }
        else {
            elist_append(elements, (void *)(-1));
            seq_numbering++;
        }
    }
}

/* Split sequence a into two sequences before site i. The postfix is
 * inserted as last sequence.
 */
void split(Genes *g, int a, int i)
{
    _split(g, a, modblocksize(i), divblocksize(i));
}

/* For each split in sequence a of g leading to a unique new set of
 * sequences, force it and insert the resulting set of sequences in an
 * EList that is returned. If event is not NULL the recombinations are
 * appended to this EList in the same order as the resulting
 * configurations are stored in the EList returned.
 */
EList *force_split(Genes *g, int a, EList *event)
{
    int i, index = 0, block = 0;
    Genes *h;
    EList *forced = elist_make();
    Event *e;

    /* No point in splitting before the first site carrying ancestral material */
    for (i = 1; i < g->length; i++) {
        if ((g->data[a].ancestral[block] & ((unsigned long)1 << index)) != 0)
            break;
        if (index == BLOCKSIZE - 1) {
            index = 0;
            block++;
        }
        else
            index++;
    }

    for (; i < g->length; i++) {
        /* Update split point */
        if (index == BLOCKSIZE - 1) {
            index = 0;
            block++;
        }
        else
            index++;
        /* We only split before site carrying ancestral material */
        if ((g->data[a].ancestral[block] & ((unsigned long)1 << index)) != 0) {
#ifdef ENABLE_VERBOSE
            if (verbose()) {
                printf("Splitting sequence %d before site %d\n",
                       a, index + mulblocksize(block));
                output_genes_indexed(g, NULL);
            }
#endif
            h = copy_genes(g);
            _split(h, a, index, block);
            elist_append(forced, h);
            if (event != NULL) {
                /* Insert corresponding event in list of events */
                e = (Event *)xmalloc(sizeof(Event));
                e->type = RECOMBINATION;
                e->event.r.seq = a;
                e->event.r.pos = index + mulblocksize(block);
                elist_append(event, e);
            }
        }
    }

    return forced;
}

/* Split sequence a before site given by index and block and
 * immediately coalesce prefix with sequence b. It is assumed that the
 * prefix is compatible with sequence b.
 */
void split_coalesceprefix(Genes *g, int a, int index, int block, int b)
{
    unsigned long filter;

    if (index != 0) {
        /* Merge in block that is split */
        filter = ((unsigned long)1 << index) - 1;
        g->data[b].type[block] |= g->data[a].type[block] & filter;
        g->data[b].ancestral[block] |= g->data[a].ancestral[block] & filter;
        g->data[a].type[block] &= ~filter;
        g->data[a].ancestral[block] &= ~filter;
    }

    /* Merge complete blocks of prefix with sequence b */
    for (block--; block >= 0; block--) {
        g->data[b].type[block] |= g->data[a].type[block];
        g->data[b].ancestral[block] |= g->data[a].ancestral[block];
        g->data[a].type[block] = g->data[a].ancestral[block] = 0;
    }
    
    if(elements != NULL) {
        elist_change(elements, b, (void *)(-1));
    }
}

/* Split sequence a before site given by index and block and
 * immediately coalesce postfix with sequence b. It is assumed that the
 * postfix is compatible with sequence b.
 */
void split_coalescepostfix(Genes *g, int a, int index, int block, int b)
{
    int blocks = divblocksize(g->length - 1) + 1;
    unsigned long filter;

    if (index != 0) {
        /* Merge in block that is split */
        filter = ((unsigned long)1 << index) - 1;
        g->data[b].type[block] |= g->data[a].type[block] & ~filter;
        g->data[b].ancestral[block] |= g->data[a].ancestral[block] & ~filter;
        g->data[a].type[block] &= filter;
        g->data[a].ancestral[block] &= filter;
        /* The block of the split has now been handled */
        block++;
    }

    /* Merge complete blocks of postfix with sequence b */
    for (; block < blocks; block++) {
        g->data[b].type[block] |= g->data[a].type[block];
        g->data[b].ancestral[block] |= g->data[a].ancestral[block];
        g->data[a].type[block] = g->data[a].ancestral[block] = 0;
    }
    
    if(elements != NULL) {
        elist_change(elements, b, (void *)(-1));
    }
}

/* Split sequence a after site given by index and block and
 * immediately coalesce postfix with sequence b. It is assumed that
 * the postfix is compatible with sequence b.
 */
void splitafter_coalescepostfix(Genes *g, int a, int index, int block, int b)
{
    int blocks = divblocksize(g->length - 1) + 1;
    unsigned long filter;

    if (index != BLOCKSIZE - 1) {
        /* Merge in block that is split */
        filter = ((unsigned long)2 << index) - 1;
        g->data[b].type[block] |= g->data[a].type[block] & ~filter;
        g->data[b].ancestral[block] |= g->data[a].ancestral[block] & ~filter;
        g->data[a].type[block] &= filter;
        g->data[a].ancestral[block] &= filter;
    }

    /* Merge complete blocks of postfix with sequence b */
    for (block++; block < blocks; block++) {
        g->data[b].type[block] |= g->data[a].type[block];
        g->data[b].ancestral[block] |= g->data[a].ancestral[block];
        g->data[a].type[block] = g->data[a].ancestral[block] = 0;
    }
    
    if(elements != NULL) {
        elist_change(elements, b, (void *)(-1));
    }
}

/* Split sequence a before site given by index and block and remove
 * prefix. It reflects the prefix being coalesced with some other
 * sequence in which it is subsumed.
 */
int split_removeprefix(Genes *g, int a, int index, int block)
{
    unsigned long filter;
    int blocks = divblocksize(g->length - 1) + 1;
    Genes *h = copy_genes(g);
    int b;
    int block2 = block;

    if (index != 0) {
        /* Remove prefix in block that is split */
        filter = ~(((unsigned long)1 << index) - 1);
        g->data[a].type[block] &= filter;
        g->data[a].ancestral[block] &= filter;
    }

    /* Remove complete blocks of prefix */
    for (block--; block >= 0; block--)
        g->data[a].type[block] = g->data[a].ancestral[block] = 0;
    
    // Now get the coalescing sequence and find index of the sequence into which it coalesces
    if (index != 0) {
        /* Remove postfix in block that is split */
        filter = ((unsigned long)1 << index) - 1;
        h->data[a].type[block2] &= filter;
        h->data[a].ancestral[block2] &= filter;
        /* The block of the split has now been handled */
        block2++;
    }
    
    /* Remove complete blocks of postfix */
    for (; block2 < blocks; block2++)
        h->data[a].type[block2] = h->data[a].ancestral[block2] = 0;
    
    b = find_safe_coalescence(h, a);
    if(elements != NULL) {
        elist_change(elements, b, (void *)(-1));
    }
    
    free_genes(h);
    
    return b;
    
}

/* Split sequence a before site given by index and block and remove
 * postfix. It reflects the postfix being coalesced with some other
 * sequence in which it is subsumed.
 */
int split_removepostfix(Genes *g, int a, int index, int block)
{
    int blocks = divblocksize(g->length - 1) + 1;
    unsigned long filter;
    Genes *h = copy_genes(g);
    int b;
    int block2 = block;

    if (index != 0) {
        /* Remove postfix in block that is split */
        filter = ((unsigned long)1 << index) - 1;
        g->data[a].type[block] &= filter;
        g->data[a].ancestral[block] &= filter;
        /* The block of the split has now been handled */
        block++;
    }

    /* Remove complete blocks of postfix */
    for (; block < blocks; block++)
        g->data[a].type[block] = g->data[a].ancestral[block] = 0;
    
    if (index != 0) {
        /* Remove prefix in block that is split */
        filter = ~(((unsigned long)1 << index) - 1);
        h->data[a].type[block2] &= filter;
        h->data[a].ancestral[block2] &= filter;
    }
    
    /* Remove complete blocks of prefix */
    for (block2--; block2 >= 0; block2--)
        h->data[a].type[block2] = h->data[a].ancestral[block2] = 0;
    
    b = find_safe_coalescence(h, a);
    if(elements != NULL) {
        elist_change(elements, b, (void *)(-1));
    }
    
    free_genes(h);
    
    return b;
}

/* Count the number of sites carrying ancestral material in
 * sequence i of g.
 */
int individual_ancestral_material(Genes *g, int i)
{
    int j, am = 0, blocks = divblocksize(g->length - 1) + 1;

    for (j = 0; j < blocks; j++)
        am += weight(g->data[i].ancestral[j]);

    return am;
}

/* Count the number of sites carrying ancestral material in the
 * sequences of g.
 */
int ancestral_material(Genes *g)
{
    int i, j, am = 0, blocks = divblocksize(g->length - 1) + 1;

    for (i = 0; i < g->n; i++)
        for (j = 0; j < blocks; j++)
            am += weight(g->data[i].ancestral[j]);

    return am;
}

/* Determine whether sequence i contains ancestral material in all sites */
int individual_all_ancestral(Genes *g, int i)
{
    int j, blocks = divblocksize(g->length - 1) + 1;

    for (j = 0; j < blocks - 1; j++)
        if (g->data[i].ancestral[j] != ~0)
            return 0;

    return (g->data[i].ancestral[blocks - 1] + 1
            == (unsigned long)1 << (modblocksize(g->length - 1) + 1));
}

/* Determine whether all sequences contain ancestral material in all sites */
int all_ancestral(Genes *g)
{
    int i;

    for (i = 0; i < g->n; i++)
        if (!individual_all_ancestral(g, i))
            return 0;
    return 1;
}

/* Count the number of sites carrying ancestral material in
 * informative sites in the sequences of g.
 */
int informative_ancestral_material(Genes *g)
{
    int i, j, am = 0, blocks = divblocksize(g->length - 1) + 1;
    unsigned long one0, two0, one1, two1, filter;

    for (i = 0; i <= blocks - 1; i++) {
        /* Compile type information from all genes for block i */
        one1 = g->data[0].type[i];
        one0 = ~g->data[0].type[i] & g->data[0].ancestral[i];
        two0 = two1 = 0;
        for (j = 1; j < g->n; j++) {
            two1 |= one1 & g->data[j].type[i];
            one1 |= g->data[j].type[i];
            filter = ~g->data[j].type[i] & g->data[j].ancestral[i];
            two0 |= one0 & filter;
            one0 |= filter;
        }
        /* two1 (two0) now has a 1 in all positions where there are at
         * least two 1s (0s) in the data.
         */
        if (gene_knownancestor)
            filter = two1;
        else
            filter = two0 & two1;

        /* Count ancestral material in informative sites for this block */
        for (j = 0; j < g->n; j++)
            am += weight(g->data[j].ancestral[i] & filter);
    }

    return am;
}

/* Determine whether there is any overlapping ancestral material for
 * the sequences in g.
 */
int ancestral_material_overlap(Genes *g)
{
    int i, j, blocks = divblocksize(g->length - 1) + 1;
    unsigned long pattern;

    /* Check block by block */
    for (i = 0; i < blocks; i++) {
        pattern = 0;
        /* ...and sequence by sequence */
        for (j = 0; j < g->n; j++) {
            if ((pattern & g->data[j].ancestral[i]) != 0)
                /* Overlap detected */
                return 1;
            pattern |= g->data[j].ancestral[i];
        }
    }

    return 0;
}

/* Determine whether sequences a and b in g share ancestral sites */
int pairwise_am_overlap(Genes *g, int a, int b)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    /* Check block by block */
    for (i = 0; i < blocks; i++)
        if ((g->data[a].ancestral[i] & g->data[b].ancestral[i]) != 0)
            /* Overlap detected */
            return 1;

    return 0;
}

/* Find maximum of the n numbers in a */
static int max(int n, int *a)
{
    int i, m;

    m = a[0];
    for (i = 1; i < n; i++)
        if (a[i] > m)
            m = a[i];

    return m;
}

/* Find the minimum number of segments we need to chop sequence a into
 * to be able to coalesce every segment with a compatible sequence. If
 * a contains a site that cannot be coalesced with any other sequence,
 * i.e. a site where all sequences have ancestral material and a is
 * unique, -1 is returned.
 */
int minimum_compatiblechops(Genes *g, int a)
{
    unsigned long filter, tmp;
    int i, j, k, chops = 0, blockmatched, blockswitch = 0,
                 *out = (int *)xcalloc(g->n, sizeof(int)),
                  blocks = divblocksize(g->length - 1) + 1;

    out[a] = -1;
    /* Run through the entire sequence */
    for (i = 0; i < blocks; i++) {
        /* Handle block i */
        filter = ~0;
        while (filter) {
            blockmatched = 0;
            for (j = 0; j < g->n; j++)
                if (!out[j]) {
                    /* Sequence j is still a contender */
                    k = lsb((g->data[a].type[i] ^ g->data[j].type[i])
                            & (g->data[a].ancestral[i] & g->data[j].ancestral[i])
                            & filter);
                    if (k >= 0)
                        /* But not anymore */
                        out[j] = k + 1;
                    else
                        /* Sequence j and a are compatible on remainder of block */
                        blockmatched = 1;
                }
                else
                    out[j] = -1;
            if (blockmatched)
                filter = 0;
            else {
                /* We could not match a with a compatible sequence to the end
                 * of the block.
                 */
                /* Update filter */
                tmp = ~(((unsigned long)1 << (max(g->n, out) - 1)) - 1);
                if (tmp == filter) {
                    /* We need to be a bit careful as the site that causes a
                     * chop can be the first site in a block. This means that
                     * even though tmp = filter, the current incompatibility
                     * does not match the previous incompatibility.
                     */
                    if (!blockswitch) {
                        /* We found a site that cannot be coalesced */
                        free(out);
                        return -1;
                    }
                    else
                        blockswitch = 0;
                }
                filter = tmp;
                /* Reset array of contenders */
                for (j = 0; j < a; j++)
                    out[j] = 0;
                for (j++; j < g->n; j++)
                    out[j] = 0;
                /* Increment number of chops */
                chops += 1;
            }
        }
        blockswitch = 1;
    }

    /* Clean up */
    free(out);

    return chops;
}

/* Find the minimum number of segment we need to chop sequence a into
 * to be able to coalesce every segment with a sequence it is subsumed
 * in. If a contains a site that cannot be coalesced with any other
 * sequence, i.e. a site where a is unique, -1 is returned.
 */
int minimum_subsumedchops(Genes *g, int a)
{
    unsigned long filter, tmp;
    int i, j, k, chops = 0, blockmatched, blockswitch = 0,
                 *out = (int *)xcalloc(g->n, sizeof(int)),
                  blocks = divblocksize(g->length - 1) + 1;

    out[a] = -1;
    /* Run through the entire sequence */
    for (i = 0; i < blocks; i++) {
        /* Handle block i */
        filter = ~0;
        while (filter) {
            blockmatched = 0;
            for (j = 0; j < g->n; j++)
                if (!out[j]) {
                    /* Sequence j is still a contender */
                    k = lsb((((g->data[a].type[i] ^ g->data[j].type[i])
                              & g->data[a].ancestral[i])
                             | (g->data[a].ancestral[i] & ~g->data[j].ancestral[i]))
                            & filter);
                    if (k >= 0)
                        /* But not anymore */
                        out[j] = k + 1;
                    else
                        /* Sequence a is subsumed in sequence j for remainder of
                         * the block.
                         */
                        blockmatched = 1;
                }
                else
                    out[j] = -1;
            if (blockmatched)
                filter = 0;
            else {
                /* We could not match a with a compatible sequence to the end
                 * of the block.
                 */
                /* Update filter */
                tmp = ~(((unsigned long)1 << (max(g->n, out) - 1)) - 1);
                if (tmp == filter) {
                    /* We need to be a bit careful as the site that causes a
                     * chop can be the first site in a block. This means that
                     * even though tmp = filter, the current incompatibility
                     * does not match the previous incompatibility.
                     */
                    if (!blockswitch) {
                        /* We found a site that cannot be coalesced */
                        free(out);
                        return -1;
                    }
                    else
                        blockswitch = 0;
                }
                filter = tmp;
                /* Reset array of contenders */
                for (j = 0; j < a; j++)
                    out[j] = 0;
                for (j++; j < g->n; j++)
                    out[j] = 0;
                /* Increment number of chops */
                chops += 1;
            }
        }
        blockswitch = 1;
    }

    /* Clean up */
    free(out);

    return chops;
}

/* Find the maximum prefix of each sequence that is subsumed in some
 * other sequence. The indeces returned are of the site immediately
 * following the subsumed prefix, i.e. if a sequence has no subsumed
 * prefix the index will be 0. If the entire sequence is subsumed in
 * another sequence the block is set to -1.
 */
Index *maximumsubsumedprefixs(Genes *g)
{
    Index *prefixs = (Index *)xmalloc(g->n * sizeof(Index));
    int a, i, j, k, prefixends, *out = (int *)xmalloc(g->n * sizeof(int)),
                                 blocks = divblocksize(g->length - 1) + 1;

    for (a = 0; a < g->n; a++) {
        /* Find prefix for sequence a */
        for (i = 0; i < g->n; i++)
            out[i] = -1;
        out[a] = -2;
        /* Run through the sequence from left to right */
        for (i = 0; i < blocks; i++) {
            /* Handle block i */
            prefixends = 1;
            for (j = 0; j < g->n; j++)
                if (out[j] == -1) {
                    /* Sequence j is still a contender */
                    k = lsb(((g->data[a].type[i] ^ g->data[j].type[i]) & g->data[a].ancestral[i])
                            | (g->data[a].ancestral[i] & ~g->data[j].ancestral[i]));
                    if (k >= 0)
                        /* But not anymore */
                        out[j] = k;
                    else
                        /* Sequence a is subsumed in sequence j for the remainder
                         * of the block.
                         */
                        prefixends = 0;
                }
                else
                    /* Register that sequence j didn't subsume sequence a all
                     * the way up to this block.
                     */
                    out[j] = -2;
            if (prefixends) {
                /* The last sequences subsuming a prefix have to give up in
                 * this block.
                 */
                prefixs[a].block = i;
                prefixs[a].index = max(g->n, out);
                break;
            }
        }
        if (!prefixends)
            /* Entire sequence is subsumed */
            prefixs[a].block = -1;
    }

    /* Clean up */
    free(out);

    return prefixs;
}

/* Find minimum of the n numbers in a - n is assumed to be larger than 0 */
static int min(int n, int *a)
{
    int i, m = a[0];

    for (i = 1; i < n; i++)
        if (a[i] < m)
            m = a[i];

    return m;
}

/* Find the maximum postfix of each sequence that is subsumed in some
 * other sequence. The indeces returned are of the first site in the
 * subsumed postfix, i.e. if a sequence has no subsumed postfix the
 * index will be equal to the sequence length. If the entire sequence
 * is subsumed in another sequence the block is set to -1.
 */
Index *maximumsubsumedpostfixs(Genes *g)
{
    Index *postfixs = (Index *)xmalloc(g->n * sizeof(Index));
    int a, i, j, k, postfixends, *out = (int *)xmalloc(g->n * sizeof(int)),
                                  blocks = divblocksize(g->length - 1) + 1;

    for (a = 0; a < g->n; a++) {
        /* Find postfix for sequence a */
        for (i = 0; i < g->n; i++)
            out[i] = 0;
        out[a] = INT_MAX;
        /* Run through the sequence from right to left */
        for (i = blocks - 1; i >= 0; i--) {
            /* Handle block i */
            postfixends = 1;
            for (j = 0; j < g->n; j++)
                if (!out[j]) {
                    /* Sequence j is still a contender */
                    k = msb(((g->data[a].type[i] ^ g->data[j].type[i])
                             & g->data[a].ancestral[i])
                            | (g->data[a].ancestral[i] & ~g->data[j].ancestral[i]));
                    if (k >= 0)
                        /* But not anymore */
                        out[j] = k + 1;
                    else
                        /* Sequence a is subsumed in sequence j for the remainder
                         * of the block.
                         */
                        postfixends = 0;
                }
                else
                    out[j] = INT_MAX;
            if (postfixends) {
                /* The last sequences subsuming a prefix have to give up in
                 * this block.
                 */
                postfixs[a].index = min(g->n, out);
                if (postfixs[a].index == BLOCKSIZE) {
                    /* The subsumed postfix does not extend into this block */
                    postfixs[a].index = 0;
                    postfixs[a].block = i + 1;
                }
                else
                    postfixs[a].block = i;
                break;
            }
        }
        if (!postfixends)
            /* Entire sequence is subsumed */
            postfixs[a].block = -1;
    }

    /* Clean up */
    free(out);

    return postfixs;
}

/* Find length of maximum prefix of all other sequences that is
 * subsumed in s.
 */
Index *maximumsubsumedprefix(Genes *g, int s)
{
    Index *prefixs = xmalloc(g->n * sizeof(Index));
    int i, j, blocks = divblocksize(g->n - 1) + 1;

    /* Run through all sequences */
    for (i = 0; i < g->n; i++) {
        /* Use -1 block to mark full compatibility */
        prefixs[i].block = -1;
        if (i != s)
            /* Look for leftmost insubsumation */
            for (j = 0; j < blocks; j++)
                if ((prefixs[i].index =
                            lsb(((g->data[s].type[j] ^ g->data[i].type[j]) &
                                 g->data[i].ancestral[j]) |
                                (g->data[i].ancestral[j] & ~g->data[s].ancestral[j]))) >= 0) {
                    prefixs[i].block = j;
                    break;
                }
    }

    return prefixs;
}

/* Find length of maximum postfix of all other sequences that
 * is subsumed in s.
 */
Index *maximumsubsumedpostfix(Genes *g, int s)
{
    Index *postfixs = xmalloc(g->n * sizeof(Index));
    int i, j, blocks = divblocksize(g->n - 1) + 1;

    /* Run through all sequences */
    for (i = 0; i < g->n; i++) {
        /* Use -1 block to mark full compatibility */
        postfixs[i].block = -1;
        if (i != s)
            /* Look for rightmost insubsumation */
            for (j = blocks - 1; j >= 0; j--)
                if ((postfixs[i].index =
                            msb(((g->data[s].type[j] ^ g->data[i].type[j]) &
                                 g->data[i].ancestral[j]) |
                                (g->data[i].ancestral[j] & ~g->data[s].ancestral[j]))) >= 0) {
                    if (postfixs[i].index == BLOCKSIZE - 1) {
                        postfixs[i].index = 0;
                        postfixs[i].block = j + 1;
                    }
                    else {
                        postfixs[i].index++;
                        postfixs[i].block = j;
                    }
                    break;
                }
    }

    return postfixs;
}

/* Find all prefixs of the sequences in g that are maximally
 * compatible with another sequence in g and perform the corresponding
 * split and coalescence. The two arrays of indeces are assumed to be
 * returned from the maximum subsumed (pre | post)fix functions on an
 * imploded gene, i.e. site a and the site to the left of b carry
 * ancestral material in s and b is smaller than the sequence
 * length. Each maximum subsumed prefix is also split off. The
 * function f is applied to the resulting HistoryFragments. It is the
 * responsibility of the calling function to free memory used for the
 * HistoryFragments.
 */
void maximal_prefix_coalesces_map(Genes *g, Index *a, Index *b,
                                  void (*f)(Genes *))
{
    int i, j, s, index, block, sindex, sblock,
        *ancestral = xmalloc(g->n * sizeof(int)),
         *out = xmalloc(g->n * sizeof(int));
    Genes *h;
    Event *e;
    LList *tmp = eventlist;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
#ifdef ENABLE_VERBOSE
    int v = verbose();

    set_verbose(0);
    if (v) {
        printf("Splitting off maximal compatible prefixes in:\n");
        output_genes_indexed(g, NULL);
    }
#endif

    /* Run through all sequences of g */
    for (s = 0; s < g->n; s++) {
        if (a[s].index + mulblocksize(a[s].block) > 0){
            /* Start by splitting off maximum subsumed prefix */
            h = copy_genes(g);
            if(elements != NULL) {
                elements = elist_make();
                elist_safeextend(elements, tmp_elements);
            }
            if(sites != NULL) {
                sites = elist_make();
                elist_safeextend(sites, tmp_sites);
            }
            if(split_removeprefix(h, s, a[s].index, a[s].block) != -1){
                if (eventlist != NULL) {
                    eventlist = MakeLList();
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = RECOMBINATION;
                    e->event.r.seq = s;
                    e->event.r.pos = a[s].index + mulblocksize(a[s].block);
                    Enqueue(eventlist, e);
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = SWAP;
                    e->event.swap.s1 = s;
                    e->event.swap.s2 = g->n;
                    Enqueue(eventlist, e);
                    e = (Event *)xmalloc(sizeof(Event));
                    e->type = COALESCENCE;
                    e->event.c.s1 = -1;
                    e->event.c.s2 = g->n;
                    Enqueue(eventlist, e);
                }
                implode_genes(h);
                f(h);
            }
            else {
                free_genes(h);
            }

    #ifdef ENABLE_VERBOSE
            if (v) {
                printf("Splitting maximum subsumed prefix off sequence %d at %d\n", s,
                    mulblocksize(a[s].block) + a[s].index);
                output_genes_indexed(h, NULL);
            }
    #endif
            /* Only consider splitting sequences where maximum subsumed prefix
            * and maximum subsumed postfix do not overlap.
            */
            if ((a[s].block < b[s].block) ||
                    ((a[s].block == b[s].block) && a[s].index < b[s].index)) {
                /* Determine remaining contenders */
                for (i = 0; i < g->n; i++) {
                    if ((i == s) ||
                            (((g->data[s].type[a[s].block] ^ g->data[i].type[a[s].block])
                            & g->data[s].ancestral[a[s].block]
                            & g->data[i].ancestral[a[s].block]
                            & (((unsigned long)2 << a[s].index) - 1))
                            != 0)) {
                        /* We will not consider prefixes between s and itself, nor
                        * between s and a sequence that has a conflict with s no
                        * later than a in the block containing a.
                        */
                        out[i] = 1;
                        continue;
                    }
                    out[i] = 0;
                    /* Nor will we consider prefixes between s and a sequence that
                    * has a conflict with s anywhere else before a.
                    */
                    for (j = 0; j < a[s].block; j++)
                        if (((g->data[s].type[j] ^ g->data[i].type[j])
                                & g->data[s].ancestral[j] & g->data[i].ancestral[j]) != 0) {
                            /* Conflict detected */
                            out[i] = 1;
                            break;
                        }
                }
                /* Convert out list to list of sequence indeces */
                j = 0;
                for (i = 0; i < g->n; i++)
                    if (!out[i])
                        out[j++] = i;
                if (j == 0)
                    /* No sequences remain */
                    continue;

                /* Determine for which of the contenders we have already
                * encountered some ancestral material.
                */
                /* First determine rightmost occurrence of ancestral material in s */
                sindex = -1;
                for (sblock = 0; sblock < a[s].block; sblock++)
                    if ((sindex = lsb(g->data[s].ancestral[sblock])) >= 0)
                        /* Found it */
                        break;
                if (sindex < 0) {
                    /* Haven't found it yet */
                    sindex = lsb(g->data[s].ancestral[sblock]);
                    /* We know that there is some ancestral before a (to end
                    * subsumation), so we can rest assured that sindex is well
                    * defined.
                    */
                    for (i = 0; i < j; i++)
                        ancestral[i] = g->data[out[i]].ancestral[sblock]
                                    & (((unsigned long)2 << a[s].index) - 1)
                                    & ~(((unsigned long)1 << sindex) - 1);
                }
                else {
                    /* Found it - first check block where maximum subsumed prefix
                    * of s is split off and block where s has leftmost ancestral
                    * material.
                    */
                    for (i = 0; i < j; i++) {
                        if ((g->data[out[i]].ancestral[a[s].block]
                                & (((unsigned long)2 << a[s].index) - 1)) != 0) {
                            /* Ancestral material encountered no later than a */
                            ancestral[i] = 1;
                            continue;
                        }
                        if ((g->data[out[i]].ancestral[sblock]
                                & ~(((unsigned long)1 << sindex) - 1)) != 0) {
                            /* Ancestral material encountered no sooner than sindex */
                            ancestral[i] = 1;
                            continue;
                        }
                        /* Then check blocks in between */
                        ancestral[i] = 0;
                        for (block = sblock + 1; block < a[s].block; block++)
                            if (g->data[out[i]].ancestral[block] != 0) {
                                /* Ancestral material encountered */
                                ancestral[i] = 1;
                                break;
                            }
                    }
                }

                /* Now determine all prefixes it makes sense to split off from s */
                index = a[s].index;
                block = a[s].block;
                do {
                    /* Update index */
                    if (index == BLOCKSIZE - 1) {
                        index = 0;
                        block++;
                    }
                    else
                        index++;

                    /* Check whether we have reached b */
                    if ((index == b[s].index) && (block == b[s].block))
                        break;

                    /* If s does not carry ancestral material at current index, we
                    * might as well postpone the split.
                    */
                    if ((((unsigned long)1 << index) & (g->data[s].ancestral[block])) != 0) {
                        /* Check remaining contenders for compatibility and maximal prefix */
                        for (i = 0; i < j;) {
                            if ((((unsigned long)1 << index) & (g->data[out[i]].ancestral[block])) == 0) {
                                /* Sequence i does not contain ancestral material in this
                                * site; split s before this site and coalesce with i.
                                */
                                if (ancestral[i]) {
                                    h = copy_genes(g);
                                    if(elements != NULL) {
                                        elements = elist_make();
                                        elist_safeextend(elements, tmp_elements);
                                    }
                                    if(sites != NULL) {
                                        sites = elist_make();
                                        elist_safeextend(sites, tmp_sites);
                                    }
                                    split_coalesceprefix(h, s, index, block, out[i]);
                                    if (eventlist != NULL) {
                                        eventlist = MakeLList();
                                        e = (Event *)xmalloc(sizeof(Event));
                                        e->type = RECOMBINATION;
                                        e->event.r.seq = s;
                                        e->event.r.pos = index + mulblocksize(block);
                                        Enqueue(eventlist, e);
                                        e = (Event *)xmalloc(sizeof(Event));
                                        e->type = SWAP;
                                        e->event.swap.s1 = s;
                                        e->event.swap.s2 = g->n;
                                        Enqueue(eventlist, e);
                                        e = (Event *)xmalloc(sizeof(Event));
                                        e->type = COALESCENCE;
                                        e->event.c.s1 = out[i];
                                        e->event.c.s2 = g->n;
                                        Enqueue(eventlist, e);
                                    }
                                    implode_genes(h);
                                    f(h);
                                    
    #ifdef ENABLE_VERBOSE
                                    if (v) {
                                        printf("Splitting off prefix at %d and coalescing with sequence %d\n", mulblocksize(block) + index, out[i]);
                                        output_genes_indexed(h, NULL);
                                    }
    #endif
                                }
                                i++;
                            }
                            else {
                                /* We now know that both s and i have ancestral material at
                                * this site.
                                */
                                if (((g->data[s].type[block] ^ g->data[out[i]].type[block])
                                        & ((unsigned long)1 << index)) != 0) {
                                    /* Sequence s and i are incompatible at this site; split
                                    * s before this site and coalesce with i.
                                    */
                                    if (ancestral[i]) {
                                        h = copy_genes(g);
                                        if(elements != NULL) {
                                            elements = elist_make();
                                            elist_safeextend(elements, tmp_elements);
                                        }
                                        if(sites != NULL) {
                                            sites = elist_make();
                                            elist_safeextend(sites, tmp_sites);
                                        }
                                        split_coalesceprefix(h, s, index, block, out[i]);
                                        if (eventlist != NULL) {
                                            eventlist = MakeLList();
                                            e = (Event *)xmalloc(sizeof(Event));
                                            e->type = RECOMBINATION;
                                            e->event.r.seq = s;
                                            e->event.r.pos = index + mulblocksize(block);
                                            Enqueue(eventlist, e);
                                            e = (Event *)xmalloc(sizeof(Event));
                                            e->type = SWAP;
                                            e->event.c.s1 = s;
                                            e->event.c.s2 = g->n;
                                            Enqueue(eventlist, e);
                                            e = (Event *)xmalloc(sizeof(Event));
                                            e->type = COALESCENCE;
                                            e->event.c.s1 = out[i];
                                            e->event.c.s2 = g->n;
                                            Enqueue(eventlist, e);
                                        }
                                        implode_genes(h);
                                        f(h);
                                        
    #ifdef ENABLE_VERBOSE
                                        if (v) {
                                            printf("Splitting off prefix at %d and coalescing with sequence %d\n", mulblocksize(block) + index, out[i]);
                                            output_genes_indexed(h, NULL);
                                        }
    #endif
                                    }
                                    /* Sequence i is no longer a contender */
                                    out[i] = out[--j];
                                    ancestral[i] = ancestral[j];
                                    continue;
                                }
                                /* Sequences s and i carry the same ancestral character at
                                * this site; postpone split but mark that ancestral
                                * material was encountered.
                                */
                                ancestral[i++] = 1;
                            }
                        }
                    }
                    else
                        /* Check for ancestral material in remaining contenders */
                        for (i = 0; i < j; i++)
                            ancestral[i] = ancestral[i] ||
                                        (g->data[out[i]].ancestral[block] & (unsigned long)1 << index) != 0;
                    /* Continue until we have eliminated all contenders */
                } while (j > 0);
            }
        }
    }

    /* Clean up */
    free(out);
    free(ancestral);
    eventlist = tmp;
    elements = tmp_elements;
    sites = tmp_sites;
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif
}

/* Find all prefixs of the sequences in g that are maximally
 * compatible with another sequence in g and perform the corresponding
 * split and coalescence. The two arrays of indeces are assumed to be
 * returned from the maximum subsumed (pre | post)fix functions on an
 * imploded gene, i.e. site a and the site to the left of b carry
 * ancestral material in s and b is smaller than the sequence
 * length. Each maximum subsumed prefix is also split off. The
 * resulting configurations are returned as a list of
 * HistoryFragments.
 */
static EList *_maximal_prefix_coalesces_list;
static void _maximal_prefix_coalesces_f(Genes *g)
{
    HistoryFragment *f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));

    f->g = g;
    f->event = eventlist;
    elist_append(_maximal_prefix_coalesces_list, f);
}
EList *maximal_prefix_coalesces(Genes *g, Index *a, Index *b)
{
    _maximal_prefix_coalesces_list = elist_make();
    maximal_prefix_coalesces_map(g, a, b, _maximal_prefix_coalesces_f);
    return _maximal_prefix_coalesces_list;
}

/* Find all postfixes of sequences in g that are maximally compatible
 * with another sequence in g and perform the corresponding split and
 * coalescence. The two arrays of indeces are assumed to be returned
 * from the maximum subsumed (pre | post)fix functions on an imploded
 * gene, i.e. site a and the site immediately to the left of b carry
 * ancestral material in s, a is larger than zero, and b is smaller
 * than the sequence length. The function f is applied to the
 * resulting HistoryFragments. It is the responsibility of the calling
 * function to free memory used for the HistoryFragments.
 */
void maximal_postfix_coalesces_map(Genes *g, Index *a, Index *b,
                                   void (*f)(Genes *))
{
    int i, j, k, s, index, block, sindex, sblock,
        *ancestral = xmalloc(g->n * sizeof(int)),
         *out = xmalloc(g->n * sizeof(int)),
          blocks = divblocksize(g->length - 1) + 1;
    Genes *h;
    Event *e;
    LList *tmp = eventlist;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
#ifdef ENABLE_VERBOSE
    int v = verbose();

    set_verbose(0);
    if (v) {
        printf("Splitting off maximal compatible postfixes in:\n");
        output_genes_indexed(g, NULL);
    }
#endif

    /* Run through all sequences of g */
    for (s = 0; s < g->n; s++) {
        /* Only consider splitting sequences where maximum subsumed prefix
         * and maximum subsumed postfix do not meet or overlap.
         */
        if (b[s].index + mulblocksize(b[s].block) < g->length) {
            if ((a[s].block < b[s].block) ||
                    ((a[s].block == b[s].block) && (a[s].index < b[s].index))) {
                /* Start by splitting off maximum subsumed postfix */
                h = copy_genes(g);
                if(elements != NULL) {
                    elements = elist_make();
                    elist_safeextend(elements, tmp_elements);
                }
                if(sites != NULL) {
                    sites = elist_make();
                    elist_safeextend(sites, tmp_sites);
                }
                if(split_removepostfix(h, s, b[s].index, b[s].block) != -1) {
                    if (eventlist != NULL) {
                        eventlist = MakeLList();
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = s;
                        e->event.r.pos = b[s].index + mulblocksize(b[s].block);
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = COALESCENCE;
                        e->event.c.s1 = -1;
                        e->event.c.s2 = g->n;
                        Enqueue(eventlist, e);
                    }
                    implode_genes(h);
                    f(h);
                }
                else {
                    free_genes(h);
                }
    #ifdef ENABLE_VERBOSE
                if (v) {
                    printf("Splitting maximum subsumed postfix off sequence %d at %d\n", s,
                        mulblocksize(a[s].block) + a[s].index);
                    output_genes_indexed(h, NULL);
                }
    #endif
                /* Determine index of site to the left of b, i.e. the first site
                * not part of a subsumed postfix for s.
                */
                if (b[s].index == 0) {
                    index = BLOCKSIZE - 1;
                    block = b[s].block - 1;
                }
                else {
                    index = b[s].index - 1;
                    block = b[s].block;
                }
                /* Only consider splitting sequences where maximum subsumed prefix
                * and maximum subsumed postfix are separated by at least two sites.
                */
                if ((a[s].block < block) ||
                        ((a[s].block == block) && a[s].index < index)) {

                    /* Determine remaining contenders */
                    for (i = 0; i < g->n; i++) {
                        if ((i == s)
                                || (((g->data[s].type[block] ^ g->data[i].type[block])
                                    & g->data[s].ancestral[block] & g->data[i].ancestral[block]
                                    & ~(((unsigned long)1 << index) - 1)) != 0)) {
                            /* We will not consider postfixes between s and itself, nor
                            * between s and a sequence that has a conflict with s to
                            * the right of the left neighbour of b in the block
                            * containing this neighbour.
                            */
                            out[i] = 1;
                            continue;
                        }
                        out[i] = 0;

                        /* Nor will we consider prefixes between s and a sequence that
                        * has a conflict with s anywhere else before a.
                        */
                        for (j = block + 1; j < blocks; j++)
                            if (((g->data[s].type[j] ^ g->data[i].type[j])
                                    & g->data[s].ancestral[j] & g->data[i].ancestral[j]) != 0) {
                                /* Conflict detected */
                                out[i] = 1;
                                break;
                            }
                    }
                    /* Convert out list to list of sequence indeces */
                    j = 0;
                    for (i = 0; i < g->n; i++)
                        if (!out[i])
                            out[j++] = i;
                    if (j == 0)
                        /* No sequences remain */
                        continue;

                    /* Determine for which of the contenders we have already
                    * encountered some ancestral material.
                    */
                    /* Start by finding rightmost occurrence of ancestral material in s */
                    sindex = -1;
                    for (sblock = blocks - 1; sblock > block; sblock--)
                        if ((sindex = msb(g->data[s].ancestral[sblock])) >= 0)
                            /* Found it */
                            break;
                    if (sindex < 0) {
                        /* Haven't found it yet */
                        sindex = msb(g->data[s].ancestral[sblock]);
                        /* We know that there is some ancestral after b (to end
                        * subsumation), so we can rest assured that sindex is well
                        * defined.
                        */
                        for (i = 0; i < j; i++)
                            if ((g->data[out[i]].ancestral[sblock]
                                    & ~(((unsigned long)1 << index) - 1) & (((unsigned long)2 << sindex) - 1)) != 0)
                                /* Ancestral material encountered no sooner than b and no
                                * later than sindex.
                                */
                                ancestral[i] = 1;
                            else
                                ancestral[i] = 0;
                    }
                    else {
                        /* Found it - first check block where maximum subsumed postfix
                        * of s is split off and block where s has rightmost ancestral
                        * material.
                        */
                        for (i = 0; i < j; i++) {
                            if ((g->data[out[i]].ancestral[block]
                                    & ~(((unsigned long)1 << index) - 1)) != 0) {
                                /* Ancestral material encountered no sooner than b */
                                ancestral[i] = 1;
                                continue;
                            }
                            if ((g->data[out[i]].ancestral[sblock]
                                    & (((unsigned long)2 << sindex) - 1)) != 0) {
                                /* Ancestral material encountered no later than sindex */
                                ancestral[i] = 1;
                                continue;
                            }
                            /* Then check blocks in between */
                            ancestral[i] = 0;
                            for (k = block + 1; k < sblock; k++)
                                if (g->data[out[i]].ancestral[k] != 0) {
                                    /* Ancestral material encountered */
                                    ancestral[i] = 1;
                                    break;
                                }
                        }
                    }

                    /* Now determine all prefixes it makes sense to split off from s */
                    do {
                        /* Update index */
                        if (index == 0) {
                            index = BLOCKSIZE - 1;
                            block--;
                        }
                        else
                            index--;
                        /* Check whether we have reached a */
                        if ((index == a[s].index) && (block == a[s].block))
                            break;

                        /* If s does not carry ancestral material at current index, we
                        * might as well postpone the split.
                        */
                        if ((((unsigned long)1 << index) & (g->data[s].ancestral[block])) != 0) {
                            /* Check remaining contenders for compatibility and
                            * maximal prefix.
                            */
                            for (i = 0; i < j;) {
                                if ((((unsigned long)1 << index) & (~g->data[out[i]].ancestral[block])) != 0) {
                                    /* Sequence i does not contain ancestral material in this
                                    * site; split s after this site and coalesce with i.
                                    */
                                    if (ancestral[i]) {
                                        h = copy_genes(g);
                                        if(elements != NULL) {
                                            elements = elist_make();
                                            elist_safeextend(elements, tmp_elements);
                                        }
                                        if(sites != NULL) {
                                            sites = elist_make();
                                            elist_safeextend(sites, tmp_sites);
                                        }
                                        splitafter_coalescepostfix(h, s, index, block, out[i]);
                                        if (eventlist != NULL) {
                                            eventlist = MakeLList();
                                            e = (Event *)xmalloc(sizeof(Event));
                                            e->type = RECOMBINATION;
                                            e->event.r.seq = s;
                                            e->event.r.pos = index + mulblocksize(block) + 1;
                                            Enqueue(eventlist, e);
                                            e = (Event *)xmalloc(sizeof(Event));
                                            e->type = COALESCENCE;
                                            e->event.c.s1 = out[i];
                                            e->event.c.s2 = g->n;
                                            Enqueue(eventlist, e);
                                        }
                                        implode_genes(h);
                                        f(h);
    #ifdef ENABLE_VERBOSE
                                        if (v) {
                                            printf("Splitting off postfix at %d and coalescing with sequence %d\n", mulblocksize(block) + index, out[i]);
                                            output_genes_indexed(h, NULL);
                                        }
    #endif
                                    }
                                    i++;
                                }
                                else {
                                    /* We now know that both s and i have ancestral material at
                                    * this site.
                                    */
                                    if (((g->data[s].type[block] ^ g->data[out[i]].type[block])
                                            & ((unsigned long)1 << index)) != 0) {
                                        /* Sequence s and i are incompatible at this site; split s
                                        * after this site and coalesce with i.
                                        */
                                        if (ancestral[i]) {
                                            h = copy_genes(g);
                                            if(elements != NULL) {
                                                elements = elist_make();
                                                elist_safeextend(elements, tmp_elements);
                                            }
                                            if(sites != NULL) {
                                                sites = elist_make();
                                                elist_safeextend(sites, tmp_sites);
                                            }
                                            splitafter_coalescepostfix(h, s, index, block, out[i]);
                                            if (eventlist != NULL) {
                                                eventlist = MakeLList();
                                                e = (Event *)xmalloc(sizeof(Event));
                                                e->type = RECOMBINATION;
                                                e->event.r.seq = s;
                                                e->event.r.pos = index + mulblocksize(block) + 1;
                                                Enqueue(eventlist, e);
                                                e = (Event *)xmalloc(sizeof(Event));
                                                e->type = COALESCENCE;
                                                e->event.c.s1 = out[i];
                                                e->event.c.s2 = g->n;
                                                Enqueue(eventlist, e);
                                            }
                                            implode_genes(h);
                                            f(h);
    #ifdef ENABLE_VERBOSE
                                            if (v) {
                                                printf("Splitting off postfix at %d and coalescing with sequence %d\n", mulblocksize(block) + index, out[i]);
                                                output_genes_indexed(h, NULL);
                                            }
    #endif
                                        }
                                        /* Sequence i is no longer a contender */
                                        out[i] = out[--j];
                                        ancestral[i] = ancestral[j];
                                    }
                                    else
                                        /* Sequences s and i carry the same ancestral character at
                                        * this site; postpone split but mark that ancestral
                                        * material was encountered.
                                        */
                                        ancestral[i++] = 1;
                                }
                            }
                        }
                        else
                            /* Check for ancestral material in remaining contenders */
                            for (i = 0; i < j; i++)
                                ancestral[i] = ancestral[i] ||
                                            (g->data[out[i]].ancestral[block] & (unsigned long)1 << index) != 0;
                        /* Continue until we have eliminated all contenders */
                    } while (j > 0);
                }
            }
        }
    }

    /* Clean up */
    free(out);
    free(ancestral);
    eventlist = tmp;
    elements = tmp_elements;
    sites = tmp_sites;
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif
}

/* Find all postfixes of sequences in g that are maximally compatible
 * with another sequence in g and perform the corresponding split and
 * coalescence. The two arrays of indeces are assumed to be returned
 * from the maximum subsumed (pre | post)fix functions on an imploded
 * gene, i.e. site a and the site immediately to the left of b carry
 * ancestral material in s, a is larger than zero, and b is smaller
 * than the sequence length.
 */
static EList *_maximal_postfix_coalesces_list;
static void _maximal_postfix_coalesces_f(Genes *g)
{
    HistoryFragment *f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));

    f->g = g;
    f->event = eventlist;
    elist_append(_maximal_postfix_coalesces_list, f);
}
EList *maximal_postfix_coalesces(Genes *g, Index *a, Index *b)
{
    _maximal_postfix_coalesces_list = elist_make();
    maximal_postfix_coalesces_map(g, a, b, _maximal_postfix_coalesces_f);
    return _maximal_postfix_coalesces_list;
}

/* Initialise subsumed to be those sequences that agree with the
 * sequence given by index and block in the site type and ancestral
 * (which are both of length blocks) are taken from. It is assumed
 * that the sequence given by index and block does contain ancestral
 * material in this site.
 */
static void initialise_subsumed(int index, int block, int blocks,
                                unsigned long *subsumed,
                                unsigned long *type,
                                unsigned long *ancestral)
{
    int i;

    if (type[block] & (unsigned long)1 << index)
        /* Central sequence has type 1 at this site */
        for (i = 0; i < blocks; i++)
            subsumed[i] = type[i];
    else
        /* Central sequence has type 0 at this site */
        for (i = 0; i < blocks; i++)
            subsumed[i] = ~type[i] & ancestral[i];

    /* The sequence given by index and block should not be included */
    subsumed[block] &= ~((unsigned long)1 << index);
}

/* Update subsumed to be those sequences that still agree with the
 * sequence given by index and block in the site type and ancestral
 * (which are both of length blocks) are taken from. It is checked
 * whether the sequence given by index and block does contain
 * ancestral material in this site. The return value is true if the
 * set of subsuming sequences does not become empty.
 */
static int update_subsumed(int index, int block, int blocks,
                           unsigned long *subsumed, unsigned long *type,
                           unsigned long *ancestral)
{
    int i, j;

    if (ancestral[block] & (unsigned long)1 << index) {
        /* Central sequence has ancestral material in current site */
        j = 0;
        if (type[block] & (unsigned long)1 << index)
            /* Central sequence has type 1 in current site */
            for (i = 0; i < blocks; i++) {
                if (subsumed[i] &= type[i])
                    j = 1;
            }
        else
            /* Central sequence has type 0 in current site */
            for (i = 0; i < blocks; i++) {
                if (subsumed[i] &= ~type[i] & ancestral[i])
                    j = 1;
            }
    }
    else
        j = 1;

    return j;
}

/* Initialise subsumed and compatible to reflect those sequences that
 * agree with and do not disagree with, respectively, the sequence
 * given by index and block in the site type and ancestral (which
 * should all be of length blocks) are taken from. It is assumed that
 * the sequence given by index and block has ancestral material in
 * this site.
 */
static void initialise_subsumedinterval(int index, int block, int blocks,
                                        unsigned long *subsumed,
                                        unsigned long *compatible,
                                        unsigned long *type,
                                        unsigned long *ancestral, Genes *g)
{
    int i;

    if (type[block] & (unsigned long)1 << index)
        /* Central sequence has type 1 at this site */
        for (i = 0; i < blocks; i++) {
            subsumed[i] = type[i];
            compatible[i] = type[i] | ~ancestral[i];
        }
    else
        /* Central sequence has type 0 at this site */
        for (i = 0; i < blocks; i++) {
            subsumed[i] = ~type[i] & ancestral[i];
            compatible[i] = ~type[i];
        }

    /* The sequence given by index and block should not be included */
    subsumed[block] &= ~((unsigned long)1 << index);
    compatible[block] &= ~((unsigned long)1 << index);
    /* Make sure to set bits not corresponding to a sequence in last
     * block to zeros.
     */
    i = modblocksize(g->n);
    if (i != 0)
        compatible[blocks - 1] &= ((unsigned long)1 << i) - 1;
}

/* If subsumed does not become empty, update subsumed and compatible
 * to reflect those sequences that agree with and do not disagree
 * with, respectively, the sequence given by index and block in the
 * site type and ancestral (which should all be of length blocks) are
 * taken from. The tmp array is used for temporary storage of the new
 * subsumed values. The return value is true if the update was
 * actually performed. It is checked whether the sequence given by
 * index and block has ancestral material in this site.
 */
static int extend_subsumedinterval(int index, int block, int blocks,
                                   unsigned long *subsumed,
                                   unsigned long *compatible,
                                   unsigned long *tmp,
                                   unsigned long *type,
                                   unsigned long *ancestral)
{
    int i, j;

    if (ancestral[block] & (unsigned long)1 << index) {
        /* Central sequence has ancestral material in this site */
        j = 0;
        if (type[block] & (unsigned long)1 << index) {
            /* Central sequence has type 1 at this site */
            /* Check whether we can extend interval to include current site */
            for (i = 0; i < blocks; i++)
                if ((tmp[i] = subsumed[i] & type[i]) != 0)
                    j = 1;
            if (j)
                /* We can extend interval to include current site */
                for (i = 0; i < blocks; i++) {
                    subsumed[i] = tmp[i];
                    compatible[i] &= type[i] | ~ancestral[i];
                }
        }
        else {
            /* Central sequence has type 0 at this site */
            /* Check whether we can extend interval to include current site */
            for (i = 0; i < blocks; i++)
                if ((tmp[i] = subsumed[i] & ~type[i] & ancestral[i]) != 0)
                    j = 1;
            if (j)
                /* We can extend interval to include current site */
                for (i = 0; i < blocks; i++) {
                    subsumed[i] = tmp[i];
                    compatible[i] &= ~type[i];
                }
        }
    }
    else
        j = 1;

    return j;
}

/* Similar to extend_subsumedinterval, but the arrays are updated even
 * if the set of subsumed sequences becomes empty.
 */
static int extend_subsumedinterval2(int index, int block, int blocks,
                                    unsigned long *subsumed,
                                    unsigned long *compatible,
                                    unsigned long *type,
                                    unsigned long *ancestral)
{
    int i, j;

    if (ancestral[block] & (unsigned long)1 << index) {
        /* Central sequence has ancestral material in this site */
        j = 0;
        if (type[block] & (unsigned long)1 << index) {
            /* Central sequence has type 1 at this site */
            /* Extend interval to include current site */
            for (i = 0; i < blocks; i++) {
                if (subsumed[i] &= type[i])
                    j = 1;
                compatible[i] &= type[i] | ~ancestral[i];
            }
        }
        else {
            /* Central sequence has type 0 at this site */
            /* Extend interval to include current site */
            for (i = 0; i < blocks; i++) {
                if (subsumed[i] &= ~type[i] & ancestral[i])
                    j = 1;
                compatible[i] &= ~type[i];
            }
        }
    }
    else
        j = 1;

    return j;
}

/* Initialise compatible2 to those sequences given by compatible that
 * do not disagree with the sequence given by index and block in the
 * site type and ancestral are taken from. All arrays should be of
 * length blocks. It is assumed that the sequence given by index and
 * block has ancestral material in the site. The return value is true
 * if compatible does not represent the empty set.
 */
static int initialise_compatibleinterval(int index, int block, int blocks,
        unsigned long *compatible,
        unsigned long *compatible2,
        unsigned long *type,
        unsigned long *ancestral)
{
    int i, j = 0;

    if (type[block] & (unsigned long)1 << index) {
        /* Central sequence has type 1 at this site */
        for (i = 0; i < blocks; i++)
            if ((compatible2[i] = compatible[i] & (type[i] | ~ancestral[i])) != 0)
                j = 1;
    }
    else {
        /* Central sequence has type 0 at this site */
        for (i = 0; i < blocks; i++)
            if ((compatible2[i] = compatible[i] & ~type[i]) != 0)
                j = 1;
    }

    return j;
}

/* Initialise maximal to those sequences that do not agree with the
 * sequence given by index and block in the site type and ancestral
 * are taken from, and that are in compatible. All arrays should be of
 * length blocks. It is checked whether the sequence given by index
 * and block has ancestral material in this site. The return value is
 * true if maximal does not indicate the empty set.
 */
static int check_maximality(int index, int block, int blocks,
                            unsigned long *maximal, unsigned long *compatible,
                            unsigned long *type, unsigned long *ancestral)
{
    int i, j = 0;

    if (ancestral[block] & (unsigned long)1 << index) {
        /* Central sequence has ancestral material in this site */
        if (type[block] & (unsigned long)1 << index) {
            /* Central sequence has type 1 in this site */
            for (i = 0; i < blocks; i++)
                if ((maximal[i] = compatible[i] & ~type[i]) != 0)
                    j = 1;
        }
        else {
            /* Central sequence has type 0 in this site */
            for (i = 0; i < blocks; i++)
                if ((maximal[i] = compatible[i] & (type[i] | ~ancestral[i])) != 0)
                    j = 1;
        }
    }

    return j;
}

/* Check whether any of the sequences in compatible are not subsumed in
 * the central sequence until this point, and remove those that are.
 */
static int check_leftflank(int blocks, unsigned long *compatible,
                           Index *prefixs, int index, int block)
{
    int i, j = 0, k;
    unsigned long u;

    for (i = 0; i < blocks; i++)
        if ((k = lsb(u = compatible[i])) >= 0) {
            u = u >> k;
            for (;;) {
                if ((prefixs[k + mulblocksize(i)].block > block) ||
                        ((prefixs[k + mulblocksize(i)].block == block)
                         && (prefixs[k + mulblocksize(i)].index >= index)))
                    /* Sequence k is subsumed - remove it from compatible */
                    compatible[i] &= ~((unsigned long)1 << k);
                else
                    /* We found one (more) that is not subsumed */
                    j = 1;
                /* If this was not the last sequence in compatible, find the
                 * next one.
                 */
                u = u >> 1;
                if (u == 0)
                    break;
                k++;
                while ((u & 1) == 0) {
                    u = u >> 1;
                    k++;
                }
            }
        }

    return j;
}

/* Check whether any of the sequences in compatible are not subsumed in
 * the central sequence after this point and remove those that are.
 */
static int check_rightflank(int blocks, unsigned long *compatible,
                            Index *postfixs, int index, int block)
{
    int i, j = 0, k;
    unsigned long u;

    for (i = 0; i < blocks; i++)
        if ((k = lsb(u = compatible[i])) >= 0) {
            u = u >> k;
            for (;;) {
                if ((postfixs[k + mulblocksize(i)].block < block) ||
                        ((postfixs[k + mulblocksize(i)].block == block)
                         && (postfixs[k + mulblocksize(i)].index <= index)))
                    /* Sequence k is subsumed - remove it from compatible */
                    compatible[i] &= ~((unsigned long)1 << k);
                else
                    /* We found one (more) that is not subsumed */
                    j = 1;
                /* If this was not the last sequence in compatible, find the
                 * next one.
                 */
                u = u >> 1;
                if (u == 0)
                    break;
                k++;
                while ((u & 1) == 0) {
                    u = u >> 1;
                    k++;
                }
            }
        }

    return j;
}

/* Split the sequence given by index and block as well as s in the two
 * places given by sindex and sblock and by k and coalesce the infix
 * part with all sequences given by maximal that do not agree with the
 * sequence given by index and block in the site type and ancestral
 * are taken from. The data sets thus obtained are inserted into
 * genes. Arrays should be of length blocks.
 */
static void perform_maximal_splits(int index, int block, int s, int blocks,
                                   int sindex, int sblock,
                                   unsigned long *maximal,
                                   unsigned long *type,
                                   unsigned long *ancestral, Genes *g, int k,
                                   void (*f)(Genes *))
{
    int i, j, eindex, eblock;
    unsigned long pattern;
    Genes *h;
    Event *e;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
#ifdef ENABLE_VERBOSE
    int v = verbose();

    set_verbose(0);
    if (v) {
        printf("Splitting off maximal compatible infixes in:\n");
        output_genes_indexed(g, NULL);
    }
#endif


    /* Check whether there are any right maximal intervals ending
     * immediately before k.
     */
    if (type[block] & ((unsigned long)1 << index)) {
        /* Central sequence has type 1 in site k */
        for (i = 0; i < blocks; i++)
            if ((pattern = maximal[i] & ~type[i]) != 0) {
                /* The sequences in this block we need to merge the current
                 * segment of the central sequence into are specified by pattern.
                 */
                eindex = modblocksize(k);
                eblock = divblocksize(k);
                j = lsb(pattern);
                pattern >>= j;
                j += mulblocksize(i);
                for (;;) {
                    /* Perform splits and coalesces with sequence j */
                    h = copy_genes(g);
                    if(elements != NULL) {
                        elements = elist_make();
                        elist_safeextend(elements, tmp_elements);
                    }
                    if(sites != NULL) {
                        sites = elist_make();
                        elist_safeextend(sites, tmp_sites);
                    }
                    _split(h, s, sindex, sblock);
                    split_coalesceprefix(h, g->n, eindex, eblock, j);
                    if (eventlist != NULL) {
                        eventlist = MakeLList();
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = s;
                        e->event.r.pos = sindex + mulblocksize(sblock);
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = g->n;
                        e->event.r.pos = k;
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = COALESCENCE;
                        e->event.c.s1 = j;
                        e->event.c.s2 = g->n;
                        Enqueue(eventlist, e);
                    }
                    implode_genes(h);
                    f(h);
#ifdef ENABLE_VERBOSE
                    if (v) {
                        printf("Splitting maximal compatible infix off of sequence %d at %d and %d and\ncoalescing with sequence %d\n", s, mulblocksize(sblock) + sindex,
                               k, j);
                        output_genes_indexed(h, NULL);
                    }
#endif
                    /* Find next sequence to coalesce with */
                    if ((pattern >>= 1) == 0)
                        /* We just handled last sequence in this block */
                        break;
                    j++;
                    while ((pattern & 1) == 0) {
                        pattern >>= 1;
                        j++;
                    }
                }
            }
    }
    else {
        /* Central sequence has type 0 in site k */
        for (i = 0; i < blocks; i++)
            if ((pattern = maximal[i] & (type[i] | ~ancestral[i]))
                    != 0) {
                /* The sequences in this block we need to merge the current
                 * segment of sequence i into is specified by pattern.
                 */
                eindex = modblocksize(k);
                eblock = divblocksize(k);
                j = lsb(pattern);
                pattern >>= j;
                j += mulblocksize(i);
                for (;;) {
                    /* Perform splits and coalesces with sequence j */
                    h = copy_genes(g);
                    if(elements != NULL) {
                        elements = elist_make();
                        elist_safeextend(elements, tmp_elements);
                    }
                    if(sites != NULL) {
                        sites = elist_make();
                        elist_safeextend(sites, tmp_sites);
                    }
                    _split(h, s, sindex, sblock);
                    split_coalesceprefix(h, g->n, eindex, eblock, j);
                    if (eventlist != NULL) {
                        eventlist = MakeLList();
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = s;
                        e->event.r.pos = sindex + mulblocksize(sblock);
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = g->n;
                        e->event.r.pos = k;
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = COALESCENCE;
                        e->event.c.s1 = j;
                        e->event.c.s2 = g->n;
                        Enqueue(eventlist, e);
                    }
                    implode_genes(h);
                    f(h);
#ifdef ENABLE_VERBOSE
                    if (v) {
                        printf("Splitting maximal compatible infix off of sequence %d at %d and %d and\ncoalescing with sequence %d\n", s, mulblocksize(sblock) + sindex,
                               k, j);
                        output_genes_indexed(h, NULL);
                    }
#endif
                    /* Find next sequence to coalesce with */
                    if ((pattern >>= 1) == 0)
                        /* We just handled last sequence in this block */
                        break;
                    j++;
                    while ((pattern & 1) == 0) {
                        pattern >>= 1;
                        j++;
                    }
                }
            }
    }
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif
    
    elements = tmp_elements;
    sites = tmp_sites;

}

/* Update compatible to include the site type and ancestral are taken
 * from, where compatible are those sequences not disagreeing with the
 * sequence given by index and block. All arrays should have length
 * blocks. It is assumed that the sequence given by index and block
 * has ancestral material in the current site. The return value is
 * true if the set of sequences indicated by compatible does not
 * become empty.
 */
static int extend_compatibleinterval(int index, int block, int blocks,
                                     unsigned long *compatible,
                                     unsigned long *type,
                                     unsigned long *ancestral)
{
    int i, j;

    /* Update compatible to reflect remaining sequences that are compatible
     * with the sequence given by index and block in the current interval.
     */
    j = 0;
    if (type[block] & ((unsigned long)1 << index)) {
        /* Central sequence has type 1 in current site */
        for (i = 0; i < blocks; i++)
            if (compatible[i] &= type[i] | ~ancestral[i])
                j = 1;
    }
    else {
        /* Central sequence has type 0 in current site */
        for (i = 0; i < blocks; i++)
            if (compatible[i] &= ~type[i])
                j = 1;
    }

    return j;
}

/* Extend interval starting at the site given by leftindex and
 * leftblock to the right as long as any sequences remain compatible,
 * and perform all maximal compatible splits of sequence i (that is
 * also given by index and block) along the way and insert the
 * resulting data sets in genes.
 */
static void find_compatibleintervals(int index, int block, int i, int blocks,
                                     int leftindex, int leftblock, int start,
                                     int end, unsigned long *compatible,
                                     Index *postfixs, Sites *s, Genes *g,
                                     void (*f)(Genes *))
{
    int c = 1;

    while (c && (++start < end)) {
        /* Check whether we have reached a point where the infix coalesce
         * is already handled by maximal_postfix_coalesces for any of the
         * sequences in compatible.
         */
        if (!check_rightflank(blocks, compatible, postfixs, modblocksize(start),
                              divblocksize(start)))
            break;
        if (s->data[start].ancestral[block] & ((unsigned long)1 << index)) {
            /* Sequence i has ancestral material in site start */
            perform_maximal_splits(index, block, i, blocks, leftindex,
                                   leftblock, compatible,
                                   s->data[start].type, s->data[start].ancestral,
                                   g, start, f);
            c = extend_compatibleinterval(index, block, blocks, compatible,
                                          s->data[start].type,
                                          s->data[start].ancestral);
        }
    }
}

/* Find all infixs of sequences in g that are maximally compatible
 * with another sequence in g and perform the corresponding splits and
 * coalesces. The two index arrays are assumed to contain return
 * values from the maximum subsumed (pre | post)fix functions on an
 * imploded gene, i.e. sites a and b carry ancestral material in s, a
 * is larger than zero, and b is smaller than the sequence length. The
 * function f is applied to all resulting HistoryFragments. It is the
 * responsibility of the calling function to free the memory used by
 * these HistoryFragments.
 */
void maximal_infix_coalesces_map(Genes *g, Index *a, Index *b,
                                 void (*f)(Genes *))
{
    int c, i, j, k, index, block, start, end, left, right, leftindex,
        leftblock, blocks = divblocksize(g->n - 1) + 1;
    Sites *s = genes2sites(g);
    unsigned long
    *subsumed = (unsigned long *)xmalloc(blocks * sizeof(unsigned long)),
     *compatible = (unsigned long *)xmalloc(blocks * sizeof(unsigned long)),
      *compatible2 = (unsigned long *)xmalloc(blocks * sizeof(unsigned long)),
       *leftmaximal = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    Genes *h;
    Index *prefixs = NULL, *postfixs = NULL;
    Event *e;
    LList *tmp = eventlist;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
#ifdef ENABLE_VERBOSE
    int v = verbose();

    set_verbose(0);
    if (v) {
        printf("Splitting off maximal compatible infixes in:\n");
        output_genes_indexed(g, NULL);
    }
#endif

    /* Run through the sequences in turn */
    for (i = 0; i < g->n; i++) {
        /* Convert i to index, block notation */
        index = modblocksize(i);
        block = divblocksize(i);
        /* Determine first and last index that is a candidate for a split */
        end = mulblocksize(b[i].block) + b[i].index;
        right = start = mulblocksize(a[i].block) + a[i].index;
        /* If we are using two cuts, the last cut should fall to the right
         * of the last cut when we are cutting a maximum subsumed prefix
         * twice.
         */
        initialise_subsumed(index, block, blocks, subsumed, s->data[right].type,
                            s->data[right].ancestral);
        c = 1;
#ifdef ENABLE_VERBOSE
        if (v > 1)
            printf("First internal subsumed infix in sequence %d is [%d,", i, right);
#endif
        while (c && (++right < end))
            c = update_subsumed(index, block, blocks, subsumed, s->data[right].type,
                                s->data[right].ancestral);
#ifdef ENABLE_VERBOSE
        if (v > 1)
            printf(" %d]\n", right - 1);
#endif

        /* The first interesting subsumed infix contains right - determine
         * this infix, and any following interesting infixes.
         */
        while (right < end) {
            /* Find maximal subsumed interval containing site right */
            /* The interval should contain the site immediately to the right
             * of the previous interval.
             */
            left = right;
            initialise_subsumedinterval(index, block, blocks, subsumed, compatible,
                                        s->data[left].type, s->data[left].ancestral,
                                        g);
            /* Start by finding left end point */
            do {
                left--;
                c = extend_subsumedinterval(index, block, blocks, subsumed, compatible,
                                            leftmaximal, s->data[left].type,
                                            s->data[left].ancestral);
            } while(c);
            /* left indicates the site immediately to the left of the
             * minimal subsumed interval that we are currently determining,
             * i.e. the first site not included in it.
             */
            /* Handle maximal compatible intervals starting between left end
             * points of previous and current minimal subsumed interval.
             */
            j = left;
            c = initialise_compatibleinterval(index, block, blocks, compatible,
                                              compatible2, s->data[left].type,
                                              s->data[left].ancestral);
            while (c && (j > start)) {
                /* Only continue as long as there are still some sequences
                 * with compatible overlap. Intervals with first split no
                 * later than start have already been handled.
                 */
                /* Update site, but remember previous site in index, block format */
                leftindex = modblocksize(j);
                leftblock = divblocksize(j);
                j--;
                /* If not already computed, compute maximum prefix subsumation
                 * in sequence i for all other sequences.
                 */
                if (prefixs == NULL)
                    prefixs = maximumsubsumedprefix(g, i);
                /* Check whether we have reached a point where the infix
                 * coalesce is already handled by maximal_prefix_coalesces for
                 * any of the sequences in compatible2.
                 */
                if (!check_leftflank(blocks, compatible2, prefixs, leftindex,
                                     leftblock))
                    break;
                /* Handle current site */
                if (check_maximality(index, block, blocks, leftmaximal, compatible2,
                                     s->data[j].type, s->data[j].ancestral)) {
                    /* There are some left maximal compatible intervals starting
                     * right after site j. Expand to the right to find all
                     * maximal intervals.
                     */
                    /* If not already computed, compute maximum postfix
                     * subsumation in sequence i for all other sequences.
                     */
                    if (postfixs == NULL)
                        postfixs = maximumsubsumedpostfix(g, i);
                    find_compatibleintervals(index, block, i, blocks, leftindex,
                                             leftblock, right, end, leftmaximal,
                                             postfixs, s, g, f);
                    c = extend_compatibleinterval(index, block, blocks, compatible2,
                                                  s->data[j].type, s->data[j].ancestral);
                }
            }

            /* Position left failed the subsumation test */
            left++;
            /* See how far we can extend the minimal subsumed interval to
             * the right.
             */
            c = 1;
            while (c && (++right < end))
                c = extend_subsumedinterval2(index, block, blocks, subsumed,
                                             compatible, s->data[right].type,
                                             s->data[right].ancestral);
            /* If subsumed interval extends all the way (in)to the maximum
             * subsumed postfix, we might as well split the sequence by
             * splitting off a maximum subsumed postfix twice.
             */
            if (right < end) {
                /* Start by splitting out subsumed infix */
                leftindex = modblocksize(left);
                leftblock = divblocksize(left);
                h = copy_genes(g);
                if(elements != NULL) {
                    elements = elist_make();
                    elist_safeextend(elements, tmp_elements);
                }
                if(sites != NULL) {
                    sites = elist_make();
                    elist_safeextend(sites, tmp_sites);
                }
                _split(h, i, leftindex, leftblock);
                if(split_removeprefix(h, g->n, modblocksize(right), divblocksize(right)) != -1){
                    if (eventlist != NULL) {
                        eventlist = MakeLList();
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = i;
                        e->event.r.pos = left;
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RECOMBINATION;
                        e->event.r.seq = g->n;
                        e->event.r.pos = right;
                        Enqueue(eventlist, e);
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = COALESCENCE;
                        e->event.c.s1 = -1;
                        e->event.c.s2 = g->n;
                        Enqueue(eventlist, e);
                    }
                    implode_genes(h);
                    f(h);
                }
                else{
                    free_genes(h);
                }
#ifdef ENABLE_VERBOSE
                if (v > 1) {
                    printf("Splitting off maximum subsumed infix of sequence %d between %d and %d\n", i, left, right);
                    set_verbose(v);
                }
#endif
                /* Now find compatible intervals starting at the same site as
                 * the subsumed infix but ending further to the right.
                 */
                elements = tmp_elements;
                sites = tmp_sites;
                for (k = 0; k < blocks; k++)
                    if (compatible[k]) {
                        if (postfixs == NULL)
                            postfixs = maximumsubsumedpostfix(g, i);
                        if(elements != NULL) {
                            elements = tmp_elements;
                        }
                        if(sites != NULL) {
                            sites = tmp_sites;
                        }
                        find_compatibleintervals(index, block, i, blocks, leftindex,
                                                 leftblock, right, end, compatible,
                                                 postfixs, s, g, f);
                        break;
                    }
                /* Prepare to look for next maximal subsumed interval */
                start = left;
#ifdef ENABLE_VERBOSE
                set_verbose(0);
#endif
            }
#ifdef ENABLE_VERBOSE
            else if (v > 1)
                printf("Last internal subsumed interval is [%d, %d]\n", left, right);
#endif
        }
        /* Clean up */
        if (prefixs != NULL) {
            free(prefixs);
            prefixs = NULL;
        }
        if (postfixs!= NULL) {
            free(postfixs);
            postfixs = NULL;
        }
    }

    /* Clean up */
    free(compatible);
    free(compatible2);
    free(leftmaximal);
    free_sites(s);
    free(subsumed);
    eventlist = tmp;
    elements = tmp_elements;
    sites = tmp_sites;
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif
    
}

/* Find all infixs of sequences in g that are maximally compatible
 * with another sequence in g and perform the corresponding splits and
 * coalesces. The two index arrays are assumed to contain return
 * values from the maximum subsumed (pre | post)fix functions on an
 * imploded gene, i.e. sites a and b carry ancestral material in s, a
 * is larger than zero, and b is smaller than the sequence length.
 */
static EList *_maximal_infix_coalesces_list;
static void _maximal_infix_coalesces_f(Genes *g)
{
    HistoryFragment *f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));

    f->g = g;
    f->event = eventlist;
    elist_append(_maximal_infix_coalesces_list, f);
}
EList *maximal_infix_coalesces(Genes *g, Index *a, Index *b)
{
    _maximal_infix_coalesces_list = elist_make();
    maximal_infix_coalesces_map(g, a, b, _maximal_infix_coalesces_f);
    return _maximal_infix_coalesces_list;
}

/* Check whether there are any incompatibilities between s1 and s2 in
 * the region between the site given by index1 and block1 and the site
 * given by b for sequence s2, whether an overlap would be right
 * maximal, and whether the prefix from s1 obtained from splitting at
 * index1, block1 overlaps any postfixes from s2 it makes senes to
 * split off. If no incompatibilities are found, the overlap is right
 * maximal, and there are overlaps between the prefix of s1 and
 * sensible postfixes of s2, initialise index2 and block2 to the
 * leftmost of index1 and block1 and the site given by b for sequence
 * s2 and return 1; if the split in s1 equals the left point of the
 * interval where it makes sense to split s2, return -1; otherwise
 * return 0.
 */
static int initialise_secondsplit(int s1, int index1, int block1, int s2,
                                  int *index2, int *block2, Index *a, Index *b,
                                  Genes *g)
{
    int i;

    /* Check whether split in s1 equals the left end point of the
     * interval where it makes sense to split s2.
     */
    if ((block1 == a[s2].block) && (index1 == a[s2].index))
        /* The prefix of s1 obtained from splitting and index1, block1
         * will not overlap any postfixes obtained from splitting s2 in
         * the interval where it makes sense to split s2.
         */
        return -1;

    /* Check whether an overlap would be right maximal */
    if ((((g->data[s1].type[block1] ^ g->data[s2].type[block1])
            | ~g->data[s2].ancestral[block1]) & ((unsigned long)1 << index1)) == 0)
        /* s1 and s2 have the same type in the site given by index1 and
         * block1, so we can extend overlap to the right.
         */
        return 0;

    /* Check for compatibility */
    if (b[s2].block < block1) {
        if (((g->data[s1].type[block1] ^ g->data[s2].type[block1])
                & (((unsigned long)1 << index1) - 1) & g->data[s1].ancestral[block1]
                & g->data[s2].ancestral[block1]) != 0)
            /* Incompatibility detected */
            return 0;
        for (i = block1 - 1; i > b[s2].block; i--)
            if (((g->data[s1].type[i] ^ g->data[s2].type[i])
                    & g->data[s1].ancestral[i] & g->data[s2].ancestral[i]) != 0)
                /* Incompatibility detected */
                return 0;
        if (((g->data[s1].type[b[s2].block] ^ g->data[s2].type[b[s2].block])
                & ~(((unsigned long)1 << b[s2].index) - 1) & g->data[s1].ancestral[b[s2].block]
                & g->data[s2].ancestral[b[s2].block]) != 0)
            /* Incompatibility detected */
            return 0;
        /* No incompatibilities - prepare to run through possible
         * splits in s2.
         */
        *index2 = b[s2].index;
        *block2 = b[s2].block;
    }
    else if ((b[s2].block == block1) && (b[s2].index < index1)) {
        if (((g->data[s1].type[block1] ^ g->data[s2].type[block1])
                & (((unsigned long)1 << index1) - 1) & ~(((unsigned long)1 << b[s2].index) - 1)
                & g->data[s1].ancestral[block1] & g->data[s2].ancestral[block1]) != 0)
            /* Incompatibility detected */
            return 0;
        /* No incompatibilities - prepare to run through possible
         * splits in s2.
         */
        *index2 = b[s2].index;
        *block2 = b[s2].block;
    }
    else {
        /* Current split in s1 falls within the interval where it
         * makes sense to split s2, so start from here.
         */
        *index2 = index1;
        *block2 = block1;
    }

    /* Make sure that the overlap entangles some ancestral material from
     * s1 and s2.
     */
    /* Find rightmost occurrence of ancestral material in s1 that is to
     * the left of the site given by index2, block2.
     */
    *index2 = msb(g->data[s1].ancestral[*block2] & (((unsigned long)1 << *index2) - 1));
    while ((*index2 == -1) && ((*block2)-- > 0))
        *index2 = msb(g->data[s1].ancestral[*block2]);
    if (*block2 < 0)
        /* No such occurrence */
        return -1;
    /* Now find rightmost occurrence of ancestral material in s2 that is
     * not to the right of the site given by index2, block2.
     */
    *index2 = msb(g->data[s2].ancestral[*block2] & (((unsigned long)2 << *index2) - 1));
    while ((*index2 == -1) && ((*block2)-- > 0))
        *index2 = msb(g->data[s2].ancestral[*block2]);
    /* Check whether we have moved the split point in s2 outside the
     * region where it makes sense to split s2.
     */
    if ((*block2 < a[s2].block) ||
            ((*block2 == a[s2].block) && (*index2 < a[s2].index)))
        return -1;
    /* Both s1 and s2 may have ancestral material in this site - if so,
     * it should be compatible.
     */
    if (((g->data[s1].type[*block2] ^ g->data[s2].type[*block2])
            & g->data[s1].ancestral[*block2] & g->data[s2].ancestral[*block2]
            & ((unsigned long)1 << *index2)) != 0)
        /* But it isn't */
        return 0;

    return 1;
}

/* Find all overlaps of sequences in g that are maximally compatible
 * and perform the corresponding splits and coalesces. The two index
 * arrays are assumed to contain return values from the maximum
 * subsumed (pre | post)fix functions on an imploded gene, i.e. sites
 * a and b carry ancestral material in s, a is larger than zero, and b
 * is smaller than the sequence length. The function f is applied to
 * the resulting HistoryFragments. It is the responsibility of the
 * calling function to free the memory used by these HistoryFragments.
 */
void maximal_overlap_coalesces_map(Genes *g, Index *a, Index *b,
                                   void (*f)(Genes *))
{
    int i, j, s1, s2, index1, block1, index2, block2,
        *in = (int *)xmalloc(g->n * sizeof(int));
    unsigned long pattern;
    Genes *h;
    Event *e;
    LList *tmp = eventlist;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
#ifdef ENABLE_VERBOSE
    int v = verbose();

    set_verbose(0);
    if (v) {
        printf("Splitting and coalescing maximal overlaps in:\n");
        output_genes_indexed(g, NULL);
    }
#endif

    /* Run through all pairs of sequences */
    for (s1 = 0; s1 < g->n; s1++)
        /* Only consider splits in sequences where the maximum subsumed
         * prefix does not overlap the maximum subsumed postfix.
         */
        if ((a[s1].block < b[s1].block)
                || ((a[s1].block == b[s1].block) && (a[s1].index < b[s1].index))) {
            /* s1 is the sequence providing the overlapping prefix */
            index1 = b[s1].index;
            block1 = b[s1].block;

            /* Create list of sequences we can ever form overlaps with */
            i = 0;
            for (s2 = 0; s2 < g->n; s2++)
                if ((s1 != s2) &&
                        /* s1 and s2 are not the same sequence */
                        ((block1 > a[s2].block) ||
                         ((block1 == a[s2].block) && (index1 > a[s2].index))) &&
                        /* Rightmost split in s1 is to the right of leftmost split in s2 */
                        ((b[s2].block > a[s2].block)
                         || ((b[s2].block == a[s2].block) && (b[s2].index > a[s2].index)))) {
                    /* The interval of sensible splits for s2 is not empty */
                    in[i++] = s2;
                }
            /* Run through possible split points for s1 */
            for(;;) {
                /* Update index1 */
                if (index1 == 0) {
                    index1 = BLOCKSIZE - 1;
                    block1--;
                }
                else
                    index1--;
                /* Check whether we have reached the end of the interval where
                 * it makes sense to split s1.
                 */
                if ((index1 == a[s1].index) && (block1 == a[s1].block))
                    break;
                /* Check whether s1 has ancestral material in this site - we
                 * only split before sites carrying ancestral material.
                 */
                if ((g->data[s1].ancestral[block1] & ((unsigned long)1 << index1)) == 0)
                    continue;
                /* Now iterate through the second sequence of the pair */
                for (s2 = 0; s2 < i;) {
                    /* s2 is the sequence providing the overlapping postfix */
                    j = initialise_secondsplit(s1, index1, block1, in[s2], &index2,
                                               &block2, a, b, g);
                    if (j == -1) {
                        /* We can form no more sensible overlaps between s1 and s2 */
                        in[s2] = in[--i];
                        continue;
                    }
                    if (j) {
                        /* s1 and s2 are compatible in the region back to the
                         * interval where it makes sense to split s2.
                         */
                        /* Create bit pattern of identities in current block */
                        pattern = (g->data[s1].type[block2] ^ g->data[in[s2]].type[block2])
                                  | ~g->data[s1].ancestral[block2];
                        for (;;) {
                            /* Update index2 */
                            if (index2 == 0) {
                                index2 = BLOCKSIZE - 1;
                                block2--;
                                pattern = (g->data[s1].type[block2]
                                           ^ g->data[in[s2]].type[block2])
                                          | ~g->data[s1].ancestral[block2];
                            }
                            else
                                index2--;
                            /* Check whether we have reached the end of the interval
                             * where it makes sense to split s2.
                             */
                            if ((block2 < a[in[s2]].block) ||
                                    ((block2 == a[in[s2]].block) && (index2 < a[in[s2]].index)))
                                break;
                            /* Check whether s2 has ancestral material in this site -
                             * we only split after sites carrying ancestral material.
                             */
                            if ((g->data[in[s2]].ancestral[block2] & ((unsigned long)1 << index2)) == 0)
                                continue;
                            /* Check whether overlap is left maximal */
                            if ((pattern & ((unsigned long)1 << index2)) == 0)
                                /* Overlap is not left maximal */
                                continue;
                            /* It makes sense to split s1 at index1, block1 and s2 at
                             * index2, block2, and s1 and s2 are compatible in the
                             * overlapping region.
                             */
                            h = copy_genes(g);
                            if(elements != NULL) {
                                elements = elist_make();
                                elist_safeextend(elements, tmp_elements);
                            }
                            if(sites != NULL) {
                                sites = elist_make();
                                elist_safeextend(sites, tmp_sites);
                            }
                            _split(h, s1, index1, block1);
                            splitafter_coalescepostfix(h, in[s2], index2, block2, s1);
                            if (eventlist != NULL) {
                                eventlist = MakeLList();
                                e = (Event *)xmalloc(sizeof(Event));
                                e->type = RECOMBINATION;
                                e->event.r.seq = s1;
                                e->event.r.pos = index1 + mulblocksize(block1);
                                Enqueue(eventlist, e);
                                e = (Event *)xmalloc(sizeof(Event));
                                e->type = RECOMBINATION;
                                e->event.r.seq = in[s2];
                                e->event.r.pos = index2 + mulblocksize(block2) + 1;
                                Enqueue(eventlist, e);
                                e = (Event *)xmalloc(sizeof(Event));
                                e->type = COALESCENCE;
                                e->event.c.s1 = s1;
                                e->event.c.s2 = g->n + 1;
                                Enqueue(eventlist, e);
                            }
                            implode_genes(h);
                            f(h);
#ifdef ENABLE_VERBOSE
                            if (v) {
                                printf("Splitting sequence %d at %d and sequence %d at %d and coalescing overlaps\n", s1, mulblocksize(block1) + index1, in[s2],
                                       mulblocksize(block2) + index2);
                                output_genes_indexed(h, NULL);
                            }
#endif
                            /* Check what caused the overlap to be maximal - if it
                             * is due to an incompatibility in the current site, we
                             * should not extend the overlap further.
                             */
                            if (g->data[s1].ancestral[block2] & ((unsigned long)1 << index2))
                                /* s1 has ancestral material in this site, so the
                                 * maximality must have been due to an
                                 * incompatibility.
                                 */
                                break;
                        }
                    }
                    /* Update sequence providing postfix of overlap */
                    s2++;
                }
            }
        }

    /* Clean up */
    free(in);
    eventlist = tmp;
    elements = tmp_elements;
    sites = tmp_sites;
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif
}

/* Find all overlaps of sequences in g that are maximally compatible
 * and perform the corresponding splits and coalesces. The two index
 * arrays are assumed to contain return values from the maximum
 * subsumed (pre | post)fix functions on an imploded gene, i.e. sites
 * a and b carry ancestral material in s, a is larger than zero, and b
 * is smaller than the sequence length. A list containing the
 * resulting HistoryFragments is returned.
 */
static EList *_maximal_overlap_coalesces_list;
static void _maximal_overlap_coalesces_f(Genes *g)
{
    HistoryFragment *f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));

    f->g = g;
    f->event = eventlist;
    elist_append(_maximal_overlap_coalesces_list, f);
}
EList *maximal_overlap_coalesces(Genes *g, Index *a, Index *b)
{
    _maximal_overlap_coalesces_list = elist_make();
    maximal_overlap_coalesces_map(g, a, b, _maximal_overlap_coalesces_f);
    return _maximal_overlap_coalesces_list;
}

typedef struct _HashGenesParameters {
    unsigned long m; /* Modulo */
    unsigned long a; /* Coefficient for number of sequences */
    unsigned long b; /* Coefficient for length of sequences */
    unsigned long c; /* Coefficient for type patterns */
    unsigned long d; /* Coefficient for ancestral patterns */
    unsigned long e; /* Coefficient for block */
} HashGenesParameters;

static HashGenesParameters *initialise_hashparameters(unsigned long m)
{
    HashGenesParameters *p = xmalloc(sizeof(HashGenesParameters));
    p->m = m;

    /* Initialise parameters */
    p->a = xrandom() % p->m;
    p->b = xrandom() % p->m;
    p->c = xrandom() % p->m;
    p->d = xrandom() % p->m;
    p->e = xrandom() % p->m;

    return p;
}

/* Compute hash value for g */
static unsigned long hash_genes(Genes *g, HashGenesParameters *p)
{
    int i, j;
    unsigned long pattern, value;

    /* Start with contributions from sequence number and length */
    value = (g->n * p->a + g->length * p->b) % p->m;
    for (i = 0; i < divblocksize(g->length - 1) + 1; i++) {
        /* Find contribution from block i */
        pattern = 0;
        for (j = 0; j < g->n; j++)
            pattern ^= p->c * (g->data[j].type[i] % p->m)
                       + p->d * (g->data[j].ancestral[i] % p->m)
                       + (g->data[j].type[i] * g->data[j].ancestral[i]) % p->m;
        value = (value * p->e + pattern) % p->m;
    }

    return value;
}

/* Compare sequences a and b in g and return -1, 0, or 1 depending on
 * whether sequence a is smaller than sequence b. The main purpose of
 * this function is to allow a sorting of sequences - the definition
 * of `smaller' does not carry information of e.g. subsumedness.
 */
int compare_sequences(Genes *g, int a, int b)
{
    int i, blocks = divblocksize(g->length - 1) + 1;

    for (i = 0; i < blocks; i++)
        if (g->data[a].type[i] < g->data[b].type[i])
            return -1;
        else if (g->data[a].type[i] > g->data[b].type[i])
            return 1;
        else if (g->data[a].ancestral[i] < g->data[b].ancestral[i])
            return -1;
        else if (g->data[a].ancestral[i] > g->data[b].ancestral[i])
            return 1;

    return 0;
}

/* Compare sites a and b in s and return -1, 0, or 1 depending on
 * whether site a is smaller than site b. The main purpose of this
 * function is to allow a sorting of sites - the definition of
 * `smaller' does not carry information of e.g. subsumedness.
 */
int compare_sites(Sites *s, int a, int b)
{
    int i, blocks = divblocksize(s->n - 1) + 1;

    for (i = 0; i < blocks; i++)
        if (s->data[a].type[i] < s->data[b].type[i])
            return -1;
        else if (s->data[a].type[i] > s->data[b].type[i])
            return 1;
        else if (s->data[a].ancestral[i] < s->data[b].ancestral[i])
            return -1;
        else if (s->data[a].ancestral[i] > s->data[b].ancestral[i])
            return 1;

    return 0;
}

/* Check whether the two genes g and h are identical */
static int compare_gene(int length, Gene *g, Gene *h)
{
    int i, blocks = divblocksize(length - 1) + 1;

    /* Check all blocks in turn */
    for (i = 0; i < blocks; i++) {
        if (g->type[i] ^ h->type[i])
            return 0;
        if (g->ancestral[i] ^ h->ancestral[i])
            return 0;
    }

    return 1;
}

/* Compare two data sets to determine whether they contain the same genes */
int compare_genes(Genes *g, Genes *h)
{
    int i, j, *companion, *count;

    if ((g->n != h->n) || (g->length != h->length))
        return 0;

    /* Rationale: first time we encounter a sequence in g we count the
     * number of times that sequence occurs in h; this count is then
     * decreased by one for every further encounter of this sequence in
     * g, and if the count reaches 0 we know that we have encountered it
     * more often in g than in h and return 0. If a sequence occurs in g
     * but not in h we can immediately return 0. Thus it is ensured that
     * no sequence occurs more often in g than in h. If a sequence
     * occurs less often in g than in h, as the two data sets contain
     * the same number of sequences some other sequence must occur more
     * often.
     */
    count = (int *)xcalloc(g->n, sizeof(int));
    companion = (int *)xcalloc(g->n, sizeof(int));
    for (i = 0; i < g->n; i++) {
        for (j = 0; j < g->n; j++)
            if (compare_gene(g->length, g->data + i, h->data + j)) {
                /* Sequence i in g and sequence j in h are identical */
                if (companion[j]) {
                    /* i is not the first occurrence of this sequence in g */
                    if (--count[companion[j] - 1] == 0) {
                        free(companion);
                        free(count);
                        return 0;
                    }
                    break;
                }
                else {
                    /* First time we encounter a sequence in g that is identical
                     * to sequence j in h.
                     */
                    companion[j] = i + 1;
                    count[i]++;
                }
            }
        if ((j >= g->n) && (count[i] == 0)) {
            /* The i'th gene of g does not have a match in h */
            free(companion);
            free(count);
            return 0;
        }
    }

    free(companion);
    free(count);
    return 1;
}

/* Determine if gene a is less than gene b in a global ordering on genes */
static int gene_less_than_blocks;
static int gene_less_than(Gene *a, Gene *b)
{
    int i;

    for (i = 0; i < gene_less_than_blocks; i++) {
        if (a->type[i] < b->type[i])
            return 1;
        if (a->type[i] > b->type[i])
            return 0;
        if (a->ancestral[i] < b->ancestral[i])
            return 1;
        if (a->ancestral[i] > b->ancestral[i])
            return 0;
    }

    /* The two genes are identical */
    return 0;
}

/* Create a more compact version of g */
PackedGenes *pack_genes(Genes *g)
{
    int i, j, k, gblock, gindex, pblock, pindex;
    unsigned int c;
    PackedGenes *p = (PackedGenes *)xmalloc(sizeof(PackedGenes));
    Genes *tmp = copy_genes(g);

    p->n = g->n;
    p->length = g->length;
    p->data = (unsigned int *)
              xcalloc((g->n * g->length - 1) / TERNARY_BLOCKSIZE + 1,
                      sizeof(unsigned int));

    /* Start by sorting the sequences */
    gene_less_than_blocks = divblocksize(g->length - 1) + 1;
    merge_sort(tmp->data, p->n, sizeof(Gene),
               (int (*)(void *, void *))gene_less_than);
    /* Now compact sequences one by one */
    pblock = pindex = 0;
    c = 1;
    for (i = 0; i < p->n; i++) {
        gblock = gindex = 0;
        for (j = 0; j < p->length; j++) {
            /* Convert type and ancestral values for site j in sequence i to
             * a ternary value.
             */
            k = ((tmp->data[i].type[gblock] >> gindex) & 1)
                + 2 - 2 * ((tmp->data[i].ancestral[gblock] >> gindex) & 1);
            p->data[pblock] += k * c;
            /* Update indeces in tmp and in p */
            if (gindex == BLOCKSIZE - 1) {
                gindex = 0;
                gblock += 1;
            }
            else
                gindex += 1;
            if (pindex == TERNARY_BLOCKSIZE - 1) {
                pindex = 0;
                pblock += 1;
                c = 1;
            }
            else {
                pindex += 1;
                c *= 3;
            }
        }
    }

    /* Clean up */
    free_genes(tmp);

    return p;
}
/* Create a less compact version of g */
Genes *unpack_genes(PackedGenes *p)
{
    int i, j, l, gindex, gblock, pindex, pblock,
        blocks = divblocksize(p->length - 1) + 1;
    unsigned int k;
    Genes *g = (Genes *)xmalloc(sizeof(Genes));

    /* Allocate structure for unpacked set of genes */
    g->n = p->n;
    g->length = p->length;
    g->data = (Gene *)xmalloc(g->n * sizeof(Gene));
    for (i = 0; i < g->n; i++) {
        g->data[i].type = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
        g->data[i].ancestral
            = (unsigned long *)xcalloc(blocks, sizeof(unsigned long));
    }

    /* Unpack p into the newly allocated structure */
    pblock = -1;
    pindex = TERNARY_BLOCKSIZE - 1;
    k = p->data[0];
    for (i = 0; i < g->n; i++) {
        gblock = gindex = 0;
        for (j = 0; j < g->length; j++) {
            /* Update index in p */
            if (pindex == TERNARY_BLOCKSIZE - 1) {
                pblock += 1;
                pindex = 0;
                k = p->data[pblock];
            }
            else {
                pindex += 1;
                k /= 3;
            }
            /* Determine type and ancestral values for site j in sequence i */
            l = k % 3;
            if (l == 1)
                g->data[i].type[gblock] |= (unsigned long)1 << gindex;
            if (l != 2)
                g->data[i].ancestral[gblock] |= (unsigned long)1 << gindex;
            /* Update index in g */
            if (gindex == BLOCKSIZE - 1) {
                gindex = 0;
                gblock += 1;
            }
            else
                gindex += 1;
        }
    }

    return g;
}

/* Deallocate a compact sequence set structure */
void free_packedgenes(PackedGenes *p)
{
    free(p->data);
    free(p);
}

typedef struct _HashPackedGenesParameters {
    unsigned long m; /* Modulo */
    unsigned long a; /* Coefficient */
} HashPackedGenesParameters;

static HashPackedGenesParameters *initialise_hashpackedparameters
(unsigned long m)
{
    static unsigned long i = 0;

    HashPackedGenesParameters *p = xmalloc(sizeof(HashPackedGenesParameters));
    p->m = m;

    /* Initialise parameters */
    p->a = (i == 0 ? (i = xrandom() % p->m) : i);

    return p;
}

/* Compute hash value for g */
static unsigned long hash_packedgenes(PackedGenes *g,
                                      HashPackedGenesParameters *p)
{
    int i;
    unsigned long value;

    /* Start with contributions from sequence number and length */
    value = (g->n * p->a % p->m + g->length) % p->m;
    /* Add contributions from the actual sequences */
    for (i = 0; i <= (g->n * g->length - 1) / TERNARY_BLOCKSIZE; i++)
        value = (value * p->a + g->data[i]) % p->m;

    return value;
}

/* Compare two packed data sets to determine whether they contain the
 * same genes.
 */
int compare_packedgenes(PackedGenes *g, PackedGenes *h)
{
    int i;

    /* First compare sequence counts and lengths */
    if ((g->n != h->n) || (g->length != h->length))
        return 0;

    /* Packaging involves sorting, so to sets are identical if and only
     * if they yield identical packages.
     */
    for (i = 0; i <= (g->n * g->length - 1) / TERNARY_BLOCKSIZE; i++)
        if (g->data[i] != h->data[i])
            return 0;

    /* No disagreements found */
    return 1;
}

HashTable *new_geneshashtable(int bits)
{
    if (bits > 20)
        bits = 20;
    if (bits < 3)
        bits = 3;
    return hashtable_new(bits, (unsigned long (*)(void *, void *))hash_genes,
                         (int (*)(void *, void *))compare_genes,
                         (void *(*)(unsigned long))initialise_hashparameters);
}

void init_geneshashtable(HashTable *table, int bits)
{
    if (bits > 20)
        bits = 20;
    if (bits < 3)
        bits = 3;
    hashtable_init(bits, table, (unsigned long (*)(void *, void *))hash_genes,
                   (int (*)(void *, void *))compare_genes,
                   (void *(*)(unsigned long))initialise_hashparameters);
}

HashTable *new_packedgeneshashtable(int bits)
{
    if (bits > 20)
        bits = 20;
    if (bits < 3)
        bits = 3;
    return hashtable_new(bits,
                         (unsigned long (*)(void *, void *))hash_packedgenes,
                         (int (*)(void *, void *))compare_packedgenes,
                         (void *(*)(unsigned long))
                         initialise_hashpackedparameters);
}

void init_packedgeneshashtable(HashTable *table, int bits)
{
    if (bits > 20)
        bits = 20;
    if (bits < 3)
        bits = 3;
    hashtable_init(bits, table,
                   (unsigned long (*)(void *, void *))hash_packedgenes,
                   (int (*)(void *, void *))compare_packedgenes,
                   (void *(*)(unsigned long))
                   initialise_hashpackedparameters);
}


/* Function to try all possible flips of sequencing errors
 */
void seqerror_flips(Genes* g, void (*f)(Genes *)) {
    int q, s, m;
    Genes *h;
    Event *e;
    LList *tmp = eventlist;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
    char c;
    
    for (q = 0; q < g->n; q++) {
        // Check that the sequence has not previously coalesced with anything
        if((int)(elist_get(tmp_elements, q)) != -1) {
            for(s = 0; s < g->length; s++) {     
                _recombinations = se_cost;
                // Get the "multiplicity" of the site (how many columns have been collapsed into it)
                m = (int)(elist_get(tmp_sites, s));
                if(m < 0) {
                    _recombinations = se_cost * (-m);
                    no_events = -m;
                }
                c = get_genes_character(g, q, s);
                // Check that the site is ancestral, if so flip and store
                if(c != 2) {
                    h = copy_genes(g);
                    if(c == 0){
                        set_genes_character(h, q, s, 1);
                    } else {
                        set_genes_character(h, q, s, 0);
                    }
                    if(elements != NULL) {
                        elements = elist_make();
                        elist_safeextend(elements, tmp_elements);
                    }
                    if(sites != NULL) {
                        sites = elist_make();
                        elist_safeextend(sites, tmp_sites);
                    }
                    if (eventlist != NULL) {
                        eventlist = MakeLList();
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = SEFLIP;
                        e->event.flip.seq = q;
                        e->event.flip.site = s;
                        Enqueue(eventlist, e);
                    }
                    implode_genes(h);
                    f(h);
                }
            }
        }
    }

    no_events = 1;
    _recombinations = se_cost;
    eventlist = tmp;
    elements = tmp_elements;
    sites = tmp_sites;
}

/* Function to try all possible flips of recurrent mutations
 */
void recmut_flips(Genes* g, void (*f)(Genes *)) {
    int q, s, m;
    Genes *h;
    Event *e;
    LList *tmp = eventlist;
    EList *tmp_elements = elements;
    EList *tmp_sites = sites;
    char c;
    
    for (q = 0; q < g->n; q++) {
        // Check that the sequence has previously coalesced with something
        if((int)(elist_get(tmp_elements, q)) == -1) {
            for(s = 0; s < g->length; s++) {     
                _recombinations = rm_cost;
                // Get the "multiplicity" of the site (how many columns have been collapsed into it)
                m = (int)(elist_get(tmp_sites, s));
                if(m < 0) {
                    _recombinations = rm_cost * (-m);
                    no_events = -m;
                }
                c = get_genes_character(g, q, s);
                // Check that the site is ancestral, if so flip and store
                if(c != 2) {
                    h = copy_genes(g);
                    if(c == 0){
                        set_genes_character(h, q, s, 1);
                    } else {
                        set_genes_character(h, q, s, 0);
                    }
                    if(elements != NULL) {
                        elements = elist_make();
                        elist_safeextend(elements, tmp_elements);
                    }
                    if(sites != NULL) {
                        sites = elist_make();
                        elist_safeextend(sites, tmp_sites);
                    }
                    if (eventlist != NULL) {
                        eventlist = MakeLList();
                        e = (Event *)xmalloc(sizeof(Event));
                        e->type = RMFLIP;
                        e->event.flip.seq = q;
                        e->event.flip.site = s;
                        Enqueue(eventlist, e);
                    }
                    implode_genes(h);
                    f(h);
                }
            }
        }
    }
    
    _recombinations = rm_cost;
    no_events = 1;
    eventlist = tmp;
    elements = tmp_elements;
    sites = tmp_sites;
}
