/*******************************************************************
 * 
 *    backtrack.c
 *  
 *    Implementation of functions for backtracking a minimum
 *    recombination history
 * 
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "backtrack.h"
#include "arg.h"
#include "common.h"
#include "gene.h"

/* Return edge indicated by edge */
static ARGEdge *getedge(ARG *arg, int edge)
{
    ARGNode *node = arg_getnode(arg, edge >> 1);
    
    if (node->type == ARGRECOMBINATION)
        if (edge & 1)
            /* Suffix edge */
            return &node->predecessor.two.suffix;
        else
            /* Prefix edge */
            return &node->predecessor.two.prefix;
        else
            return &node->predecessor.one;
}

/* Create an evolutionary history from the sequence of events
 * resulting in a minimum number of recombinations. The history is
 * returned as an ancestral recombination graph (ARG). If output is
 * not NULL, the history is printed to output.
 */
ARG *eventlist2history(AnnotatedGenes *a, FILE *output)
{
    int i, j, k, l, n = a->g->n, next_seq = a->g->n, *edges, n_se = 0, n_rm = 0, n_re = 0;
    LList *positions, *sequences, *tmp;
    LListCounter *lcounter, *lpos, *lseq, *lc;
    Genes *g = NULL, *h = copy_genes(a->g);
    Event *e, *f;
    char *pfix, *s, *t, **sequence, c;
    ARG *arg = NULL;
    ARGNode *node;
    ARGEdge *edge;
    PackedGenes *p;
    #ifdef DEBUG
    Genes *old;
    #endif
    
    if ((eventlist != NULL) && (Length(eventlist) > 0)){
        /* Determine number of nodes in ARG */
        lcounter = MakeCounter(eventlist, FIRST);
        while ((e = (Event *)Next(lcounter)) != NULL)
            if (e->type == RECOMBINATION)
                n++;
        /* Create initial ARG information */
        arg = arg_new();
        sequence = genes2string(a->g);
        sequence = (char **)xrealloc(sequence, n * sizeof(char *));
        for (i = a->g->n; i < n; i++)
            sequence[i] = (char *)xmalloc((a->g->length + 1) * sizeof(char));
        edges = (int *)xmalloc(n * sizeof(int));
        /* Create nodes corresponding to the sampled sequences */
        if (a->sequences != NULL){
            /* Sampled sequences are labelled */
            lc = MakeCounter(a->sequences, FIRST);
            for (i = 0; i < a->g->n; i++){
                s = (char *)Next(lc);
                if (s != NULL){
                    t = (char *)xmalloc((strlen(s) + 1) * sizeof(char));
                    strcpy(t, s);
                }
                else
                    t = NULL;
                s = (char *)xmalloc((a->g->length + 1) * sizeof(char));
                strcpy(s, sequence[i]);
                edges[i] = arg_addnode(arg, ARGSAMPLE, t, s) << 1;
            }
            DestroyCounter(lc);
        }
        else
            /* Sampled sequences are not labelled - use index as label */
            for (i = 0; i < a->g->n; i++){
                s = i2a(i + 1);
                t = (char *)xmalloc((a->g->length + 1) * sizeof(char));
                strcpy(t, sequence[i]);
                edges[i] = arg_addnode(arg, ARGSAMPLE, s, t) << 1;
            }
            for (; i < n; i++)
                edges[i] = -1;
            
            /* Construct llists for mapping from sequence and site numbers in
             * working gene set copy to sequence and site numbers in original
             * gene set.
             */
            sequences = MakeLList();
            positions = MakeLList();
            for (i = 0; i < h->length; i++)
                if (segregating_site(h, i)){
                    tmp = MakeLList();
                    Push(tmp, (void *)i);
                    Enqueue(positions, (void *)tmp);
                }
                for (i = 0; i < h->n; i++)
                    Enqueue(sequences, (void *)i);
            lpos = MakeCounter(positions, FIRST);
            lseq = MakeCounter(sequences, FIRST);
            lc = MakeCounter(positions, FIRST);
            InitCounter(lcounter, eventlist, FIRST);
            
            /* Go through the events recorded */
            while ((e = (Event *)Next(lcounter)) != NULL){
                #ifdef DEBUG
                old = copy_genes(h);
                #endif
                /* If a reduced state was remembered from a REMOVE in the
                 * previous round, delete it unless this event is also a
                 * REMOVE.
                 */
                if ((g != NULL) && (e->type != REMOVE)){
                    free_genes(g);
                    g = NULL;
                }
                
                switch(e->type){
                    case SWAP:
                        /* Sanity check - all SWAPs should be handled under
                         * RECOMBINATION events.
                         */
                        fprintf(stderr, "[SWAP:] Error in backtracking - please email error report\n");
                        exit(1);
                        break;
                    case COLLAPSE:
                        /* Next event was a collapse of neighbouring siamese twins;
                         * merge the two corresponding position lists.
                         */
                        SetCounter(lpos, e->event.collapse + 1);
                        tmp = (LList *)RemoveMoveLeft(positions, lpos);
                        Append((LList *)Current(lpos), tmp);
                        break;
                    case SUBSTITUTION:
                        /* Next event is a mutation of a singleton */
                        SetCounter(lpos, e->event.s.site);
                        tmp = RemoveMoveRight(positions, lpos);
                        /* Make sure to carry out the mutation in all original
                         * positions that has been collapsed into this position in the
                         * working copy.
                         */
                        while (Length(tmp) > 0){
                            i = (int)Pop(tmp);
                            if ((j = mutate(h, i, -1)) < 0){
                                /* Sanity check - a mutation should be possible in this site */
                                fprintf(stderr, "[SUBSTITUTION:] Error in backtracking - please email error report\n");
                                exit(1);
                            }
                            j = (int)SetCounter(lseq, j);
                            if (output != NULL){
                                /* Print mutation to output */
                                if (a->positions != NULL)
                                    s = (char *)GetByIndex(a->positions, i);
                                else
                                    s = i2a(i + 1);
                                t = NULL;
                                if ((a->sequences != NULL) && (Length(a->sequences) > j))
                                    t = (char *)GetByIndex(a->sequences, j);
                                if (t != NULL) {
                                    if(howverbose != -1) {
                                        fprintf(output, "Mutation of site %s in sequence %s\n", s, t);
                                    }
                                }
                                else {
                                    if(howverbose != -1) {
                                        fprintf(output, "Mutation of site %s in sequence %d\n", s, j + 1);
                                    }
                                }
                                if (a->positions == NULL)
                                        free(s);
                            }
                            /* Update ARG */
                            if (sequence[j][i] == '0')
                                sequence[j][i] = '1';
                            else
                                sequence[j][i] = '0';
                            edge = getedge(arg, edges[j]);
                            Enqueue(edge->mutations, (void *)(i));
                        }
//                         output_genes(h, output, "New:\n");
                        /* Position has now been removed from working copy */
                        DestroyLList(tmp);
                        break;
                    case SEFLIP:
                        i = (int)SetCounter(lseq, e->event.flip.seq); //sequence
                        tmp = SetCounter(lpos, e->event.flip.site); //get all the collapsed sites
                        if(Length(tmp) > 1) {
                            if(howverbose != -1 && output != NULL) {
                                fprintf(output, "---->Stretch of sequencing errors spanning %d sites:\n", Length(tmp));
                                if(howverbose > 0 && output != stdout && output != NULL) {
                                    printf("---->Stretch of sequencing errors spanning %d sites:\n", Length(tmp));
                                }
                            }
                        }
                        for(k = 0; k < Length(tmp); k++) {
                            n_se++;
                            j = (int)GetByIndex(tmp, k);
                            if (output != NULL){
                                /* Print mutation to output */
                                if (a->positions != NULL)
                                    s = (char *)GetByIndex(a->positions, j);
                                else
                                    s = i2a(j + 1);
                                t = NULL;
                                if ((a->sequences != NULL) && (Length(a->sequences) > i))
                                    t = (char *)GetByIndex(a->sequences, i);
                                if (t != NULL) {
                                    if(howverbose != -1) {
                                        fprintf(output, "---->");
                                        fprintf(output, "Sequencing error at site %s in sequence %s\n", s, t);
                                        if(howverbose > 0 && output != stdout) {
                                            printf("Sequencing error at site %s in sequence %s\n", s, t);
                                        }
                                    }
                                }
                                else {
                                    if(howverbose != -1) {
                                        fprintf(output, "---->");
                                        fprintf(output, "Sequencing error at site %s in sequence %d\n", s, i + 1);
                                        if(howverbose > 0 && output != stdout) {
                                            printf("Sequencing error at site %s in sequence %d\n", s, i + 1);
                                        }
                                    }
                                }
                                if (a->positions == NULL)
                                    free(s);
                            }
                            /* Update ARG */
                            if (sequence[i][j] == '0') {
                                sequence[i][j] = '1';
                                set_genes_character(h, e->event.flip.seq, j, 1);
                            }
                            else if (sequence[i][j] == '1'){
                                sequence[i][j] = '0';
                                set_genes_character(h, e->event.flip.seq, j, 0);
                            }
                            edge = getedge(arg, edges[i]);
                            if(j == 0) {
                                Enqueue(edge->mutations, (void *)(-INT_MAX));
                            }
                            else {
                                Enqueue(edge->mutations, (void *)(-j));
                            }
                        }
                        break;
                    case RMFLIP:
                        i = (int)SetCounter(lseq, e->event.flip.seq); //sequence
                        tmp = SetCounter(lpos, e->event.flip.site); //get all the collapsed sites
                        if(Length(tmp) > 1) {
                            if(howverbose != -1 && output != NULL) {
                                fprintf(output, "---->Stretch of recurrent mutations spanning %d sites:\n", Length(tmp));
                                if(howverbose > 0 && output != stdout && output != NULL) {
                                    printf("---->Stretch of recurrent mutations spanning %d sites:\n", Length(tmp));
                                }
                            }
                        }
                        for(k = 0; k < Length(tmp); k++) {
                            n_rm++;
                            j = (int)GetByIndex(tmp, k);
                            if (output != NULL){
                                /* Print mutation to output */
                                if (a->positions != NULL)
                                    s = (char *)GetByIndex(a->positions, j);
                                else
                                    s = i2a(j + 1);
                                t = NULL;
                                if ((a->sequences != NULL) && (Length(a->sequences) > i))
                                    t = (char *)GetByIndex(a->sequences, i);
                                if (t != NULL) {
                                    if(howverbose != -1) {
                                        fprintf(output, "---->");
                                        fprintf(output, "Recurrent mutation at site %s in sequence %s\n", s, t);
                                        if(howverbose > 0 && output != stdout && output != NULL) {
                                            printf("Recurrent mutation at site %s in sequence %s\n", s, t);
                                        }
                                    }
                                }
                                else {
                                    if(howverbose != -1) {
                                        fprintf(output, "---->");
                                        fprintf(output, "Recurrent mutation at site %s in sequence %d\n", s, i + 1);
                                        if(howverbose > 0 && output != stdout && output != NULL) {
                                            printf("Recurrent mutation at site %s in sequence %d\n", s, i + 1);
                                        }
                                    }
                                }
                                if (a->positions == NULL)
                                    free(s);
                            }
                            /* Update ARG */
                            if (sequence[i][j] == '0') {
                                sequence[i][j] = '1';
                                set_genes_character(h, e->event.flip.seq, j, 1);
                            }
                            else if (sequence[i][j] == '1') {
                                sequence[i][j] = '0';
                                set_genes_character(h, e->event.flip.seq, j, 0);
                            }
                            edge = getedge(arg, edges[i]);
                            if(j == 0) {
                                Enqueue(edge->mutations, (void *)(-INT_MAX));
                            }
                            else {
                                Enqueue(edge->mutations, (void *)(-j));
                            }
                        }
                        break;
                    case COALESCENCE:
                        /* Next event is a coalescence of two compatible sequences */
                        if (e->event.c.s1 == -1){
                            /* When coalescing two sequences where one is subsumed in
                             * the other we did not always know the sequence subsuming
                             * the other; find this, remembering that a removal of
                             * siamese twins may alter the picture.
                             */
                            g = copy_genes(h);
                            tmp = eventlist;
                            eventlist = NULL;
                            remove_nonsegregating(g);
                            if (g->length > 0){
                                remove_siamesetwins(g);
                                e->event.c.s1 = find_safe_coalescence(g, e->event.c.s2);
                            }
                            else
                                /* We should never get in this situation - if all
                                 * informative columns are removed, the remaining events
                                 * should be a sequences of coalescences with the first
                                 * sequence and the next remaining sequence.
                                 */
                                e->event.c.s1 = (e->event.c.s2 == 0 ? 1 : 0);
                            eventlist = tmp;
                            free_genes(g);
                            g = NULL;
                        }
                        i = (int)SetCounter(lseq, e->event.c.s1);
                        j = (int)SetCounter(lseq, e->event.c.s2);
                        if (output != NULL){
                            /* Print coalescence to output */
                            s = NULL;
                            t = NULL;
                            if (a->sequences != NULL){
                                if (Length(a->sequences) > i)
                                    s = (char *)GetByIndex(a->sequences, i);
                                if (Length(a->sequences) > j)
                                    t = (char *)GetByIndex(a->sequences, j);
                            }
                            if (s != NULL) {
                                if(howverbose != -1) {
                                    fprintf(output, "Coalescing sequences %s", s);
                                }
                            }
                            else {
                                if(howverbose != -1) {
                                    fprintf(output, "Coalescing sequences %d", i + 1);
                                }
                            }
                            if (t != NULL) {
                                if(howverbose != -1) {
                                    fprintf(output, " and %s\n", t);
                                }
                            }
                            else {
                                if(howverbose != -1) {
                                    fprintf(output, " and %d\n", j + 1);
                                }
                            }
                        }
                        /* Update ARG */
                        for (k = 0; k < a->g->length; k++)
                            if (sequence[i][k] == 'x')
                                sequence[i][k] = sequence[j][k];
                        s = (char *)xmalloc((a->g->length + 1) * sizeof(char));
                        strcpy(s, sequence[i]);
                        k = arg_addnode(arg, ARGCOALESCENCE, NULL, s);
                        edge = getedge(arg, edges[i]);
                        edge->target = k;
                        edge = getedge(arg, edges[j]);
                        edge->target = k;
                        edges[i] = k << 1;
                        edges[j] = -1;
                        /* Remove j from llist of sequences */
                        ChangeCurrent(lseq, Bottom(sequences));
                        /* Move counter in sequences one to the left to make sure that
                         * it will not be referring to the element we are about to
                         * remove from sequences.
                         */
                        Prev(lseq);
                        Dequeue(sequences);
                        coalesce(h, e->event.c.s1, e->event.c.s2);
//                         output_genes(h, output, "New:\n");
                        break;
                    case REMOVE:
                            /* Next event is a coalescence of two sequences where one is
                             * subsumed in the other. Compared to the COALESCENCE event,
                             * this event differs in the ordering of the remaining
                             * sequences after the coalescence; after a COALESCENCE event
                             * the last sequence is moved to the place of the sequence
                             * disappearing by the coalescence, while after a REMOVE event
                             * the entire block of sequences following the disappearing
                             * sequence is moved one place up.
                             */
                            if (g == NULL){
                                /* We do not have a reduced state from previous round,
                                 * i.e. this is the first REMOVE in a REMOVE sequence so it
                                 * was immediately preceeded by a removal of uninformative
                                 * sites and collapse of Siamese twins. Recreate this
                                 * picture. If the previous event was a REMOVE, it may have
                                 * created new Siamese twins whose collapse may allow new
                                 * safe coalesces that can confuse the backtrack, so in this
                                 * case no further reductions should be made.
                                 */
                                g = copy_genes(h);
                                tmp = eventlist;
                                eventlist = NULL;
                                remove_nonsegregating(g);
                                if (g->length > 0)
                                    remove_siamesetwins(g);
                            }
                            /* Determine sequence the removed sequence is subsumed in */
                            if (g->length > 0)
                                i = find_safe_coalescence(g, e->event.remove);
                            else
                                /* We should never get in this situation - if all
                                 * informative columns are removed, the remaining events
                                 * should be a sequences of coalescences with the first
                                 * sequence and the next remaining sequence.
                                 */
                                i = (e->event.remove == 0 ? 1 : 0);
                            eventlist = tmp;
                            j = (int)SetCounter(lseq, i);
                            k = (int)SetCounter(lseq, e->event.remove);
                            if (output != NULL){
                                s = NULL;
                                t = NULL;
                                if (a->sequences != NULL){
                                    if (Length(a->sequences) > j)
                                        s = (char *)GetByIndex(a->sequences, j);
                                    if (Length(a->sequences) > k)
                                        t = (char *)GetByIndex(a->sequences, k);
                                }
                                if (s != NULL) {
                                    if(howverbose != -1) {
                                        fprintf(output, "Coalescing sequences %s", s);
                                    }
                                }
                                else {
                                    if(howverbose != -1) {
                                        fprintf(output, "Coalescing sequences %d", j + 1);
                                    }
                                }
                                if (t != NULL) {
                                    if(howverbose != -1) {
                                        fprintf(output, " and %s\n", t);
                                    }
                                }
                                else {
                                    if(howverbose != -1) {
                                        fprintf(output, " and %d\n", k + 1);
                                    }
                                }
                            }
                            /* Now remove subsumed sequence and compact sequence set */
                            coalesce(h, i, e->event.remove);
                            if (e->event.remove < h->n - 1)
                                memmove(&(h->data[e->event.remove]), &(h->data[e->event.remove + 1]),
                                        (h->n - e->event.remove) * sizeof(Gene));
                                coalesce(g, i, e->event.remove); /* Also remove it in reduced copy */
                                if (e->event.remove < g->n - 1)
                                    memmove(&(g->data[e->event.remove]), &(g->data[e->event.remove + 1]),
                                            (g->n - e->event.remove) * sizeof(Gene));
                                    RemoveMoveLeft(sequences, lseq);
                                /* Update ARG */
                                s = (char *)xmalloc((a->g->length + 1) * sizeof(char));
                            strcpy(s, sequence[j]);
                            i = arg_addnode(arg, ARGCOALESCENCE, NULL, s);
                            edge = getedge(arg, edges[j]);
                            edge->target = i;
                            edge = getedge(arg, edges[k]);
                            edge->target = i;
                            edges[j] = i << 1;
                            edges[k] = -1;
//                             output_genes(h, output, "New:\n");
                            break;
                        case RECOMBINATION:
                                n_re++;
                                /* Next event is a recombination, splitting one sequence into
                                 * two sequences; the prefix will retain the place of the old
                                 * sequence while the postfix is inserted at the end of the
                                 * list of existing sequences.
                                 */
                                i = (int)Top(SetCounter(lpos, e->event.r.pos));
                                j = (int)SetCounter(lseq, e->event.r.seq);
                                split(h, e->event.r.seq, i);
                                Enqueue(sequences, (void *)next_seq);
                                /* Update ARG */
                                s = (char *)xmalloc((a->g->length + 1) * sizeof(char));
                                strcpy(s, sequence[j]);
                                strcpy(sequence[next_seq], sequence[j]);
                                l = arg_addnode(arg, ARGRECOMBINATION, NULL, s, i);
                                edge = getedge(arg, edges[j]);
                                edge->target = l;
                                /* Determine whether a SWAP event immediately following will
                                 * change whether the prefix or postfix retains the place of
                                 * the old sequence.
                                 */
                                f = Next(lcounter);
                                if ((f != NULL) && (f->type == SWAP)){
                                    pfix = "pre";
                                    swap_genes(h, f->event.swap.s1, f->event.swap.s2);
                                    /* Update sequences for ARG reconstruction */
                                    for (k = 0; k < i; k++)
                                        sequence[j][k] = 'x';
                                    for (; k < a->g->length; k++)
                                        sequence[next_seq][k] = 'x';
                                    /* Update ARG */
                                    edges[j] = (l << 1) + 1;
                                    edges[next_seq++] = l << 1;
                                }
                                else{
                                    pfix = "suf";
                                    Prev(lcounter);
                                    /* Update sequences for ARG reconstruction */
                                    for (k = 0; k < i; k++)
                                        sequence[next_seq][k] = 'x';
                                    for (; k < a->g->length; k++)
                                        sequence[j][k] = 'x';
                                    /* Update ARG */
                                    edges[j] = l << 1;
                                    edges[next_seq++] = (l << 1) + 1;
                                }
                                /* Determine site label */
                                if (a->positions != NULL)
                                    s = (char *)GetByIndex(a->positions, i - 1);
                                else
                                    s = i2a(i);
                                /* Label recombination node with recombination site */
                                k = strlen(s);
                                t = (char *)xmalloc((k + 2) * sizeof(char));
                                strcpy(t, s);
                                t[k] = '-';
                                t[k + 1] = '\0';
                                node = arg_getnode(arg, l);
                                node->label = t;
                                if (output != NULL){
                                    /* Determine sequence label */
                                    t = NULL;
                                    if ((a->sequences != NULL) && (Length(a->sequences) > j))
                                        t = (char *)GetByIndex(a->sequences, j);
                                    if (t != NULL) {
                                        if(howverbose != -1) {
                                            fprintf(output, "---->");
                                            fprintf(output, "Recombination in sequence %s after site %s; %sfix is new sequence %d\n", t, s, pfix, next_seq);
                                        }
                                        if(howverbose > 0 && output != stdout) {
                                            printf("Recombination in sequence %s after site %s; %sfix is new sequence %d\n", t, s, pfix, next_seq);
                                        }
                                    }
                                    else {
                                        if(howverbose != -1) {
                                            fprintf(output, "---->");
                                            fprintf(output, "Recombination in sequence %d after site %s; %sfix is new sequence %d\n",j + 1, s, pfix, next_seq);
                                        }
                                        if(howverbose > 0 && output != stdout) {
                                            printf("Recombination in sequence %d after site %s; %sfix is new sequence %d\n",j + 1, s, pfix, next_seq);
                                        }
                                    }
                                }
                                /* Clean up */
                                if (a->positions == NULL)
                                    free(s);
                                break;
                                case LOOKUP:
                                    if (output != NULL){
                                        fprintf(output, "Ancestral state");
                                        output_genes_indexed(h, output);
                                        fprintf(output,
                                                "found in lookup table; it requires %d recombinations\n",
                                                e->event.lookup);
                                    }
                                    break;
                }
                #ifdef DEBUG
                /* Sanity check - did we see this ancestral state in the forward pass? */
                if ((ancestral_state_trace != NULL) && (e->type != RECOMBINATION)){
                    tmp = eventlist;
                    eventlist = NULL;
                    g = copy_genes(h);
                    implode_genes(g);
                    eventlist = tmp;
                    if (!no_recombinations_required(g)){
                        p = pack_genes(g);
                        if (!hashtable_lookup(p, ancestral_state_trace, NULL)){
                            fprintf(stderr, "Warning - did not encounter state\n\n");
                            output_genes_indexed(h, stderr);
                            fprintf(stderr, "\nreached from\n\n");
                            output_genes_indexed(old, stderr);
                            fprintf(stderr, "\nin forward pass.\n");
                        }
                        free_packedgenes(p);
                    }
                    free_genes(g);
                    g = NULL;
                }
                free_genes(old);
                #endif
            }
            if(howverbose != -1 && output != stdout && output != NULL) {
                fprintf(output, "Total: %d sequencing errors, %d recurrent mutations, %d recombinations.\n", n_se, n_rm, n_re);
            }

            /* Free copy of genes, if present */
            if (g != NULL)
                free_genes(g);
            /* Finalise ARG */
            for (i = a->g->n; i < arg->n; i++){
                node = arg_getnode(arg, i);
                if (node->type == ARGRECOMBINATION){
                    if (node->predecessor.two.prefix.target == -1)
                        /* Unterminated recombination edge */
                        for (j = 0; j < n; j++)
                            if (edges[j] >> 1 == i){
                                node->predecessor.two.prefix.target
                                = arg_addnode(arg, ARGANCESTOR, NULL, sequence[j]);
                                sequence[j] = NULL;
                                break;
                            }
                            if (node->predecessor.two.suffix.target == -1)
                                /* Unterminated recombination edge */
                                for (j = 0; j < n; j++)
                                    if (edges[j] >> 1 == i){
                                        node->predecessor.two.suffix.target
                                        = arg_addnode(arg, ARGANCESTOR, NULL, sequence[j]);
                                        sequence[j] = NULL;
                                        break;
                                    }
                }
            }
            for (i = 0; i < n; i++)
                if (edges[i] != -1){
                    edge = getedge(arg, edges[i]);
                    if (Length(edge->mutations) > 0)
                        /* Mutations have appeared on this edge since its
                         * introduction - terminate it with a ANCESTOR node.
                         */
                        edge->target = arg_addnode(arg, ARGANCESTOR, NULL, sequence[i]);
                    else
                        if (sequence[i] != NULL)
                            free(sequence[i]);
                }
                else
                    if (sequence[i] != NULL)
                        free(sequence[i]);
                    arg_finalise(arg);
                free(sequence);
                free(edges);
                
                if (Length(positions) > 0){
                    /* Sanity check */
                    fprintf(stderr, "Possible error in backtrack\n");
                    while (Length(positions) > 0)
                        DestroyLList(Pop(positions));
                }
                DestroyLList(positions);
                DestroyLList(sequences);
                DestroyCounter(lcounter);
                DestroyCounter(lseq);
                DestroyCounter(lpos);
                DestroyCounter(lc);
    }
    else if (output != NULL)
        fprintf(output, "No valid history found\n");
    
    free_genes(h);
    
    return arg;
}
