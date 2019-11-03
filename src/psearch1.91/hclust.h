/****************************************************************
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *   Please direct all questions related to this program or bug reports or
 *updates to Yevgeniy Gelfand (ygelfand@bu.edu)
 *
 ****************************************************************/

#ifndef HCLUST_H
#define HCLUST_H

#include <math.h>

#include "../libs/easylife/easylife.h"
#include "bitwise LCS multiple word.h"
#include "bitwise edit distance alignment multiple word no end penalty.h"
#include "doublehash.h"
#include "profile.h"

//#pragma pack(push)  /* push current alignment to stack */
//#pragma pack(1)     /* set alignment to 1 byte boundary */

typedef struct {
    short int pmin;
    char      dir;
    PROFILE * prof;
} SEED_HIT;

//#pragma pack(pop)

typedef struct tagSTATS {
    __int64 stTuplesProcessed; // various statistics
    __int64 stCandidatesConsidered;
    __int64 stLCSalignments;
    __int64 stProfileAlignments;
    __int64 stProfileAlignmentsSuccess;
    __int64 stPairsConfirmedAfterFlankAlignments;
    __int64 CellsProcessed;
} BLASTSTATS;

BLASTSTATS blaststats;

int Trange_Profile_Range[MAXPROFILESIZE + 1];

#define TRUNC_ROUND_CEIL \
    ( 1.0 ) // truncate = 0.0, round = 0.5, ciel = 1.0 (default)

#define MAX_ARRAY_TRANGE ( 4 )
#define MAX_LIST_TRANGE ( 14 )
#define MAX_N_ONES ( 29 )

EASY_ARRAY *
  *Easy_Array_Tuples[MAX_LIST_TRANGE + 1]; // range 1 - MAX_ARRAY_TRANGE

SEED_HIT **Array_Tuples[MAX_LIST_TRANGE + 1]; // range 1 - MAX_ARRAY_TRANGE
EASY_ARRAY *
  *Array_Tuples_Index[MAX_LIST_TRANGE + 1]; // range 1 - MAX_ARRAY_TRANGE
unsigned int
  *Array_Tuples_Length[MAX_LIST_TRANGE + 1];      // range 1 - MAX_ARRAY_TRANGE
unsigned int Trange_N_Codes[MAX_LIST_TRANGE + 1]; // range 1 - MAX_ARRAY_TRANGE

GDHASH *Hash_Tuples[MAX_LIST_TRANGE +
                    1]; // range (MAX_ARRAY_TRANGE+1) - MAX_LIST_TRANGE

int   Trange_Tuple_Size[MAX_LIST_TRANGE + 1]; // range 1 - MAX_LIST_TRANGE
int   Trange_Ones[MAX_LIST_TRANGE + 1];       // range 1 - MAX_LIST_TRANGE
char *Trange_Seed[MAX_LIST_TRANGE + 1];       // range 1 - MAX_LIST_TRANGE

long int *Index, *Complement;
char *    Complementascii = NULL;

typedef struct {
    char *seed;
    int   length; // seed length
    int   num_ones;
    int   hasX;
    // int max_code_num;
    int  num_patterns; // 2 with X, 1 no X
    int *foffset[2];   // forward offset arrays - 0 unshifted, 1 shifted
    int *roffset[2];   // reverse offset arrays
} SEED_STRUCT;

SEED_STRUCT seedstruct;

/************************************************************************************************************************/
int _HCLUST_process_sequence( char *Sequence, SEED_STRUCT *pSS, int TRANGE,
  int MaxTupleSize_SN, PROFILE *prof, int rc, int pmin ) {
    int          Length          = prof->seq->length;
    int          tuplesprocessed = 0, where;
    int          g = 0, badcharindex = -1;
    int          offset;
    unsigned int seedcode, seedcode2;
    int          i;
    int *        SeedSeenHash = NULL;
    int          code1Bad = 0, code2Bad = 0;
    EASY_ARRAY * artemp;
    SEED_HIT *   sh;

    // printf(stderr, "\ntrange: %d sequence: %s\n",TRANGE,Sequence);

    /* start processing */
    for ( i = 0; i < ( Length ); i++ ) {

        if ( ( strchr( "ACGT", Sequence[i] ) ==
               NULL ) ) { /* not one of A,C,G,T */
            badcharindex = i;
            /* find first good string of mintupsize characters */
            g = 0;

            while ( ( g < MaxTupleSize_SN ) && ( i < ( Length - 1 ) ) ) {
                i++;

                /* adjust the running max and min centers */
                // for(k=1;k<=NT_P;k++)
                // Update_Runningmaxmincenters3(i,k);
                if ( strchr( "ACGT", Sequence[i] ) == NULL ) {
                    badcharindex = i;
                    g            = 0;
                } else
                    g++;
            }

            if ( g < MaxTupleSize_SN )
                break; /* i=Length and minimum tuple not found */
        }

        // fprintf(stderr,"\npSS->roffset[0]: %d",pSS->roffset[0]);
        // fprintf(stderr,"\npSS->roffset[1]: %d",pSS->roffset[1]);

        // fprintf(stderr,"\n\npSS->foffset[0]: %d",pSS->foffset[0]);
        // fprintf(stderr,"\npSS->foffset[1]: %d",pSS->foffset[1]);

        // fprintf(stderr,"\n\npSS->num_patterns: %d",pSS->num_patterns);

        // fprintf(stderr,"\n\npIndex: %d",Index);

        // printf("\n\npIndex['A']: %d",Index['A']);

        // exit(1);

        if ( ( i - badcharindex ) >=
             MaxTupleSize_SN ) { /* check for valid tuple */

            /* produce code */
            {
                {
                    seedcode = seedcode2 = 0;

                    for ( offset = 0; offset < pSS->num_ones; offset++ ) {

                        // if (0 == rc)
                        //  where = i + pSS->roffset[0][pSS->num_ones - offset
                        //  -1];
                        // else
                        where = i + pSS->foffset[0][offset];

                        // fprintf(stderr,"\nOFFSET: %d, value: %d, i: %d,
                        // where:
                        // %d\n\n",offset,pSS->foffset[0][offset],i,where);

                        if ( where < 0 ) {
                            doCriticalErrorAndQuit( "WHERE is <0!!" );
                        }

                        if ( where >= ( Length + MaxTupleSize_SN ) ) {
                            doCriticalErrorAndQuit( "where too large!" );
                        }

                        if ( !( Sequence[where] == 'A' ||
                                Sequence[where] == 'C' ||
                                Sequence[where] == 'G' ||
                                Sequence[where] == 'T' ) ) {
                            doCriticalErrorAndQuit(
                              "Unknown character in sequence! (%c) seqlen: %d, "
                              "trange: %d, MaxTupleSize_SN: %d, where: %d, rc: "
                              "%d, sequence: %s\n",
                              Sequence[where], Length, TRANGE, MaxTupleSize_SN,
                              where, rc, Sequence );
                        }

                        if ( offset <= 14 )
                            seedcode = seedcode * 4 + Index[Sequence[where]];
                        else if ( offset <= 29 )
                            seedcode2 = seedcode2 * 4 + Index[Sequence[where]];
                        else
                            doCriticalErrorAndQuit(
                              "Seed too large! Aborting!" );
                    }
                }
                // Tuplecode[numpatterns]=code;
                // Seedcode[numpatterns]=seedcode;
                // Seedcode2[numpatterns]=seedcode2;
            }

            /* create an entry for seed hit, use -location to indicated reverse
             * complemented hit, distinction between shifted and unshifted is
             * not needed */
            sh       = (SEED_HIT *) scalloc( 1, sizeof( SEED_HIT ) );
            sh->pmin = pmin;
            sh->dir  = rc;
            sh->prof = prof;

            if ( pSS->num_ones <= 13 && TRANGE <= MAX_ARRAY_TRANGE &&
                 TRANGE >= 1 ) {

                if ( NULL ==
                     ( artemp = Easy_Array_Tuples[TRANGE][seedcode] ) ) {
                    Easy_Array_Tuples[TRANGE][seedcode] = artemp =
                      EasyArrayCreate( 2, NULL, free );
                }

                EasyArrayInsert( artemp, (void *) sh );

            } else if ( pSS->num_ones <= MAX_N_ONES &&
                        TRANGE > MAX_ARRAY_TRANGE &&
                        TRANGE <= MAX_LIST_TRANGE ) {

                if ( NULL ==
                     ( artemp = (EASY_ARRAY *) GetDoubleHashItem(
                         Hash_Tuples[TRANGE], seedcode, seedcode2 ) ) ) {
                    artemp = EasyArrayCreate( 2, NULL, free );
                    SetDoubleHashItem( Hash_Tuples[TRANGE], seedcode, seedcode2,
                      (void *) artemp, NULL, 0 );
                }

                EasyArrayInsert( artemp, (void *) sh );

            } else {
                doCriticalErrorAndQuit( "Seed with more than %d 1s or TRANGE "
                                        "higher than 7, not yet implemented!",
                  MAX_N_ONES );
            }

            /*
                sh->location=i;   // i is index of current tuple
                sh->seq=pseq;
                sh->tsize=MaxTupleSize_SN;


                if
               (NULL==(rc=GetTripleHashItem(SHash,seedcode,seedcode2,seedcode3)))
               { wordList=EasyListCreate(NULL,free);
                        SetTripleHashItem(SHash,seedcode,seedcode2,seedcode3,wordList);
                } else {
                        wordList = rc;
                }


                EasyListInsertHead(wordList,sh);

            */

            tuplesprocessed++;
        }
    }

    return tuplesprocessed;
}

/************************************************************************************************************************/
int _HCLUST_find_candidates( char *Sequence, FASTASEQUENCE *pseq,
  SEED_STRUCT *pSS, int TRANGE, int MaxTupleSize_SN, int rc, int pmin,
  EASY_ARRAY *candidates ) {

    int Length = pseq->length;
    // int tuplesprocessed=0;
    int          where;
    int          g = 0, badcharindex = -1;
    int          offset;
    unsigned int seedcode, seedcode2;
    size_t       length;
    int          i;
    int *        SeedSeenHash = NULL;
    int          code1Bad = 0, code2Bad = 0;
    EASY_ARRAY * iatemp;
    SEED_HIT *   shtemp;

    /* start processing */
    for ( i = 0; i < ( Length ); i++ ) {

        if ( ( strchr( "ACGT", Sequence[i] ) ==
               NULL ) ) { /* not one of A,C,G,T */
            badcharindex = i;
            /* find first good string of mintupsize characters */
            g = 0;

            while ( ( g < MaxTupleSize_SN ) && ( i < ( Length - 1 ) ) ) {
                i++;

                if ( strchr( "ACGT", Sequence[i] ) == NULL ) {
                    badcharindex = i;
                    g            = 0;
                } else
                    g++;
            }

            if ( g < MaxTupleSize_SN )
                break; /* i=Length and minimum tuple not found */
        }

        if ( ( i - badcharindex ) >=
             MaxTupleSize_SN ) { /* check for valid tuple */

            SEED_HIT *start_hit = NULL, *end_hit = NULL, *temp_hit;
            GDHITEM * gdhit;
            size_t    uj;

            /* produce code */
            {
                {
                    seedcode = seedcode2 = 0;

                    for ( offset = 0; offset < pSS->num_ones; offset++ ) {

                        // if (0 == rc)
                        //  where = i + pSS->roffset[0][pSS->num_ones - offset
                        //  -1];
                        // else
                        where = i + pSS->foffset[0][offset];

                        if ( where < 0 ) {
                            doCriticalErrorAndQuit( "WHERE is <0!!" );
                        }

                        if ( where >= ( Length + MaxTupleSize_SN ) ) {
                            doCriticalErrorAndQuit( "where too large!" );
                        }

                        if ( !( Sequence[where] == 'A' ||
                                Sequence[where] == 'C' ||
                                Sequence[where] == 'G' ||
                                Sequence[where] == 'T' ) ) {
                            doCriticalErrorAndQuit(
                              "Unknown character in sequence! (%c) seqlen: %d, "
                              "trange: %d, MaxTupleSize_SN: %d, where: %d, rc: "
                              "%d, sequence: %s\n",
                              Sequence[where], Length, TRANGE, MaxTupleSize_SN,
                              where, rc, Sequence );
                        }

                        if ( offset <= 14 )
                            seedcode = seedcode * 4 + Index[Sequence[where]];
                        else if ( offset <= 29 )
                            seedcode2 = seedcode2 * 4 + Index[Sequence[where]];
                        else
                            doCriticalErrorAndQuit(
                              "Seed too large! Aborting!" );
                    }
                }
            }

            /* create an entry for seed hit, use -location to indicated reverse
             * complemented hit, distinction between shifted and unshifted is
             * not needed */
            if ( pSS->num_ones <= 13 && TRANGE <= MAX_ARRAY_TRANGE &&
                 TRANGE >= 1 ) {

                if ( NULL !=
                     ( iatemp = Array_Tuples_Index[TRANGE][seedcode] ) ) {

                    for ( uj = 0; uj < iatemp->size; uj++ ) {

                        temp_hit = (SEED_HIT *) ( iatemp->array[uj] );

                        // fprintf(stderr,"\n\tGOT HERE!, TRANGE: %d,
                        // temp_hit->pmin: %u",TRANGE,temp_hit->pmin);

                        if ( !start_hit &&
                             temp_hit->pmin >=
                               ( pmin -
                                 (int) ( pmin * PATLEN_SIZE_ERR_FRACTION +
                                         TRUNC_ROUND_CEIL ) ) ) {
                            start_hit = temp_hit;
                        }

                        if ( !end_hit &&
                             temp_hit->pmin >
                               ( pmin +
                                 (int) ( pmin * PATLEN_SIZE_ERR_FRACTION +
                                         TRUNC_ROUND_CEIL ) ) ) {
                            end_hit = temp_hit;
                        }

                        if ( start_hit && end_hit ) {
                            break;
                        }
                    }

                    //                    if (!start_hit) { start_hit = &
                    //                    Array_Tuples[TRANGE][seedcode][ 0 ]; }
                    if ( !end_hit ) {
                        end_hit =
                          &Array_Tuples[TRANGE][seedcode]
                                       [Array_Tuples_Length[TRANGE][seedcode] -
                                         1];
                    }

                    // fprintf(stderr,"\nGOT HERE!!, start_hit: %llu, end_hit:
                    // %llu",start_hit,end_hit); exit(1);

                    if ( start_hit && end_hit && start_hit <= end_hit ) {

                        for ( ;; start_hit++ ) {
                            EasyArrayInsert( candidates, (void *) start_hit );

                            if ( start_hit == end_hit )
                                break;
                        }
                    }
                }

            } else if ( pSS->num_ones <= MAX_N_ONES &&
                        TRANGE > MAX_ARRAY_TRANGE &&
                        TRANGE <= MAX_LIST_TRANGE ) {

                if ( NULL !=
                     ( gdhit = (GDHITEM *) GetDoubleHash(
                         Hash_Tuples[TRANGE], seedcode, seedcode2 ) ) ) {

                    iatemp = (EASY_ARRAY *) gdhit->index;
                    shtemp = (SEED_HIT *) gdhit->data;
                    length = gdhit->length;

                    for ( uj = 0; uj < iatemp->size; uj++ ) {

                        temp_hit = (SEED_HIT *) iatemp->array[uj];

                        // temp_hit = (SEED_HIT*)EasyArrayItem(iatemp,uj);
                        if ( !start_hit &&
                             temp_hit->pmin >=
                               ( pmin -
                                 (int) ( pmin * PATLEN_SIZE_ERR_FRACTION +
                                         TRUNC_ROUND_CEIL ) ) ) {
                            start_hit = temp_hit;
                        }

                        if ( !end_hit &&
                             temp_hit->pmin >
                               ( pmin +
                                 (int) ( pmin * PATLEN_SIZE_ERR_FRACTION +
                                         TRUNC_ROUND_CEIL ) ) ) {
                            end_hit = temp_hit;
                        }

                        if ( start_hit && end_hit ) {
                            break;
                        }
                    }

                    if ( !end_hit ) {
                        end_hit = &shtemp[length - 1];
                    }

                    //         start_hit = & shtemp[ 0 ];
                    //         end_hit = & shtemp[ length - 1 ];

                    if ( start_hit && end_hit && start_hit <= end_hit ) {
                        for ( ;; start_hit++ ) {
                            // printf("\nstart_hit: %llu end_hit: %llu length:
                            // %d",start_hit,end_hit,length);
                            EasyArrayInsert( candidates, (void *) start_hit );

                            if ( start_hit == end_hit )
                                break;
                        }
                    }
                }

                // EasyArrayInsert(  artemp , (void*) sh);

            } else {
                doCriticalErrorAndQuit( "Seed with more than %d 1s or TRANGE "
                                        "higher than 7, not yet implemented!",
                  MAX_N_ONES );
            }

            // tuplesprocessed++;
        }
    }

    return 0;
}

/************************************************************************************************************************/
int _HCLUST_align_candidates( PROFILE *readprof, int rc, int pmin,
  EASY_ARRAY *candidates, unsigned char *dt1, CLUSTERBASE *cb ) {

    PROFPAIR *   refpair;
    SEED_HIT *   pshit;
    PAP *        pap1;
    EASY_NODE *  nof1;
    PROFILE *    refprof;
    double       bestsim, sizeerror;
    int          oldkey;
    int          concensusLCS, shortest, minlen;
    char         olddir;
    unsigned int uj;
    unsigned int written = 0;

    // go through all candidates
    refprof = NULL;

    for ( uj = 0; uj < candidates->size; uj++ ) {

        pshit = (SEED_HIT *) candidates->array[uj];

        blaststats.stCandidatesConsidered++;

        // skip same ref same dir, sorted
        if ( pshit->prof == refprof ) {
            continue;
        }

        refprof = pshit->prof;

        // use lcs on concensus sequences firest
#ifdef CONCENSUS_LCS_CHECK

        blaststats.stLCSalignments++;

        concensusLCS =
          LCS_multiple_word( refprof->seq->sequence, readprof->seq->sequence,
            refprof->seq->length, readprof->seq->length );

        shortest = min( readprof->seq->length, refprof->seq->length );
        minlen =
          (int) ( shortest * LCS_CUTOFF + .5 ); // + .5 in there would round and
                                                // would filter more stuff out
        minlen = min( minlen, ( shortest - 1 ) );

        if ( concensusLCS < minlen ) {
            continue;
        }

#endif

        // align profiles
        if ( min( refprof->patlen, readprof->patlen ) >= 30 ) {

            if ( refprof->proflen > readprof->proflen ) {
                sizeerror =
                  ( refprof->proflen / (double) readprof->proflen - 1.0 ) * 100;
            } else {
                sizeerror =
                  ( readprof->proflen / (double) refprof->proflen - 1.0 ) * 100;
            }

            // narrowband sometimes returns error when patterns are too
            // different (100% diff) in size so use regular
            if ( sizeerror > 100.0 ) {
                pap1 = GetFixedProfileAP( refprof, readprof, dt1 );
            } else {
                pap1 = GetFixedProfileAPNarrowband(
                  refprof, readprof, dt1, TARGET_HOMOLOGY_INT );
            }

            bestsim = pap1->similarity;
            FreePAP( pap1 );

        } else {

            pap1    = GetFixedProfileAP( refprof, readprof, dt1 );
            bestsim = pap1->similarity;
            FreePAP( pap1 );
        }

        blaststats.stProfileAlignments++;

        /* DEBUG
        if (refprof->key == 175388148 || refprof->key == -175388148 ||
        refprof->key == 176296442 || refprof->key == -176296442) { if
        (readprof->key == 175388148 || readprof->key == -175388148 ||
        readprof->key == 176296442 || readprof->key == -176296442) {
        printf("\nDEBUG %d%s vs %d%s,\tprofscore = %.2lf\n",refprof->key,
        refprof->rcflag ? "RC" : "", readprof->key,readprof->rcflag ? "RC" :
        "",bestsim);
        }}
        */

        // add to final result
        if ( bestsim >= TARGET_HOMOLOGY ) {

            EASY_NODE *nrotf1, *nrotf2;
            PROFILE *  prof1ROT, *prof1rcROT, *prof2ROT, *prof2rcROT;
            PROFPAIR * profpair1, *profpair2;
            int        flanksPasses, lerr, rerr;
            PROFILE *  read_ptr, *ref_ptr;

            blaststats.stProfileAlignmentsSuccess++;

            /* cyclic rotations, check all flanks, 1.87 change */
            for ( nrotf1 = refprof->rotlist->head; nrotf1 != NULL;
                  nrotf1 = nrotf1->next ) {
                for ( nrotf2 = readprof->rotlist->head; nrotf2 != NULL;
                      nrotf2 = nrotf2->next ) {
                    // fprintf(fpo,"%d,%d,%lf\n",prof1->key, prof2->key,
                    // bestsim);

                    int newdir;

                    profpair1 = (PROFPAIR *) EasyListItem( nrotf1 );
                    profpair2 = (PROFPAIR *) EasyListItem( nrotf2 );

                    prof1ROT   = profpair1->prof;
                    prof1rcROT = profpair1->profrc;
                    prof2ROT   = profpair2->prof;
                    prof2rcROT = profpair2->profrc;

                    // printf("\nbestsim: %.2l, fref ->: %d <-: %d   read -> :
                    // %d <-:
                    // %d",bestsim,prof1ROT->key,prof1rcROT->key,prof2ROT->key,prof2rcROT->key);
                    // printf("\nREF: %s|%s RC:
                    // %s|%s",prof1ROT->left,prof1ROT->right,prof1rcROT->left,prof1rcROT->right);
                    // printf("\nREAD: %s|%s RC:
                    // %s|%s",prof2ROT->left,prof2ROT->right,prof2rcROT->left,prof2rcROT->right);
                    // printf("\n");

                    // set newdir to read orientation
                    newdir = rc;

                    // flip dir if ref is reversed
                    if ( 1 == pshit->dir )
                        newdir = !newdir;

                    // flip dir if redundant ref is reversed
                    if ( profpair1->prof->dir )
                        newdir = !newdir;

                    // flip dir if redundant read is reversed
                    if ( profpair2->prof->dir )
                        newdir = !newdir;

                    read_ptr = prof2ROT;

                    if ( 0 == newdir ) {
                        ref_ptr = prof1ROT;
                    } else {
                        ref_ptr = prof1rcROT;
                    }

                    // printf("\nAligning flanks ");

                    /* refs vs reads */
                    if ( OPTION != 'R' ) {

                        int MAXEL1, MAXER1, LARGESERRORALLOWED;

                        if ( 0 != MAXERRORS ) {
                            MAXEL1 = MAXERRORS;
                            MAXER1 = MAXERRORS;
                        } else {
                            MAXEL1 =
                              min( 8, (int) ( 0.4 * read_ptr->leftlen + .01 ) );
                            MAXER1 = min(
                              8, (int) ( 0.4 * read_ptr->rightlen + .01 ) );
                        }

                        /* reflen is almost always 60, readlen is <=50. If
                         * reflen is from ends of chromosome it could be
                         * shorter.Thats why readlen will be shortened to reflen
                         */
                        lerr = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->left, read_ptr->left, ref_ptr->leftlen,
                          ( ref_ptr->leftlen > read_ptr->leftlen )
                            ? read_ptr->leftlen
                            : ref_ptr->leftlen );
                        rerr = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->right, read_ptr->right, ref_ptr->rightlen,
                          ( ref_ptr->rightlen > read_ptr->rightlen )
                            ? read_ptr->rightlen
                            : ref_ptr->rightlen );

                        // fprintf(stderr,"\nref: %d read: %d (lflank: %d,
                        // rflank: %d) ===> lerr: %d (MAXEL1=%d), rerr: %d
                        // (MAXER1:
                        // %d)\n\n",ref_ptr->key,read_ptr->key,read_ptr->leftlen,read_ptr->rightlen,lerr,MAXEL1,rerr,MAXER1);

                        if ( lerr == -1 || rerr == -1 ) {
                            doCriticalErrorAndQuit(
                              "Edit_Distance_multiple_word_NoEndPenaltySeq1 "
                              "returned -1. Aborting!" );
                        }

                        // flanksPasses = (max(lerr,rerr) <= MAXERRORS);
                        flanksPasses = ( lerr <= MAXEL1 && rerr <= MAXER1 );

                        /* refs vs refs */
                    } else {

                        int MAXEL1, MAXER1;

                        /* Dec 5, 2014: Yozen asked to  make ref-ref scoring
                         * just like regular pipeline runs */
                        MAXEL1 =
                          min( 8, (int) ( 0.4 * read_ptr->leftlen + .01 ) );
                        MAXER1 =
                          min( 8, (int) ( 0.4 * read_ptr->rightlen + .01 ) );

                        lerr = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->left, read_ptr->left, ref_ptr->leftlen,
                          read_ptr->leftlen );
                        rerr = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                          ref_ptr->right, read_ptr->right, ref_ptr->rightlen,
                          read_ptr->rightlen );

                        if ( lerr == -1 || rerr == -1 ) {
                            doCriticalErrorAndQuit(
                              "Edit_Distance_multiple_word_NoEndPenaltySeq1 "
                              "returned -1. Aborting!" );
                        }

                        // flanksPasses = (min(lerr,rerr) <= MAXERRORS);
                        // TESTING BOTH FLANKS FOR PARAMTER SELECTION
                        flanksPasses = ( lerr <= MAXEL1 && rerr <= MAXER1 );

                        /* there is a possibility that reversing these might
                         * have different result, lets try it, 1.90 testing for
                         * mapping */
                        if ( !flanksPasses ) {

                            lerr = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                              read_ptr->left, ref_ptr->left, read_ptr->leftlen,
                              ref_ptr->leftlen );
                            rerr = Edit_Distance_multiple_word_NoEndPenaltySeq1(
                              read_ptr->right, ref_ptr->right,
                              read_ptr->rightlen, ref_ptr->rightlen );

                            if ( lerr == -1 || rerr == -1 ) {
                                doCriticalErrorAndQuit(
                                  "Edit_Distance_multiple_word_"
                                  "NoEndPenaltySeq1 returned -1. Aborting!" );
                            }

                            // flanksPasses = (min(lerr,rerr) <= MAXERRORS);
                            // TESTING BOTH FLANKS FOR PARAMTER SELECTION
                            flanksPasses = ( lerr <= MAXEL1 && rerr <= MAXER1 );
                        }

                    } // end of R option

                    /* if flanks failed, try the other way (for palindromes).
                    This is more liberal for pipeline, but not sure if this
                    should be in final version

                    TURNS OUT THIS IS NOT NEEDED, because palindromes would be
                    find on reverse seeds too, very unlikely for this to not be
                    found, even impossible for true palindromes

                    if (!flanksPasses) {

                            if (0 == newdir) {
                                    ref_ptr = prof1rcROT;
                            } else {
                                    ref_ptr = prof1ROT;
                            }

                            newdir=!newdir;

                            // reflen is almost always 60, readlen is <=50. If
                    reflen is from ends of chromosome it could be shorter.Thats
                    why readlen will be shortened to reflen

                            lerr =
                    Edit_Distance_multiple_word_NoEndPenaltySeq1(ref_ptr->left,read_ptr->left,
                    ref_ptr->leftlen,(ref_ptr->leftlen>read_ptr->leftlen)?read_ptr->leftlen:ref_ptr->leftlen);
                            rerr =
                    Edit_Distance_multiple_word_NoEndPenaltySeq1(ref_ptr->right,read_ptr->right,
                    ref_ptr->rightlen,(ref_ptr->rightlen>read_ptr->rightlen)?read_ptr->rightlen:ref_ptr->rightlen);

                            //fprintf(stderr,"\nlerr: %d, rerr:
                    %d\n\n",lerr,rerr);

                            if (lerr==-1 || rerr==-1)
                            {
                                 doCriticalErrorAndQuit("Edit_Distance_multiple_word_NoEndPenaltySeq1
                    returned -1. Aborting!");
                            }


                            flanksPasses = (max(lerr,rerr) <= MAXERRORS);
                    }
                    */

                    /* DEBUG
                    if (refprof->key == 175388148 || refprof->key == -175388148
                    || refprof->key == 176296442 || refprof->key == -176296442)
                    { if (readprof->key == 175388148 || readprof->key ==
                    -175388148 || readprof->key == 176296442 || readprof->key ==
                    -176296442) { printf("\nDEBUG %d%s vs %d%s,\tprofscore =
                    %.2lf lerr=%d, rerr=%d\n",refprof->key, refprof->rcflag ?
                    "RC" : "", readprof->key,readprof->rcflag ? "RC" :
                    "",bestsim,lerr,rerr);
                    }}
                    */

                    /* all green */
                    if ( flanksPasses ) {

                        // printf("\nAdding link!");

                        ClusterBaseAddLink(
                          cb, prof1ROT->key, prof2ROT->key, newdir );
                        blaststats.stPairsConfirmedAfterFlankAlignments++;

                        //  write to map file (new after  pscearch1.9)
                        if ( MAPFILE ) {
                            fprintf( MAPFILE, "%d%c=>%d:%.2lf:%d:%d\n",
                              prof1ROT->key, ( 0 == newdir ) ? '\'' : '\"',
                              prof2ROT->key, bestsim, lerr, rerr );
                        }

                        /* mark connections for ref-to-ref */
                        if ( OPTION == 'R' ) {

                            CLUSTERKEY *pk1, *pk2;

                            pk1 = cl_findkey( cb, prof1ROT->key );
                            pk2 = cl_findkey( cb, prof2ROT->key );

                            /*  exclude itself */
                            if ( abs( prof1ROT->key ) ==
                                 abs( prof2ROT->key ) ) {

                            } else {
                                if ( lerr <= MAXERRORS ) {
                                    pk1->leftcon = 1;
                                    pk2->leftcon = 1;
                                }

                                if ( rerr <= MAXERRORS ) {
                                    pk1->rightcon = 1;
                                    pk2->rightcon = 1;
                                }
                            }
                        }
                    }
                }
            } // end of cyclic rotation

        } // end if meeting homology

    } // end of going though candidate list

    return 0;
}

/************************************************************************************************************************/
char *seqdupPlusTuple( int length, char *in, int tsize ) {

    char *res;

    int start, stop, i;

    if ( tsize > length )
        doCriticalErrorAndQuit(
          "Tuple size (%d) is larger than sequence(%d). Aborting!", tsize,
          length );

    res = smalloc( ( ( length + tsize + 2 ) * sizeof( char ) ) );
    strcpy( res, in );
    memcpy( res + length, in, tsize );
    res[length + tsize] = 0;

    return res;
}

/************************************************************************************************************************/
void init_complement( void )

{
    int i;

    /* complement has 256 entries so that finding the entries for A, C, G and T
     */
    /* require no calculation */
    Complement = (long int *) scalloc( 256, sizeof( long int ) );

    // memory_stats_print("\n init_complement: requesting", 256 * sizeof(long
    // int) );
    if ( Complement == NULL ) {
        doCriticalErrorAndQuit( "\nInit_complement: Out of memory!" );
        exit( 0 );
    }

    for ( i = 0; i < 256; i++ )
        Complement[i] = -1;

    /*  */
    {
        Complement['A'] = 3;
        Complement['C'] = 2;
        Complement['G'] = 1;
        Complement['T'] = 0;
    }
}

/************************************************************************************************************************/
void init_index( void )

{
    int i;

    /* index has 256 entries so that finding the entries for A, C, G and T */
    /* require no calculation */
    Index = (long int *) scalloc( 256, sizeof( long int ) );

    // memory_stats_print("\n init_index: requesting", 256 * sizeof(long int) );
    if ( Index == NULL ) {
        doCriticalErrorAndQuit( "\nInit_index: Out of memory!" );
        exit( 0 );
    }

    for ( i = 0; i < 256; i++ )
        Index[i] = -1;

    Index['A'] = 0;
    Index['C'] = 1;
    Index['G'] = 2;
    Index['T'] = 3;
}

/************************************************************************************************************************/
void init_complement_ascii( void )

{
    /* complement has 256 entries so that finding the entries for A, C, G and T
     * which are complement ascii values */
    /* require no calculation */
    int i;
    Complementascii = (char *) scalloc( 256, sizeof( long int ) );

    // memory_stats_print("\n init_complement_ascii: requesting", 256 *
    // sizeof(long int) );
    if ( Complementascii == NULL ) {
        doCriticalErrorAndQuit( "\ninit_complement_ascii: Out of memory!" );
        exit( 0 );
    }

    for ( i = 0; i < 256; i++ )
        Complementascii[i] = i;

    /* */
    {
        Complementascii['A'] = 'T';
        Complementascii['C'] = 'G';
        Complementascii['G'] = 'C';
        Complementascii['T'] = 'A';
        Complementascii['a'] = 't';
        Complementascii['c'] = 'g';
        Complementascii['g'] = 'c';
        Complementascii['t'] = 'a';
    }
}

/************************************************************************************************************************/
char *GetReverseComplement( char *original ) {
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    if ( NULL == original )
        return NULL;

    /* safety check */
    if ( Complementascii == NULL ) {
        init_complement_ascii();
    }

    /* find out how long the string is */
    length = strlen( original );

    /* allocate the new string */
    buffer = (char *) malloc( sizeof( char ) * ( length + 1 ) );

    if ( buffer == NULL ) {
        return NULL;
    }

    /* reverse complement pattern */
    sourceptr  = original;
    destinptr  = buffer + length; /* position at termination */
    *destinptr = '\0';            /* terminate string */
    destinptr--;

    while ( *sourceptr != '\0' ) {
        ( *destinptr ) = Complementascii[( *sourceptr )];
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return buffer;
}

/***********************************************************************/
char *GetComplement( char *original ) {
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    /* safety check */
    if ( Complementascii == NULL ) {
        init_complement_ascii();
    }

    /* find out how long the string is */
    length = strlen( original );

    /* allocate the new string */
    buffer = (char *) malloc( sizeof( char ) * ( length + 1 ) );

    if ( buffer == NULL ) {
        return NULL;
    }

    /* complement pattern */
    sourceptr = original;
    destinptr = buffer;

    while ( *sourceptr != '\0' ) {
        ( *destinptr ) = Complementascii[( *sourceptr )];
        destinptr++;
        sourceptr++;
    }

    *destinptr = '\0'; /* terminate string */

    /* return buffer containg reverse complement */
    return buffer;
}

/************************************************************************************************************************/
char *GetReverse( char *original ) {
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    /* find out how long the string is */
    length = strlen( original );

    /* allocate the new string */
    buffer = (char *) malloc( sizeof( char ) * ( length + 1 ) );

    if ( buffer == NULL ) {
        return NULL;
    }

    /* reverse pattern */
    sourceptr  = original;
    destinptr  = buffer + length; /* position at termination */
    *destinptr = '\0';            /* terminate string */
    destinptr--;

    while ( *sourceptr != '\0' ) {
        ( *destinptr ) = *sourceptr;
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return buffer;
}

/************************************************************************************************************************/
void retrieve_seed_info( char *seed, SEED_STRUCT *pSS ) {
    int i, j, k, num, length, hasX, offset, value, shift;

    if ( NULL == pSS )
        return;

    pSS->seed   = seed;
    length      = strlen( seed );
    pSS->length = length;

    // check seed for bad characters
    for ( i = 0; i < length; i++ ) {
        if ( strchr( "1xX*", seed[i] ) == NULL ) {
            doCriticalErrorAndQuit( "\nSeed has invalid characters!\n\n" );
            exit( 0 );
        }
    }

    // find number of ones in seed
    num  = 0;
    hasX = 0;

    for ( i = 0; i < length; i++ ) {
        if ( seed[i] == '1' )
            num++;

        if ( seed[i] == 'X' )
            hasX = 1;

        if ( seed[i] == 'x' )
            hasX = 1;
    }

    if ( hasX )
        pSS->hasX = 1;

    pSS->num_patterns = ( hasX ) ? 2 : 1;
    pSS->num_ones     = num;

    if ( num > MAX_N_ONES ) {
        doCriticalErrorAndQuit(
          "\n\nSorry, number of ones in the seed cannot exceed %d.",
          MAX_N_ONES );
        exit( 0 );
    }

    // pSS->max_code_num = four_to_the[num] - 1;

    // initialize offset arrays
    for ( i = 0; i <= 1; i++ ) {
        pSS->foffset[i] = (int *) scalloc( num, sizeof( int ) );
        pSS->roffset[i] = (int *) scalloc( num, sizeof( int ) );
    }

    // calculate offset values
    // forward
    shift = 0;

    for ( k = 0; k <= hasX; k++ ) {
        j = num - 1;
        i = length;

        while ( i > 0 ) {
            if ( k == 1 )
                offset = i - length + shift;
            else
                offset = i - length;

            if ( seed[i - 1] == '1' ) {
                pSS->foffset[k][j] = offset;
                j--;
            }

            if ( k == 1 && ( seed[i - 1] == 'X' || seed[i - 1] == 'x' ) ) {
                pSS->foffset[k][j] = offset;
                j--;
                i--;
                shift = 1;
            }

            i--;
        }
    }

    // reverse
    for ( k = 0; k <= hasX; k++ ) {
        j     = num - 1;
        value = -pSS->foffset[k][0];

        for ( i = 0; i < num; i++ ) {
            pSS->roffset[k][j] = -( pSS->foffset[k][i] + value );
            j--;
        }
    }
}

/************************************************************************************************************************/
void free_seed_info( SEED_STRUCT *pSS ) {

    int i;

    for ( i = 0; i <= 1; i++ ) {
        if ( pSS->foffset[i] )
            sfree( pSS->foffset[i] );

        pSS->foffset[i] = NULL;

        if ( pSS->roffset[i] )
            sfree( pSS->roffset[i] );

        pSS->roffset[i] = NULL;
    }
}

#endif
