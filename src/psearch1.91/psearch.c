/*

this is a restructured version of proclu, optimized for speed for refs vs reads comparison,
no cyclic alignment

  1.92 - passing 0 for maxerrors will now make the program pick one based on length of the flank

  1.91 - fixed seeds from not starting at 0 position, should add results
       - MAX_FLANK_CONSIDERED is now passed to psearch.exe

  1.90 - added MAPFILE written out with -m option

  1.89 - added commmand line param to do ref-to-ref (only one flank has to match)

  1.88 - first stable version

  restractred2 - seed++ and shortened sequence (no longer includes rotated part)

*/

#include <stdio.h>

#define CONCENSUS_LCS_CHECK (1)
#define PATLEN_SIZE_ERR (10.0)
#define PATLEN_SIZE_ERR_FRACTION (0.10)
#define LCS_CUTOFF (.85)


#ifdef _WIN_32_YES
#include <windows.h>
#endif

#ifndef _WIN_32_YES
#define LARGE_INTEGER time_t
#define I64d ld
#endif

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


int REFLEN = 50; // for ref-vs-ref
int  MAXFLANKCONSIDERED = 1000000;  // for ref-vs-read
char OPTION = 0;  // for ref-vs-ref
char OPTION2 = 0; // v1.9, for mapfile

static time_t startTime;

#include "../libs/easylife/easylife.h"

#include "doublehash.h"
#include "profile.h"


#include "cluster.h"

double  TARGET_HOMOLOGY = 0.0;
int     TARGET_HOMOLOGY_INT = 0;

// this is for flanks, maximum distance in sum of flank alignments
int  MAXERRORS = 0;

typedef struct {
    PROFILE *prof;
    PROFILE *profrc;
} PROFPAIR;

GSHASH *profileHash = NULL;
FILE *MAPFILE = NULL; // written in version 1.9

#include "hclust.h"

int LoadRotated(FILE *fp);

//#include "narrowbandDistanceAlignment.h"
#include "bitwise edit distance alignment multiple word no end penalty.h"
#include "bitwise LCS single word.h"
#include "bitwise LCS multiple word.h"


int doEdgesBruteForce(FILE *fpi, FILE *fpi2, FILE *edgesIn, unsigned char *dt1, int fixed, char *outputfile);

/*******************************************************************************************/
int patcmp( const void *item1, const void *item2 )
{

    short int d1, d2;

    d1 = ((SEED_HIT *)item1)->pmin;
    d2 = ((SEED_HIT *)item2)->pmin;

    if (d1 > d2)
        return 1;
    else if (d1 < d2)
        return -1;

    return 0;
}

/*******************************************************************************************/
int candsort( const void *item1, const void *item2 )
{

    PROFILE *d1, *d2;

    d1 = ((SEED_HIT *)item1)->prof;
    d2 = ((SEED_HIT *)item2)->prof;

    if (d1 > d2)
        return 1;
    else if (d1 < d2)
        return -1;

    return 0;
}

/*******************************************************************************************/
CLUSTERBASE *doSearchSimilarities(FILE *fpi, FILE *fpi2, FILE *fpirot, FILE *fpi2rot, unsigned char *dt1)
{


    EASY_LIST       *profileList = NULL, *profileList2 = NULL, *tempList = NULL;
    EASY_ARRAY      *artemp, *iatemp;
    EASY_NODE       *nof1, *nof2;
    SEED_HIT        *shtemp;
    GDHASH      *gdtemp;
    EASY_NODE       *iter;
    GDHITEM     *gdrack, *curr;
    PROFPAIR        *profpair, *profpair2;
    CLUSTERBASE     *cb;
    PROFILE     *prof1, *prof2, *prof1rc, *prof2rc;
    int         off1, i, j, pmin, TRANGE, NEWRANGE, lowerpat, higherpat, readcount, readpercent;
    size_t    ui, uj;
    char        *src, *sequence;


    // initializations
    seedstruct.foffset[0] = NULL;
    seedstruct.foffset[1] = NULL;
    seedstruct.roffset[0] = NULL;
    seedstruct.roffset[1] = NULL;

    init_complement_ascii();
    init_index();
    init_complement();

    blaststats.stTuplesProcessed = 0;
    blaststats.CellsProcessed = 0;

    blaststats.stCandidatesConsidered = 0;
    blaststats.stLCSalignments = 0;
    blaststats.stProfileAlignments = 0;
    blaststats.stProfileAlignmentsSuccess = 0;
    blaststats.stPairsConfirmedAfterFlankAlignments = 0;


    // easy array of tuples, it automatically reallocs when grows bigger
    Easy_Array_Tuples[0] = NULL;
    Trange_Tuple_Size[0] = 0;
    Trange_N_Codes[0] = 0;
    Trange_Ones[0] = 0;
    Trange_Seed[0] = NULL;

    for (TRANGE = 1; TRANGE <= MAX_LIST_TRANGE; TRANGE++) {

        /* SEED++ */
        if (TRANGE == 14) { Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11**111*1*11");  }
        else if (TRANGE == 13) { Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11**111");  }
        else if (TRANGE == 12) { Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11");  }
        else if (TRANGE == 11) { Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE == 10) { Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE == 9) { Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }

        else if (TRANGE == 8) { Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE == 7) { Trange_Seed[TRANGE] = strdup("111*11111*111*111*11");  }
        else if (TRANGE == 6) { Trange_Seed[TRANGE] = strdup("111111*111111*11");  }
        else if (TRANGE == 5) { Trange_Seed[TRANGE] = strdup("111111*11111*11");  }      // kb1 lost 1
        else if (TRANGE == 4) { Trange_Seed[TRANGE] = strdup("11111*1111");  }       // kb1 lost 1
        else if (TRANGE == 3) { Trange_Seed[TRANGE] = strdup("11111*111");  }
        // else if (TRANGE==2) { Trange_Seed[TRANGE] = strdup("11111*11");  }         // kb1 lost 8
        else if (TRANGE == 2) { Trange_Seed[TRANGE] = strdup("1111*11");  }          // added 1 one
        // else if (TRANGE==1) { Trange_Seed[TRANGE] = strdup("111*11");  } // 111*1
        else if (TRANGE == 1) { Trange_Seed[TRANGE] = strdup("111*1");  } // because its too short sometimes
        else { doCriticalErrorAndQuit("Unsupported TRANGE(%d)", TRANGE); };

        /* SEED+
             if (TRANGE==14){ Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11**111*1*11");  }
        else if (TRANGE==13){ Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11**111");  }
        else if (TRANGE==12){ Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11");  }
        else if (TRANGE==11){ Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE==10){ Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE==9) { Trange_Seed[TRANGE] = strdup("111*11111*111*111*11");  }

        else if (TRANGE==8) { Trange_Seed[TRANGE] = strdup("111*11111*111*111*11");  }
        else if (TRANGE==7) { Trange_Seed[TRANGE] = strdup("111111*111111*11");  }
        else if (TRANGE==6) { Trange_Seed[TRANGE] = strdup("111111*11111*11");  }
        else if (TRANGE==5) { Trange_Seed[TRANGE] = strdup("11111*1111");  }
        else if (TRANGE==4) { Trange_Seed[TRANGE] = strdup("11111*111");  }
        else if (TRANGE==3) { Trange_Seed[TRANGE] = strdup("11111*11");  }
        else if (TRANGE==2) { Trange_Seed[TRANGE] = strdup("111*11");  }
        else if (TRANGE==1) { Trange_Seed[TRANGE] = strdup("111*1");  } // 111*1
        else { doCriticalErrorAndQuit("Unsupported TRANGE(%d)",TRANGE); };
        */

        /* default
             if (TRANGE==14){ Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11**111*1*11");  }
        else if (TRANGE==13){ Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11**111");  }
        else if (TRANGE==12){ Trange_Seed[TRANGE] = strdup("1**11***1111**11*11**1*11*11*111**11");  }
        else if (TRANGE==11){ Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE==10){ Trange_Seed[TRANGE] = strdup("11*111*11111**11**111*11*11");  }
        else if (TRANGE==9) { Trange_Seed[TRANGE] = strdup("111*11111*111*111*11");  }
        else if (TRANGE==8) { Trange_Seed[TRANGE] = strdup("111111*111111*11");  }
        else if (TRANGE==7) { Trange_Seed[TRANGE] = strdup("111111*11111*11");  }
        else if (TRANGE==6) { Trange_Seed[TRANGE] = strdup("11111*1111");  }
        else if (TRANGE==5) { Trange_Seed[TRANGE] = strdup("11111*111");  }
        else if (TRANGE==4) { Trange_Seed[TRANGE] = strdup("11111*11");  }
        else if (TRANGE==3) { Trange_Seed[TRANGE] = strdup("111*11");  }
        else if (TRANGE==2) { Trange_Seed[TRANGE] = strdup("111*1");  }
        else if (TRANGE==1) { Trange_Seed[TRANGE] = strdup("11*1");  } // 111*1
        else { doCriticalErrorAndQuit("Unsupported TRANGE(%d)",TRANGE); };
        */

        Trange_Ones[TRANGE] = 0;

        for (src = Trange_Seed[TRANGE]; *src != '\0'; src++) { if (*src == '1') { Trange_Ones[TRANGE]++;}  }

        Trange_Tuple_Size[TRANGE] = strlen(Trange_Seed[TRANGE]);

        // for arrays only
        if (TRANGE <= MAX_ARRAY_TRANGE) {

            Trange_N_Codes[TRANGE] = pow(4, Trange_Ones[TRANGE]);
            Easy_Array_Tuples[TRANGE]   = (EASY_ARRAY **) scalloc( Trange_N_Codes[TRANGE], sizeof(EASY_ARRAY **) );
            Hash_Tuples[TRANGE]    = NULL;

            // for hashes only
        }
        else {

            Trange_N_Codes[TRANGE] = 0;
            Easy_Array_Tuples[TRANGE]   = NULL;
            Hash_Tuples[TRANGE]    = (GDHASH *) CreateDoubleHash(1000000);

        }
    }

    // structure to quickly find TRANGE from patsize
    for (pmin = 0; pmin <= MAXPROFILESIZE; pmin++) {

        if (pmin >= 1111)     { Trange_Profile_Range[pmin] = 14;  }
        else if (pmin >= 715) { Trange_Profile_Range[pmin] = 13;  }
        else if (pmin >= 474) { Trange_Profile_Range[pmin] = 12;  }
        else if (pmin >= 303) { Trange_Profile_Range[pmin] = 11;   }
        else if (pmin >= 201) { Trange_Profile_Range[pmin] = 10;  }

        else if (pmin >= 134) { Trange_Profile_Range[pmin] = 9;  }
        else if (pmin >= 89)  { Trange_Profile_Range[pmin] = 8;  }
        else if (pmin >= 59)  { Trange_Profile_Range[pmin] = 7;  }
        else if (pmin >= 39)  { Trange_Profile_Range[pmin] = 6;  }
        else if (pmin >= 30)  { Trange_Profile_Range[pmin] = 5;  }
        else if (pmin >= 25)  { Trange_Profile_Range[pmin] = 4;  }
        else if (pmin >= 17)  { Trange_Profile_Range[pmin] = 3;  }
        else if (pmin >= 10)  { Trange_Profile_Range[pmin] = 2;  }
        else if (pmin >= 6)   { Trange_Profile_Range[pmin] = 1; }
        else Trange_Profile_Range[pmin] = 0;

    }


    // create a clusterbase to find connections
    cb = ClusterBaseCreate();


    // create profile lists
    profileList = EasyListCreate(NULL, NULL);
    profileList2 = EasyListCreate(NULL, NULL);


    // read ref profiles into list
    fprintf(stderr, "\nReading reference profiles..."); fflush(stderr);
    off1 = 0;

    while (1) {

        prof1 = NULL; prof1rc = NULL;
        prof1 = ReadProfileWithRC(fpi, off1, &prof1rc, 0);

        if (NULL == prof1 || NULL == prof1rc) break;

        EasyListInsertTail(profileList, prof1);
        EasyListInsertTail(profileList, prof1rc);

        off1 = prof1->nextoffset;
    }

    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));
    fprintf(stderr, "\nLoaded %zu profiles from reference file!", profileList->size / 2); fflush(stderr);


    // read read profiles into list
    fprintf(stderr, "\nReading read profiles..."); fflush(stderr);
    off1 = 0;

    while (1) {

        prof1 = NULL; prof1rc = NULL;
        prof1 = ReadProfileWithRC(fpi2, off1, &prof1rc, 1);

        if (NULL == prof1 || NULL == prof1rc) break;

        EasyListInsertTail(profileList2, prof1);
        EasyListInsertTail(profileList2, prof1rc);

        off1 = prof1->nextoffset;
    }

    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));
    fprintf(stderr, "\nLoaded %zu profiles from reads file!", profileList2->size / 2); fflush(stderr);


    // creating profile lookup hash
    fprintf(stderr, "\nCreating profile lookup hash..."); fflush(stderr);
    {
        unsigned int pread = 0;

        profileHash = CreateSingleHash((profileList->size + profileList2->size));

        for (nof1 = profileList->head; nof1 != NULL; ) {
            profpair = scalloc(1, sizeof(PROFPAIR));
            profpair->prof = (PROFILE *)EasyListItem(nof1);
            profpair->profrc = (PROFILE *)EasyListItem(nof1->next);

            profpair->prof->rotlist = profpair->profrc->rotlist = EasyListCreate(NULL, NULL);
            EasyListInsertTail(profpair->prof->rotlist, profpair);

            if (NULL != GetSingleHashItem(profileHash, profpair->prof->key)) {
                printf("\nERROR: duplicate index %d!", profpair->prof->key);
                exit(1);
            }

            SetSingleHashItem(profileHash, profpair->prof->key, profpair );

            //if ((pread%1000)==0) printf("\nLoaded %u profiles",pread);

            nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
            pread++;
        }

        for (nof1 = profileList2->head; nof1 != NULL; ) {

            profpair = scalloc(1, sizeof(PROFPAIR));
            profpair->prof = (PROFILE *)EasyListItem(nof1);
            profpair->profrc = (PROFILE *)EasyListItem(nof1->next);

            profpair->prof->rotlist = profpair->profrc->rotlist = EasyListCreate(NULL, NULL);
            EasyListInsertTail(profpair->prof->rotlist, profpair);

            if (NULL != GetSingleHashItem(profileHash, profpair->prof->key)) {
                printf("\nERROR: duplicate index %d!", profpair->prof->key);
                exit(1);
            }

            SetSingleHashItem(profileHash, profpair->prof->key, profpair );

            //if ((pread%1000)==0) printf("\nLoaded %u profiles",pread);

            nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
            pread++;
        }
    }
    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));


    // redundcy files
    if (fpirot) {
        fprintf(stderr, "\nUsing redundancy file to speed up alignments..."); fflush(stderr);

        if (0 != LoadRotated(fpirot)) doCriticalErrorAndQuit("Error reading rotated refs file!");

        // if not masterrotated, remove from lists
        tempList = EasyListCreate(NULL, NULL);

        for (nof1 = profileList->head; nof1 != NULL; ) {
            prof1 = (PROFILE *)EasyListItem(nof1);
            prof1rc = (PROFILE *)EasyListItem(nof1->next);

            if (1 == prof1->rotationmaster) {
                EasyListInsertTail(tempList, prof1);
                EasyListInsertTail(tempList, prof1rc);
            }

            nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
        }

        EasyListDestroy(profileList);
        profileList = tempList;

        fprintf(stderr, "new ref list size %zu!", profileList->size / 2);
    }

    if (fpi2rot) {
        fprintf(stderr, "\nUsing redundancy file to speed up alignments..."); fflush(stderr);

        if (0 != LoadRotated(fpi2rot)) doCriticalErrorAndQuit("Error reading rotated reads file!");

        // if not masterrotated, remove from lists
        tempList = EasyListCreate(NULL, NULL);

        for (nof1 = profileList2->head; nof1 != NULL; ) {
            prof1 = (PROFILE *)EasyListItem(nof1);
            prof1rc = (PROFILE *)EasyListItem(nof1->next);

            if (1 == prof1->rotationmaster) {
                EasyListInsertTail(tempList, prof1);
                EasyListInsertTail(tempList, prof1rc);
            }

            nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
        }

        EasyListDestroy(profileList2);
        profileList2 = tempList;

        fprintf(stderr, "new read list size %zu!", profileList2->size / 2);
    }


    // seed references
    fprintf(stderr, "\nSeeding reference profiles' concensus sequences..."); fflush(stderr);

    for (nof1 = profileList->head; nof1 != NULL; ) {

        prof1 = (PROFILE *)EasyListItem(nof1);
        prof1rc = (PROFILE *)EasyListItem(nof1->next);

        pmin = min(prof1->patlen, prof1rc->patlen);

        if (pmin < 7 || pmin > MAXPROFILESIZE) { doCriticalErrorAndQuit("pmin(%d) must be in 7-%d range!", pmin, MAXPROFILESIZE); }

        TRANGE = Trange_Profile_Range[pmin];

        if (TRANGE > 0) {

            // current range
            free_seed_info(&seedstruct);
            retrieve_seed_info(Trange_Seed[TRANGE], &seedstruct);

            if (seedstruct.hasX)
                doCriticalErrorAndQuit("Seeds with X are not allowed in this version. Aborting!");


            sequence = seqdupPlusTuple(prof1->seq->length, prof1->seq->sequence, Trange_Tuple_Size[TRANGE]);
            blaststats.stTuplesProcessed += _HCLUST_process_sequence(sequence, &seedstruct, TRANGE, Trange_Tuple_Size[TRANGE], prof1, 0, pmin);
            sfree(sequence);

            sequence = seqdupPlusTuple(prof1rc->seq->length, prof1rc->seq->sequence, Trange_Tuple_Size[TRANGE]);
            blaststats.stTuplesProcessed += _HCLUST_process_sequence(sequence, &seedstruct, TRANGE, Trange_Tuple_Size[TRANGE], prof1rc, 1, pmin);
            sfree(sequence);

            // lower range
            lowerpat = pmin -  (int)(pmin * PATLEN_SIZE_ERR_FRACTION + TRUNC_ROUND_CEIL);

            if (lowerpat >= 7 && Trange_Profile_Range[lowerpat] <= (TRANGE - 1)) {

                NEWRANGE = TRANGE - 1;

                free_seed_info(&seedstruct);
                retrieve_seed_info(Trange_Seed[NEWRANGE], &seedstruct);

                if (seedstruct.hasX)
                    doCriticalErrorAndQuit("Seeds with X are not allowed in this version. Aborting!");


                sequence = seqdupPlusTuple(prof1->seq->length, prof1->seq->sequence, Trange_Tuple_Size[NEWRANGE]);
                blaststats.stTuplesProcessed += _HCLUST_process_sequence(sequence, &seedstruct, NEWRANGE, Trange_Tuple_Size[NEWRANGE], prof1, 0, pmin);
                sfree(sequence);

                sequence = seqdupPlusTuple(prof1rc->seq->length, prof1rc->seq->sequence, Trange_Tuple_Size[NEWRANGE]);
                blaststats.stTuplesProcessed += _HCLUST_process_sequence(sequence, &seedstruct, NEWRANGE, Trange_Tuple_Size[NEWRANGE], prof1rc, 1, pmin);
                sfree(sequence);
            }


            // higher range
            higherpat = pmin +  (int)(pmin * PATLEN_SIZE_ERR_FRACTION + TRUNC_ROUND_CEIL);

            if (Trange_Profile_Range[higherpat] >= (TRANGE + 1)) {

                NEWRANGE = TRANGE + 1;

                free_seed_info(&seedstruct);
                retrieve_seed_info(Trange_Seed[NEWRANGE], &seedstruct);

                if (seedstruct.hasX)
                    doCriticalErrorAndQuit("Seeds with X are not allowed in this version. Aborting!");


                sequence = seqdupPlusTuple(prof1->seq->length, prof1->seq->sequence, Trange_Tuple_Size[NEWRANGE]);
                blaststats.stTuplesProcessed += _HCLUST_process_sequence(sequence, &seedstruct, NEWRANGE, Trange_Tuple_Size[NEWRANGE], prof1, 0, pmin);
                sfree(sequence);

                sequence = seqdupPlusTuple(prof1rc->seq->length, prof1rc->seq->sequence, Trange_Tuple_Size[NEWRANGE]);
                blaststats.stTuplesProcessed += _HCLUST_process_sequence(sequence, &seedstruct, NEWRANGE, Trange_Tuple_Size[NEWRANGE], prof1rc, 1, pmin);
                sfree(sequence);
            }

        }

        nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
    }

    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));


    // sort
    fprintf(stderr, "\nSorting ref seed arrays..."); fflush(stderr);

    for (TRANGE = 1; TRANGE <= MAX_LIST_TRANGE; TRANGE++) {

        if (TRANGE <= MAX_ARRAY_TRANGE) {

            for (ui = 0; ui < Trange_N_Codes[TRANGE]; ui++) {
                if ((artemp = Easy_Array_Tuples[TRANGE][ui])) {
                    EasyArrayQuickSort(artemp, patcmp);
                }
            }

        }
        else {

            for (iter = Hash_Tuples[TRANGE]->walker->head; iter != NULL; iter = iter->next) {
                ui = *(size_t*)EasyListItem(iter);
                gdrack = Hash_Tuples[TRANGE]->rack[ui];

                while (gdrack) {
                    curr = gdrack;
                    gdrack = gdrack->next;
                    artemp = (EASY_ARRAY *)(curr->data);
                    EasyArrayQuickSort(artemp, patcmp);
                }
            }

        }
    }

    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));


    // Optimizing ref seed arrays and creating ref patsize lookups
    fprintf(stderr, "\nOptimizing ref seed arrays and creating ref patsize lookups..."); fflush(stderr);
    Array_Tuples[0] = NULL;
    Array_Tuples_Length[0] = 0;
    Array_Tuples_Index[0] = NULL;

    for (TRANGE = 1; TRANGE <= MAX_LIST_TRANGE; TRANGE++) {


        if (TRANGE <= MAX_ARRAY_TRANGE) {

            int prev_pat;

            Array_Tuples[TRANGE] = (SEED_HIT **)scalloc(Trange_N_Codes[TRANGE], sizeof(*(Array_Tuples[TRANGE])));
            Array_Tuples_Length[TRANGE] = scalloc(Trange_N_Codes[TRANGE], sizeof(*(Array_Tuples_Length[TRANGE])));
            Array_Tuples_Index[TRANGE]   = (EASY_ARRAY **) scalloc( Trange_N_Codes[TRANGE], sizeof(EASY_ARRAY **) );

            for (ui = 0; ui < Trange_N_Codes[TRANGE]; ui++) {
                if ((artemp = Easy_Array_Tuples[TRANGE][ui])) {
                    Array_Tuples[TRANGE][ui] = (shtemp = (SEED_HIT *)scalloc(artemp->size, sizeof(SEED_HIT)));
                    Array_Tuples_Index[TRANGE][ui] = (iatemp = EasyArrayCreate(2, NULL, NULL));

                    for (uj = 0; uj < artemp->size; uj++) {
                        shtemp[uj] = *(SEED_HIT *)artemp->array[uj];

                        if (uj == 0 || shtemp[uj].pmin != prev_pat) {
                            prev_pat = shtemp[uj].pmin;
                            EasyArrayInsert(iatemp, &shtemp[uj]);
                        }
                    }

                    Array_Tuples_Length[TRANGE][ui] = artemp->size;
                    EasyArrayDestroy(artemp);
                }
            }

            sfree(Easy_Array_Tuples[TRANGE]);
            Easy_Array_Tuples[TRANGE] = NULL;

        }
        else {

            int prev_pat;

            Array_Tuples[TRANGE] = NULL;
            Array_Tuples_Length[TRANGE] = 0;
            Array_Tuples_Index[TRANGE] = NULL;

            for (iter = Hash_Tuples[TRANGE]->walker->head; iter != NULL; iter = iter->next) {
                ui = *(size_t*)EasyListItem(iter);

                gdrack = Hash_Tuples[TRANGE]->rack[ui];

                while (gdrack) {
                    curr = gdrack;
                    gdrack = gdrack->next;
                    artemp = (EASY_ARRAY *)(curr->data);


                    (shtemp = (SEED_HIT *)scalloc(artemp->size, sizeof(SEED_HIT)));
                    iatemp = EasyArrayCreate(2, NULL, NULL);

                    for (uj = 0; uj < artemp->size; uj++) {
                        shtemp[uj] = *(SEED_HIT *)artemp->array[uj];

                        if (uj == 0 || shtemp[uj].pmin != prev_pat) {
                            prev_pat = shtemp[uj].pmin;
                            EasyArrayInsert(iatemp, &shtemp[uj]);
                        }
                    }

                    curr->length = artemp->size;
                    EasyArrayDestroy(artemp);
                    curr->data = (void *)shtemp;
                    curr->index = (void *)iatemp;
                }
            }


        }
    }

    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));


    // Rehashing ref hashes to optimal sizes
    fprintf(stderr, "\nRehashing ref hashes to optimal sizes..."); fflush(stderr);

    for (TRANGE = MAX_ARRAY_TRANGE + 1; TRANGE <= MAX_LIST_TRANGE; TRANGE++) {

        gdtemp = (GDHASH *) CreateDoubleHash(Hash_Tuples[TRANGE]->nitems * 2);

        for (iter = Hash_Tuples[TRANGE]->walker->head; iter != NULL; iter = iter->next) {
            ui = *(size_t*)EasyListItem(iter);
            gdrack = Hash_Tuples[TRANGE]->rack[ui];

            while (gdrack) {
                curr = gdrack;
                gdrack = gdrack->next;
                SetDoubleHashItem(gdtemp, curr->key1, curr->key2, curr->data, curr->index, curr->length);
            }
        }

        DestroyDoubleHash(Hash_Tuples[TRANGE], NULL);
        Hash_Tuples[TRANGE] = gdtemp;
    }

    fprintf(stderr, "(total time: %.1lf secs)", (double)(time(NULL) - startTime));


    // scanning reads
    fprintf(stderr, "\nScanning reads..."); fflush(stderr);
    readcount = 0; readpercent = 0;

    for (nof1 = profileList2->head; nof1 != NULL; ) {

        prof1 = (PROFILE *)EasyListItem(nof1);
        prof1rc = (PROFILE *)EasyListItem(nof1->next);

        pmin = min(prof1->patlen, prof1rc->patlen);

        if (pmin < 7 || pmin > MAXPROFILESIZE) { doCriticalErrorAndQuit("pmin(%d) must be in 7-%d range!", pmin, MAXPROFILESIZE); }

        TRANGE = Trange_Profile_Range[pmin];

        if (TRANGE > 0) {

            EASY_ARRAY *candidates1 = EasyArrayCreate(100, NULL, NULL);
            EASY_ARRAY *candidates2 = EasyArrayCreate(100, NULL, NULL);

            free_seed_info(&seedstruct);
            retrieve_seed_info(Trange_Seed[TRANGE], &seedstruct);

            if (seedstruct.hasX)
                doCriticalErrorAndQuit("Seeds with X are not allowed in this version. Aborting!");

            //fprintf(stderr,"\nreadcount: %d => candidates1: %d candidates2: %d (key: %d, pmin: %d, TRANGE: %d, read: %s, readrc: %s)",readcount/2,candidates1->size,candidates2->size,prof1->key,pmin,TRANGE,prof1->seq->sequence,prof1rc->seq->sequence); fflush(stderr);

            sequence = seqdupPlusTuple(prof1->seq->length, prof1->seq->sequence, Trange_Tuple_Size[TRANGE]);
            _HCLUST_find_candidates(sequence, prof1->seq, &seedstruct, TRANGE, Trange_Tuple_Size[TRANGE], 0, pmin, candidates1);
            sfree(sequence);

            sequence = seqdupPlusTuple(prof1rc->seq->length, prof1rc->seq->sequence, Trange_Tuple_Size[TRANGE]);
            _HCLUST_find_candidates(sequence, prof1rc->seq, &seedstruct, TRANGE, Trange_Tuple_Size[TRANGE], 1, pmin, candidates2);
            sfree(sequence);

            EasyArrayQuickSort(candidates1, candsort);
            EasyArrayQuickSort(candidates2, candsort);


            _HCLUST_align_candidates(prof1, 0, pmin, candidates1, dt1, cb);
            _HCLUST_align_candidates(prof1rc, 1, pmin, candidates2, dt1, cb);

            EasyArrayDestroy(candidates1);
            EasyArrayDestroy(candidates2);
        }

        nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
        readcount += 2;

        if ( (int) (100 * readcount / (double) profileList2->size) >= readpercent)  { fprintf(stderr, "\n%d%% (ProfileAlignments: %llu, total time: %.1lf secs)", readpercent, blaststats.stProfileAlignments, (double)(time(NULL) - startTime)); fflush(stderr); readpercent++;  }
    }


    fprintf(stderr, "\nscan finished!\n\t"
            "CandidatesConsidered: %llu\n\t"
            "LCSalignments: %llu\n\t"
            "ProfileAlignments: %llu\n\t"
            "ProfileAlignmentsSuccess: %llu\n\t"
            "pairs confirmed after flank LCS alignments: %llu\n\t"
            "total clustered: %llu\n\t"
            "total time: %.1lf secs\n\t"
            "writing results ...",

            blaststats.stCandidatesConsidered,
            blaststats.stLCSalignments,
            blaststats.stProfileAlignments,
            blaststats.stProfileAlignmentsSuccess,
            blaststats.stPairsConfirmedAfterFlankAlignments,
            cb->trcount,
            (double)(time(NULL) - startTime));

    fflush(stderr);

    return cb;
}

/*******************************************************************************************/
int main(int argc, char **argv)
{
    FILE      *fpiRefs, *fpirot, *fpiReads, *fpi2rot, *fpo;
    int       cutoff2;
    unsigned char *dt1;
    CLUSTERBASE   *cb;
    char      inputfile[1000] = "", inputfile2[1000] = "", *edgesIn = NULL;
    char      buffer[1000] = "", distfile1[1000] = "", outputfile[1000] = "", mapfile[1000] = "", edgesfile[1000] = "";

//printf("\nseed_hit size: %d",sizeof(SEED_HIT));
//exit(1);

    // verify parameter count
    if (argc < 7) {
        printf("\nPSEARCH v1.91 - Finds clusters from list of profiles.");
        printf("\nauthors: Yevgeniy Gefland, Alfredo Rodriguez, Gary Benson");
        printf("\n\nUsage:  %s REFPROFILES READPROFILES DISTANCEFILE CUTOFF(70-100) MAXERRORS MAXFLANKCONSIDERED [EDGESINFILE] [OPTION]\n\n", argv[0]);
        printf("\tOPTIONS: \n");
        printf("\t\t-r ref-to-ref alignments, only one flank has to match, 50 is default length but may be followed by integer flank length, ex -r 50 \n");
        printf("\t\t-m produce map file  \n");
        printf("\n");
        exit(1);
    }

    // process command line
    strcpy(inputfile, argv[1]);
    strcpy(inputfile2, argv[2]);
    strcpy(distfile1, argv[3]);
    sprintf(outputfile, "%s.clu", inputfile2);
    sprintf(mapfile, "%s.map", inputfile2);
    sprintf(edgesfile, "%s.edges", inputfile2);
    sscanf(argv[4], "%d", &cutoff2);

    if (0 == cutoff2 || cutoff2 < 70 || cutoff2 > 100) {
        printf("\nERROR: Please enter similarity cutoff (70-100)!");
        exit(1);
    }

    TARGET_HOMOLOGY = cutoff2 / 100.0;
    MAXERRORS = atoi(argv[5]);
    MAXFLANKCONSIDERED = min ( 1000000, atoi(argv[6]) ); // Shortening of flanks for for ref-to-read. Set this to anything and use -r [reflen] for ref-to-ref flank shortening.


    // is option the next to last parameter?
    if (argc >= 8 &&  (0 == strcmp(argv[7], "-r") || 0 == strcmp(argv[7], "-R"))) {

        OPTION = 'R';

        if (argc >= 9 && 0 != atoi(argv[8])) { REFLEN = atoi(argv[8]);}

        // if not, is it the EDGESIN file?
    }
    else {

        edgesIn = ( argc >= 8 ) ? argv[7] : NULL;

        if (edgesIn && 0 != strcmp(edgesIn + strlen(edgesIn) - 7, ".edgein")) {
            printf("\nERROR: Please make sure edgesin file has extension .edgein!");
            exit(1);
        }
    }

    // -r
    if (argc >= 9 &&  (0 == strcmp(argv[8], "-r") || 0 == strcmp(argv[8], "-R"))) {

        OPTION = 'R';

        if (argc >= 10 && 0 != atoi(argv[9])) { REFLEN = atoi(argv[9]);}

    }

    if (argc >= 10 &&  (0 == strcmp(argv[9], "-r") || 0 == strcmp(argv[9], "-R"))) {

        OPTION = 'R';

        if (argc >= 11 && 0 != atoi(argv[10])) { REFLEN = atoi(argv[10]);}

    }

    // -m
    if (argc >= 8 &&  (0 == strcmp(argv[7], "-m") || 0 == strcmp(argv[7], "-M"))) {

        OPTION2 = 'M';

    }

    if (argc >= 9 &&  (0 == strcmp(argv[8], "-m") || 0 == strcmp(argv[8], "-M"))) {

        OPTION2 = 'M';

    }

    if (argc >= 10 &&  (0 == strcmp(argv[9], "-m") || 0 == strcmp(argv[9], "-M"))) {

        OPTION2 = 'M';

    }


    // open map file
    if (OPTION2 == 'M') {
        MAPFILE = fopen(mapfile, "w");

        if (MAPFILE == NULL) {
            fprintf(stderr, "\nERROR: Unable to open mapfile file '%s'", mapfile);
            exit(1);
        }
    }



    /* 0 for maxerror only allowed for ref-read, not implemented for other mode */
    if (MAXERRORS != 0 && OPTION == 'R') {
        fprintf(stderr, "\nERROR: maxerror must be 0 because it will be dynamically generated based on flank length.");
        exit(1);
    }


    // set parameters
    startTime = time(NULL);
    TARGET_HOMOLOGY_INT = (int)(TARGET_HOMOLOGY * 100);

    // open the ref input file for reading
    fpiRefs = fopen(inputfile, "r");

    if (fpiRefs == NULL) {
        fprintf(stderr, "\nERROR: Unable to open input file '%s'", inputfile);
        exit(1);
    }

    // open the reads input file for reading
    fpiReads = fopen(inputfile2, "r");

    if (fpiReads == NULL) {
        fprintf(stderr, "\nERROR: Unable to open input file '%s'", inputfile2);
        exit(1);
    }


    // load distance table
    dt1 = LoadDistanceTable(distfile1);

    if (dt1 == NULL) {
        fprintf(stderr, "\nERROR: Unable to load distance table '%s'", distfile1);
        exit(1);
    }


    // just doing edges?
    if (edgesIn) {

        FILE *fpeIn = fopen(edgesIn, "r");

        if (NULL == fpeIn) {
            doCriticalErrorAndQuit("Unable to open input file '%s'!", edgesIn);
        }

        doEdgesBruteForce(fpiRefs, fpiReads, fpeIn, dt1, 1, edgesfile);
        fclose(fpeIn);

        exit(0);
    }


    // open redundancy index
    sprintf(buffer, "%s.rotindex", inputfile);
    fpirot = fopen(buffer, "r");
    sprintf(buffer, "%s.rotindex", inputfile2);
    fpi2rot = fopen(buffer, "r");


    cb = doSearchSimilarities(fpiRefs, fpiReads, fpirot, fpi2rot, dt1);


    fclose(fpiRefs);
    fclose(fpiReads);

    if (fpirot) fclose(fpirot);

    if (fpi2rot) fclose(fpi2rot);

    if (MAPFILE) fclose(MAPFILE);

    // print clusters to output file
    fpo = fopen(outputfile, "w");

    if (fpo == NULL) {
        fprintf(stderr, "\nERROR: Unable to create output file '%s'", outputfile);
        exit(1);
    }

    ClusterBasePrint(cb, fpo);
    fclose(fpo);
    fprintf(stderr, "done (total time: %.1lf secs)\n\n", (double)(time(NULL) - startTime)); fflush(stderr);

    return 0;
}

/*******************************************************************************************/
int LoadRotated(FILE *fp)
{


    int repeatcount, direction, i = 1, repeatkey;
    char symbol, firstdir, newdir;
    PROFPAIR *profpair, *temppair;

    if (NULL == fp) return 0; /* no rotated file is fine */

    /* loop for every cluster */
    while (1) {
        repeatcount  = 0;
        firstdir = 0;
        profpair = NULL; /* indicates master repeat */

        /* read the next symbol (ignore white space) */
        i = fscanf(fp, " %c", &symbol);

        if (i != 1) break;

        /* put character back */
        ungetc(symbol, fp);

        //printf("\n");

        /* loop for every member */
        while (1) {
            /* read number */
            i = fscanf(fp, "%d", &repeatkey);

            if (i != 1) {
                doCriticalErrorAndQuit("Invalid format detected in the redundancy file!");
            }

            repeatcount++;

            //printf(" %d",repeatkey);

            /* read symbol */
            i = fscanf(fp, " %c", &symbol);

            if (i != 1 || (symbol != '\'' && symbol != '"')) {
                doCriticalErrorAndQuit("\n\nUnexpected value qualifier detected in the redundancy file!");
            }

            /* assign direction */
            if (symbol == '\'')
                direction = 0;
            else
                direction = 1;

            /* get the master repeat from the hash, create a list on it */
            if (1 == repeatcount) {
                if (NULL == (profpair = GetSingleHashItem(profileHash, repeatkey )))
                    doCriticalErrorAndQuit("\n\nError looking up key %d in profileHash!", repeatkey);

                firstdir = direction;

                profpair->prof->dir = 0;
                profpair->profrc->dir = 0;

            }
            else {

                if (NULL == (temppair = GetSingleHashItem(profileHash, repeatkey )))
                    doCriticalErrorAndQuit("\n\nError looking up key %d in profileHash!", repeatkey);

                EasyListInsertTail(profpair->prof->rotlist, temppair);
                temppair->prof->rotationmaster = 0;
                temppair->profrc->rotationmaster = 0;

                newdir = (firstdir) ?  (!direction) : (direction) ;
                temppair->prof->dir = newdir;
                temppair->profrc->dir = newdir;
            }


            /* read whitespace character */
            symbol = ' ';

            while (symbol == ' ' || symbol == '\t') {
                symbol = 0;
                i = fscanf(fp, "%c", &symbol);

                if (i != 1) { break; }
            }

            /* if no more characters then break */
            if (i != 1 || symbol == 13 || symbol == 10) break;

            /* put character back */
            ungetc(symbol, fp);

            /* if character was a digit continue otherwise break */
            if ((symbol >= '0' && symbol <= '9') || symbol == '-') continue;
            else break;
        }

        if (i != 1) break;

    }

    return 0;

}

/*******************************************************************************************/
void doCriticalErrorAndQuit(const char *format, ... )
{

    va_list                         argp;


    fprintf(stderr, "\n\nERROR!!!!:\n");

    if (format == NULL) exit(0);;

    va_start(argp, format);
    vfprintf(stderr, format, argp);
    fflush(stderr);
    va_end(argp);

    exit(1);
}

/*******************************************************************************************/
int doEdgesBruteForce(FILE *fpi, FILE *fpi2, FILE *edgesIn, unsigned char *dt1, int fixed, char *outputfile)
{

    char buffer[1010] = "";
    double bestsim, totsim = 0;
    __int64 count = 0;
    EASY_LIST *profileList = NULL, *profileList2 = NULL;
    EASY_NODE *nof1, *nof2;
    FILE *fpo;
    int dir, off1, id1, id2;
    PROFILE *prof1, *prof2, *prof1rc, *prof2rc;
    PAP *pap1, *pap2, *bestpap;
    GSHASH *profileHash = NULL;
    PROFPAIR *profpair, *profpair2;

    // open output file
    fpo = fopen(outputfile, "w");

    if (fpo == NULL) {
        fprintf(stderr, "\nERROR: Unable to create output file '%s'", outputfile);
        exit(1);
    }

    // create lists
    profileList = EasyListCreate(NULL, NULL);
    profileList2 = EasyListCreate(NULL, NULL);

    // read 1st file profiles into list
    off1 = 0;

    while (1) {

        prof1 = NULL; prof1rc = NULL;
        prof1 = ReadProfileWithRC(fpi, off1, &prof1rc, 0);

        if (NULL == prof1 || NULL == prof1rc) break;

        EasyListInsertTail(profileList, prof1);
        EasyListInsertTail(profileList, prof1rc);

        off1 = prof1->nextoffset;
    }

    fprintf(stdout, "\nLoaded %zu profiles from 1st file", profileList->size / 2);

    // read 2nd file profiles into list
    off1 = 0;

    while (1) {

        prof1 = NULL; prof1rc = NULL;
        prof1 = ReadProfileWithRC(fpi2, off1, &prof1rc, 0);

        if (NULL == prof1 || NULL == prof1rc) break;

        EasyListInsertTail(profileList2, prof1);
        EasyListInsertTail(profileList2, prof1rc);

        off1 = prof1->nextoffset;
    }

    fprintf(stdout, "\nLoaded %zu profiles from 2nd file", profileList2->size / 2);



    // crate profile hash
    profileHash = CreateSingleHash((profileList->size + profileList2->size));

    for (nof1 = profileList->head; nof1 != NULL; ) {
        profpair = scalloc(1, sizeof(PROFPAIR));
        profpair->prof = (PROFILE *)EasyListItem(nof1);
        profpair->profrc = (PROFILE *)EasyListItem(nof1->next);


        if (NULL != GetSingleHashItem(profileHash, profpair->prof->key)) {
            fprintf(stderr, "\nERROR: duplicate index!");
            exit(1);
        }

        SetSingleHashItem(profileHash, profpair->prof->key, profpair );
        //printf("\nSetting hash id: %d",profpair->prof->key);
        //if ((pread%1000)==0) printf("\nLoaded %u profiles",pread);

        nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
    }

    for (nof1 = profileList2->head; nof1 != NULL; ) {

        profpair = scalloc(1, sizeof(PROFPAIR));
        profpair->prof = (PROFILE *)EasyListItem(nof1);
        profpair->profrc = (PROFILE *)EasyListItem(nof1->next);

        if (NULL != GetSingleHashItem(profileHash, profpair->prof->key)) {
            fprintf(stderr, "\nERROR: duplicate index!");
            exit(1);
        }

        //printf("\nSetting hash id: %d",profpair->prof->key);
        SetSingleHashItem(profileHash, profpair->prof->key, profpair );

        //if ((pread%1000)==0) printf("\nLoaded %u profiles",pread);

        nof1 = nof1->next; if (nof1 != NULL) {nof1 = nof1->next;}
    }

    // for each edge do an alignment
    while ( fgets(buffer, 1000, edgesIn)) {
        if (2 == sscanf(buffer, "%d,%d", &id1, &id2)) {
            profpair = GetSingleHashItem(profileHash, id1);
            profpair2 = GetSingleHashItem(profileHash, id2);

            if (NULL == profpair) {
                fprintf(stderr, "\nERROR: could not lookup profile id: %d!", id1);
                exit(1);
            }

            if (NULL == profpair2) {
                fprintf(stderr, "\nERROR: could not lookup profile id: %d!", id2);
                exit(1);
            }

            bestsim = 0.0;
            {

                pap1 = GetFixedProfileAP(profpair->prof, profpair2->prof, dt1);
                pap2 = GetFixedProfileAP(profpair->profrc, profpair2->prof, dt1);

                if (pap1->similarity >= pap2->similarity) {
                    bestpap = pap1;
                    dir = 0;
                }
                else {
                    bestpap = pap2;
                    dir = 1;
                }

                bestsim = bestpap->similarity;

                FreePAP(pap1);
                FreePAP(pap2);

                // now doing all 4 ways
                pap1 = GetFixedProfileAP(profpair->prof, profpair2->profrc, dt1);
                pap2 = GetFixedProfileAP(profpair->profrc, profpair2->profrc, dt1);

                if (pap1->similarity >= pap2->similarity) {
                    bestpap = pap1;
                }
                else {
                    bestpap = pap2;
                }

                bestsim = max(bestsim, bestpap->similarity);

                FreePAP(pap1);
                FreePAP(pap2);



            }

            count++;
            totsim += bestsim;
            fprintf(fpo, "%d%c,%d,%.4lf\n", profpair->prof->key, (dir == 0) ? '\'' : '\"', profpair2->prof->key,  bestsim);

        }
    }

    fprintf(fpo, "#ave: %.4lf\n", totsim / count);
    fprintf(stdout, "\nAlignment complete\n");

}


