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

/****************************************************************
 *   patupdt.h   :   Given a tandem repeat this routine recomputes
 *                   the consensus to try to find a better scoring
 *                   one.
 *
 *   Developed by:   Alfredo Rodriguez
 *   Last Update :   3/25/2002
 *
 *****************************************************************/

#ifndef PATUPDT_H
#define PATUPDT_H

#include "align.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    short int A;
    short int C;
    short int G;
    short int T;
    short int dash;
    short int other;
} COMPOSIT;

typedef struct tagREPEAT {
    int            updates;
    int            concensussize;
    int            firstindex;
    int            lastindex;
    int            percentmatch;
    int            percentindels;
    int            score;
    int            countmatch;
    int            countindel;
    int            countmismatch;
    unsigned char *pattern;
    unsigned char *subsequence;
    int            proflen;
    COMPOSIT *     profile;
    int            gappedlength;
} REPEAT;

typedef struct {
    short int A;
    short int C;
    short int G;
    short int T;
    short int dash;
    short int other;
    char      consensus;
} COMPOSITION;

/*************************************************************
 * This routine returns one to indicate an update has occurred
 * return minus one to indicate memory error and returns 0
 * if patterns is already optimal.
 **************************************************************/
int ComputeBetterPattern( REPEAT *repeat, int alpha, int beta, int *sm ) {
    int           update, tryagain, i, count, newlength, match, indel, mismatch;
    WDPALIGNPAIR *wdpap, *newwdpap;
    GAPPEDALIGNPAIR *gap;
    COMPOSITION *    comps, *pcomp;
    char *           ptr, *newpattern, *ptrpat;
    COMPOSIT *       pprof;

    /* try to update pattern until it doesn't get any better */
    repeat->updates = update = 0;
    tryagain                 = 1;
    while ( tryagain ) {

        /* get a WDP alignment */
        wdpap = GetWDPGlobalAlignPair( (char *) repeat->subsequence,
          (char *) repeat->pattern, repeat->lastindex - repeat->firstindex + 1,
          repeat->concensussize, alpha, beta, sm );
        if ( wdpap == NULL ) {
            return -1;
        }

        /* recalculation of %match and %indel to fix data from TRF */
        ptr      = wdpap->sequenceside;
        ptrpat   = wdpap->patternside;
        match    = 0;
        indel    = 0;
        mismatch = 0;
        while ( *ptr != '\0' ) {
            if ( *ptr == *ptrpat )
                match++;
            else {
                if ( *ptr != '-' && *ptrpat != '-' )
                    mismatch++;
                else
                    indel++;
            }
            ptr++;
            ptrpat++;
        }

        repeat->percentmatch  = match * 100 / wdpap->length;
        repeat->percentindels = indel * 100 / wdpap->length;
        repeat->score         = wdpap->score;
        repeat->countmatch    = match;
        repeat->countmismatch = mismatch;
        repeat->countindel    = indel;

        /* get a gapped alignment */
        gap = GetGappedAlignPair( wdpap );

        /* allocate a composition row as zeroes */
        comps = calloc( gap->gappedpatternlength, sizeof( COMPOSITION ) );
        if ( NULL == comps )
            return -1;

        for ( i = 0, ptr = gap->sequenceside; *ptr != '\0'; ptr++ ) {
            /* increase composition for given column */
            switch ( *ptr ) {
            case 'A':
                comps[i].A++;
                break;
            case 'C':
                comps[i].C++;
                break;
            case 'G':
                comps[i].G++;
                break;
            case 'T':
                comps[i].T++;
                break;
            case '-':
                comps[i].dash++;
                break;
            default:
                comps[i].other++;
                break;
            }
            /* move index */
            i++;
            if ( i == gap->gappedpatternlength )
                i = 0;
        }
        /* calculate the concensus for each column */
        for ( i = 0; i < gap->gappedpatternlength; i++ ) {
            count = -1;
            if ( comps[i].A > count ) {
                count              = comps[i].A;
                comps[i].consensus = 'A';
            }
            if ( comps[i].C > count ) {
                count              = comps[i].C;
                comps[i].consensus = 'C';
            }
            if ( comps[i].G > count ) {
                count              = comps[i].G;
                comps[i].consensus = 'G';
            }
            if ( comps[i].T > count ) {
                count              = comps[i].T;
                comps[i].consensus = 'T';
            }
            if ( comps[i].dash > count ) {
                count              = comps[i].dash;
                comps[i].consensus = '-';
            }
            if ( comps[i].other > count ) {
                count              = comps[i].other;
                comps[i].consensus = 'N';
            }
        }

        /* free the repeat profile if it existed */
        if ( repeat->profile != NULL ) {
            free( repeat->profile );
            repeat->profile = NULL;
        }

        /* allocate new profile */
        repeat->proflen      = gap->gappedpatternlength;
        repeat->gappedlength = gap->gappedlength;
        repeat->profile =
          (COMPOSIT *) malloc( sizeof( COMPOSIT ) * repeat->proflen );
        if ( NULL == repeat->profile )
            return -1;

        /* fill in profile with compositions from above */
        pprof = repeat->profile;
        pcomp = comps;
        for ( i = 0; i < repeat->proflen; i++ ) {
            memcpy( pprof, pcomp, sizeof( COMPOSIT ) );
            pprof++;
            pcomp++;
        }

        /* count how many have something other than a dash at that position */
        for ( i = 0, newlength = 0; i < gap->gappedpatternlength; i++ ) {
            if ( comps[i].consensus != '-' )
                newlength++;
        }

        /* allocate and fill in new pattern */
        newpattern = (char *) malloc( newlength + 1 );
        for ( i = 0, ptr = newpattern; i < gap->gappedpatternlength; i++ ) {
            if ( comps[i].consensus != '-' ) {
                *ptr = comps[i].consensus;
                ptr++;
            }
        }
        *ptr = '\0';

        /* only proceed if pattern obtained is not the original pattern */
        if ( strcmp( newpattern, (char *) repeat->pattern ) ) {
            /* get an alignment */
            newwdpap = GetWDPGlobalAlignPair( (char *) repeat->subsequence,
              newpattern, repeat->lastindex - repeat->firstindex + 1, newlength,
              alpha, beta, sm );

            if ( NULL == newwdpap )
                return -1;

            /* if score is higher update repeat */
            if ( newwdpap->score > wdpap->score ) {
                /* swap old pattern with new */
                ptr             = newpattern;
                newpattern      = (char *) repeat->pattern;
                repeat->pattern = (unsigned char *) ptr;

                /* fix length and score */
                repeat->concensussize = newlength;
                repeat->score         = newwdpap->score;
                repeat->updates++;

                update = 1;
            } else /* if score is not better don't try again */
            {
                tryagain = 0;
            }
            FreeWDPAlignPair( newwdpap );
        } else /* don't try again if pattern is the same */
        {
            tryagain = 0;
        }

        /* free all the memory allocated */
        FreeWDPAlignPair( wdpap );
        FreeGappedAlignPair( gap );
        free( comps );
        free( newpattern );
    }

    return update; /* 1 if repeat has been adjusted, 0 otherwise */
}

#endif
