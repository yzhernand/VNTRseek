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
*   Please direct all questions related to this program or bug reports or updates to
*   Yevgeniy Gelfand (ygelfand@bu.edu)
*
****************************************************************/

/*************************************************************************
*
*   align.h :   Quick Alignment Routines
*               Copyright (c) 2000, Dr. Gary Benson
*               All rights reserved.
*
*   Project Development : Alfredo Rodriguez
*   Library Version     : 1.5L
*   Last Revision Date  : 11/13/2002
*
*   These routines are light-weight versions of those found in align.h.
*   There is also only a subset of the routines implemented here. The
*   main reason for having this library is to minimize memory requirements
*   and execution time by eliminating some of the generality in the
*   routines. Having smaller memory requirements means that larger
*   alignments are possible.
*   A maximum length alignment in align.h requires 12x2,000x100,000 ~
*   2.24GB wich is not possible under 32bit Windows. The same alignment
*   requires only 5x2,000x100,000 ~ 954MB when using this version.
*   As an added benefit the code for wraparound alignment is simpler
*   an executes faster.
*   
*   NOTE: all DNA data passed in must be uppercased. (Gelfand 04/20/2004)
*************************************************************************/

#ifndef ALIGN_H
#define ALIGN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>

#define MAX_ALPHABETSIZ 25
int ALPHABETSIZ = MAX_ALPHABETSIZ;              /* the alphabet size for compositional alignment */



/*************************************************************************
*                      -type definitions-
*************************************************************************/

typedef struct
{
    char*   sequence1side;
    char*   sequence2side;
    int     sequence1start;
    int     sequence1end;
    int     sequence2start;
    int     sequence2end;
    int     score;
    int     length;
}   BASICALIGNPAIR;


typedef struct
{
    char*   sequenceside;
    char*   patternside;
    int     sequencestart;
    int     sequenceend;
    int     patternstart;
    int     patternend;
    int     score;
    int     length;
    int     patternlength;
}   WDPALIGNPAIR;

typedef struct
{
    char*   patternside;
    char*   sequenceside;
    int     gappedlength;
    int     gappedpatternlength;
}   GAPPEDALIGNPAIR;

typedef struct
{
    char*   pattern1side;
    char*   pattern2side;
    int     pattern1start;
    int     pattern1end;
    int     pattern2start;
    int     pattern2end;
    int     score;
    int     length;
}   CYCLICALIGNPAIR;


/* added by Gary Benson 1/2003 */
typedef struct
{
    char*   sequence1side;
    char*   sequence2side;
	char*   matchboundaries;
    int     sequence1start;
    int     sequence1end;
    int     sequence2start;
    int     sequence2end;
    int     score;
    int     length;
}   COMPOSITIONALIGNPAIR;

typedef struct
{
	int C[MAX_ALPHABETSIZ];
} COMPOSITIONVECTOR;

/* #pragma pack is MSVC++ specific */

#pragma pack(push,1)

typedef struct
{
    int		 score;
    char     direction;
}   SD;

#pragma pack(pop)



/************************************************************************
*                       - function prototypes -
*************************************************************************/
BASICALIGNPAIR* GetBasicLocalAlignPair (char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix,SD **Sret);

BASICALIGNPAIR* GetBasicGlobalAlignPair(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix,SD **Sret);

BASICALIGNPAIR* GetBasicPatternGlobalTextLocalAlignPair( char* sequence, char* pattern,
                                    int sequencelength, int patternlength,
                                    int alpha, int indel, int* submatrix);


WDPALIGNPAIR* GetWDPLocalAlignPair( char* sequence, char* pattern,
                                    int sequencelength, int patternlength,
                                    int alpha, int indel, int* submatrix);


WDPALIGNPAIR*    GetWDPGlobalAlignPair(char* sequence, char* pattern,
                                    int sequencelength, int patternlength,
                                    int alpha, int indel, int* submatrix);

GAPPEDALIGNPAIR*   GetGappedAlignPair(WDPALIGNPAIR* wdpap);

CYCLICALIGNPAIR*    GetCyclicDiffAlignPair (char* pat1, char* pat2,
                                            int len1, int len2);


/* added by Gary Benson 1/2003 */
COMPOSITIONALIGNPAIR* GetCompositionGlobalAlignPair(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix, int* M);

/* added by Gary Benson 1/2003 */
COMPOSITIONALIGNPAIR* GetCompositionLocalAlignPair(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix, int* M);


/* added by Gary Benson 1/2003 */
int* GetCompositionSubstringMatchLengths(char* seq1, char* seq2, int len1, int len2,
																				 int limit, int ALPHABETSIZE);

void FreeBasicAlignPair(BASICALIGNPAIR* ap);

void FreeWDPAlignPair(WDPALIGNPAIR* ap);

void FreeGappedAlignPair(GAPPEDALIGNPAIR* ap);

void FreeCyclicAlignPair(CYCLICALIGNPAIR* ap);

/* added by Gary Benson 1/2003 */
void FreeCompositionAlignPair(COMPOSITIONALIGNPAIR* ap);

char*   GetReverseComplement( char* original);

int*    CreateSubstitutionMatrix(int match,int mismatch);

int*    CreateSubstitutionMatrixGT(int match,int mismatch,int gtmatch);

int     GetCharsNotDash(char*);

char*   GetLeftRotatedString(char* original, int count);

void    PrintGappedAlignPair(FILE *fp, GAPPEDALIGNPAIR* gappedalignpair);


/************************************************************************
*                   - Macros, Constants, and Globals -
*************************************************************************/

/* contants, and macros */
#define ZERO            0
#define START           0
#define RIGHT           1
#define DOWN            2
#define DIAGONAL        3
#define WRAPRIGHT       4
#define WRAPDIAGONAL    5
#define LONGDIAGONAL    6
#define DIRECTIONNONE	99
#define MINIMUM(a,b) (a<=b?a:b)
#define MAXIMUM(a,b) (a>=b?a:b)
#define min3switch(a,b,c) ((a<=b)?((a<=c)?1:3):((b<=c)?2:3))
#define max3switch(a,b,c) ((a>=b)?((a>=c)?1:3):((b>=c)?2:3))
#define max4switch(a,b,c,d) ((a>=b)?((a>=c)?((a>=d)?1:4):((c>=d)?3:4)):\
                            ((b>=c)?((b>=d)?2:4):((c>=d)?3:4)))

/* maxpatsize must include +1. some routines need one extra element */
#define MAXPATSIZE      2500
#define MAXSEQSIZE      120000
#define VLN             999999



/***********************************************************************
*                   - Function Definitions -
************************************************************************/


int*    CreateSubstitutionMatrix(int match,int mismatch)
{
    int *sm,*currint;
    int  i,j;

    /* allocate memory for matrix */
    sm = (int*) malloc(26*26*sizeof(int));
    if(sm==NULL) return NULL;

    /* generate 26x26 into matrix */
    for(i=1,currint=sm;i<=26;i++)
    {
        for(j=1;j<=26;j++,currint++)
        {
            *currint=mismatch;
			if(i==j && (i==('A'-'A'+1) || i==('C'-'A'+1) || i==('G'-'A'+1) || i==('T'-'A'+1))) {
				*currint=match;
			}
        }
    }

    /* return matrix */
    return sm;
}

/***********************************************************************/

int*    CreateSubstitutionMatrixGT(int match,int mismatch, int gtmatch)
{
    int *sm,*currint;
    int  i,j;

    /* allocate memory for matrix */
    sm = (int*) malloc(26*26*sizeof(int));
    if(sm==NULL) return NULL;

    /* generate 26x26 into matrix */
    for(i=1,currint=sm;i<=26;i++)
    {
        for(j=1;j<=26;j++,currint++)
        {
            if(i==j) *currint=match;
            else *currint=mismatch;
        }
    }

    /* add the GT match */
    /* Note that the second character is complemented!!! */
    sm[26*('G'-'A')+('A'-'A')] = gtmatch;
    sm[26*('T'-'A')+('C'-'A')] = gtmatch;

    /* return matrix */
    return sm;
}

/***********************************************************************/
BASICALIGNPAIR* GetBasicPatternGlobalTextLocalAlignPair( char* seq1, char* seq2,
                                    int len1, int len2,
									int alpha, int indel, int* submatrix) 
{

	int     col, row,length;
    int     eb;             /* equals e+indel */
    int     fb;             /* equals f+indel */
    int     sa;             /* equals s+match or mismatch */
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */


    BASICALIGNPAIR* ap;
    char *seq1pos,*seq2pos,*seq1sidepos,*seq2sidepos;
    int bestscore,bestscorerow;

#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

	/* set up to find bestscore in last column */
	bestscore=0;
	bestscorerow=0;

    /* compute SDD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top, pattern - global */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = col*indel;
        Sp->direction = RIGHT;
    }

    /* Down the left side, sequence - local */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = 0;
        Sp->direction = START;
    }

    /* body of matrix */
    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=len2;col++)
        {
            /* E */
            eb = Scm1p->score+indel;

            /* F */
            fb = Srm1p->score+indel; 

            /* S */
            sa = Srcm1p->score+
                 submatrix[26*(seq1[row-1]-'A')+(seq2[col-1]-'A')];
						
            switch(max3switch(sa,eb,fb))
            {
                case 1:
                    Sp->score = sa;
                    Sp->direction = DIAGONAL;
                    break;
                case 2:
                    Sp->score = eb;
                    Sp->direction = RIGHT;
                    break;
                case 3:
                    Sp->score = fb;
                    Sp->direction = DOWN;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }

		/* store bestscore in final column */
		if((Sp-1)->score>bestscore)
		{
			bestscore=(Sp-1)->score;
			bestscorerow=row;
		}
    }
    /* get the length of the alignment */
    Sp = S+(len2+1)*bestscorerow+len2; /* set at bestscore element in final column */

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
    /*     length ++; */
        switch(Sp->direction)
        {
            case START:
                 break;
            case DIAGONAL:
                 Sp = Sp -(len2+2);
				 length++;
                 break;
            case RIGHT:
                 Sp = Sp - 1;
				 length++;
                 break;
            case DOWN:
                 Sp = Sp - (len2+1);
				 length++;
                 break;
        }
    }

    /*  allocate memory     */
    ap = (BASICALIGNPAIR*) malloc(sizeof(COMPOSITIONALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->sequence1side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence1side==NULL) { free(ap); return NULL;}

    ap->sequence2side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence2side==NULL)
    {
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;

    Sp = S+(len2+1)*bestscorerow+len2; /* set at bestscore element in final column */

    ap->score  = Sp->score;
    ap->sequence1end = bestscorerow;
    ap->sequence2end = len2;
    ap->sequence2start = 1;

    seq1pos = seq1+bestscorerow-1;    /* position at bestscore character */
    seq2pos = seq2+len2-1;						/* position at end */
    seq1sidepos = ap->sequence1side+length;  /* position at termination */
    seq2sidepos = ap->sequence2side+length;
    *seq1sidepos=*seq2sidepos='\0'; /* terminate strings */
    seq1sidepos--;  /* move down to last char */
    seq2sidepos--;


    /* Sp has ben set to bestscore element above */
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *seq1sidepos = *seq1pos;
                 *seq2sidepos = *seq2pos;
                 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
				 Sp = Sp -(len2+2);
                 break;

            case RIGHT:
                 *seq1sidepos='-';
                 *seq2sidepos=*seq2pos;
				 seq1sidepos--;
                 seq2sidepos--;
                 seq2pos--;
                 Sp = Sp - 1;
                 break;
            case DOWN:
                 *seq1sidepos=*seq1pos;
                 *seq2sidepos='-';
                 seq1sidepos--;
                 seq1pos--;
                 seq2sidepos--;
                 Sp = Sp - (len2+1);
                 break;
        }
    }

		/* count characters in alignment to find start for sequence1 */
		length=0;
		for(row=0;row<ap->length;row++)
			if(ap->sequence1side[row]!='-') length++; 
    ap->sequence1start = ap->sequence1end-length+1;

    
		/* free alignment matrix */
		free(S);
    


    return ap;

}

/***********************************************************************/

BASICALIGNPAIR* GetBasicLocalAlignPair (char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix,SD ** Sret)
{
    int     col, row;
    int     eb;             /* equals e+beta */
    int     fb;             /* equals f+beta */
    int     sa;             /* equals s+alpha+beta */
    SD* S=NULL;
    int  Slen = 0;

    SD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */

    BASICALIGNPAIR* ap;
    int length,maxrow,maxcol,minrow,mincol,maxscore;
    char *seq1pos,*seq2pos,*seq1sidepos,*seq2sidepos;

#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* compute alignment Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* init top row to zero and start for local alignment */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = 0;
        Sp->direction = START;
    }

    /* init left column to zero and start for local alignment */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = 0;
        Sp->direction = START;
    }

    /* body of matrix */
    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=len2;col++)
        {
            /* E */
            eb = Scm1p->score+indel;

            /* F */
            fb = Srm1p->score+indel;

            /* S */
            sa = Srcm1p->score+
                 submatrix[26*(seq1[row-1]-'A')+(seq2[col-1]-'A')];
            switch(max3switch(sa,eb,fb))
            {
                case 1:
                    Sp->score = sa;
                    Sp->direction = DIAGONAL;
                    break;
                case 2:
                    Sp->score = eb;
                    Sp->direction = RIGHT;
                    break;
                case 3:
                    Sp->score = fb;
                    Sp->direction = DOWN;
                    break;
            }

            /* since it's local, never go below zero */
            if(Sp->score<0)
            {
                Sp->score = 0;
                Sp->direction = START;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }
    }

    /* find the highest scoring position for local trace back */
    maxscore = maxrow = maxcol = 0;
    Sp = S;
    for(row=0;row<=len1;row++)
    {
        for(col=0;col<=len2;col++)
        {
            if(Sp->score>=maxscore) /* choose longer if same score */
            {
                maxscore = Sp->score;
                maxrow = row;
                maxcol = col;
            }
            Sp++;
        }
    }
    Sp = S+maxrow*(len2+1)+maxcol; /* set at maximum row and column */

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
        length ++;
        switch(Sp->direction)
        {
            case START:
                 break;
            case DIAGONAL:
                 Sp = Sp -(len2+2);
                 break;
            case RIGHT:
                 Sp = Sp - 1;
                 break;
            case DOWN:
                 Sp = Sp - (len2+1);
                 break;
        }
    }

    /*  allocate memory     */
    ap = (BASICALIGNPAIR*) malloc(sizeof(BASICALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->sequence1side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence1side==NULL) { free(ap); return NULL;}

    ap->sequence2side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence2side==NULL)
    {
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;
    ap->score  = maxscore;
    ap->sequence1end = maxrow;
    ap->sequence2end = maxcol;
    seq1pos = seq1+(maxrow-1);    /* position at maximum score position */
    seq2pos = seq2+(maxcol-1);
    seq1sidepos = ap->sequence1side+length;  /* position at termination */
    seq2sidepos = ap->sequence2side+length;
    *seq1sidepos=*seq2sidepos='\0'; /* terminate strings */
    seq1sidepos--;  /* move down to last char */
    seq2sidepos--;

    Sp = S+maxrow*(len2+1)+maxcol; /* set at maximum row and column */
    minrow = maxrow; /* will updated in trace back */
    mincol = maxcol;
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *seq1sidepos = *seq1pos;
                 *seq2sidepos = *seq2pos;
                 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
                 Sp = Sp -(len2+2);
                 minrow--;mincol--;
                 break;
            case RIGHT:
                 *seq1sidepos='-';
                 *seq2sidepos=*seq2pos;
                 seq1sidepos--;
                 seq2sidepos--;
                 seq2pos--;
                 Sp = Sp - 1;
                 mincol--;
                 break;
            case DOWN:
                 *seq1sidepos=*seq1pos;
                 *seq2sidepos='-';
                 seq1sidepos--;
                 seq1pos--;
                 seq2sidepos--;
                 Sp = Sp - (len2+1);
                 minrow--;
                 break;
        }
    }
    ap->sequence1start = minrow+1;
    ap->sequence2start = mincol+1;

    /* free alignment matrix */
 	if (Sret)
		*Sret=S;
	else
		free(S);

    return ap;
}

/***********************************************************************/
BASICALIGNPAIR* GetBasicGlobalAlignPair(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix,SD ** Sret)
{
    int     col, row,length;
    int     eb;             /* equals e+indel */
    int     fb;             /* equals f+indel */
    int     sa;             /* equals s+indel */
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */


    BASICALIGNPAIR* ap;
    char *seq1pos,*seq2pos,*seq1sidepos,*seq2sidepos;


#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

	
    /* compute SDD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = col*indel;
        Sp->direction = RIGHT;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = row*indel;
        Sp->direction = DOWN;
    }

    /* body of matrix */
    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=len2;col++)
        {
            /* E */
            eb = Scm1p->score+indel;

            /* F */
            fb = Srm1p->score+indel;

            /* S */
            sa = Srcm1p->score+
                 submatrix[26*(seq1[row-1]-'A')+(seq2[col-1]-'A')];
            switch(max3switch(sa,eb,fb))
            {
                case 1:
                    Sp->score = sa;
                    Sp->direction = DIAGONAL;
                    break;
                case 2:
                    Sp->score = eb;
                    Sp->direction = RIGHT;
                    break;
                case 3:
                    Sp->score = fb;
                    Sp->direction = DOWN;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }
    }

    /* get the length of the alignment */
    Sp = S+(len2+1)*len1+len2;      /* set at last element */

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
        length ++;
        switch(Sp->direction)
        {
            case START:
                 break;
            case DIAGONAL:
                 Sp = Sp -(len2+2);
                 break;
            case RIGHT:
                 Sp = Sp - 1;
                 break;
            case DOWN:
                 Sp = Sp - (len2+1);
                 break;
        }
    }

    /*  allocate memory     */
    ap = (BASICALIGNPAIR*) malloc(sizeof(BASICALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->sequence1side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence1side==NULL) { free(ap); return NULL;}

    ap->sequence2side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence2side==NULL)
    {
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;

    Sp = S+(len2+1)*len1+len2;      /* set at last element */
    ap->score  = Sp->score;
    ap->sequence1end = len1;
    ap->sequence2end = len2;
    ap->sequence1start = 1;
    ap->sequence2start = 1;

    seq1pos = seq1+len1-1;    /* position last character */
    seq2pos = seq2+len2-1;
    seq1sidepos = ap->sequence1side+length;  /* position at termination */
    seq2sidepos = ap->sequence2side+length;
    *seq1sidepos=*seq2sidepos='\0'; /* terminate strings */
    seq1sidepos--;  /* move down to last char */
    seq2sidepos--;

    /* Sp has ben set to last element above */
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *seq1sidepos = *seq1pos;
                 *seq2sidepos = *seq2pos;
                 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
                 Sp = Sp -(len2+2);
                 break;
            case RIGHT:
                 *seq1sidepos='-';
                 *seq2sidepos=*seq2pos;
                 seq1sidepos--;
                 seq2sidepos--;
                 seq2pos--;
                 Sp = Sp - 1;
                 break;
            case DOWN:
                 *seq1sidepos=*seq1pos;
                 *seq2sidepos='-';
                 seq1sidepos--;
                 seq1pos--;
                 seq2sidepos--;
                 Sp = Sp - (len2+1);
                 break;
        }
    }

    /* free alignment matrix */
 	if (Sret)
		*Sret=S;
	else
		free(S);
    


    return ap;
}

/***********************************************************************/
WDPALIGNPAIR* GetWDPLocalAlignPair( char* sequence, char* pattern,
                                    int sequencelength, int patternlength,
                                    int alpha, int indel, int* submatrix)
{
    int     col, row;
    int     ei;             /* equals e+indel */
    int     fi;             /* equals f+indel */
    int     sa;             /* equals s+match */
    int     dirhold;
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp,*Srm1p,*Scm1p,*Srcm1p;
    WDPALIGNPAIR* ap;
    int length,maxrow,maxcol,minrow,mincol,maxscore;
    char *sequencepos,*patternpos,*sequencesidepos,*patternsidepos;


#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(sequencelength+1)*(__int64)(patternlength+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (sequencelength+1)*(patternlength+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* compute SDD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* init top row to zero and start for local alignment */
    for(col=1; col<=patternlength; col++)
    {
        Sp++;
        Sp->score = 0;
        Sp->direction = START;
    }

    /* init left column to zero and start for local alignment */
    Sp=S;
    for(row=1; row<=sequencelength; row++)
    {
        Sp+=(patternlength+1);
        Sp->score = 0;
        Sp->direction = START;
    }

    /* body of matrix */
    for(row=1;row<=sequencelength;row++)
    {
        /* first pass of wrap around */
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(patternlength+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-patternlength;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=patternlength;col++)
        {
            /* compute score from the left */
            ei = Scm1p->score+indel;

            /* compute score from the top */
            fi = Srm1p->score+indel;

            /* S */
            if(col==1&&(Sp-2)->score>Srcm1p->score)
            {
                sa = (Sp-2)->score+
                     submatrix[26*(sequence[row-1]-'A')+(pattern[col-1]-'A')];
                Sp->direction = WRAPDIAGONAL;
            }
            else
            {
                sa = Srcm1p->score+
                     submatrix[26*(sequence[row-1]-'A')+(pattern[col-1]-'A')];
                Sp->direction = DIAGONAL;
            }
            switch(max4switch(0,sa,ei,fi))
            {
                case 1:
                    Sp->score = 0;
                    Sp->direction = START;
                    break;
                case 2:
                    Sp->score = sa; /* direction was set above */
                    break;
                case 3:
                    Sp->score = ei;
                    Sp->direction = RIGHT;
                    break;
                case 4:
                    Sp->score = fi;
                    Sp->direction = DOWN;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }

        /* second pass of wrap around */
        Sp      =   S+1+row*(patternlength+1);
        Scm1p   =   Sp-1;

        for(col=1;col<=patternlength;col++)
        {
            /* compute score from the left */
            if(col==1)
            {
                ei = (Scm1p+patternlength)->score+indel;
                dirhold = WRAPRIGHT;
            }
            else
            {
                ei = Scm1p->score+indel;
                dirhold = RIGHT;
            }

            if(ei>Sp->score)
            {
                Sp->score = ei;
                Sp->direction = dirhold;
            }
            else break;

            /* update pointers */
            Sp++;
            Scm1p++;
        }
    }

    /* find the highest scoring position for local trace back */
    maxscore = maxrow = maxcol = 0;
    Sp = S;


    for(row=0;row<=sequencelength;row++)
    {
        for(col=0;col<=patternlength;col++)
        {
            if(Sp->score>=maxscore) /* choose longer if same score */
            {
                maxscore = Sp->score;
                maxrow = row;
                maxcol = col;
            }
            Sp++;
        }
    }
    Sp = S+maxrow*(patternlength+1)+maxcol; /* set at maximum row and column */

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
        length += 1;
        switch(Sp->direction)
        {
            case DIAGONAL:
                 Sp = Sp -(patternlength+2);
                 break;
            case WRAPDIAGONAL:
                 Sp = Sp-2;
                 break;
            case RIGHT:
                 Sp = Sp-1;
                 break;
            case WRAPRIGHT:
                 Sp = Sp + patternlength - 1;
                 break;
            case DOWN:
                 Sp = Sp - (patternlength+1);
                 break;
        }
    }

    /*  allocate memory     */
    ap = (WDPALIGNPAIR*) malloc(sizeof(WDPALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->sequenceside = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequenceside==NULL) { free(ap); return NULL;}

    ap->patternside = (char*) malloc((length+1)*sizeof(char));
    if(ap->patternside==NULL)
    {
        free(ap->sequenceside);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;
    ap->score  = maxscore;
    ap->sequenceend = maxrow;
    ap->patternend = maxcol;
    ap->patternlength = patternlength;
    sequencepos = sequence+(maxrow-1);    /* position at maximum score position */
    patternpos = pattern+(maxcol-1);
    sequencesidepos = ap->sequenceside+length;  /* position at termination */
    patternsidepos = ap->patternside+length;
    *sequencesidepos=*patternsidepos='\0'; /* terminate strings */
    sequencesidepos--;  /* move down to last char */
    patternsidepos--;

    Sp=S+maxrow*(patternlength+1)+maxcol; /* set at maximum row and column */
    minrow = maxrow; /* will updated in trace back */
    mincol = maxcol;
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *sequencesidepos = *sequencepos;
                 *patternsidepos = *patternpos;
                 sequencepos--;patternpos--;
                 sequencesidepos--;patternsidepos--;
                 Sp = Sp -(patternlength+2);
                 minrow--;mincol--;
                 break;
            case WRAPDIAGONAL:
                 *sequencesidepos = *sequencepos;
                 *patternsidepos = *patternpos;
                 sequencepos--;patternpos+=(patternlength-1);
                 sequencesidepos--;patternsidepos--;
                 Sp = Sp-2;
                 minrow--;mincol+=(patternlength-1);
                 break;
            case RIGHT:
                 *sequencesidepos='-';
                 *patternsidepos=*patternpos;
                 sequencesidepos--;
                 patternsidepos--;
                 patternpos--;
                 Sp = Sp - 1;
                 mincol--;
                 break;
            case WRAPRIGHT:
                 *sequencesidepos='-';
                 *patternsidepos=*patternpos;
                 sequencesidepos--;
                 patternsidepos--;
                 patternpos--;
                 if(patternpos<pattern)
                 {
                     patternpos = pattern+(patternlength-1);
                 }
                 Sp = Sp + patternlength - 1;
                 mincol--;
                 if(mincol<1) mincol += patternlength;
                 break;
            case DOWN:
                 *sequencesidepos=*sequencepos;
                 *patternsidepos='-';
                 sequencesidepos--;
                 sequencepos--;
                 patternsidepos--;
                 Sp = Sp - (patternlength+1);
                 minrow--;
                 break;
        }
    }
    ap->sequencestart = minrow+1;
    ap->patternstart = mincol+1;


    /* free alignment matrix */
    free(S);
    
    return ap;
}

/***********************************************************************/
WDPALIGNPAIR* GetWDPGlobalAlignPair( char* sequence, char* pattern,
                                    int sequencelength, int patternlength,
                                    int alpha, int indel, int* submatrix)
{
    int     col, row;
    int     ei;             /* equals e+indel */
    int     fi;             /* equals f+indel */
    int     sa;             /* equals s+match */
    int     dirhold;
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp,*Srm1p,*Scm1p,*Srcm1p;
    WDPALIGNPAIR* ap;
    int length,maxcol,mincol,maxscore;
    char *sequencepos,*patternpos,*sequencesidepos,*patternsidepos;

#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(sequencelength+1)*(__int64)(patternlength+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif


    /* allocate the alignment matrix */
    Slen = (sequencelength+1)*(patternlength+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* compute SD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top */
    for(col=1; col<=patternlength; col++)
    {
        Sp++;
        Sp->score = 0;
        Sp->direction = START;
    }

    /* down the left side */
    Sp=S;
    for(row=1; row<=sequencelength; row++)
    {
        Sp+=(patternlength+1);
        Sp->score = row*indel;
        Sp->direction = DOWN;
    }

    /* body of matrix */
    for(row=1;row<=sequencelength;row++)
    {
        /* first pass of wrap around */
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(patternlength+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-patternlength;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=patternlength;col++)
        {
            /* compute score from left */
            ei = Scm1p->score+indel;

            /* compute score from the top */
            fi = Srm1p->score+indel;

            /* compute score from the diagonal */
            if(col==1&&(Sp-2)->score>Srcm1p->score)
            {
                sa = (Sp-2)->score+
                     submatrix[26*(sequence[row-1]-'A')+(pattern[col-1]-'A')];
                dirhold = WRAPDIAGONAL;
            }
            else
            {
                sa = Srcm1p->score+
                     submatrix[26*(sequence[row-1]-'A')+(pattern[col-1]-'A')];
                dirhold = DIAGONAL;
            }

            /* decide the optimal approach */
            switch(max3switch(sa,ei,fi))
            {
                case 1:
                    Sp->score = sa;
                    Sp->direction = dirhold;
                    break;
                case 2:
                    Sp->score = ei;
                    Sp->direction = RIGHT;
                    break;
                case 3:
                    Sp->score = fi;
                    Sp->direction = DOWN;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }


        /* second pass of wrap around */
        Sp      =   S+1+row*(patternlength+1);
        Scm1p   =   Sp-1;

        for(col=1;col<=patternlength;col++)
        {
            /* compute score from the left */
            if(col==1)
            {
                ei = (Scm1p+patternlength)->score+indel;
                dirhold = WRAPRIGHT;
            }
            else
            {
                ei = Scm1p->score+indel;
                dirhold = RIGHT;
            }

            if(ei>Sp->score)
            {
                Sp->score = ei;
                Sp->direction = dirhold;
            }
            else break;

            /* update pointers */
            Sp++;
            Scm1p++;
        }
    }

    /* find highest score of last row for global WDP trace back */
    Sp = S+sequencelength*(patternlength+1);
    maxscore = Sp->score;
    maxcol = 0;
    for(col=0;col<=patternlength;col++)
    {
        if(Sp->score>=maxscore) /* choose longer if same score */
        {
            maxscore = Sp->score;
            maxcol = col;
        }
        Sp++;
    }
    Sp = S+sequencelength*(patternlength+1)+maxcol; /* at maxcol,last row */

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
        length++;
        switch(Sp->direction)
        {
            case DIAGONAL:
                 Sp = Sp -(patternlength+2);
                 break;
            case WRAPDIAGONAL:
                 Sp = Sp-2;
                 break;
            case RIGHT:
                 Sp = Sp-1;
                 break;
            case WRAPRIGHT:
                 Sp = Sp + patternlength - 1;
                 break;
            case DOWN:
                 Sp = Sp - (patternlength+1);
                 break;
        }
    }

    /*  allocate memory     */
    ap = (WDPALIGNPAIR*) malloc(sizeof(WDPALIGNPAIR));
    if(ap==NULL)
    {
        return NULL;
    }

    ap->sequenceside = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequenceside==NULL)
    {
        free(ap);
        return NULL;
    }

    ap->patternside = (char*) malloc((length+1)*sizeof(char));
    if(ap->patternside==NULL)
    {
        free(ap->sequenceside);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;
    ap->score  = maxscore;
    ap->sequenceend = sequencelength;
    ap->patternend = maxcol;
    ap->patternlength = patternlength;
    sequencepos = sequence+(sequencelength-1);/* last row */
    patternpos = pattern+(maxcol-1); /* maximum column */
    sequencesidepos = ap->sequenceside+length;  /* position at termination */
    patternsidepos = ap->patternside+length;
    *sequencesidepos=*patternsidepos='\0'; /* terminate strings */
    sequencesidepos--;  /* move down to last char */
    patternsidepos--;

    Sp=S+sequencelength*(patternlength+1)+maxcol; /* last row max column */
    mincol = maxcol;
	while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *sequencesidepos = *sequencepos;
                 *patternsidepos = *patternpos;
                 sequencepos--;patternpos--;
                 sequencesidepos--;patternsidepos--;
                 Sp = Sp -(patternlength+2);
                 mincol--;
                 break;
            case WRAPDIAGONAL:
                 *sequencesidepos = *sequencepos;
                 *patternsidepos = *patternpos;
                 sequencepos--;patternpos+=(patternlength-1);
                 sequencesidepos--;patternsidepos--;
                 Sp = Sp-2;
                 mincol+=(patternlength-1);
                 break;
            case RIGHT:
                 *sequencesidepos='-';
                 *patternsidepos=*patternpos;
                 sequencesidepos--;
                 patternsidepos--;
                 patternpos--;
                 Sp = Sp - 1;
                 mincol-=1;
                 break;
            case WRAPRIGHT:
                 *sequencesidepos='-';
                 *patternsidepos=*patternpos;
                 sequencesidepos--;
                 patternsidepos--;
                 patternpos--;
                 if(patternpos<pattern)
                 {
                     patternpos = pattern+(patternlength-1);
                 }
                 Sp = Sp + patternlength - 1;
                 mincol-=1;
                 if(mincol<1) mincol += patternlength;
                 break;
            case DOWN:
                 *sequencesidepos=*sequencepos;
                 *patternsidepos='-';
                 sequencesidepos--;
                 sequencepos--;
                 patternsidepos--;
                 Sp = Sp - (patternlength+1);
                 break;
        }
    }
    ap->sequencestart = 1;
    ap->patternstart = mincol+1;

    /* free alignment matrix */
    free(S);

    return ap;
}

/***************************************************************/
int* GetCompositionSubstringMatchLengths(char* seq1, char* seq2, int len1, int len2,
																				 int limit, int ALPHABETSIZE)
																				 /*
																				 Returns a two dimensional array (created and used as a one dimensional array) 
																				 M of size (len1+1)x(len2+1) which holds the length of the shortest 
																				 composition matching substrings (with zero < length <= limit) 
																				 ending at seq1[i] and seq2[j] or zero if no such matching substrings exist.
																				 Assumes the definition of a global ALPHABETSIZE.
																				 */
{
	
	int Mlen=0;
	COMPOSITIONVECTOR *Diffvectors, *dfi, *dfic;
  int *M, *Radsort1, *Radsort2, *Counts;
  int g,k,size,m,n,d,xstart,ystart,length;
	int minval,maxval,val,*rsi,*trs;
	int current,previous,match,i,j,substring_length;
  char *X, *Y;
  int *Index;	
	
#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif
	
  /* allocate the substring composition match array */
	Mlen = (len1+1)*(len2+1);
	M = (int*) malloc(sizeof(int)*Mlen);
	if(M==NULL)
	{
		return NULL; /* in case memory allocation fails */
	}
	
  size=min(len1,len2);
  Diffvectors=(COMPOSITIONVECTOR *)malloc(sizeof(COMPOSITIONVECTOR)*(size+1));
	if(Diffvectors==NULL)
	{
		return NULL; /* in case memory allocation fails */
	}
  Radsort1=(int *)malloc(sizeof(int)*(size+1));
	if(Radsort1==NULL)
	{
		return NULL; /* in case memory allocation fails */
	}
  Radsort2=(int *)malloc(sizeof(int)*(size+1));
	if(Radsort2==NULL)
	{
		return NULL; /* in case memory allocation fails */
	}
  Counts=(int *)malloc(sizeof(int)*((2*size)+1));
	if(Counts==NULL)
	{
		return NULL; /* in case memory allocation fails */
	}
	Index=(int *)calloc(256,sizeof(int));
  if(Index==NULL){return NULL;}
  if(ALPHABETSIZE == 5)
  {
  Index['A']=0;
  Index['C']=1;
  Index['G']=2;
  Index['T']=3;
  Index['X']=4;
  }
  else
  {
  Index['A']=0;
  Index['B']=1;
  Index['C']=2;
  Index['D']=3;
  Index['E']=4;
  Index['F']=5;
  Index['G']=6;
  Index['H']=7;
  Index['I']=8;
  Index['J']=9;
  Index['K']=10;
  Index['L']=11;
  Index['M']=12;
  Index['N']=13;
  Index['O']=14;
  Index['P']=15;
  Index['Q']=16;
  Index['R']=17;
  Index['S']=18;
  Index['T']=19;
  Index['U']=20;
  Index['V']=21;
  Index['W']=22;
  Index['X']=23;
  Index['Y']=24;
  }
	/******/
	m=len1;
	n=len2;
	/* compute once for each diagonal of M */
	for(d=-m+1;d<=n-1;d++)
	{
		/* calculate starting locations in X and Y (strings start at 
		index zero) and the length of diagonal d */
		
		if(d<0)
		{
			xstart=0;
			ystart=-d;
			length=min(n,m+d);
		}
		else
		{
			xstart=d;
			ystart=0;
			length=min(n-d,m);
		}
    
		/* compute substring differences for diagonal d and store in Diffvectors */
		
		
		dfi=dfic=Diffvectors;
		Y=seq1+ystart;
		X=seq2+xstart;
		
		/* zeroth vector is all zeros */
		for(k=0;k<ALPHABETSIZE;k++)
			dfi->C[k]=0;
		
		/* step through the substrings computing the difference X-Y at each length */
		for(g=0;g<length;g++)
		{
			dfic++;
			/* current composition is equal to the previous composition */
			for(k=0;k<ALPHABETSIZE;k++)
			{
				dfic->C[k]=dfi->C[k];
			}
			/* and then adjusted for +X[g] and -Y[g] */
			dfic->C[Index[X[g]]]++;
			dfic->C[Index[Y[g]]]--;
			dfi++;
		}
		
	
		/* radixsort the composition vectors in Diffvectors */
		/* sort into Radsort1 using Radsort2 for swapping */
		/* length+1 is number of vectors to sort */
		/*********************/
		
		
		/* initialize Radsort1 to the order of the vectors in Diffvectors */
		for(g=0;g<=length;g++)
			Radsort1[g]=g;
		
		/* sort on each component starting with last */  
		for(k=ALPHABETSIZE-1;k>=0;k--)
		{
			
			/* get counts for counting sort */
			/* first zero counts */

			for(g=-length;g<=length;g++)
				Counts[g+length]=0;        /* +length so 0<= index <= 2*length */
			
			/* get counts in each component of vector */
			
			dfi=Diffvectors;
			minval=0;
			maxval=0;
			for(g=0;g<=length;g++)
			{
				val=dfi->C[k];
				Counts[val+length]++;
				if(val>maxval) maxval=val;
				if(val<minval) minval=val;
				dfi++;
			}
			
			/* make counts cumulative, use minval and maxval */
			
			for(g=minval+1;g<=maxval;g++)
				Counts[g+length]+=Counts[g+length-1];
			
				/* go through Diffvectors (last to first) and save new positions
			in Radsort2 */
			
			rsi=Radsort1+length;
			for(g=0;g<=length;g++)
			{
				val=Diffvectors[(*rsi)].C[k];
				Radsort2[Counts[val+length]-1]=(*rsi); /* subtract 1 because Radsort is zero based */
				Counts[val+length]--;
				rsi--;
			}
			
			/* swap Radsort2 into Radsort1 */
			trs=Radsort1;
			Radsort1=Radsort2;
			Radsort2=trs;
		}
		
		
		/**********************************/
		
		/* fill M */
		/* first fill zero row and zero column */
		for(g=0;g<=len2;g++)
			*(M+g)=0;
		for(g=1;g<=len1;g++)
			*(M+g*(len2+1))=0;
		
		/* find matching compositions in the sorted Diffvectors */
		/* length+1 is the number of vectors in the sorted list */
		
		rsi=Radsort1;
		previous=current=(*rsi);  /* previous not used the first time through the 
		loop */
		
		/* get location in M to save matching substring lengths */
		
		for(g=0;g<=length;g++)
		{
			/* first vector always has no previous one to look at so no match */
			if(g==0)
			{
				substring_length=0;
			}
			else
			{
				/* test if current and previous vectors match */
				/* but only if the matching substring would be no longer than limit */
				if((substring_length=current-previous)>limit) 
					substring_length=0; /* store a default length of zero */
				else
				{
					match=1;
					for(k=0;k<ALPHABETSIZE;k++)
					{
						if(Diffvectors[current].C[k]!=Diffvectors[previous].C[k])
						{
							match=0;
							break;
						}
					}
					
					/* if they match, compute matching substring length */
					if(match)
						substring_length=current-previous;
					else /* store a default length of zero */
						substring_length=0;
					/* previous=current; */ /* uncomment for longest substring match */
				}
				/* calculate substring ending indices */
			}
			if(d>=0)
			{
				j=current+d;
				i=current;
			}
			else
			{
				j=current;
				i=current-d;
			}
			
			/* store substring length at M[i][j] */
			*(M+i*(len2+1)+j)=substring_length;
			
			/* update pointers */
			previous=current;  /* uncomment for shortest substring */
			rsi++;
			current=(*rsi);
			/**********************************/
			/******/
		}
	}
	free(Index);
	free(Diffvectors);
	free(Radsort1);
	free(Radsort2);
	free(Counts);
	return(M);
}

/***********************************************************************/
COMPOSITIONALIGNPAIR* GetCompositionGlobalAlignPair(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix, int* M)
/*
Produces a composition alignment for two sequences.  The aligned pair has an 
extra string (matchboundaries) which marks the boundaries of the shortest composition matching substrings.
M[i][j] is a matrix (len1+1)x(len2+1) which holds the length of the shortest 
composition matching substrings (with zero< length < limit) ending at seq1[i] and seq2[j] 
or zero if no such matching substrings exist.
The limit is the longest of these short substring that can be included
in the alignment. 
M is computed with the function GetCompositionSubstringMatchLengths which includes 
the limit as a parameter.  This function must be passed M as a parameter.
In the alignment, a substring composition match is scored as the length of the 
substring times a matchscore which is assumed to be a constant.  Here it is set 
to the score of A vs A.
*/
{
    int     col, row,length;
    int     eb;             /* equals e+indel */
    int     fb;             /* equals f+indel */
    int     sa;             /* equals s+indel */
    int     stemp;          /* Gary Benson 10/13/03 */
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */


    COMPOSITIONALIGNPAIR* ap;
    char *seq1pos,*seq2pos,*seq1sidepos,*seq2sidepos;

    int substringlen, matchscore, satype, i;
		char *mbpos;
		int *Mp;

		/* matchscore for substring composition match is the score of A vs A */
		matchscore=submatrix[26*('A'-'A')+('A'-'A')];

#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

	
    /* compute SDD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = col*indel;
        Sp->direction = RIGHT;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = row*indel;
        Sp->direction = DOWN;
    }

    /* body of matrix */
    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;
				Mp			=		M+1+row*(len2+1);

        for(col=1;col<=len2;col++)
        {
            /* E */
            eb = Scm1p->score+indel;

            /* F */
            fb = Srm1p->score+indel;

            /* S */
						
						/* changed S computation, Gary Benson, 10/13/03 to allow comparison of mismatch
						and composition match greater than length 1 */
						
						/* calculate match or mismatch */
						
						sa = Srcm1p->score+submatrix[26*(seq1[row-1]-'A')+(seq2[col-1]-'A')];
						satype=DIAGONAL;
						
						/* calculate composition match and test against mismatch */
						if((substringlen=*Mp)>1)
						{
							stemp=(Sp-substringlen*(len2+2))->score+substringlen*matchscore;
                            if(stemp>sa) {
                                sa=stemp;
								satype=LONGDIAGONAL;
                            }
						}
						
		    switch(max3switch(sa,eb,fb))			
           {
                case 1:
                    Sp->score = sa;
                    Sp->direction = satype; /* DIAGONAL	or LONGDIAGONAL */
                    break;
                case 2:
                    Sp->score = eb;
                    Sp->direction = RIGHT;
                    break;
                case 3:
                    Sp->score = fb;
                    Sp->direction = DOWN;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
						Mp++;
        }
    }

    /* get the length of the alignment */
    Sp = S+(len2+1)*len1+len2; /* set at last element */
		Mp = M+(len2+1)*len1+len2;

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
    /*     length ++; */
        switch(Sp->direction)
        {
            case START:
                 break;
            case DIAGONAL:
                 Sp = Sp -(len2+2);
								 Mp = Mp -(len2+2);
								 length++;
                 break;
						case LONGDIAGONAL:
							   substringlen=*Mp;
								 Sp=Sp-substringlen*(len2+2);
								 Mp=Mp-substringlen*(len2+2);
								 length+=substringlen;
								 break;
            case RIGHT:
                 Sp = Sp - 1;
								 Mp = Mp - 1;
								 length++;
                 break;
            case DOWN:
                 Sp = Sp - (len2+1);
								 Mp = Mp - (len2+1);
								 length++;
                 break;
        }
    }

    /*  allocate memory     */
    ap = (COMPOSITIONALIGNPAIR*) malloc(sizeof(COMPOSITIONALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->sequence1side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence1side==NULL) { free(ap); return NULL;}

    ap->sequence2side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence2side==NULL)
    {
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }
    ap->matchboundaries = (char*) malloc((length+1)*sizeof(char));
    if(ap->matchboundaries==NULL)
    {
			  free(ap->sequence2side);
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;

    Sp = S+(len2+1)*len1+len2;      /* set at last element */
    Mp = M+(len2+1)*len1+len2;      /* set at last element */

    ap->score  = Sp->score;
    ap->sequence1end = len1;
    ap->sequence2end = len2;
    ap->sequence1start = 1;
    ap->sequence2start = 1;

    seq1pos = seq1+len1-1;    /* position last character */
    seq2pos = seq2+len2-1;
    seq1sidepos = ap->sequence1side+length;  /* position at termination */
    seq2sidepos = ap->sequence2side+length;
    *seq1sidepos=*seq2sidepos='\0'; /* terminate strings */
    seq1sidepos--;  /* move down to last char */
    seq2sidepos--;

    mbpos=ap->matchboundaries+length; /* matchboundaries string */
		*mbpos='\0';
		mbpos--;


    /* Sp has ben set to last element above */
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *seq1sidepos = *seq1pos;
                 *seq2sidepos = *seq2pos;
								 if(*seq1sidepos==*seq2sidepos) *mbpos='|';
								 else *mbpos=' ';
                 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
								 mbpos--;
					       Sp = Sp -(len2+2);
                 Mp = Mp -(len2+2);
                 break;
            case LONGDIAGONAL:
							   substringlen=*Mp;
								 for(i=1;i<=substringlen;i++)
								 {
									 *seq1sidepos = *seq1pos;
									 *seq2sidepos = *seq2pos;
									 if(i==1) *mbpos='>';
									 else if (i==substringlen) *mbpos='<';
									 else *mbpos='-';
									 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
									 mbpos--;
									 Sp = Sp -(len2+2);
                   Mp = Mp -(len2+2);
                 }
								 break;
            case RIGHT:
                 *seq1sidepos='-';
                 *seq2sidepos=*seq2pos;
                 *mbpos=' ';
								 seq1sidepos--;
                 seq2sidepos--;
                 seq2pos--;
								 mbpos--;
                 Sp = Sp - 1;
                 Mp = Mp - 1;
                 break;
            case DOWN:
                 *seq1sidepos=*seq1pos;
                 *seq2sidepos='-';
								 *mbpos=' ';
                 seq1sidepos--;
                 seq1pos--;
                 seq2sidepos--;
								 mbpos--;
                 Sp = Sp - (len2+1);
                 Mp = Mp - (len2+1);
                 break;
        }
    }

    /* free alignment matrix */
    free(S);
    

    return ap;
}

/***********************************************************************/
/***********************************************************************/
COMPOSITIONALIGNPAIR* GetCompositionLocalAlignPair(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix, int* M)
/*
Produces a local composition alignment for two sequences.  The aligned pair has an 
extra string (matchboundaries) which marks the boundaries of the shortest composition matching substrings.
M[i][j] is a matrix (len1+1)x(len2+1) which holds the length of the shortest 
composition matching substrings (with zero< length < limit) ending at seq1[i] and seq2[j] 
or zero if no such matching substrings exist.
The limit is the longest of these short substring that can be included
in the alignment. 
M is computed with the function GetCompositionSubstringMatchLengths which includes 
the limit as a parameter.  This function must be passed M as a parameter.
In the alignment, a substring composition match is scored as the length of the 
substring times a matchscore which is assumed to be a constant.  Here it is set 
to the score of A vs A.
*/
{
    int     col, row,length;
    int     eb;             /* equals e+indel */
    int     fb;             /* equals f+indel */
    int     sa;             /* equals s+indel */
    int			stemp;          /* Gary Benson 10/13/03 */
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */


    COMPOSITIONALIGNPAIR* ap;
    char *seq1pos,*seq2pos,*seq1sidepos,*seq2sidepos;

    int substringlen, matchscore, satype, i;
		char *mbpos;
		int *Mp;
		int bestscore,bestscorerow,bestscorecol;


		/* matchscore for substring composition match is the score of A vs A */
		matchscore=submatrix[26*('A'-'A')+('A'-'A')];

#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

	
    /* compute SDD Matrix */
		/* set up for finding best score */
		bestscore=0;
		bestscorerow=0;
		bestscorecol=0;

    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = 0;
        Sp->direction = START;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = 0;
        Sp->direction = START;
    }

    /* body of matrix */
    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;
				Mp			=		M+1+row*(len2+1);

        for(col=1;col<=len2;col++)
        {
            /* E */
            eb = Scm1p->score+indel;

            /* F */
            fb = Srm1p->score+indel;

            /* S */
						
						/* changed S computation, Gary Benson, 10/13/03 to allow comparison of mismatch
						and composition match greater than length 1 */
						
						/* calculate match or mismatch */
						
						sa = Srcm1p->score+submatrix[26*(seq1[row-1]-'A')+(seq2[col-1]-'A')];
						satype=DIAGONAL;
						
						/* calculate composition match and test against mismatch */
						if((substringlen=*Mp)>1)
						{
							stemp=(Sp-substringlen*(len2+2))->score+substringlen*matchscore;
                            if(stemp>sa) {
                                sa=stemp;
								satype=LONGDIAGONAL;
                            }
						}
	
            switch(max4switch(0,sa,eb,fb))
            {
								case 1:
										Sp->score=0;
										Sp->direction=START;
										break;
                case 2:
                    Sp->score = sa;
                    Sp->direction = satype; /* DIAGONAL	or LONGDIAGONAL */
                    break;
                case 3:
                    Sp->score = eb;
                    Sp->direction = RIGHT;
                    break;
                case 4:
                    Sp->score = fb;
                    Sp->direction = DOWN;
                    break;
            }
						if(Sp->score>bestscore)
						{
							bestscore=Sp->score;
							bestscorerow=row;
							bestscorecol=col;
						}

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
						Mp++;
        }
    }

    /* get the length of the alignment */
    Sp = S+bestscorerow*(len2+1)+bestscorecol; /* set at bestscore element */
		Mp = M+bestscorerow*(len2+1)+bestscorecol;

    /* find out how long the alignment is */
    length = 0;
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case START:
                 break;
            case DIAGONAL:
                 Sp = Sp -(len2+2);
								 Mp = Mp -(len2+2);
								 length++;
                 break;
						case LONGDIAGONAL:
							   substringlen=*Mp;
								 Sp=Sp-substringlen*(len2+2);
								 Mp=Mp-substringlen*(len2+2);
								 length+=substringlen;
								 break;
            case RIGHT:
                 Sp = Sp - 1;
								 Mp = Mp - 1;
								 length++;
                 break;
            case DOWN:
                 Sp = Sp - (len2+1);
								 Mp = Mp - (len2+1);
								 length++;
                 break;
        }
    }

    /*  allocate memory     */
    ap = (COMPOSITIONALIGNPAIR*) malloc(sizeof(COMPOSITIONALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->sequence1side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence1side==NULL) { free(ap); return NULL;}

    ap->sequence2side = (char*) malloc((length+1)*sizeof(char));
    if(ap->sequence2side==NULL)
    {
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }
    ap->matchboundaries = (char*) malloc((length+1)*sizeof(char));
    if(ap->matchboundaries==NULL)
    {
			  free(ap->sequence2side);
        free(ap->sequence1side);
        free(ap);
        return NULL;
    }

    /***************************
    * trace back and fill in ap
    ****************************/
    ap->length = length;

    Sp = S+bestscorerow*(len2+1)+bestscorecol;      /* set at bestscore element */
    Mp = M+bestscorerow*(len2+1)+bestscorecol;      /* set at bestscore element */

    ap->score  = Sp->score;
    ap->sequence1end = bestscorerow; /* was len1 */
    ap->sequence2end = bestscorecol; /* was len2 */

    seq1pos = seq1+bestscorerow-1;    /* position last character */
    seq2pos = seq2+bestscorecol-1;
    seq1sidepos = ap->sequence1side+length;  /* position at termination */
    seq2sidepos = ap->sequence2side+length;
    *seq1sidepos=*seq2sidepos='\0'; /* terminate strings */
    seq1sidepos--;  /* move down to last char */
    seq2sidepos--;

    mbpos=ap->matchboundaries+length; /* matchboundaries string */
		*mbpos='\0';
		mbpos--;


    /* Sp has been set to bestscore element above */
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                 *seq1sidepos = *seq1pos;
                 *seq2sidepos = *seq2pos;
								 if(*seq1sidepos==*seq2sidepos) *mbpos='|';
								 else *mbpos=' ';
                 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
								 mbpos--;
					       Sp = Sp -(len2+2);
                 Mp = Mp -(len2+2);
                 break;
            case LONGDIAGONAL:
							   substringlen=*Mp;
								 for(i=1;i<=substringlen;i++)
								 {
									 *seq1sidepos = *seq1pos;
									 *seq2sidepos = *seq2pos;
									 if(i==1) *mbpos='>';
									 else if (i==substringlen) *mbpos='<';
									 else *mbpos='-';
									 seq1pos--;seq2pos--;seq1sidepos--;seq2sidepos--;
									 mbpos--;
									 Sp = Sp -(len2+2);
                   Mp = Mp -(len2+2);
                 }
								 break;
            case RIGHT:
                 *seq1sidepos='-';
                 *seq2sidepos=*seq2pos;
                 *mbpos=' ';
								 seq1sidepos--;
                 seq2sidepos--;
                 seq2pos--;
								 mbpos--;
                 Sp = Sp - 1;
                 Mp = Mp - 1;
                 break;
            case DOWN:
                 *seq1sidepos=*seq1pos;
                 *seq2sidepos='-';
								 *mbpos=' ';
                 seq1sidepos--;
                 seq1pos--;
                 seq2sidepos--;
								 mbpos--;
                 Sp = Sp - (len2+1);
                 Mp = Mp - (len2+1);
                 break;
        }
    }

		/* count characters in alignment to find start */
		length=0;
		for(row=0;row<ap->length;row++)
			if(ap->sequence1side[row]!='-') length++; 
    ap->sequence1start = ap->sequence1end-length+1;
		
		length=0;
		for(col=0;col<ap->length;col++)
			if(ap->sequence2side[col]!='-') length++; 
    ap->sequence2start = ap->sequence2end-length+1;
		
    /* free alignment matrix */
    free(S);
    

    return ap;
}





/***********************************************************************/

void PrintGappedAlignPair(FILE *fp, GAPPEDALIGNPAIR* gapap)
{
    int i,gpl;
    char *lowerp, *upperp;

    /* print a copy of the pattern */
    gpl = gapap->gappedpatternlength;
    lowerp = gapap->patternside;
    fprintf(fp,"\n\nDifference Format:\n\n    ");
    for(i=0;i<gpl&&*lowerp!='\0';i++,lowerp++)
    {
        fprintf(fp,"%c",*lowerp);
    }

    upperp = gapap->sequenceside;
    lowerp = gapap->patternside;
    while(*upperp!='\0')
    {
        fprintf(fp,"\n    ");
        /* print differeces or spaces */
        for(i=0;i<gpl;i++)
        {
            if(*upperp==*lowerp) fprintf(fp,".");
            else fprintf(fp,"%c",*upperp);
            upperp++;
            lowerp++;
            if(*upperp=='\0')break;
        }
    }

    /* print another copy of the pattern */
    lowerp = gapap->patternside;
    fprintf(fp,"\n    ");
    for(i=0;i<gpl&&*lowerp!='\0';i++,lowerp++)
    {
        fprintf(fp,"%c",*lowerp);
    }

    return;
}


/***********************************************************************/

void FreeBasicAlignPair(BASICALIGNPAIR* ap)
{
    free(ap->sequence1side);
    free(ap->sequence2side);
    free(ap);
    return;
}

/***********************************************************************/

void FreeWDPAlignPair(WDPALIGNPAIR* ap)
{
    free(ap->sequenceside);
    free(ap->patternside);
    free(ap);
    return;
}


/***********************************************************************/

void FreeGappedAlignPair(GAPPEDALIGNPAIR* ap)
{
    free(ap->sequenceside);
    free(ap->patternside);
    free(ap);
    return;
}


/***********************************************************************/

void FreeCyclicAlignPair(CYCLICALIGNPAIR* ap)
{
    free(ap->pattern1side);
    free(ap->pattern2side);
    free(ap);
    return;
}


/***********************************************************************/


void FreeCompositionAlignPair(COMPOSITIONALIGNPAIR* ap)
{
    free(ap->sequence1side);
    free(ap->sequence2side);
		free(ap->matchboundaries);
    free(ap);
    return;
}

/***********************************************************************/
int     GetCharsNotDash(char* string)
{
    static int count;
    static char* ptr;

    for(count=0,ptr=string;*ptr!='\0';ptr++)
    {
        if(*ptr!='-') count++;
    }
    return count;
}

/***********************************************************************/

char*   GetLeftRotatedString(char* original, int count)
{
    int length,i;
    char* new;
    char* source;
    char* destination;

    length = strlen(original);
    new = (char*) malloc((length+1)*sizeof(char));

    count = count%length; /* correct in case count is too big */

    /* copy the first count characters to end of new string */
    source=original;
    destination=new+(length-count);
    for(i=0;i<count;i++)
    {
        *destination = *source;
        source++;
        destination++;
    }
    *destination = '\0';

    /* copy the rest to the front of the new string */
    destination=new;
    for(i=length-count;i>0;i--)
    {
        *destination = *source;
        source++;
        destination++;
    }

    return new;
}


/***********************************************************************/
typedef struct
{
    int left;
    int right;
} LIMIT;

typedef struct
{
    int     done;
    int     score;
    LIMIT*  limits;
} SPINEOFLIMITS;

/* LOGDEPTH is used by ControStack and limits Reserve */

#define LOGDEPTH 14 /* ceil(log2(MAXPATTERNSIZE))+3 */

/* global variables */
SPINEOFLIMITS L[MAXPATSIZE];
char   doublepat[MAXPATSIZE*2];
LIMIT  outerleft[MAXPATSIZE];
LIMIT  outerright[MAXPATSIZE];

struct
{
    int next;
    LIMIT*  limits[LOGDEPTH];
} Reserve;

struct 
{
    int height;
    int left[LOGDEPTH];
    int right[LOGDEPTH];
} ControlStack;    /* used to get correct order for the alignments */


void CyclicAPToLimit(int start, int slength, int plength, CYCLICALIGNPAIR* pa)
{
    int row,col,i,length;
    char *ss, *ps;

    length=pa->length;
    ss=pa->pattern1side;
    ps=pa->pattern2side;
	  
    row=0;
	col=start;
	L[start].limits[row].right=col;
	
	for(i=0;i<length;i++) 
	{
        if((ss[i]!='-')&&(ps[i]!='-'))
        {
            row++;
            col++;
        }
        else if (ss[i]=='-') {col++;}
        else /* (ps[i]=='-') */ {row++;}
        L[start].limits[row].right=col;
	}
	
	row=slength;
	col=start+plength; 
	
	for(i=length-1;i>=0;i--) 
	{
        L[start].limits[row].left=col;
        if((ss[i]!='-')&&(ps[i]!='-'))
        {
            row--;
            col--;
        }
        else if (ss[i]=='-') {col--;}
        else /* (ps[i]=='-') */ {row--;}
	}
	
	L[start].limits[row].left=col;
}



CYCLICALIGNPAIR* GetBoundedDiffAlignment(char* pattern1, char* pattern2,
                        int pattern1len, int start,
                        int end,LIMIT* mins, LIMIT* maxs,SD* S)
{
    int min,max;
    int r,c,i,width;
    SD  *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */
    int sa,f,e;
    CYCLICALIGNPAIR *ap;
    char *pat1sidepos,*pat2sidepos,*pattern1pos,*pattern2pos;

    /****************************************
    * Do scores and directions for row zero
    *****************************************/
    r=0;
    Sp = S; /* place at first position */
    Sp->score = 0; S->direction = START;Sp++;
    max = MINIMUM(maxs[r].right,end);
    for(c=start+1,i=1;c<=max;c++,i++)
    {
        Sp->direction = RIGHT;
        Sp->score = i;
		Sp++;
    }
    if(max<MINIMUM(maxs[r+1].right,end)) /* need to pad with VLN at end */
    {
        for(i=MINIMUM(maxs[r+1].right,end);i>max;i--)
        {
             Sp->score=VLN;
             Sp++;
        }
    }

    /******************************************
    * Do scores and directions for start column
    *******************************************/
    Sp = S; /* reset to first element */
    width = end-start+1;
    for(r=1;r<=pattern1len;r++)
    {
        Sp+= width;
        Sp->score = r;
        Sp->direction =DOWN;
        min = MAXIMUM(mins[r+1].left,start); /* min of next row */
        if(min>start) break;  /* break if left is cut-off */
    }

    /**********************************************
    * Do scores and directions in center of matrix
    ***********************************************/
    for(r=1;r<=pattern1len;r++)
    {
        Sp  = S+(width*r);
        min = MAXIMUM(mins[r].left,start);
        max = MINIMUM(maxs[r].right,end);
        if(min>start) /* it's a left cut-off row */
        {
            Sp+= (min-start); /* position on first valid column */
            (Sp-1)->score=VLN;  /* pad with VLN if left cut-off */
			c=min;
        }
        else
		{
			Sp++; /* otherwise just move to column start+1 */
			c=min+1;
		}
        /* set the other pointers */
        Scm1p = Sp-1;
        Srm1p = Sp-width;
        Srcm1p= Scm1p-width;
        for(;c<=max;c++)
        {
            e = Scm1p->score+1;
            f = Srm1p->score+1;
            sa= Srcm1p->score+(pattern1[r-1]==pattern2[c-1]?0:1);
            switch(min3switch(sa,f,e))
            {
                case 1:
                    Sp->score = sa;
                    Sp->direction = DIAGONAL;
                    break;
                case 2:
                    Sp->score = f;
                    Sp->direction = DOWN;
                    break;
                case 3:
                    Sp->score = e;
                    Sp->direction = RIGHT;
                    break;
            }



            /* increase pointers */
            Sp++;Scm1p++;Srm1p++;Srcm1p++;
        }
        if(r<pattern1len&&max<MINIMUM(maxs[r+1].right,end))/* need to pad with VLN at end */
        {
            for(i=MINIMUM(maxs[r+1].right,end);i>max;i--)
            {
                 Sp->score=VLN;
                 Sp++;
            }
        }
    }

    /**********************************************
    * Allocate alignment an do trace back
    ***********************************************/
    ap = (CYCLICALIGNPAIR*) malloc(sizeof(CYCLICALIGNPAIR));
    if(ap==NULL) return NULL;

    ap->pattern2start = start+1; /* start is zero based */
    Sp=S+(pattern1len+1)*(width)-1; /* position in last column of last row */
    ap->score = Sp->score;
    /*figure out the length of the alignment */
    ap->length = 0;
    for(;;)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                Sp -=(width+1);
                break;
            case RIGHT:
                Sp--;
                break;
            case DOWN:
                Sp -= width;
                break;
            default:
                goto exit1;
        }
        ap->length++;
    }
    exit1:; /* goes here in case START */

    /* allocate memory for the alignment */
    ap->pattern1side  = (char*) malloc((ap->length+1)*sizeof(char));
    ap->pattern2side   = (char*) malloc((ap->length+1)*sizeof(char));
    if(ap->pattern1side==NULL||ap->pattern2side==NULL)
    {
        free(ap->pattern1side);
        free(ap->pattern2side);
        free(ap);
        return NULL;
    }

    /* trace back and fill in alignment */
    pat1sidepos  = ap->pattern1side+ap->length;
    pat2sidepos  = ap->pattern2side+ap->length;
    *pat1sidepos=*pat2sidepos='\0';
    pat1sidepos--; pat2sidepos--;
    pattern2pos = pattern2+(end-1); 
    pattern1pos = pattern1+(pattern1len-1);
    Sp=S+(pattern1len+1)*(width)-1; /* position in last column of last row */
    for(;;)
    {
        switch(Sp->direction)
        {
            case DIAGONAL:
                *pat1sidepos=*pattern1pos; *pat2sidepos=*pattern2pos;
                pat1sidepos--;pattern1pos--;pat2sidepos--;pattern2pos--;
                Sp -=(width+1);
                break;
            case RIGHT:
                *pat1sidepos = '-';
                *pat2sidepos = *pattern2pos;
                pat1sidepos--;pat2sidepos--;pattern2pos--;
                Sp--;
                break;
            case DOWN:
                *pat2sidepos='-';
                *pat1sidepos=*pattern1pos;
                pat2sidepos--;pat1sidepos--;pattern1pos--;
                Sp -= width;
                break;
            default:
                goto exit2;
        }
    }
    exit2:; /* goes here in case START */

    /* return pointer to alignment */
    return ap;

}



CYCLICALIGNPAIR* GetCyclicDiffAlignPair( char* pattern1, char* pattern2,
                                         int length1, int length2)
{
    static int gotreserve=0;
    int i,m,left,right;
    CYCLICALIGNPAIR *bestsofar,*ap;
    SD* S=NULL;
    int  Slen = 0;


    /* allocate the alignment matrix */
    Slen = (length1+1)*(length2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }


    if(!gotreserve)
    {
        /* allocate LOGDEPTH row limits */
        for(i=0; i<LOGDEPTH; i++)
        {
            Reserve.limits[i] = (LIMIT*) malloc(sizeof(LIMIT)*(MAXPATSIZE+1));
            if(Reserve.limits[i]==NULL)
            {
                free(S);
                return NULL;
            }
        }
        Reserve.next=-1;
        gotreserve=1;
    }


    ControlStack.height=-1; /* reset the control stack */

    /* add limits */
    Reserve.next++;
    L[0].limits=Reserve.limits[Reserve.next];
    Reserve.limits[Reserve.next]=NULL;

    /* add limits */
    Reserve.next++;
    L[length2].limits=Reserve.limits[Reserve.next];
    Reserve.limits[Reserve.next]=NULL;

    /* create the boundaries used in first alignment (rectangular limits)*/
    for(i=0;i<=length1;i++)
	{
		L[0].limits[i].right=0;
		L[0].limits[i].left=0;
        L[length2].limits[i].right=length2;
        L[length2].limits[i].left=length2;
	}

    /* duplicate the top pattern string */
    strcpy(doublepat,pattern2);
    strcat(doublepat,pattern2);

	/* get first alignment */
    bestsofar = GetBoundedDiffAlignment(pattern1, doublepat, length1, 0, length2,
                                    L[0].limits, L[length2].limits,S);

	/* copy alignment to L[0].limits */
    CyclicAPToLimit(0, length1, length2, bestsofar);

    /* copy L[0].limits to L[length2].limits */
    for(m=0;m<=length1;m++)
    {
        L[length2].limits[m].right = L[0].limits[m].right+length2;
        L[length2].limits[m].left = L[0].limits[m].left+length2;
    }

    /* push control stack */
    ControlStack.height++;
    ControlStack.left[ControlStack.height]=0;
    ControlStack.right[ControlStack.height]=length2;

    while(ControlStack.height!=-1) /* while stack is not empty */
	{
        /* pop control stack */
        left=ControlStack.left[ControlStack.height];
        right=ControlStack.right[ControlStack.height];
        ControlStack.height--;


        if(right-left>1)
        {
            m=(left+right)/2;
            /* the next is the routine for the alignment */
            ap = GetBoundedDiffAlignment(pattern1, doublepat, length1, m, m+length2,
                                        L[left].limits, L[right].limits,S);

            /* add limits */
            Reserve.next++;
            L[m].limits=Reserve.limits[Reserve.next];
            Reserve.limits[Reserve.next]=NULL;

            /* copy alignment to L[m].limits */
            CyclicAPToLimit(m,length1,length2,ap);
            if(ap->score<bestsofar->score)
            {
                FreeCyclicAlignPair(bestsofar);
                bestsofar = ap;
            }
            else
            {
                FreeCyclicAlignPair(ap);
            }

            /* push control stack */
            ControlStack.height++;
            ControlStack.left[ControlStack.height]=m;
            ControlStack.right[ControlStack.height]=right;

            /* push control stack */
            ControlStack.height++;
            ControlStack.left[ControlStack.height]=left;
            ControlStack.right[ControlStack.height]=m;
        }
        else
        {
            /* remove limits */
            Reserve.limits[Reserve.next]=L[left].limits;
            L[left].limits=NULL;
            Reserve.next--;
        }
    }
    /* remove limits */
    Reserve.limits[Reserve.next]=L[length2].limits;
    L[length2].limits=NULL;
    Reserve.next--;

    /* assign alignment starting positions */
    bestsofar->pattern1start = 1;  /* pattern1 is fixed, pattern2 is rotated */
    bestsofar->pattern1end   = length1;
    /* pattern2start was assigned inside of GetBoundedDiffAlignment() */
    if(bestsofar->pattern2start==1) bestsofar->pattern2end = length2;
    else bestsofar->pattern2end = bestsofar->pattern2start-1;


    /* free alignment matrix */
    free(S);

	return bestsofar;
}

/********************************	init_complement_ascii				****************************/

char* init_complement_ascii(void)

{
  /* complement has 256 entries so that finding the entries for A, C, G and T which are complement ascii values */
  /* require no calculation */
  int i;
  char *Complementascii=(char *)calloc(256,sizeof(long int));
  
  if(Complementascii==NULL){ return NULL; }
  
  for (i=0; i<256; i++)
	Complementascii[i]=i;

  Complementascii['A']='T';
  Complementascii['C']='G';
  Complementascii['G']='C';
  Complementascii['T']='A';
  Complementascii['a']='t';
  Complementascii['c']='g';
  Complementascii['g']='c';
  Complementascii['t']='a';

  return Complementascii;
}

/***********************************************************************/
char*   GetReverseComplement(char* original)
{
    static char *Complementascii = NULL;            /* precomputed array used for the reverse complement function */
    char *sourceptr,*destinptr;
    int length;
    char *buffer;



    /* safety check */
	if (Complementascii==NULL) {
        Complementascii = init_complement_ascii();
        if ( NULL == Complementascii ) return NULL;
	}

    /* find out how long the string is */
    length = strlen(original);

    /* allocate the new string */
    buffer = (char*) malloc(sizeof(char)*(length+1));
    if(buffer==NULL)
    {
		return NULL;
    }

    /* reverse complement pattern */
    sourceptr = original;
    destinptr = buffer+length; /* position at termination */
    *destinptr = '\0'; /* terminate string */
    destinptr--;
    while(*sourceptr!='\0')
    {
		(*destinptr)=Complementascii[(*sourceptr)];
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return buffer;
}

/***********************************************************************/
char*   GetComplement(char* original)
{
    char *sourceptr,*destinptr;
    int length;
    char *buffer;


    /* find out how long the string is */
    length = strlen(original);

    /* allocate the new string */
    buffer = (char*) malloc(sizeof(char)*(length+1));
    if(buffer==NULL)
    {
		return NULL;
    }

    /* complement pattern */
    sourceptr = original;
    destinptr = buffer+length; /* position at termination */
    *destinptr = '\0'; /* terminate string */
    destinptr--;
    while(*sourceptr!='\0')
    {
		(*destinptr)=(*sourceptr);
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return buffer;
}


/***********************************************************************/

GAPPEDALIGNPAIR*     GetGappedAlignPair(WDPALIGNPAIR* wdpap)
{
    GAPPEDALIGNPAIR* gapap;
    int *counts,*cntpos,*maxcntpos;
    char *patsrc,*seqsrc,*patdst,*seqdst;
    int i,length;

    /* allocate "blank" memory for the counts */
    counts = (int*) calloc(wdpap->patternlength,sizeof(int));

    /* find out maximum number of gaps before each each position */
    cntpos = counts;
    maxcntpos = counts+wdpap->patternlength-1;
    for(i=0,patsrc = wdpap->patternside ;*patsrc!='\0';patsrc++)
    {
        if(*patsrc=='-') i++;
        else
        {
            if(i>*cntpos) *cntpos=i;
            cntpos++; if(cntpos>maxcntpos) cntpos = counts;
            i=0;
        }
    }
    if(i>*cntpos) *cntpos=i; /* check for dashes at the end */

    /* calculate the total gapped length of alignment */
    for(length=0,cntpos=counts,patsrc=wdpap->patternside;*patsrc!='\0';)
    {
        if(*patsrc!='-')
        {
            length +=(1+*cntpos);
            cntpos++; if(cntpos>maxcntpos) cntpos = counts;
            patsrc++;
        }
        else do patsrc++; while(*patsrc=='-');
    }
    patsrc--;   /* back up to count gaps at end of string */
    while(*patsrc=='-'){length++;patsrc--;} /* add gaps at end to length */

    /* allocate memory */
    gapap = (GAPPEDALIGNPAIR*) malloc(sizeof(GAPPEDALIGNPAIR));
    if(gapap==NULL) return NULL;
    gapap->patternside = (char*)  malloc(sizeof(char)*(length+1));
    gapap->sequenceside = (char*) malloc(sizeof(char)*(length+1));
    if(gapap->patternside==NULL||gapap->sequenceside==NULL)
    {
        free(gapap->patternside);free(gapap->sequenceside);free(gapap);
        return NULL;
    }
    gapap->gappedlength = length;
    gapap->gappedpatternlength = wdpap->patternlength;
    for(cntpos=counts;cntpos<=maxcntpos;cntpos++)
    {
        gapap->gappedpatternlength += *cntpos;
    }

    /* fill in alignment */
    patsrc = wdpap->patternside; patdst = gapap->patternside;
    seqsrc = wdpap->sequenceside; seqdst = gapap->sequenceside;
    for(i=0,cntpos=counts;*patsrc!='\0';patdst++,patsrc++,seqdst++,seqsrc++)
    {
        if(*patsrc=='-') i++; /* i counts number of dashes actually found */
        else
        {
            for(;i<*cntpos;i++) {*patdst=*seqdst='-'; patdst++; seqdst++;}
            cntpos++;if(cntpos>maxcntpos) cntpos = counts;
            i = 0;
        }
        *patdst = *patsrc; *seqdst = *seqsrc;
    }
    *patdst = *seqdst = '\0'; /* terminate */

    free(counts); /* free buffer holding counts */
    return gapap;
}

/*****************************************************************************************************/
void printBasicGlobalMatrix(char* seq1, char* seq2,int len1, int len2, SD *S) {

    int     col, row, maxrow, val, maxv;
	SD		*Sp;

	printf("<BR><BR><PRE>\n");

	// top row indices
	maxv = max(len1,len2) * 10;
	for(row=1,maxrow=1;row<=maxv;row*=10) {
		maxrow*=10;
	} // maxrow now has largest dicimal bucket
	for(row=maxrow;row>=1;row/=10) {
		printf("      ");
		for(col=1;col<=len2;col++)
        {
			val = (len2-col) / row % 10;
			if (val || (len2-col)>row || (col==len2 && row==1))
				printf("%d ",val );
			else
				printf("  ");
		}
		printf("\n");
	}
	

	printf("--- - ");
    for(col=1;col<=len2;col++)
    {
		printf("%c ",seq2[col-1]);
	}
	printf("\n");


    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);

		printf("%3d %c ",len1-row, seq1[row-1]);

        for(col=1;col<=len2;col++)
        {
			if (Sp->score==2)
				printf("<font color=blue>%d</font> ",Sp->score);
			else
				printf("%d ",Sp->score);
            /* update pointers */
            Sp++;
        }

		printf("\n");

    }


	printf("</PRE><BR><BR>\n");


}


/*****************************************************************************************************/
SD* GetBasicGlobalMatrix(char* seq1, char* seq2, int len1,
                            int len2, int alpha, int indel, int* submatrix)
{
    int     col, row;
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp;

#ifdef _WIN_32_YES
    {
	__int64 SlenCheck; //to check for overflow

	/* this check was put here because of a problem of overflow when calculating the Slen
	and malloc() was just using the smaller number that came out as a result. A very hard and
	pesky error. Gelfand. */
	SlenCheck=(__int64)(len1+1)*(__int64)(len2+1)*(__int64)(sizeof(SD));
	if (SlenCheck>(__int64)INT_MAX) {
		return NULL;
	}
    }
#endif

    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

	
    /* compute SDD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = col*indel;
        Sp->direction = DIRECTIONNONE;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = row*indel;
        Sp->direction = DIRECTIONNONE;
    }


    for(row=1;row<=len1;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);



        for(col=1;col<=len2;col++)
        {
			Sp->score = submatrix[26*(seq1[row-1]-'A')+(seq2[col-1]-'A')];
			Sp->direction = DIRECTIONNONE;



            /* update pointers */
            Sp++;
        }

    }


    return S;
}

typedef struct {
	int start1,start2,end1, end2;
	BASICALIGNPAIR *ap;
} SEED_HIT;

   
#endif
