#ifndef _ALN_NB_C
#define _ALN_NB_C

/***********************************************************************/
/***********************************************************************/
/*                                                                     */
/*                          Narrowband Alignment                       */
/*                                                                     */
/***********************************************************************/
/***********************************************************************/


#include <stdio.h>

#include "aln.h"


/***********************************************************************/
/***********************************************************************/

BASICALIGNPAIR* align_nb( char* seq1, char* seq2, int len1, int len2, int *submatrix, int alpha, int indelstart, int indelextend, int maxgaps) 
{

    int     col, row,length, gaps, mism, from, to;
    int     eb;             /* equals e+indel */
    int     fb;             /* equals f+indel */
    int     sa;             /* equals s+match or mismatch */
    SD* S=NULL;
    int  Slen = 0;
    SD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */


    BASICALIGNPAIR* ap;
    char *seq1pos,*seq2pos,*seq1sidepos,*seq2sidepos;
    int bestscore,bestscorerow;


    /* allocate the alignment matrix */
    Slen = (len1+1)*(len2+1);
    S = (SD*) malloc(sizeof(SD)*Slen);
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* set up to find bestscore in last column */ 
    bestscore=-1000000;
    bestscorerow=0;

    /* compute SDD Matrix */
    /* first element */
    Sp=S;
    Sp->score = 0;
    Sp->direction = START;

    /* across the top, pattern */
    for(col=1; col<=len2; col++)
    {
        Sp++;
        Sp->score = indelstart + (col-1)*indelextend;
        Sp->direction = RIGHT;
    }

    /* Down the left side, sequence */
    Sp=S;
    for(row=1; row<=len1; row++)
    {
        Sp+=(len2+1);
        Sp->score = indelstart + (row-1)*indelextend;
        Sp->direction = DOWN;
    }

    /* body of matrix */
    for(row=1;row<=len1;row++)
    {

        /* band */
        from = row - maxgaps;
        to = row + maxgaps;
        from = max(from ,1);
        to = min(to ,len2);

        /* set all pointer around column one of current row */
        Sp      =   S + from + row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;

        /* add bound value on left */
        if (from>1) {
            (Sp-1)->direction = DOWN;
            (Sp-1)->score = -1000000;          
        }

        for(col=from;col<=to;col++)
        {

            /* E */
            eb = Scm1p->score + ((Scm1p->direction >= RIGHT) ? indelextend : indelstart);

            /* F */
            fb = Srm1p->score + ((Srm1p->direction >= RIGHT) ? indelextend : indelstart);

            /* S */
            sa = Srcm1p->score+
                 match(seq1[row-1],seq2[col-1]);
						
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

        /* add bound values on right for next row */
        if (to<len2) {
            Sp->direction = RIGHT;
            Sp->score = -1000000;
        }

        if((Sp-1)->score>bestscore && to>=len2)
        {
          bestscore=(Sp-1)->score;
          bestscorerow=row;
        }

    }

    /* get the length of the alignment */
    Sp = S+(len2+1)*bestscorerow+len2; /* set at bestscore element in final column */



/* if we are not printing alignments, let's quickly finish up */
#ifndef PRINT_ALIGNMENTS
{
    char *src1,*src2;

    ap = (BASICALIGNPAIR*) malloc(sizeof(BASICALIGNPAIR));
    if(ap==NULL) return NULL;
    
    src1 = seq1 + bestscorerow - 1;
    src2 = seq2 + len2 - 1;

    /* Sp has ben set to bestscore element above */
    ap->score  = Sp->score;

    length = 0; mism = 0; gaps = 0;
    while(Sp->direction!=START)
    {
        switch(Sp->direction)
        {

            case DIAGONAL:
                 length++;
                 if (alpha != match(*src1,*src2)) mism++;
                 //printf("\n%d",match(*src1,*src2));
                 Sp = Sp -(len2+2);
                 src1--;
                 src2--;
                 break;

            case RIGHT:
                 length++; 
                 gaps++;
                 Sp = Sp - 1;
                 src2--;
                 break;
            case DOWN:
                 length++; 
                 gaps++;
                 Sp = Sp - (len2+1);
                 src1--;
                 break;
        }
    }

    ap->length = length;
    ap->mism = mism;
    ap->gaps = gaps;

    /* free alignment matrix */
    free(S);

    return ap;   
}
#else



    /***********************************************************************/
    /***********************************************************************/
    /**  Code below only gets executed if PRINT_ALIGNMENTS is not defined **/
    /***********************************************************************/
    /***********************************************************************/



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

    /* count characters in alignment to find start for sequence1 and indels,mismatches */
    length=mism=gaps=0;
    for(row=0;row<ap->length;row++) {
      if (ap->sequence1side[row]!='-') length++; 

      if (ap->sequence1side[row]=='-' || ap->sequence2side[row]=='-') 
        gaps++;
      else if (alpha != match(ap->sequence1side[row],ap->sequence2side[row] ))
        mism++;
    }
                
    ap->sequence1start = ap->sequence1end-length+1;
    ap->mism = mism;
    ap->gaps = gaps;

    


    {
      fprintf(stderr,"\n\nTP: (%d-%d)", ap->sequence1start, ap->sequence1end);
      fprintf(stderr,"\n>>: ");
      for(row=0;row<ap->length;row++)       
        fprintf(stderr,"%c",ap->sequence1side[row]);

      fprintf(stderr,"\n    ");
      for(row=0;row<ap->length;row++)       
        fprintf(stderr,"%c", (alpha != match(ap->sequence1side[row],ap->sequence2side[row])) ? ' ' : '*');

      fprintf(stderr,"\n>>: ");  
      for(col=0;col<ap->length;col++)       
        fprintf(stderr,"%c",ap->sequence2side[col]);
      fprintf(stderr,"\nRD: (%d-%d)", ap->sequence2start, ap->sequence2end);

      fprintf(stderr,"bestrow: %d  bestscore: %d  mism: %d gaps: %d\n",bestscorerow, bestscore, mism, gaps);


      fprintf(stderr,"\n");
    }

    /* print patrix 
    fprintf(stderr,"\n");
    Sp=S;
    for(row=0;row<=len1;row++)
    {
        fprintf(stderr,"\n");
        for(col=0;col<=len2;col++) {
          fprintf(stderr,"\t%d",Sp->score);
          Sp++;
        }
    }
    fprintf(stderr,"\n");
    */

    /* free alignment matrix */
    free(S);


    return ap;

#endif

}

/***********************************************************************/

#endif

