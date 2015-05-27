#ifndef _ALN_C
#define _ALN_C

/***********************************************************************/
/***********************************************************************/
/*                                                                     */
/*                          Alignment                                  */
/*                                                                     */
/***********************************************************************/
/***********************************************************************/


#include <stdio.h>

int *SM = NULL;
#define match(a, b) (submatrix[256*(a)+(b)])


/* contants, and macros */
#define START           0
#define DIAGONAL        3
#define LONGDIAGONAL    6
#define RIGHT           8
#define DOWN            9
#define MINIMUM(a,b) (a<=b?a:b)
#define MAXIMUM(a,b) (a>=b?a:b)
#define min3switch(a,b,c) ((a<=b)?((a<=c)?1:3):((b<=c)?2:3))
#define max3switch(a,b,c) ((a>=b)?((a>=c)?1:3):((b>=c)?2:3))
#define max4switch(a,b,c,d) ((a>=b)?((a>=c)?((a>=d)?1:4):((c>=d)?3:4)):\
                            ((b>=c)?((b>=d)?2:4):((c>=d)?3:4)))

/* structures */
typedef struct
{

#ifdef PRINT_ALIGNMENTS
    char*   sequence1side;
    char*   sequence2side;
    int     sequence1start;
    int     sequence1end;
    int     sequence2start;
    int     sequence2end;
#endif

    int     score;
    int     length;
    int     mism;
    int     gaps;
}   BASICALIGNPAIR;


typedef struct
{
    int          score;
    char     direction;
}   SD;


/***********************************************************************/
/***********************************************************************/

int *init_sm(int match,int mismatch)

{
  int i,j,*currint;
 
  if (NULL==SM) {
    if (NULL==(SM=(int *)scalloc(256*256,sizeof(int)))) { return NULL; }
  }

  /* generate 256x256 into matrix */
  for(i=0,currint=SM;i<=255;i++)
  {
    for(j=0;j<=255;j++,currint++)
    {
        *currint=mismatch;
    }
  }

  SM['A'*256+'A']=match;
  SM['A'*256+'a']=match;
  SM['a'*256+'A']=match;
  SM['a'*256+'a']=match;

  SM['C'*256+'C']=match;
  SM['C'*256+'c']=match;
  SM['c'*256+'C']=match;
  SM['c'*256+'c']=match;

  SM['G'*256+'G']=match;
  SM['G'*256+'g']=match;
  SM['g'*256+'G']=match;
  SM['g'*256+'g']=match;

  SM['T'*256+'T']=match;
  SM['T'*256+'t']=match;
  SM['t'*256+'T']=match;
  SM['t'*256+'t']=match;

  return SM;
}

/***********************************************************************/
/***********************************************************************/

BASICALIGNPAIR* align( char* seq1, char* seq2, int len1, int len2, int *submatrix, int alpha, int indelstart, int indelextend) 
{

    int     col, row,length, gaps, mism;
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
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(len2+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-len2;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=len2;col++)
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

	/* store bestscore in final column */
	if((Sp-1)->score>bestscore)
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
    
    /* free alignment matrix */
    free(S);


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


      fprintf(stderr,"\n");
    }


    return ap;

#endif

}

/***********************************************************************/
void FreeBasicAlignPair(BASICALIGNPAIR* ap)
{
#ifdef PRINT_ALIGNMENTS
    free(ap->sequence1side);
    free(ap->sequence2side);
#endif
    free(ap);
    return;
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
char*   GetReverse(char* original)
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

    /* reverse pattern */
    sourceptr = original;
    destinptr = buffer+length; /* position at termination */
    *destinptr = '\0'; /* terminate string */
    destinptr--;
    while(*sourceptr!='\0')
    {
		(*destinptr)=*sourceptr;
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return buffer;
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

#endif

