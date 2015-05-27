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

/*///////////////////////////////////////////////////////////////////

    PROFILE.H:  This file unifies all code needed to deal with
                profiles and profile alignment as it pertains
                to tandem and inverted repeats.

                It includes routines for:
                - Reading and writing profiles (see below)
                - Convert to and from standard composition
                - Profile alignment
                - Reading distance table from file


    Copyright (c) Gary Benson and Alfredo Rodriguez, 2004
    Last Revision Date April 8, 2004



    PROFILE STORAGE FORMAT EXPLANATION
    The file format used to store profiles is text based.
    Each profile is stored in its own line and in the following
	format:

    KEY PATLEN COPYNUMBER PROFLEN LEB36PROFILE

    For example, a profile record might look like:

    4321 6 6.67 7 E0XA5BG13QRK8M

    Here the key of the profile in an existing database is 4321,
    the size of the consensus pattern is 6, there are 6.67 copies
	of the pattern, the length of the profile is 7, and there are
	14 LEB36 digits representing the seven compositionsin the
	profile. Exactly two LEB36 digits per standard composition.

    LEB36 stands for little endian base-36 or alphanumeric:

    Decimal     Base36 (2-digit)    LEB36
    15          0F                  F0


*////////////////////////////////////////////////////////////////////


#ifndef PROFILE_H
#define PROFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

//////////////////////// BASIC DATA TYPES ///////////////////////////

typedef struct
{
    int key;
    int patlen;
    float copynum;
    int proflen;
    int rcflag;
    int* indices;
    int a,c,g,t,n,acgtCount;
    int nextoffset;
} PROFILE;

typedef struct
{
    int     start;
    int     distance;
    double  similarity;
    int     length;
    int     rcflag;
    int*    prof1side;
    int*    prof2side;
} PAP;



//////////////////////// INTERFACE ROUTINES   ///////////////////////

// Reads profile form file. Use -1 as offset if reading sequentially
// otherwise use nextoffset member of previous profile to ensure that
// a profile following a specific one is read. Returns NULL at the
// end of the file.
PROFILE*    ReadProfile(FILE* fp, int offset); //NO LONGER USED
PROFILE*    ReadProfileWithRC(FILE* fp, int offset, PROFILE** profrcret);

// Writes the profile to the current position of text file.
void        WriteProfile(FILE* fp, PROFILE* prof); //NO LONGER USED
void        WriteProfileWithRC(FILE* fp, PROFILE* prof, PROFILE* profrc, char *left, char *right);

// Frees memory associate with a profile
void        FreeProfile(PROFILE* prof);

// copy profile structure and data
PROFILE*    CopyProfile(PROFILE* prof);

// Creates a reverse complement of a profile
PROFILE* GetReverseComplementProfile(PROFILE* original);

// Generates a random profile
PROFILE* CreateRandomProfile(int minsize, int maxsize);



// Gets the id of the standard composition that is closest. The counts
// parameter is an array of five ints constaining counts for A,C,G,T,
// and - or insertions in that order.
int         GetCompositionId(int* counts);

// Gets the counts for a given standard composition id. Standard
// counts allways add up to 10.
void        GetCompositionCounts(int id, int* counts);

// Gets the Id of the standard composition that has the complement
// counts of the given one. If the id given has counts 2,4,3,0,1
// then it returns the id with counts 0,3,4,2,1.
int         GetComplementId(int id);

// Returns the majority character for a given composition. A more
// sophisticated implementation should return FASTA mixture symbols.
char     GetCompositionSymbol(int id);

// Prints symbols for array of compositions to file or stdout.
void    WriteSymbols(FILE* fp, int* indices, int length);



// Loads a distance table from file. The file is in binary format
// with the integers 777,10,and 5 at the beginning followed by
// array of distances as upper triagular including diagonal using
// unsigned characters. The size of the file is 12 + 501,501 =
// 501,513 bytes. Returns NULL if invalid file.
unsigned char* LoadDistanceTable(char* filename);

// Frees memory. Same as free(disttable);
void FreeDistanceTable(unsigned char* disttable);

// Aligns p1 against the cyclic permuation of p2 that produces the lowest distance
PAP* GetBestProfileAP(PROFILE* p1, PROFILE* p2, unsigned char* submatrix);

// Align p1 against current rotation of p2 (not cyclic)
PAP* GetFixedProfileAP(PROFILE* p1, PROFILE* p2, unsigned char* submatrix);


// Frees memory associated with a profile alignemnt struct and its members
void FreePAP(PAP* ap);


// Rescores a profile alignment with another substitution matrix.
// New escore updates the distance member of the PAP structure
void RescorePAP(PAP *ap, unsigned char* submatrix);









////////////////////////  IMPLEMENTATION  ///////////////////////////
/********************************************************************
*********************************************************************
********************************************************************/


// globals used for computing compositions
int cmp_ids[11][11][11][11];
int cmp_counts[1001][4];
char cmp_symbol[1001];
int cmp_init = 0;

// routines used for converting to/from little endian base-36
void LEB36ToInt(char* leb36,int* dec);
void IntToLEB36(int* dec, char* leb36);


void        FreeProfile(PROFILE* prof)
{
    if(prof!=NULL)
    {
        free(prof->indices);
        free(prof);
    }
    return;
}


PROFILE* CopyProfile(PROFILE* prof)
{
    int i;
    PROFILE* result;

    if(prof==NULL)
        return NULL;

    result = (PROFILE*) malloc(sizeof(PROFILE));
    if(result==NULL) return NULL;
    result->indices = (int*) malloc(sizeof(int)*(prof->proflen));
    if(result->indices==NULL)
    {
        free(result);
        return NULL;
    }

    // set members
    result->key = prof->key;
    result->patlen = prof->patlen;
    result->copynum = prof->copynum;
    result->proflen = prof->proflen;
    result->rcflag = prof->rcflag;
    result->nextoffset = prof->nextoffset;
    for(i=0;i<(result->proflen);i++) {
        result->indices[i] = prof->indices[i];
    }
    result->acgtCount =  prof->acgtCount;

    return result;
}


PROFILE*    ReadProfileWithRC(FILE* fp, int offset, PROFILE** profrcret)
{
    int c=0,i,key,patlen,proflen,proflenrc,number;
    float copynum;
    char digits[4];
    PROFILE *prof, *profrc;

    // go to the offset
	if(offset>=0) fseek(fp,offset,SEEK_SET);

    // read key, pattern length and profile length
    c = fscanf(fp,"%d %d %f %d %d",&key,&patlen,&copynum,&proflen,&proflenrc);
    if(c!=5) {
        return NULL;
    }

    // allocate a profile object and the corresponding array of indices
    prof = (PROFILE*) malloc(sizeof(PROFILE));
    if(prof==NULL) return NULL;
    prof->indices = (int*) malloc(sizeof(int)*proflen);
    if(prof->indices==NULL)
    {
        free(prof);
        printf("\nload error!!!2\n");
        return NULL;
    }
    // set members
    prof->key = key;
    prof->patlen  = patlen;
    prof->copynum  = copynum;
    prof->proflen  = proflen;
    prof->rcflag = 0;


    if (proflenrc>0) {
        profrc = (PROFILE*) malloc(sizeof(PROFILE));
        if (profrc) profrc->indices = (int*) malloc(sizeof(int)*proflenrc);
        if(profrc==NULL || profrc->indices==NULL)
        {
            free(prof->indices);
            free(prof);
            free(profrc);
            printf("\nload error!!!3\n");
            return NULL;
        }
        // set members
        profrc->key = key;
        profrc->patlen =  patlen;
        profrc->copynum =  copynum;
        profrc->proflen  = proflenrc;
        profrc->rcflag = 1;
    }

    for(i=0;i<proflen;i++)
    {
        // get next two digits
        c = fscanf(fp,"%2s",digits);
        if(c!=1)
        {
            free(prof->indices);
            free(prof);
            free(profrc->indices);
            free(profrc);
            printf("\nload error!!!4\n");
            return NULL;
        }

        // translate to decimal and save
        LEB36ToInt(digits,&number);
        prof->indices[i] = number;
    }

    // read the space
    if (proflenrc<=0) {

        prof->nextoffset = ftell(fp);

        *profrcret = NULL;
        return prof;
    }

    // read the rc profile
    for(i=0;i<proflenrc;i++)
    {
        // get next two digits
        c = fscanf(fp,"%2s",digits);
        if(c!=1)
        {
            free(prof->indices);
            free(prof);
            free(profrc->indices);
            free(profrc);
            return NULL;
        }

        // translate to decimal and save
        LEB36ToInt(digits,&number);
        profrc->indices[i] = number;
    }

    // read counts
    c = fscanf(fp," %d %d %d %d",&(prof->a),&(prof->c),&(prof->g),&(prof->t));
    if(c!=4) {
        return NULL;
    }
	prof->acgtCount = prof->a + prof->c + prof->g + prof->t;
	profrc->a=prof->t;
	profrc->t=prof->a;
	profrc->c=prof->g;
	profrc->g=prof->c;
	profrc->acgtCount = prof->acgtCount;

    if (profrcret) 
        *profrcret = profrc;

    prof->nextoffset = ftell(fp);

    return prof;
}

PROFILE*    ReadProfile(FILE* fp, int offset)
{
    int c,i,key,patlen,proflen,number;
    float copynum;
    char digits[4];
    PROFILE* prof;

    // go to the offset
	if(offset>=0) fseek(fp,offset,SEEK_SET);

    // read key, pattern length and profile length
    c = fscanf(fp,"%d %d %f %d",&key,&patlen,&copynum,&proflen);
    if(c!=4) return NULL;

    // allocate a profile object and the corresponding array of indices
    prof = (PROFILE*) malloc(sizeof(PROFILE));
    if(prof==NULL) return NULL;
    prof->indices = (int*) malloc(sizeof(int)*proflen);
    if(prof->indices==NULL)
    {
        free(prof);
        return NULL;
    }

    // set members
    prof->key = key;
    prof->patlen = patlen;
    prof->copynum = copynum;
    prof->proflen = proflen;
    prof->rcflag = 0;

    for(i=0;i<proflen;i++)
    {
        // get next two digits
        c = fscanf(fp,"%2s",digits);
        if(c!=1)
        {
            free(prof->indices);
            free(prof);
            return NULL;
        }

        // translate to decimal and save
        LEB36ToInt(digits,&number);
        prof->indices[i] = number;
    }

    prof->nextoffset = ftell(fp);

    return prof;
}


void  WriteProfile(FILE* fp, PROFILE* prof)
{
    int i,number,proflen;
    char digits[4];

    fprintf(fp,"%d %d %.2f %d ",prof->key,prof->patlen,prof->copynum,prof->proflen);
    proflen = prof->proflen;
    for(i=0;i<proflen;i++)
    {
        number = prof->indices[i];
        IntToLEB36(&number,digits);
        if(digits[1]=='\0')
        {
            digits[1] = '0';
            digits[2] = '\0';
        }
        fprintf(fp,"%s",digits);
    }
    fprintf(fp,"\n");

    return;
}

void  WriteProfileWithRC(FILE* fp, PROFILE* prof, PROFILE* profrc, char *left, char *right)
{
    int i,number,proflen;
    char digits[4];

	/* write repeat data */
    fprintf(fp,"%d %d %.2f %d %d ",prof->key,prof->patlen,prof->copynum,prof->proflen,profrc ? profrc->proflen : -1);
    proflen = prof->proflen;

	/* write profile */
    for(i=0;i<proflen;i++)
    {
        number = prof->indices[i];
        IntToLEB36(&number,digits);
        if(digits[1]=='\0')
        {
            digits[1] = '0';
            digits[2] = '\0';
        }
        fprintf(fp,"%s",digits);
    }

    if (NULL == profrc) { fprintf(fp,"\n"); return; }

    fprintf(fp," ");

	/* write profile RC */
    proflen = profrc->proflen;
    for(i=0;i<proflen;i++)
    {
        number = profrc->indices[i];
        IntToLEB36(&number,digits);
        if(digits[1]=='\0')
        {
            digits[1] = '0';
            digits[2] = '\0';
        }
        fprintf(fp,"%s",digits);
    }
    
    /* write a,c,g,t counts (from source sequence) for clustering hints */
    fprintf(fp," %d %d %d %d",prof->a,prof->c,prof->g,prof->t);

    /* write flanks */
    if (left&&right) { fprintf(fp," %s|%s",left,right); } 

	fprintf(fp,"\n");

    return;
}

void LEB36ToInt(char* leb36,int* dec)
{
    int i,mult,num;
    char *p;

    // step through leb36 digits
    num = 0;
    mult = 1;
    for(p=leb36;*p!='\0';p++)
    {
        if(*p<65) i = (*p)&0xF;
        else i = ((*p)&0x1F)+9;

        num += (i*mult);
        mult *= 36;
    }
    *dec = num; // assign before returning

    return;
}


void IntToLEB36(int* dec, char* leb36)
{
    int i,div;
    char *p,c;

    if(*dec==0)
    {
        leb36[0] = '0';
        leb36[1] = '\0';
    }

    div = *dec;
    p = leb36;
    while(div>0)
    {
        i = div%36;
        (i<10)?(c=i+48):(c=i+55);
        *p = c;
        div /= 36;
        p++;
    }
    if(p==leb36) *(p++) = '0';
    *p = '\0';

    return;
}


void InitializeCompositions(void)
{
    int a,c,g,t,dash,id=0,ms,mc;

    for(a=0;a<=10;a++)
    for(c=0;c<=10;c++)
    for(g=0;g<=10;g++)
    for(t=0;t<=10;t++)
    {
        dash = 10-a-c-g-t;
        if(dash>=0)
        {
            cmp_ids[a][c][g][t]=id;
            cmp_counts[id][0]=a;
            cmp_counts[id][1]=c;
            cmp_counts[id][2]=g;
            cmp_counts[id][3]=t;

            ms = 'A';
            mc = a;
            if(c>mc){mc=c;ms='C';};
            if(g>mc){mc=g;ms='G';};
            if(t>mc){mc=t;ms='T';};
            if(dash>mc){mc=dash;ms='-';}
            cmp_symbol[id] = ms;

            id++;
        }
    }

    cmp_init = 1;

    return;
}


int GetCompositionId(int* counts)
{
    int i,total=0;
    int aceil,cceil,gceil,tceil;
    int a,c,g,t,d;
    double an,cn,gn,tn,dn;
    double difa,difc,difg,dift,difd;
    int besta,bestc,bestg,bestt,bestd;
    double distance,bestdist=10000.0;

    if(!cmp_init) InitializeCompositions();

    for(i=0;i<5;i++) total+= counts[i];
    if(total==0) return 602; // id of perfect mixture
    if(total==10) return cmp_ids[counts[0]][counts[1]][counts[2]][counts[3]];

    // normalize counts as doubles;
    an = counts[0]*10.0/total;
    cn = counts[1]*10.0/total;
    gn = counts[2]*10.0/total;
    tn = counts[3]*10.0/total;
    dn = 10.0-an-cn-gn-tn;
    if(dn<0.0) dn=0.0;

    // figure out which count is closer in euclidean space
    aceil = counts[0]*10/total+1;
    cceil = counts[1]*10/total+1;
    gceil = counts[2]*10/total+1;
    tceil = counts[3]*10/total+1;
    for(a=aceil-1;a<=aceil;a++)
    for(c=cceil-1;c<=cceil;c++)
    for(g=gceil-1;g<=gceil;g++)
    for(t=tceil-1;t<=tceil;t++)
    {
        d = 10-a-c-g-t;
        if(d>=0) // valid composition
        {
            difa = an - a;
            difc = cn - c;
            difg = gn - g;
            dift = tn - t;
            difd = dn - d;

            distance = difa*difa+difc*difc+difg*difg+dift*dift+difd*difd;
            if(distance<bestdist)
            {
                bestdist = distance;
                besta = a;
                bestc = c;
                bestg = g;
                bestt = t;
                bestd = d;
            }
        }
    }

    return cmp_ids[besta][bestc][bestg][bestt];
}


void GetCompositionCounts(int id, int* counts)
{
    int i;

    if(!cmp_init) InitializeCompositions();

    if(id<0||id>1000)
    {
        for(i=0;i<5;i++) counts[i]=2; // perfect mixture
        return;
    }
    for(i=0;i<4;i++)
    {
        counts[i] = cmp_counts[id][i];
    }
    counts[4] = 10-counts[0]-counts[1]-counts[2]-counts[3];

    return;
}


int GetComplementId(int id)
{
    int a,c,g,t;

    if(!cmp_init) InitializeCompositions();

    if(id<1) return 0; // protect against gap or newgap

    a = cmp_counts[id][3];
    c = cmp_counts[id][2];
    g = cmp_counts[id][1];
    t = cmp_counts[id][0];

    return cmp_ids[a][c][g][t];
}


char    GetCompositionSymbol(int id)
{
    if(!cmp_init) InitializeCompositions();

    if(id<0) return '~'; // account for new gap index
    return cmp_symbol[id];
}


void    WriteSymbols(FILE* fp, int* indices, int length)
{
    int i,index;
    char symbol;

    for(i=0;i<length;i++)
    {
        index = indices[i];
        if(index>1000||index<0) index = -1;
        symbol = GetCompositionSymbol(index);
        fprintf(fp,"%c",symbol);
    }

    return;
}
void    WriteSymbols_stdout(int* indices, int length)
{
    int i,index;
    char symbol;

    for(i=0;i<length;i++)
    {
        index = indices[i];
        if(index>1000||index<0) index = -1;
        symbol = GetCompositionSymbol(index);
        fprintf(stdout,"%c",symbol);
    }

    return;
}



unsigned char* LoadDistanceTable(char* filename)
{
    int header[3],i;
    FILE* fp;
    unsigned char* buffer;

    // open the file for reading
    fp = fopen(filename,"rb");
    if(fp==NULL) return NULL;

    // read header
    i = fread(header,sizeof(int),3,fp);
    if(i!=3)
    {
        fclose(fp);
        return NULL;
    }

    // verify header
    if(header[0]!=777||header[1]!=10||header[2]!=5)
    {
        fclose(fp);
        return NULL;
    }

    // allocate buffer
    buffer = (unsigned char*) malloc(501501);
    if(buffer==NULL)
    {
        fclose(fp);
        return NULL;
    }

    // read data
    i = fread(buffer,1,501501,fp);
    if(i!=501501)
    {
        fclose(fp);
        return NULL;
    }

    // close and return;
    fclose(fp);
    return buffer;
}


void FreeDistanceTable(unsigned char* disttable)
{
    free(disttable);
    return;
}


#define PA_START           0
#define PA_DIAGONAL        1
#define PA_DOWN            2
#define PA_RIGHT           3
#define MAXPROFILESIZE     5000
#define PA_VLN             600000000
#define GAPINDEX           0
#define NEWGAPINDEX        -1

/* macros */
#define MINIMUM(a,b) (a<=b?a:b)
#define MAXIMUM(a,b) (a>=b?a:b)
#define min3switch(a,b,c) ((a<=b)?((a<=c)?1:3):((b<=c)?2:3))
#define EVALDISTANCE(r,c,sm) ((c)>=(r)?sm[1000*(r)-((r)*(r)-(r))/2+(c)]:\
                                    sm[1000*(c)-((c)*(c)-(c))/2+(r)])

#define EVALDISTANCENEW(r,rindex,c,sm) ((c)>=(r)?sm[(rindex)+(c)]:\
                                    sm[1000*(c)-((c)*(c)-(c))/2+(r)])


struct tagPALIMIT
{
    int left;
    int right;
};
typedef struct tagPALIMIT PALIMIT;

struct tagSPINEOFPALIMITS
{
    int     done;
    int     distance;
    PALIMIT*  limits;
};
typedef struct tagSPINEOFPALIMITS SPINEOFPALIMITS;

struct tagPSD
{
    int distance;
    int direction; 
};
typedef struct tagPSD PSD;


/* global variables */
SPINEOFPALIMITS paL[MAXPROFILESIZE];
PSD*    paS;
int    padoublepat[MAXPROFILESIZE*2];
int pa_need_init = 1;


/* paLOGDEPTH is used by Controlstack and limits paReserve */

#define paLOGDEPTH 16 /* ceil(log2(MAXPATTERNSIZE))+3 */

/* paReserve is a stack of limit structures */

struct tagRESERVEOFPALIMITS
{
    int next;
    PALIMIT*  limits[paLOGDEPTH];
};
typedef struct tagRESERVEOFPALIMITS RESERVE;

RESERVE paReserve;

/* addlimits and removelimits handle limit structures hanging off
   paL[] spine */

void addlimits(int a)
{
 paReserve.next++;
 paL[a].limits=paReserve.limits[paReserve.next];
 paReserve.limits[paReserve.next]=NULL;
}

 
void removelimits(int a)
{
 paReserve.limits[paReserve.next]=paL[a].limits;
 paL[a].limits=NULL;
 paReserve.next--;
}

/* Controlstack is used to get the correct order for the alignments */

struct tagCONTROLSTACK {
  int height;
  int left[paLOGDEPTH];
  int right[paLOGDEPTH];
};

typedef struct tagCONTROLSTACK CONTROLSTACK;

CONTROLSTACK Controlstack;

/* Controlstack functions */

int Controlstackempty()
{return(Controlstack.height==-1);}
 
void newControlstack() {Controlstack.height=-1;}


void pushControlstack(int a, int b)
{
   Controlstack.height++;
   Controlstack.left[Controlstack.height]=a;
   Controlstack.right[Controlstack.height]=b;
}

void popControlstack(int *a, int *b)
{
   *a=Controlstack.left[Controlstack.height];
   *b=Controlstack.right[Controlstack.height];
   Controlstack.height--;
}

/* changed this procedure to remove allocation for paL[], added allocation
   for paReserve, removed paouterleft */

int PA_InitializeGlobals()
{
    int i;
    
    /* allocate the maximum DS space that can possibly be needed */
    paS = (PSD*) malloc((MAXPROFILESIZE+1)*(MAXPROFILESIZE+1)*sizeof(PSD));
    if(paS==NULL) return 1;

    /* allocate paLOGDEPTH row limits */
    for(i=0; i<paLOGDEPTH; i++)
    {
        paReserve.limits[i] = (PALIMIT*) malloc(sizeof(PALIMIT)*(MAXPROFILESIZE+1));
        if(paReserve.limits[i]==NULL) return 1;
    }
    paReserve.next=-1;

    pa_need_init =0;

	return 0;
}

// creates the boundaries that are used in the first alignment
void init_rectangular_limits(int plength,int slength)
{
	static int i;

	for(i=0;i<=slength;i++)
	{
        paL[0].limits[i].right=0;
        paL[0].limits[i].left=0;
        paL[plength].limits[i].right=plength;
        paL[plength].limits[i].left=plength;
	}
}

void PAPToLimit(int start, int slength, int plength, PAP *pa)

{
    int row,col,i,length;
    int *ss, *ps;

    length=pa->length;
    ss=pa->prof1side;
    ps=pa->prof2side;
	  
	row=0;
	col=start;
    paL[start].limits[row].right=col;
	
	for(i=0;i<length;i++) 
	{
      if((ss[i]!=NEWGAPINDEX)&&(ps[i]!=NEWGAPINDEX))
	  {
	    row++;
	    col++;
	  }
      else if (ss[i]==NEWGAPINDEX) {col++;}
      else /* (ps[i]==NEWGAPINDEX) */ {row++;}
      paL[start].limits[row].right=col;
	}
	
	row=slength;
	col=start+plength; 
	
	for(i=length-1;i>=0;i--) 
	{
      paL[start].limits[row].left=col;
      if((ss[i]!=NEWGAPINDEX)&&(ps[i]!=NEWGAPINDEX))
	  {
	    row--;
	    col--;
	  }
      else if (ss[i]==NEWGAPINDEX) {col--;}
      else /* (ps[i]==NEWGAPINDEX) */ {row--;}
	}
	
    paL[start].limits[row].left=col;
}



PAP* GetBoundedProfileAP(int* left, int* top, int leftlen, int start,
                        int end,PALIMIT* mins, PALIMIT* maxs,
                        unsigned char* sm)
{
    int min,max,rindex2,leftRminus1;
    int r,c,i,width;
    PSD  *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */
    int sa,f,e;
    PAP *ap;
    int *seqpos,*patpos,*leftpos,*toppos;

    /****************************************
    * Do distances and directions for row zero
    *****************************************/
    r=0;
    Sp = paS; /* place at first position */
    Sp->distance = 0; paS->direction = PA_START;Sp++;
    max = MINIMUM(maxs[r].right,end);
    for(c=start+1,i=1;c<=max;c++,i++)
    {
        Sp->direction = PA_RIGHT;
        /* Sp->distance = i; */
        Sp->distance = EVALDISTANCENEW(0,0,top[c-1],sm)+(Sp-1)->distance;
		Sp++;
    }
    if(max<MINIMUM(maxs[r+1].right,end)) /* need to pad with PA_VLN at end */
    {
        for(i=MINIMUM(maxs[r+1].right,end);i>max;i--)
        {
             Sp->distance=PA_VLN;
             Sp++;
        }
    }

    /******************************************
    * Do distances and directions for start column
    *******************************************/
    Sp = paS; /* reset to first element */
    width = end-start+1;
    for(r=1;r<=leftlen;r++)
    {
        Sp+= width;
        /* Sp->distance = r; */
        Sp->distance = EVALDISTANCENEW(0,0,left[r-1],sm)+(Sp-width)->distance;
        Sp->direction =PA_DOWN;
        min = MAXIMUM(mins[r+1].left,start); /* min of next row */
        if(min>start) break;  /* break if left is cut-off */
    }

    /**********************************************
    * Do distances and directions in center of matrix
    ***********************************************/
	//rindex = 1000*(GAPINDEX)-((GAPINDEX)*(GAPINDEX)-(GAPINDEX))/2;

    for(r=1;r<=leftlen;r++)
    {
        Sp  = paS+(width*r);
        min = MAXIMUM(mins[r].left,start);
        max = MINIMUM(maxs[r].right,end);

        if(min>start) /* it's a left cut-off row */
        {
            Sp+= (min-start); /* position on first valid column */
            (Sp-1)->distance=PA_VLN;  /* pad with PA_VLN if left cut-off */
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

		leftRminus1 = left[r-1];
		rindex2 = 1000*(leftRminus1)-((leftRminus1)*(leftRminus1)-(leftRminus1))/2;

        for(;c<=max;c++)
        {
            /*
            e = Scm1p->distance+1;
            f = Srm1p->distance+1;
            sa= Srcm1p->distance+(left[r-1]==top[c-1]?0:1);
            */

			/*
            e = Scm1p->distance+EVALDISTANCE(GAPINDEX,top[c-1],sm);
            f = Srm1p->distance+EVALDISTANCE(GAPINDEX,left[r-1],sm);
            sa= Srcm1p->distance+EVALDISTANCE(left[r-1],top[c-1],sm);
			*/
			e = Scm1p->distance+EVALDISTANCENEW(0,0,top[c-1],sm);
            f = Srm1p->distance+EVALDISTANCENEW(0,0,leftRminus1,sm);
            sa= Srcm1p->distance+EVALDISTANCENEW(leftRminus1,rindex2,top[c-1],sm);

            switch(min3switch(sa,f,e))
            {
                case 1:
                    Sp->distance = sa;
                    Sp->direction = PA_DIAGONAL;
                    break;
                case 2:
                    Sp->distance = f;
                    Sp->direction = PA_DOWN;
                    break;
                case 3:
                    Sp->distance = e;
                    Sp->direction = PA_RIGHT;
                    break;
            }

            /* increase pointers */
            Sp++;Scm1p++;Srm1p++;Srcm1p++;
        }


        if(r<leftlen&&max<MINIMUM(maxs[r+1].right,end))/* need to pad with PA_VLN at end */
        {
            for(i=MINIMUM(maxs[r+1].right,end);i>max;i--)
            {
                 Sp->distance=PA_VLN;
                 Sp++;
            }
        }
    }

    /**********************************************
    * Allocate alignment an do trace back
    ***********************************************/
    ap = (PAP*) malloc(sizeof(PAP));
    if(ap==NULL) return NULL;

    ap->start = start;
    Sp=paS+(leftlen+1)*(width)-1; /* position in last column of last row */
    ap->distance = Sp->distance;
    /*figure out the length of the alignment */
    ap->length = 0;
    for(;;)
    {
        switch(Sp->direction)
        {
            case PA_DIAGONAL:
                Sp -=(width+1);
                break;
            case PA_RIGHT:
                Sp--;
                break;
            case PA_DOWN:
                Sp -= width;
                break;
            default:
                goto exit1;
        }
        ap->length++;
    }
    exit1:; /* goes here in case PA_START */

    /* allocate memory for the alignment */
    ap->prof1side  = (int*) malloc((ap->length)*sizeof(int));
    ap->prof2side   = (int*) malloc((ap->length)*sizeof(int));
    if(ap->prof1side==NULL||ap->prof2side==NULL)
    {
        free(ap->prof1side);
        free(ap->prof2side);
        free(ap);
        return NULL;
    }

    /* trace back and fill in alignment */
    seqpos  = ap->prof1side+ap->length-1;
    patpos  = ap->prof2side+ap->length-1;
    toppos = top+(end-1); 
    leftpos = left+(leftlen-1);
    Sp=paS+(leftlen+1)*(width)-1; /* position in last column of last row */
    for(;;)
    {
        switch(Sp->direction)
        {
            case PA_DIAGONAL:
                *seqpos=*leftpos; *patpos=*toppos;
                seqpos--;leftpos--;patpos--;toppos--;
                Sp -=(width+1);
                break;
            case PA_RIGHT:
                *seqpos = NEWGAPINDEX;
                *patpos = *toppos;
                seqpos--;patpos--;toppos--;
                Sp--;
                break;
            case PA_DOWN:
                *patpos=NEWGAPINDEX;
                *seqpos=*leftpos;
                patpos--;seqpos--;leftpos--;
                Sp -= width;
                break;
            default:
                goto exit2;
        }
    }
    exit2:; /* goes here in case PA_START */

    /* return pointer to alignment */
    return ap;
}

void DoubleProfile(int *p,int len, int* dest)
{
    int i,run;
    int *ps,*pd;

    run = len*2;
    ps = p;
    pd = dest;
    for(i=0;i<run;i++)
    {
        if(i==len) ps = p;
        *pd = *ps;
        ps++;
        pd++;
    }
    return;
}


PAP* GetBestProfileAP(PROFILE* p1, PROFILE* p2, unsigned char* submatrix)
{
    int m,left,right;
    PAP *bestsofar,*ap;
    int *s,*p,slength,plength;


    /* make sure memory for workspace and limits has been allocated */
    if(pa_need_init)
    {
        if(PA_InitializeGlobals()) return NULL;
        pa_need_init = 0;
    }

    /* set local variables */
    s = p1->indices;
    slength = p1->proflen;
    p = p2->indices;
    plength = p2->proflen;

	newControlstack();
	addlimits(0);
	addlimits(plength);
	init_rectangular_limits(plength,slength);

    /* duplicate the top profile into global padoublepat */
    DoubleProfile(p,plength,padoublepat);

	/* get first alignment */
    bestsofar = GetBoundedProfileAP(s, padoublepat, slength, 0, plength,
                                    paL[0].limits, paL[plength].limits,
                                    submatrix);

    /* copy alignment to paL[0].limits */
    PAPToLimit(0, slength, plength, bestsofar);

    /* copy paL[0].limits to paL[plength].limits */
    for(m=0;m<=slength;m++)
    {
        paL[plength].limits[m].right = paL[0].limits[m].right+plength;
        paL[plength].limits[m].left = paL[0].limits[m].left+plength;
    }

	pushControlstack(0,plength);

	while(!Controlstackempty())
	{
	  popControlstack(&left,&right);

	  if(right-left>1)
	  {
	    m=(left+right)/2;


	    /* the next is the routine for the alignment */
        ap = GetBoundedProfileAP(s, padoublepat, slength, m, m+plength,
                                    paL[left].limits, paL[right].limits,
                                    submatrix);

	    addlimits(m);
        /* copy alignment to paL[m].limits */
        PAPToLimit(m,slength,plength,ap);
        if(ap->distance<bestsofar->distance)
		//if (((ap->distance)/(255.0*ap->length)) < ((bestsofar->distance)/(255.0*bestsofar->length)) )
        {
            FreePAP(bestsofar);
            bestsofar = ap;
        }
        else
        {
            FreePAP(ap);
        }
	    pushControlstack(m,right);
	    pushControlstack(left,m);
	  }
      else removelimits(left);
	}

    removelimits(plength);

    // set rc flag
    if(p1->rcflag==p2->rcflag) bestsofar->rcflag = 0;
    else bestsofar->rcflag = 1;

    // set similarity member
    bestsofar->similarity = (1.0-(bestsofar->distance)/(255.0*bestsofar->length));

	return bestsofar;
}

void printSP(PSD *S, int leftlen, int toplen) {

	char c;
	int col,row;
    PSD *Sp; /* temporary pointers */

	Sp=S;

	/* top 
    for(col=1; col<=toplen; col++)
    {
        Sp++;
		if (PA_RIGHT ==  Sp->direction) c='R';
		if (PA_DOWN ==  Sp->direction) c='D';
		if (PA_DIAGONAL ==  Sp->direction) c='D';
		if (PA_START ==  Sp->direction) c='S';
		printf("%d%c ",Sp->direction,c);
    }
*/

    /* body of matrix */
	printf("\n\n");
    for(row=0;row<=leftlen;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+row*(toplen+1);

        for(col=0;col<=toplen;col++)
        {
			if (PA_RIGHT ==  Sp->direction) c='R';
			if (PA_DOWN ==  Sp->direction) c='D';
			if (PA_DIAGONAL ==  Sp->direction) c='\\';
			if (PA_START ==  Sp->direction) c='S';
			printf("%d-%c ",Sp->distance,c);

            /* update pointers */
            Sp++;
        }
		printf("\n");
    }

}

#define MAGIC_CENTER ((row<halfleft)?(startradius + gradvector * row):(startradius + gradvector * (leftlen-row)))

PAP* GetFixedProfileAPNarrowbandOld(PROFILE* p1, PROFILE* p2, unsigned char* sm, int homology) {
    
	int     col, row,length,lstart,rend,center,toplenreached,halfleft,oldrend;
    int     e;             /* equals e+indel */
    int     f;             /* equals f+indel */
    int     sa;             /* equals s+indel */
    PSD* S=NULL;
    PSD *Sp, *Srm1p, *Scm1p, *Srcm1p,*SrcmRp,*SL; /* temporary pointers */
    PAP* ap;
    int leftlen,toplen;
    int *left,*top;
	double svector,gradvector,diagonal;
	//double angle = acos(50.0 / ((200.0-homology-5)/2) ); /* to figure out the middle radius to allow 60% similarity, 100 diagonal is 140 sum of 2 other triangle sides */
	double angle = acos(((homology-5)/2.0) / 50.0 ); /* to figure out the middle radius to allow 60% similarity, 100 diagonal is 140 sum of 2 other triangle sides */
	int startradius,midradius;

    /* set lengths to be used */
	startradius=6;
    leftlen = p1->proflen;
    toplen = p2->proflen;
	halfleft = leftlen / 2;
	svector = (double)toplen/(double)leftlen;
	diagonal = sqrt(pow(toplen,2) + pow(leftlen,2)) ;
	midradius = (int) ( tan(angle) * (diagonal / 2.0) );
	midradius = max(midradius,startradius);

	//printf("midradius: %d\n",midradius); exit(1);
	
	//printf("%d-%d...",p1->key,p2->key); fflush(stdout);

	gradvector = (midradius - startradius) / (double)(halfleft);

    /* set shorthand for array of profile indices */
    left = p1->indices;
    top = p2->indices;
    
    /* allocate the alignment matrix */
    S = (PSD*) malloc(sizeof(PSD)*(leftlen+1)*(toplen+1));
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* compute Alignment Matrix */
    /* first element */
    Sp=S;
    Sp->distance  = 0;
    Sp->direction = PA_START;


    /* across the top */
    for(col=1; col<=toplen; col++)
    {
		center = (int)(0 * svector); 
		rend = center + startradius + 3;			rend=min(toplen,(rend));
        Sp++;
        Sp->direction = PA_RIGHT;
		if (col>(center + startradius))
			Sp->distance = 100000000;
		else
			Sp->distance = EVALDISTANCE(GAPINDEX,top[col-1],sm)+(Sp-1)->distance;
		if (rend==col) break;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=leftlen; row++)
    {
		center = (int)(row * svector); 
		lstart = center - (int)MAGIC_CENTER;		lstart=max(1,(lstart));
        Srm1p = Sp;
        Sp+=(toplen+1);
        Sp->direction = PA_DOWN;
        Sp->distance = EVALDISTANCE(GAPINDEX,left[row-1],sm)+Srm1p->distance;
		if (lstart>3) break;
    }

    /* body of matrix */
	toplenreached = 0;
	oldrend = startradius;
    for(row=1;row<=leftlen;row++)
    {


		/* calculate narrowband bounds */
		center = (int)(row * svector); 
		lstart = center - (int)MAGIC_CENTER;			lstart=max(1,(lstart));
		rend = center + (int)MAGIC_CENTER;			rend=min(toplen,(rend));

		//printf("row: %d center: %d svector: %.2lf\n",row,center,svector);

        /* set all pointer around column one of current row */
        Sp      =   S+lstart+row*(toplen+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-toplen;
        Srcm1p  =   Srm1p-1;
		SrcmRp	=	S+(row-1)*(toplen+1) + oldrend+1;

		if (lstart>1) {
			SL = Scm1p;
			for (e=0; e<5; e++) {
				if (SL < (S+row*(toplen+1))) break;
				SL->direction = PA_DOWN;
				SL->distance = 100000000;
				SL--;
			}
		}
		if (!toplenreached) {
			if (row>1) {
				for (e=0; e<5; e++) {
					if (SrcmRp >= (S+row*(toplen+1))) break;
					SrcmRp->direction = PA_RIGHT;
					SrcmRp->distance = 100000000;
					SrcmRp++; 
				}
			}
			if (rend==toplen) toplenreached=1;
		}

        for(col=lstart;col<=rend;col++)
        {
            e = Scm1p->distance+EVALDISTANCE(GAPINDEX,top[col-1],sm);
            f = Srm1p->distance+EVALDISTANCE(GAPINDEX,left[row-1],sm);
            sa = Srcm1p->distance+EVALDISTANCE(left[row-1],top[col-1],sm);

            switch(min3switch(sa,f,e))
            {
                case 1:
                    Sp->distance = sa;
                    Sp->direction = PA_DIAGONAL;
                    break;
                case 2:
                    Sp->distance = f;
                    Sp->direction = PA_DOWN;
                    break;
                case 3:
                    Sp->distance = e;
                    Sp->direction = PA_RIGHT;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }

		oldrend = rend;
    }
    /* set at last element */
    Sp = S+(toplen+1)*leftlen+toplen;


    /* allocate alignment structure */
    ap = (PAP*) malloc(sizeof(PAP));
    if(ap==NULL)
    {
        free(S);
        return NULL;
    }
    ap->start = 0;
    ap->distance = Sp->distance;

    /* get the length of the alignment */
    length = 0;
    while(Sp->direction!=PA_START)
    {
        length ++;
        switch(Sp->direction)
        {
            case PA_START:
                 break;
            case PA_DIAGONAL:
                 Sp = Sp -(toplen+2);
                 break;
            case PA_RIGHT:
                 Sp = Sp - 1;
                 break;
            case PA_DOWN:
                 Sp = Sp - (toplen+1);
                 break;
        }
		if (length>5000) {
			printf("ERROR! Cannot find start!\n");
			printSP(S, leftlen, toplen);
			printf("\n\n%d-vs-%d\n",leftlen,toplen);
			exit(1);
		}
    }
    ap->length = length;
	ap->similarity = (1.0-(ap->distance)/(255.0*ap->length));

	ap->prof1side = NULL;
	ap->prof2side = NULL;
    /* allocate memory for the alignment 
    ap->prof1side  = (int*) malloc((ap->length)*sizeof(int));
    ap->prof2side   = (int*) malloc((ap->length)*sizeof(int));
    if(ap->prof1side==NULL||ap->prof2side==NULL)
    {
        free(ap->prof1side);
        free(ap->prof2side);
        free(ap);
        free(S);
        return NULL;
    }
	*/

    /* trace back and fill in alignment 
    prof1sidepos  = ap->prof1side+ap->length-1;
    prof2sidepos  = ap->prof2side+ap->length-1;
    toppos = top+toplen-1; 
    leftpos = left+leftlen-1;
    Sp = S+(toplen+1)*leftlen+toplen;
    while(Sp->direction!=PA_START)
    {
        switch(Sp->direction)
        {
            case PA_START:
                break;
            case PA_DIAGONAL:
                *prof1sidepos=*leftpos; *prof2sidepos=*toppos;
                prof1sidepos--;leftpos--;prof2sidepos--;toppos--;
                Sp -=(toplen+2);
                break;
            case PA_RIGHT:
                *prof1sidepos = NEWGAPINDEX;
                *prof2sidepos = *toppos;
                prof1sidepos--;prof2sidepos--;toppos--;
                Sp--;
                break;
            case PA_DOWN:
                *prof2sidepos=NEWGAPINDEX;
                *prof1sidepos=*leftpos;
                prof2sidepos--;prof1sidepos--;leftpos--;
                Sp -= (toplen+1);
                break;
        }
    }*/

//116958278-117526956
//if (p1->key==116958278 && p2->key==117526956)
/*
{
	printSP(S, leftlen, toplen);
	printf("\n\n%d-vs-%d Distance: %d Length: %d sim: %.2lf\n",leftlen,toplen,ap->distance,ap->length, ap->similarity);
	exit(1);
}
*/

	//printf("|\n"); fflush(stdout);

    /* free alignment matrix */
    free(S);


    /* return pointer to alignment */
    return ap;
}

PAP* GetFixedProfileAPNarrowband(PROFILE* p1, PROFILE* p2, unsigned char* sm, int homology) {
    
    int     col, row,length,lstart,rend,center,toplenreached,halfleft,oldrend;
    int     e;             /* equals e+indel */
    int     f;             /* equals f+indel */
    int     sa;             /* equals s+indel */
    PSD* S=NULL;
    PSD *Sp, *Srm1p, *Scm1p, *Srcm1p,*SrcmRp,*SL; /* temporary pointers */
    PAP* ap;
    int leftlen,toplen;
    int *left,*top;
    double svector,diagonal,virtdiagonal,angle;
    double topangle, bandradius, midradius, minlen;
    int startradius,bandradius_int;


	// set lengths to be used 
    leftlen = p1->proflen;
	toplen = p2->proflen;
    diagonal = sqrt(pow(toplen,2) + pow(leftlen,2)) ;    // actual diagonal lenght
    minlen = min(toplen,leftlen);
    virtdiagonal = sqrt(pow(minlen,2) * 2) ;             // modified virtdiagonal to accountfor matrix not being a perfect square matrix
    angle = acos(((homology)/2.0) / 50.0 );              // to figure out the middle radiusto allow 60% similarity, 100 diagonal is 140 sum of 2 other triangle sides 


    startradius=6;


    topangle = atan( (double)leftlen / (double)toplen);
    halfleft = leftlen / 2;
    svector = (double)toplen/(double)leftlen;
    midradius =  tan(angle) * (virtdiagonal / 2.0);
    bandradius =  midradius / sin(topangle);
    bandradius_int = (int) max(bandradius,startradius) + 1; 

	//printf("\n\nleftlen:%d toptlen: %d midradius: %.2lf  bandradius: %d \n",leftlen, toplen, midradius, bandradius_int); exit(1);
	
	//printf("%d-%d...",p1->key,p2->key); fflush(stdout);

	//gradvector = (midradius - startradius) / (double)(halfleft);

    /* set shorthand for array of profile indices */
    left = p1->indices;
    top = p2->indices;
    
    /* allocate the alignment matrix */
    S = (PSD*) malloc(sizeof(PSD)*(leftlen+1)*(toplen+1));
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* compute Alignment Matrix */
    /* first element */
    Sp=S;
    Sp->distance  = 0;
    Sp->direction = PA_START;


    /* across the top */
    for(col=1; col<=toplen; col++)
    {
		center = (int)(0 * svector); 
		rend = center + bandradius_int + 3;			rend=min(toplen,(rend));
        Sp++;
        Sp->direction = PA_RIGHT;
		if (col>(center + bandradius_int))
			Sp->distance = 100000000;
		else
			Sp->distance = EVALDISTANCE(GAPINDEX,top[col-1],sm)+(Sp-1)->distance;
		if (rend==col) break;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=leftlen; row++)
    {
		center = (int)(row * svector); 
		lstart = center - bandradius_int;		lstart=max(1,(lstart));
        Srm1p = Sp;
        Sp+=(toplen+1);
        Sp->direction = PA_DOWN;
        Sp->distance = EVALDISTANCE(GAPINDEX,left[row-1],sm)+Srm1p->distance;
		if (lstart>3) break;
    }

    /* body of matrix */
	toplenreached = 0;
	oldrend = startradius;
    for(row=1;row<=leftlen;row++)
    {


		/* calculate narrowband bounds */
		center = (int)(row * svector); 
		lstart = center - bandradius_int;			lstart=max(1,(lstart));
		rend = center + bandradius_int;			rend=min(toplen,(rend));

		//printf("row: %d center: %d svector: %.2lf\n",row,center,svector);

        /* set all pointer around column one of current row */
        Sp      =   S+lstart+row*(toplen+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-toplen;
        Srcm1p  =   Srm1p-1;
		SrcmRp	=	S+(row-1)*(toplen+1) + oldrend+1;

		if (lstart>1) {
			SL = Scm1p;
			for (e=0; e<3; e++) {
				if (SL < (S+row*(toplen+1))) break;
				SL->direction = PA_DOWN;
				SL->distance = 100000000;
				SL--;
			}
		}
		if (!toplenreached) {
			if (row>1) {
				for (e=0; e<3; e++) {
					if (SrcmRp >= (S+row*(toplen+1))) break;
					SrcmRp->direction = PA_RIGHT;
					SrcmRp->distance = 100000000;
					SrcmRp++; 
				}
			}
			if (rend==toplen) toplenreached=1;
		}

        for(col=lstart;col<=rend;col++)
        {
            e = Scm1p->distance+EVALDISTANCE(GAPINDEX,top[col-1],sm);
            f = Srm1p->distance+EVALDISTANCE(GAPINDEX,left[row-1],sm);
            sa = Srcm1p->distance+EVALDISTANCE(left[row-1],top[col-1],sm);

            switch(min3switch(sa,f,e))
            {
                case 1:
                    Sp->distance = sa;
                    Sp->direction = PA_DIAGONAL;
                    break;
                case 2:
                    Sp->distance = f;
                    Sp->direction = PA_DOWN;
                    break;
                case 3:
                    Sp->distance = e;
                    Sp->direction = PA_RIGHT;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }

		oldrend = rend;
    }
    /* set at last element */
    Sp = S+(toplen+1)*leftlen+toplen;


    /* allocate alignment structure */
    ap = (PAP*) malloc(sizeof(PAP));
    if(ap==NULL)
    {
        free(S);
        return NULL;
    }
    ap->start = 0;
    ap->distance = Sp->distance;

    /* get the length of the alignment */
    length = 0;
    while(Sp->direction!=PA_START)
    {
        length ++;
        switch(Sp->direction)
        {
            case PA_START:
                 break;
            case PA_DIAGONAL:
                 Sp = Sp -(toplen+2);
                 break;
            case PA_RIGHT:
                 Sp = Sp - 1;
                 break;
            case PA_DOWN:
                 Sp = Sp - (toplen+1);
                 break;
        }
		if (length>5000) {
			printf("ERROR! Cannot find start!\n");
			printSP(S, leftlen, toplen);
			printf("\n\n%d-vs-%d\n",leftlen,toplen);
			exit(1);
		}
    }
    ap->length = length;
	ap->similarity = (1.0-(ap->distance)/(255.0*ap->length));

	ap->prof1side = NULL;
	ap->prof2side = NULL;
    /* allocate memory for the alignment 
    ap->prof1side  = (int*) malloc((ap->length)*sizeof(int));
    ap->prof2side   = (int*) malloc((ap->length)*sizeof(int));
    if(ap->prof1side==NULL||ap->prof2side==NULL)
    {
        free(ap->prof1side);
        free(ap->prof2side);
        free(ap);
        free(S);
        return NULL;
    }
	*/

    /* trace back and fill in alignment 
    prof1sidepos  = ap->prof1side+ap->length-1;
    prof2sidepos  = ap->prof2side+ap->length-1;
    toppos = top+toplen-1; 
    leftpos = left+leftlen-1;
    Sp = S+(toplen+1)*leftlen+toplen;
    while(Sp->direction!=PA_START)
    {
        switch(Sp->direction)
        {
            case PA_START:
                break;
            case PA_DIAGONAL:
                *prof1sidepos=*leftpos; *prof2sidepos=*toppos;
                prof1sidepos--;leftpos--;prof2sidepos--;toppos--;
                Sp -=(toplen+2);
                break;
            case PA_RIGHT:
                *prof1sidepos = NEWGAPINDEX;
                *prof2sidepos = *toppos;
                prof1sidepos--;prof2sidepos--;toppos--;
                Sp--;
                break;
            case PA_DOWN:
                *prof2sidepos=NEWGAPINDEX;
                *prof1sidepos=*leftpos;
                prof2sidepos--;prof1sidepos--;leftpos--;
                Sp -= (toplen+1);
                break;
        }
    }*/

//116958278-117526956
//if (p1->key==116958278 && p2->key==117526956)
/*
{
	printSP(S, leftlen, toplen);
	printf("\n\n%d-vs-%d Distance: %d Length: %d sim: %.2lf\n",leftlen,toplen,ap->distance,ap->length, ap->similarity);
	exit(1);
}
*/

	//printf("|\n"); fflush(stdout);

    /* free alignment matrix */
    free(S);


    /* return pointer to alignment */
    return ap;
}

PAP* GetFixedProfileAP(PROFILE* p1, PROFILE* p2, unsigned char* sm)
{
    int     col, row,length;
    int     e;             /* equals e+indel */
    int     f;             /* equals f+indel */
    int     sa;             /* equals s+indel */
    PSD* S=NULL;
    PSD *Sp, *Srm1p, *Scm1p, *Srcm1p; /* temporary pointers */
    PAP* ap;
    int *prof1sidepos,*prof2sidepos,*toppos,*leftpos;
    int leftlen,toplen;
    int *left,*top;


    /* set lengths to be used */
    leftlen = p1->proflen;
    toplen = p2->proflen;

    /* set shorthand for array of profile indices */
    left = p1->indices;
    top = p2->indices;

    
    /* allocate the alignment matrix */
    S = (PSD*) malloc(sizeof(PSD)*(leftlen+1)*(toplen+1));
    if(S==NULL)
    {
        return NULL; /* in case memory allocation fails */
    }

    /* compute Alignment Matrix */
    /* first element */
    Sp=S;
    Sp->distance  = 0;
    Sp->direction = PA_START;

    /* across the top */
    for(col=1; col<=toplen; col++)
    {
        Sp++;
        Sp->direction = PA_RIGHT;            
        Sp->distance = EVALDISTANCE(GAPINDEX,top[col-1],sm)+(Sp-1)->distance;
    }

    /* Down the left side */
    Sp=S;
    for(row=1; row<=leftlen; row++)
    {
        Srm1p = Sp;
        Sp+=(toplen+1);
        Sp->direction = PA_DOWN;
        Sp->distance = EVALDISTANCE(GAPINDEX,left[row-1],sm)+Srm1p->distance;
    }

    /* body of matrix */
    for(row=1;row<=leftlen;row++)
    {
        /* set all pointer around column one of current row */
        Sp      =   S+1+row*(toplen+1);
        Scm1p   =   Sp-1;
        Srm1p   =   Scm1p-toplen;
        Srcm1p  =   Srm1p-1;

        for(col=1;col<=toplen;col++)
        {
            e = Scm1p->distance+EVALDISTANCE(GAPINDEX,top[col-1],sm);
            f = Srm1p->distance+EVALDISTANCE(GAPINDEX,left[row-1],sm);
            sa = Srcm1p->distance+EVALDISTANCE(left[row-1],top[col-1],sm);

            switch(min3switch(sa,f,e))
            {
                case 1:
                    Sp->distance = sa;
                    Sp->direction = PA_DIAGONAL;
                    break;
                case 2:
                    Sp->distance = f;
                    Sp->direction = PA_DOWN;
                    break;
                case 3:
                    Sp->distance = e;
                    Sp->direction = PA_RIGHT;
                    break;
            }

            /* update pointers */
            Sp++;
            Srm1p++;
            Scm1p++;
            Srcm1p++;
        }
    }
    /* set at last element */
    Sp = S+(toplen+1)*leftlen+toplen;


    /* allocate alignment structure */
    ap = (PAP*) malloc(sizeof(PAP));
    if(ap==NULL)
    {
        free(S);
        return NULL;
    }
    ap->start = 0;
    ap->distance = Sp->distance;

    /* get the length of the alignment */
    length = 0;
    while(Sp->direction!=PA_START)
    {
        length ++;
        switch(Sp->direction)
        {
            case PA_START:
                 break;
            case PA_DIAGONAL:
                 Sp = Sp -(toplen+2);
                 break;
            case PA_RIGHT:
                 Sp = Sp - 1;
                 break;
            case PA_DOWN:
                 Sp = Sp - (toplen+1);
                 break;
        }
    }
    ap->length = length;


    /* allocate memory for the alignment */
    ap->prof1side  = (int*) malloc((ap->length)*sizeof(int));
    ap->prof2side   = (int*) malloc((ap->length)*sizeof(int));
    if(ap->prof1side==NULL||ap->prof2side==NULL)
    {
        free(ap->prof1side);
        free(ap->prof2side);
        free(ap);
        free(S);
        return NULL;
    }

    /* trace back and fill in alignment */
    prof1sidepos  = ap->prof1side+ap->length-1;
    prof2sidepos  = ap->prof2side+ap->length-1;
    toppos = top+toplen-1; 
    leftpos = left+leftlen-1;
    Sp = S+(toplen+1)*leftlen+toplen;
    while(Sp->direction!=PA_START)
    {
        switch(Sp->direction)
        {
            case PA_START:
                break;
            case PA_DIAGONAL:
                *prof1sidepos=*leftpos; *prof2sidepos=*toppos;
                prof1sidepos--;leftpos--;prof2sidepos--;toppos--;
                Sp -=(toplen+2);
                break;
            case PA_RIGHT:
                *prof1sidepos = NEWGAPINDEX;
                *prof2sidepos = *toppos;
                prof1sidepos--;prof2sidepos--;toppos--;
                Sp--;
                break;
            case PA_DOWN:
                *prof2sidepos=NEWGAPINDEX;
                *prof1sidepos=*leftpos;
                prof2sidepos--;prof1sidepos--;leftpos--;
                Sp -= (toplen+1);
                break;
        }
    }


    /* free alignment matrix */
    free(S);


    /* return pointer to alignment */
    return ap;
}




PROFILE* GetReverseComplementProfile(PROFILE* original)
{
    PROFILE* newprof;
    int i, *sourceptr,*destinptr;

    /* allocate the new PROFILE */
    newprof = (PROFILE*) malloc(sizeof(PROFILE));
    if(newprof==NULL) return NULL;
    newprof->key = original->key;
    newprof->patlen = original->patlen;
    newprof->copynum = original->copynum;
    newprof->proflen = original->proflen;
    newprof->nextoffset = original->nextoffset;
    newprof->rcflag = !original->rcflag;

    newprof->indices = (int*) malloc(sizeof(int)*newprof->proflen);
    if(newprof->indices==NULL)
    {
        free(newprof);
        return NULL;
    }

    /* reverse complement the profile */
    sourceptr = original->indices;
    destinptr = newprof->indices+newprof->proflen-1; /* place on last index */
    for(i=0;i<newprof->proflen;i++)
    {
        *destinptr = GetComplementId(*sourceptr);

        /* update pointers */
        destinptr--;
        sourceptr++;
    }

    /* return buffer containg reverse complement */
    return newprof;
}


#define p_RANDRANGE(MIN,MAX) (rand()%(MAX-MIN+1)+MIN)


PROFILE* CreateRandomProfile(int minsize, int maxsize)
{
    PROFILE* prof;
    int i,size,*p;
    static int key=1;


    // allocate a profile stucture
    prof = (PROFILE*) malloc(sizeof(PROFILE));
    if(prof==NULL) return NULL;

    // generate a random size
    if(minsize<1) minsize = 1;
    if(maxsize<1) maxsize = 1;
    if(minsize>maxsize){i=minsize;minsize=maxsize;maxsize=i;}
    size = p_RANDRANGE(minsize,maxsize);

    // set members 
    prof->key = key++;
    prof->copynum = 10.0f;
    prof->patlen = 0;
    prof->proflen = size;
    prof->nextoffset = 0;
    prof->rcflag = 0;

    // allocate the indices
    prof->indices = (int*) malloc(sizeof(int)*size);
    if(prof->indices==NULL)
    {
        free(prof);
        return NULL;
    }

    // generate the indices
    for(i=0,p=prof->indices;i<size;i++,p++)
    {
        *p = p_RANDRANGE(0,1000);
        if(GetCompositionSymbol(*p)!='-') prof->patlen++;
    }


    return prof;
}




void FreePAP(PAP* ap)
{
    free(ap->prof2side);
    free(ap->prof1side);
	free(ap);
	return;
}


void RescorePAP(PAP *ap, unsigned char* submatrix)
{
    int i,distance,c,r;

    /* for every position add the distance */
    for(i=distance=0;i<ap->length;i++)
    {
        c = ap->prof2side[i];
        r = ap->prof1side[i];
        if(c<0) c=0;
        if(r<0) r=0;
        distance += EVALDISTANCE(r,c,submatrix);
    }

    ap->distance = distance;
    ap->similarity = (1.0-(ap->distance)/(255.0*ap->length));

    return;
}



#endif
