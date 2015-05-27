/*
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
*   ************************************************************
*
*   Please direct all questions related to this program or bug reports or updates to
*   Yevgeniy Gelfand (ygelfand@bu.edu)
*
*   Authors: Yevgeniy Gelfand
*/


/* CHANGES

  Apr 4 - changed minrepresentation to use reverse compliment from the file instead of
          RC of the forward concensus
 
        - reading profile now exits correctly on error

*/

//#define RECORDS_PER_FILE (5000)
#define RECORDS_PER_FILE (100000)
#define OUTPREFIX "reads"

/***************************************************************
    redund.c    :   Program that takes a file with a list of
                    profiles and remove redundancy based on alphabetic rotation of profile

	            Usage: 
		 
 		    redund.exe INPUTFILE OUTPUTFILE [-s]
 			OR
 		    redund.exe INPUTFOLDER OUTPUTFILE [-n]		    


		    This also produces an index file called "outfile".removed that has an
		    index of the preserved repeat followed by ids of removed repeats on each line.
		    Preserved repeats that have no removed redundant repeats ARE NOT in that
		    file.

		    If -s switch is provided, this simply sorts the file. You will need to call this
		    again to actually remove redundancy.

		    Use -n switch to make the program output a single file (not broken up in multiples.)

    VERSION 	:   1.00
		
*/

#include <sys/types.h>
#include <dirent.h>
#include <malloc.h>
#include <strings.h>


//#define _WIN_32_YES

#ifndef _WIN_32_YES
#define __int64 long long
#define LARGE_INTEGER time_t
#define I64d ld
#endif


#define max3(a,b,c) (((a)>=(b))?(((a)>=(c))?(a):(c)):(((b)>=(c))?(b):(c)))
#define min3(a,b,c) (((a)<=(b))?(((a)<=(c))?(a):(c)):(((b)<=(c))?(b):(c)))
#define max(a,b) (((a)>=(b))?(a):(b))
#define min(a,b) (((a)<=(b))?(a):(b))


/***************************************************************/
//#define EASY_LIFE_DEBUG

//#define PROCLU_DEBUG


#include <stdio.h>

#ifdef _WIN_32_YES
#include <windows.h>
#endif

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

char *Complementascii = NULL;

#include "easylife.h"
#include "easylife.c"
#include "profile.h"


/*******************************************************************************************/

typedef struct {
    char            *inputfile;
    char            *outputfile;
    char            *outputfile2;
    char            *d_name;
    FILE *in;
    PROFILE *prof, *profrc;
    int *minRepresentation;
    int minrlen;
    int dir;
} FITEM_STRUCT;

typedef struct {
    PROFILE *prof, *profrc;
    int *minRepresentation;
    int minrlen;
} PITEM_STRUCT;


/*******************************************************************************************/
int* intdup(int a) { 
	
    int *b = (int*)scalloc(1,sizeof(int));
    if(b==NULL)
    {
        printf("\nERROR: Insuficient memory to create an integer\n\n");
        exit(1);
    } 
    *b = a;

    return b;
}

/*******************************************************************************************/
int pindcmp(int *p1, int *p2, int len1, int len2){

  int i,len;

  if (len1 > len2)
        return 1;
  else if (len1 < len2)
        return -1;

  len = min(len1,len2);

  for (i=0; i<len; i++) {
   if (p1[i] > p2[i]) return 1;
   else if (p1[i] < p2[i]) return -1;
  }

  return 0;
}

/*******************************************************************************************/
int name_cmp( const void *item1, const void *item2 ) {

    FITEM_STRUCT *fiptr1,*fiptr2;

    fiptr1 = ((FITEM_STRUCT*)item1);
    fiptr2 = ((FITEM_STRUCT*)item2);

    return strcmp(fiptr1->d_name,fiptr2->d_name);
}

/*******************************************************************************************/
int arsize_and_min_rep_cmp( const void *item1, const void *item2 ) {

    FITEM_STRUCT *fiptr1,*fiptr2;

    fiptr1 = ((FITEM_STRUCT*)item1);
    fiptr2 = ((FITEM_STRUCT*)item2);

    // 0th - NULL?
    if (NULL == fiptr1->minRepresentation && NULL == fiptr2->minRepresentation) return 0;
    if (NULL == fiptr1->minRepresentation) return 1;
    if (NULL == fiptr2->minRepresentation) return -1;

/*
    // 1st - arlen
    if (fiptr1->prof->acgtCount > fiptr2->prof->acgtCount)
        return 1;
    else if (fiptr1->prof->acgtCount < fiptr2->prof->acgtCount)
        return -1;
*/
    // 2nd
    return pindcmp(fiptr1->minRepresentation,fiptr2->minRepresentation,fiptr1->minrlen,fiptr2->minrlen);
}

/*******************************************************************************************/
int arsize_and_min_rep_cmp_pitem( const void *item1, const void *item2 ) {

    PITEM_STRUCT *fiptr1,*fiptr2;

    fiptr1 = ((PITEM_STRUCT*)item1);
    fiptr2 = ((PITEM_STRUCT*)item2);

    // 0th - NULL?
    if (NULL == fiptr1->minRepresentation && NULL == fiptr2->minRepresentation) return 0;
    if (NULL == fiptr1->minRepresentation) return 1;
    if (NULL == fiptr2->minRepresentation) return -1;

/*
    // 1st - arlen
    if (fiptr1->prof->acgtCount > fiptr2->prof->acgtCount)
        return 1;
    else if (fiptr1->prof->acgtCount < fiptr2->prof->acgtCount)
        return -1;
*/

    // 2nd
    return pindcmp(fiptr1->minRepresentation,fiptr2->minRepresentation,fiptr1->minrlen,fiptr2->minrlen);
}

/*******************************************************************************************/
int * pintdup(int *src,int size) {
 int i,*res;

 res = smalloc(sizeof(int) * size);
 for (i=0; i<size; i++) {
  res[i]=src[i];
 } 

 return res;
}

/*******************************************************************************************/
int ftouch(char *in) {

	FILE *fp;

	if ((fp=fopen(in,"r")) == NULL) return 0;

	fclose(fp); return 1;
}

int RC = 0;

/*******************************************************************************************/
int* MinimumRepresentation(int *indptr, int len, int *indptrrc, int lenrc, int *minrlen, int identical_only) {

 int i,j;
 int *src,*srcrc,*tar,*sourceptr,*destinptr,tmp,*indices;

 RC = 0;

 if (NULL == indptr || NULL == indptrrc) return;

 src = pintdup(indptr,len);
 srcrc = pintdup(indptrrc,lenrc);

 tar = pintdup(indptr,len);


/* added on feb 6, 2013 */
if (identical_only) {

   // minimum of the forward and reverse profiles
   if (pindcmp(srcrc,tar,lenrc,len)<0) { free(tar); tar = pintdup(srcrc,lenrc); RC=1; }

} else {


 // FORWARD - test each rotation and find the minimum one
 for (i=1; i<len; i++) {

   //printf("\n\nstart: %s",src);

   // rotate by 1
   tmp=src[0];
   for (j=1; j<len; j++) {
    src[j-1]=src[j];
   }
   src[len-1]=tmp;

   //printf("\nend: %s comp: %d\n\n",src,strcmp(src,tar));
   //exit(1);

   if (pindcmp(src,tar,len,len)<0) { free(tar); tar = pintdup(src,len);}

 }


 // REVERSE - test each rotation and find the minimum one
 for (i=1; i<lenrc; i++) {

   //printf("\n\nstart: %s",src);

   // rotate by 1
   tmp=srcrc[0];
   for (j=1; j<len; j++) {
    srcrc[j-1]=srcrc[j];
   }
   srcrc[lenrc-1]=tmp;

   //printf("\nend: %s comp: %d\n\n",src,strcmp(src,tar));
   //exit(1);

   if (pindcmp(srcrc,tar,lenrc,len)<0) { free(tar); tar = pintdup(srcrc,lenrc); RC=1; }
 }

}

 // save rotated
 free(src);
 free(srcrc);
 
 *minrlen = (0==RC) ? len : lenrc;

 return tar;
}


/*******************************************************************************************/
int main(int argc, char** argv)
{
    int SINGLE_OUTFILE,SORT_ONLY,IDENTICAL_ONLY;
    FILE   *fpto,*fpto2;
    char    bigtempbuf[2000],inputfile[1000], outputfile[1000], outputfile2[1000];
    int i,filescreated=1;;
    time_t startTime;
    FITEM_STRUCT *fiptr,*fiptr2,*lastwrite;
    PITEM_STRUCT *piptr;
    EASY_ARRAY *FARRAY = NULL;
    struct dirent *de=NULL;
    DIR *d=NULL;
    long long int nwritten,nread;

    // verify parameter count
    if(argc<3)
    {
        printf("\nREDUND v1.00 - removes redundancy based on rotated profiles.");
        printf("\nauthors: Yevgeniy Gefland");
        printf("\n\nUsage:\n\n%s INPUTFILE OUTPUTFILE [-s]\n\tOR\n%s INPUTFOLDER OUTPUFFILE\n",argv[0],argv[0]);
	printf("   -s options will sort the file (only for single file input) \n\n\n");
	printf("   -n options will output single file (will break in multiple otherwise) \n\n\n");
	printf("   -i options will only remove identical profiles, without rotating \n\n\n");

        exit(1);
    }

    init_complement_ascii();

    // parse parameters
    strcpy(inputfile,argv[1]);

    IDENTICAL_ONLY = 0;
    if (argc>=4 && (0==strcmp("-I",argv[3]) || 0==strcmp("-i",argv[3]))) {
     IDENTICAL_ONLY = 1;
    }
    if (argc>=5 && (0==strcmp("-I",argv[4]) || 0==strcmp("-i",argv[4]))) {
     IDENTICAL_ONLY = 1;
    }
    if (argc>=6 && (0==strcmp("-I",argv[5]) || 0==strcmp("-i",argv[5]))) {
     IDENTICAL_ONLY = 1;
    }


    SINGLE_OUTFILE = 0;
    if (argc>=4 && (0==strcmp("-N",argv[3]) || 0==strcmp("-n",argv[3]))) {
     SINGLE_OUTFILE = 1;
    }
    if (argc>=5 && (0==strcmp("-N",argv[4]) || 0==strcmp("-n",argv[4]))) {
     SINGLE_OUTFILE = 1;
    }
    if (argc>=6 && (0==strcmp("-N",argv[5]) || 0==strcmp("-n",argv[5]))) {
     SINGLE_OUTFILE = 1;
    }

    SORT_ONLY = 0;
    if (argc>=4 && (0==strcmp("-S",argv[3]) || 0==strcmp("-s",argv[3]))) {
     SORT_ONLY = 1;
     SINGLE_OUTFILE = 1;
    }
    if (argc>=5 && (0==strcmp("-S",argv[4]) || 0==strcmp("-s",argv[4]))) {
     SORT_ONLY = 1;
     SINGLE_OUTFILE = 1;
    }
    if (argc>=6 && (0==strcmp("-S",argv[5]) || 0==strcmp("-s",argv[5]))) {
     SORT_ONLY = 1;
     SINGLE_OUTFILE = 1;
    }

    if (SINGLE_OUTFILE) {
      sprintf(outputfile,"%s",argv[2]);
      sprintf(outputfile2,"%s.rotindex",argv[2]);
    } else {
      sprintf(outputfile,"%s.1",argv[2]);
      sprintf(outputfile2,"%s.rotindex.1",argv[2]);
    }

	

    // create a file array
    FARRAY = EasyArrayCreate(1000, NULL, NULL);


    // list file(s)
    d=opendir(inputfile);
    if(d != NULL) {

      while(de = readdir(d)) {
	
       if (strlen(de->d_name)>17 && 0==strcasecmp(".leb36.renumbered",de->d_name+strlen(de->d_name)-17)) {
        fiptr = smalloc(sizeof(FITEM_STRUCT));

        bigtempbuf[0]=0;
	strcat(bigtempbuf,inputfile);
	if ('/' != bigtempbuf[strlen(bigtempbuf)-1]) strcat(bigtempbuf,"/");
	strcat(bigtempbuf,de->d_name);
       
        fiptr->inputfile = strdup(bigtempbuf);
        fiptr->d_name = strdup(de->d_name);
        EasyArrayInsert(FARRAY,fiptr);
       }
      }
      closedir(d);

    } else {

     if (ftouch(inputfile)) {

       // -s options to sort single file
       if (SORT_ONLY) {

	     EASY_LIST *profileList;
	     EASY_NODE *nof1;
	     FILE *fpi;

	
      	     fpi = fopen(inputfile,"r");
	     if(NULL == fpi) {
	        printf("\nERROR: Unable to open input file '%s'\n\n",inputfile);
	        exit(1);
	     }
	     
		             
	     profileList = EasyListCreate(NULL,NULL);
	     while(1){

		piptr=smalloc(sizeof(PITEM_STRUCT));

	        piptr->prof = ReadProfileWithRC(fpi, &(piptr->profrc));
		piptr->minRepresentation = NULL;
        	if (NULL==piptr->prof || NULL==piptr->profrc) {
		 free(piptr);
		 break;
		}
			

	        // find minimum representation and put in list
		piptr->minRepresentation = MinimumRepresentation(piptr->prof->indices,piptr->prof->proflen,piptr->profrc->indices,piptr->profrc->proflen,&(piptr->minrlen),IDENTICAL_ONLY); 
		if (NULL == piptr->minRepresentation) {
	          printf("\nERROR: minrepresentation is NULL!\n\n");
	          exit(1);			
		}
	        EasyListInsertTail(profileList,piptr);

	    }

	    // open the output file for writing
	    fpto = fopen(outputfile,"w");
	    if(fpto==NULL)
	    {
	        printf("\nERROR: Unable to open output file '%s'\n\n",outputfile);
	        exit(1);
	    }


	    // sort and output
	    fclose(fpi);
	    EasyListQuickSort(profileList,arsize_and_min_rep_cmp_pitem);
	    for (nof1 = profileList->head; nof1!=NULL; nof1=nof1->next) {
	     piptr = (PITEM_STRUCT*)EasyListItem(nof1);
	     WriteProfileWithRC(fpto,piptr->prof,piptr->profrc);
            }
	    fclose(fpto);
	    printf("\n\nredund.exe: Input File sorted by minimum representation. Now call this program again without -s switch to remove redundancy. \n\n"); fflush(stdout);

	    exit(0);
       }

       fiptr = smalloc(sizeof(FITEM_STRUCT));
       fiptr->inputfile = strdup(inputfile);
       fiptr->d_name = strdup(inputfile); // this won't be used in this case
       EasyArrayInsert(FARRAY,fiptr);

     } else {
        printf("\nERROR: Unable to read directory or file '%s'\n\n",inputfile);
        exit(1);
     }
    }


    // first order by file name, because readdir is not sorted by default
    EasyArrayQuickSort(FARRAY,name_cmp);

    // open input file(s) for reading
    printf("Opening files...\n");
    if(0==FARRAY->size)
    {
        printf("\nERROR: No input files found in the directory.\n\n");
        exit(1);
    }
    for (i=0; i<FARRAY->size; i++) {
      fiptr = (FITEM_STRUCT*)EasyArrayItem(FARRAY,i);
      fiptr->in = fopen(fiptr->inputfile,"r");
      if(NULL == fiptr->in) {
        printf("\nERROR: Unable to open input file '%s'.\n\n",fiptr->inputfile);
        exit(1);
      }
      printf("\t%s\n",fiptr->inputfile);
    }


    // open the output file for writing
    fpto = fopen(outputfile,"w");
    if(fpto==NULL)
    {
        printf("\nERROR: Unable to open output file '%s'.\n\n",outputfile);
        exit(1);
    }


    // open the index file for writing
    fpto2 = fopen(outputfile2,"w");
    if(fpto2==NULL)
    {
        printf("\nERROR: Unable to open index file '%s'\n\n",outputfile2);
        exit(1);
    }


    // read a profile for each one
    for (i=0; i<FARRAY->size; i++) {
      fiptr = (FITEM_STRUCT*)EasyArrayItem(FARRAY,i);
      fiptr->prof = ReadProfileWithRC(fiptr->in, &(fiptr->profrc));
      if (NULL!=fiptr->prof) {
       fiptr->minRepresentation = MinimumRepresentation(fiptr->prof->indices,fiptr->prof->proflen,fiptr->profrc->indices,fiptr->profrc->proflen,&(fiptr->minrlen),IDENTICAL_ONLY); 
       if (NULL == fiptr->minRepresentation) {
          printf("\nERROR: minrepresentation is NULL!\n\n");
          exit(1);			
       }

       fiptr->dir = RC;
      } else {
       fiptr->minRepresentation = NULL;
      }
    }


    // this is to make this deterministic (ties in sorting)

    // order list by arraysize, minimum representations (all NULLS go to end)
    startTime = time(NULL);
    EasyArrayQuickSort(FARRAY,arsize_and_min_rep_cmp);

// HERE

    // keep writing out, loading and reordering the file list
    nread=0; nwritten=0;
    lastwrite = NULL;
    while (1) {
     fiptr = (FITEM_STRUCT*)EasyArrayItem(FARRAY,0);

     // we are done
     if (NULL == fiptr->prof) break;

     // by now repeats are sorted by NOT NULL, arlen and mininimum representation of the pattern which is stored in fiptr->minRepresentation
     nread++;

     if (NULL!=lastwrite && 0 == arsize_and_min_rep_cmp(lastwrite,fiptr) && 

    // (ADDED apr 18, 2013)
   (    (0==pindcmp(fiptr->prof->indices,lastwrite->prof->indices,fiptr->prof->proflen,lastwrite->prof->proflen) 
     && 0==pindcmp(fiptr->profrc->indices,lastwrite->profrc->indices,fiptr->profrc->proflen,lastwrite->profrc->proflen))
    ||  (0==pindcmp(fiptr->prof->indices,lastwrite->profrc->indices,fiptr->prof->proflen,lastwrite->profrc->proflen)
     && 0==pindcmp(fiptr->profrc->indices,lastwrite->prof->indices,fiptr->profrc->proflen,lastwrite->prof->proflen))
    ) 
         ) {

        fprintf(fpto2," %d%c",fiptr->prof->key, (0==fiptr->dir) ? '\'' : '\"');

	// in proclu 1.87 we need all profiles written
	WriteProfileWithRC(fpto,fiptr->prof,fiptr->profrc);

     } else {

       nwritten++;

       // write out 
       WriteProfileWithRC(fpto,fiptr->prof,fiptr->profrc);

       if (NULL != lastwrite) { fprintf(fpto2,"\n"); }
       fprintf(fpto2,"%d%c",fiptr->prof->key, (0==fiptr->dir) ? '\'' : '\"');  // add preserved index entry

       if (NULL != lastwrite) {
        FreeProfile(lastwrite->prof);
        FreeProfile(lastwrite->profrc);
        free(lastwrite->minRepresentation);
        free(lastwrite);    
       }

       lastwrite = smalloc(sizeof(FITEM_STRUCT));
       lastwrite->inputfile = fiptr->inputfile;
       lastwrite->outputfile = fiptr->outputfile;
       lastwrite->outputfile2 = fiptr->outputfile2;
       lastwrite->in = fiptr->in;
       lastwrite->prof = CopyProfile(fiptr->prof);
       lastwrite->profrc = CopyProfile(fiptr->profrc);
       lastwrite->minRepresentation = pintdup(fiptr->minRepresentation,fiptr->minrlen);
       lastwrite->minrlen = fiptr->minrlen;

       // open new file?
       if (!SINGLE_OUTFILE && (nwritten%(RECORDS_PER_FILE))==0) {
	 fclose(fpto);
	 fclose(fpto2);
 
	 filescreated++;

	 // open the output file for writing
	 sprintf(bigtempbuf,"%s.%d",argv[2],filescreated);
	 fpto = fopen(bigtempbuf,"w");
	 if(fpto==NULL)
	 {
	     printf("\nERROR: Unable to open output file '%s'\n\n",bigtempbuf);
	     exit(1);
	 }

	 // open the index file for writing
	 sprintf(bigtempbuf,"%s.rotindex.%d",argv[2],filescreated);
	 fpto2 = fopen(bigtempbuf,"w");
	 if(fpto2==NULL)
	 {
	     printf("\nERROR: Unable to open output file '%s'\n\n",bigtempbuf);
	     exit(1);
	 }
       }

     }

     // load new profile
     fiptr2 = smalloc(sizeof(FITEM_STRUCT)); 
     fiptr2->inputfile = fiptr->inputfile;
     fiptr2->outputfile = fiptr->outputfile; 
     fiptr2->outputfile2 = fiptr->outputfile2;
     fiptr2->in = fiptr->in;
     fiptr2->prof = ReadProfileWithRC(fiptr2->in, &(fiptr2->profrc));
     if (NULL!=fiptr2->prof) {
       fiptr2->minRepresentation = MinimumRepresentation(fiptr2->prof->indices,fiptr2->prof->proflen,fiptr2->profrc->indices,fiptr2->profrc->proflen,&(fiptr2->minrlen),IDENTICAL_ONLY); 
       fiptr2->dir = RC;
     } else {
       fiptr2->minRepresentation = NULL;
     }

     // make sure new minrep is equal or larger than old minrep
     if (arsize_and_min_rep_cmp(fiptr2,fiptr)<0) {
        printf("\nERROR: file '%s' is not sorted by minimum representation!\n\n",fiptr2->inputfile);
        
 	printf("old: \n");
        WriteProfileWithRC(stdout,fiptr->prof,fiptr->profrc);

 	printf("prof: ");  
	for(i=0;i<fiptr->prof->proflen;i++) {
	  	printf(" %d",fiptr->prof->indices[i]);           
	}
        printf("  minrep: ");
        for(i=0;i<fiptr->minrlen;i++) {
                printf(" %d",fiptr->minRepresentation[i]);

        }
 	printf("\n\n");  


 	printf("new: \n");  
        WriteProfileWithRC(stdout,fiptr2->prof,fiptr2->profrc);

 	printf("prof: ");  
	for(i=0;i<fiptr2->prof->proflen;i++) {
	  	printf(" %d",fiptr2->prof->indices[i]);                    
	}
 	printf("  minrep: ");  
	for(i=0;i<fiptr2->minrlen;i++) {
	  	printf(" %d",fiptr2->minRepresentation[i]);  
         
	}
 	printf("\n\n");  

        exit(1);
     }

     // free old one and set new
     FreeProfile(fiptr->prof);
     FreeProfile(fiptr->profrc);
     free(fiptr->minRepresentation);
     free(fiptr);
     FARRAY->array[0] = fiptr2;

     // reorder
     EasyArrayQuickSort(FARRAY,arsize_and_min_rep_cmp);

    }

  
    printf("\n\n%llu profiles read, %llu profiles marked nonredundant. (time: %d seconds)\n\n",nread,nwritten,time(NULL) - startTime); fflush(stdout);


    return 0;
}


/*******************************************************************************************/
void doCriticalErrorAndQuit(const char *format, ... ) { 
	
	va_list				argp;
	

	printf("\n\nERROR!!!!:\n");
	if (format==NULL) exit(0);;
	
	va_start(argp, format);
	vprintf(format, argp); 
	fflush(stdout);
	va_end(argp);

    exit(1);
}


