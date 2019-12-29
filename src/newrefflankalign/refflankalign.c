/****************************************************************************************************

   refflankalign.c - creates .MAP files (using narrowband alignment)

   Usage: ./refflankalign outdir xmldir inputfile MAXERRORS MAXTHREADS
PATLEN_SIZE_ERR REFLEN


****************************************************************************************************/

#include <stdio.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define max( a, b ) ( ( ( a ) >= ( b ) ) ? ( a ) : ( b ) )
#define min( a, b ) ( ( ( a ) <= ( b ) ) ? ( a ) : ( b ) )

#include "../libs/easylife/easylife.h"

#include "narrowbandDistanceAlignment.h"
//#include "aln.h"
//#include "aln_nb.h"

//#define PER_READ_STATS
//#define PRINT_ALIGNMENTS

#define MATCH_SCORE ( 1 )
#define MISM_PEN ( -3 )
#define GAP_OPEN_PEN ( -11 )
#define GAP_EXT_PEN ( -4 )

int MAXERRORS  = ( 0 );
int MAXTHREADS = ( 1 );
int REFLEN     = ( 1000000 );

#define MAXLINESIZE ( 10000000 )

int PATLEN_SIZE_ERR = 0;

char *Complementascii =
  NULL; /* precomputed array used for the reverse complement function */
int ThreadCounter = 0;

/********************************       init_complement_ascii
 * ****************************/
char *init_complement_ascii( void )

{
    /* complement has 256 entries so that finding the entries for A, C, G and T
     * which are complement ascii values */
    /* require no calculation */
    int   i;
    char *ca = (char *) calloc( 256, sizeof( int ) );

    if ( ca == NULL ) {
        return NULL;
    }

    for ( i = 0; i < 256; i++ )
        ca[i] = i;

    ca['A'] = 'T';
    ca['C'] = 'G';
    ca['G'] = 'C';
    ca['T'] = 'A';
    ca['a'] = 't';
    ca['c'] = 'g';
    ca['g'] = 'c';
    ca['t'] = 'a';

    return ca;
}

/*******************************************************************************************/
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

/***********************************************************************/
char *GetComplement( char *original ) {
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    /* safety check */
    if ( Complementascii == NULL ) {
        Complementascii = init_complement_ascii();
        if ( NULL == Complementascii )
            return NULL;
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

/***********************************************************************/
char *GetReverseComplement( char *original ) {
    static char *Complementascii =
      NULL; /* precomputed array used for the reverse complement function */
    char *sourceptr, *destinptr;
    int   length;
    char *buffer;

    /* safety check */
    if ( Complementascii == NULL ) {
        Complementascii = init_complement_ascii();
        if ( NULL == Complementascii )
            return NULL;
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

int RC = 0; // this will be set by MinimumRepresentation if the reverse
            // complement is the minimum

/*******************************************************************************************/
void MinimumRepresentation( char **indptr ) {

    int   i, j, len;
    char *src, *srcrc, *tar, *sourceptr, *destinptr, tmp, *indices;

    indices = *indptr;

    src = indices;
    if ( NULL == src )
        return;
    len = strlen( src );

    tar   = strdup( src );
    srcrc = GetReverseComplement( src );
    if ( NULL == tar || NULL == srcrc )
        return;

    RC = 0;

    // FORWARD - test each rotation and find the minimum one
    for ( i = 1; i < len; i++ ) {

        // printf("\n\nstart: %s",src);

        // rotate by 1
        tmp = src[0];
        for ( j = 1; j < len; j++ ) {
            src[j - 1] = src[j];
        }
        src[len - 1] = tmp;

        // printf("\nend: %s comp: %d\n\n",src,strcmp(src,tar));
        // exit(1);

        if ( strcmp( src, tar ) < 0 ) {
            free( tar );
            tar = strdup( src );
        }
    }

    // REVERSE - test each rotation and find the minimum one
    for ( i = 1; i < len; i++ ) {

        // printf("\n\nstart: %s",src);

        // rotate by 1
        tmp = srcrc[0];
        for ( j = 1; j < len; j++ ) {
            srcrc[j - 1] = srcrc[j];
        }
        srcrc[len - 1] = tmp;

        // printf("\nend: %s comp: %d\n\n",src,strcmp(src,tar));
        // exit(1);

        if ( strcmp( srcrc, tar ) < 0 ) {
            free( tar );
            tar = strdup( srcrc );
            RC  = 1;
        }
    }

    // save rotated
    free( src );
    free( srcrc );
    indices = tar;

    *indptr = indices;
}

/*******************************************************************************************/
void doCriticalErrorAndQuit( const char *format, ... ) {

    va_list argp;

    printf( "\n\nERROR!!!!:\n" );
    if ( format == NULL )
        exit( 0 );
    ;

    va_start( argp, format );
    vprintf( format, argp );
    fflush( stdout );
    va_end( argp );

    exit( 1 );
}

typedef struct {
    long long int id;
    int           patsize, leftlen, rightlen, count, scoresum, comparisons;
    char          leftcon, rightcon;
    char *        left;
    char *        right;
    char *        leftcp;
    char *        rightcp;
    char *        pattern;
    char *        sequence;
} FLANK;

/*******************************************************************************************/
void flankDestroy( void *in ) {

    FLANK *g = (FLANK *) in;

    free( g->left );
    free( g->right );
    free( g->leftcp );
    free( g->rightcp );
    free( g->pattern );
    free( g->sequence );
    free( g );
}

/*******************************************************************************************/
int __patsizeCmp( const void *item1, const void *item2 ) {

    int d1, d2;

    d1 = ( (FLANK *) item1 )->patsize;
    d2 = ( (FLANK *) item2 )->patsize;

    if ( d1 > d2 )
        return 1;
    else if ( d1 < d2 )
        return -1;

    return 0;
}

/*******************************************************************************************/
int main( int argc, char *argv[] ) {

    char *outdir, *xmldir, *filename, *buffer, *src1, *src2, *src3, *src4,
      *src5, *dst, tempbuf[2000];
    FILE *        rd_fp, *RFP;
    long long int rid, read_count;
    long long int processed = 0;
    char *        hasdata   = NULL;
    int           a1, a2, status, i, *PARR;
    int           distcount, undistcount, leftcount, rightcount, anycount;

    /* cmd arguments */
    if ( argc < 8 )
        doCriticalErrorAndQuit(
          "\n\nRefFlankAlign - Please use ./Refflankalign  outdir xmldir "
          "inputfile MAXERRORS MAXTHREADS PATLEN_SIZE_ERR REFLEN\n\n" );

    outdir   = argv[1];
    xmldir   = argv[2];
    filename = argv[3];

    if ( argc >= 5 ) {
        MAXERRORS = atoi( argv[4] );
    }

    if ( argc >= 6 ) {
        MAXTHREADS = atoi( argv[5] );
    }

    if ( MAXTHREADS < 1 ) {
        doCriticalErrorAndQuit( "\n\nRefFlankAlign - maxthreads should be more "
                                "than or equal to one!\n\n" );
    }

    if ( argc >= 7 ) {
        PATLEN_SIZE_ERR = atoi( argv[6] );
    }

    if ( PATLEN_SIZE_ERR < 1 ) {
        doCriticalErrorAndQuit(
          "\n\nRefFlankAlign - PATLEN_SIZE_ERR must be at least 1%%!\n\n" );
    }

    if ( argc >= 8 ) {
        REFLEN = atoi( argv[7] );
    }

    if ( REFLEN < 1 ) {
        doCriticalErrorAndQuit(
          "\n\nRefFlankAlign - REFLEN must be at least 1%%!\n\n" );
    }

    PARR = (int *) smalloc( sizeof( int ) * MAXTHREADS );

    /* alignment options */
    // if (NULL==(sm=init_sm(MATCH_SCORE,MISM_PEN)))
    //   doCriticalErrorAndQuit("Memory error. Aborting!\n\n");

    if ( Complementascii == NULL ) {
        Complementascii = init_complement_ascii();
        if ( NULL == Complementascii )
            return 1;
    }

    /* open representatives file */
    sprintf( tempbuf, "%s/representatives.txt", outdir );
    RFP = fopen( tempbuf, "w" );

    if ( !RFP )
        doCriticalErrorAndQuit( "\n\nRefFlankAlign - can't open output "
                                "representative file (%s). Aborting!\n\n",
          tempbuf );

    /* parse file one cluster at a time */
    buffer = smalloc( sizeof( char ) * ( MAXLINESIZE + 1 ) );

    if ( filename[0] == '-' )
        rd_fp = stdin;
    else
        rd_fp = fopen( filename, "r" );
    if ( !rd_fp )
        doCriticalErrorAndQuit(
          "\n\nRefFlankAlign - Could not open input file. Aborting!\n\n" );

    while ( ( hasdata = fgets( buffer, MAXLINESIZE, rd_fp ) ) != NULL ) {

        if ( buffer[0] == '@' ) {

            /* read cluster info */
            a1 = a2 = 0;
            sscanf( buffer, "@(%d_%d)", &a1, &a2 );
            if ( a1 <= 0 || a2 <= 0 ) {
                printf( "Can't read cluster number." );
                exit( 1 );
            }
            break;
        }
    };

    if ( hasdata )
        while ( 1 ) {

            EASY_LIST *ref_list;
            EASY_NODE *nd1, *nd2;
            int start, end, length, i, largestCount, scoresum, comparisons;
            long long int totlensum;
            FLANK *       templ_flank_ptr, *largestPtr;

            ref_list   = EasyListCreate( NULL, flankDestroy );
            read_count = 0;
            /* read cluster contents one line at a time */
            while ( 1 ) {

                rid = strtol( buffer, NULL, 10 );

                if ( rid < 0 ) {

                    // printf("rid: %ld\n",rid);

                    src1 = buffer;
                    while ( *src1 != ',' && *src1 != '\0' ) {
                        src1++;
                    }
                    if ( *src1 == ',' ) {
                        *src1 = '\0';
                        src1++;
                    }

                    src2 = src1;
                    while ( *src2 != ',' && *src2 != '\0' ) {
                        *src2 = toupper( *src2 );
                        if ( NULL == strchr( "ACGT", *src2 ) ) {
                            *src2 = 'N';
                        }
                        src2++;
                    }
                    if ( *src2 == ',' ) {
                        *src2 = '\0';
                        src2++;
                    }

                    src3 = src2;
                    while ( *src3 != ',' && *src3 != '\0' ) {
                        src3++;
                    }
                    if ( *src3 == ',' ) {
                        *src3 = '\0';
                        src3++;
                    }

                    src4 = src3;
                    while ( *src4 != ',' && *src4 != '\0' ) {
                        *src4 = toupper( *src4 );
                        if ( NULL == strchr( "ACGT", *src4 ) ) {
                            *src4 = 'N';
                        }
                        src4++;
                    }
                    if ( *src4 == ',' ) {
                        *src4 = '\0';
                        src4++;
                    }

                    src5 = src4 + strlen( src4 ) - 1;
                    while ( *src5 == '\n' || *src5 == '\r' ) {
                        *src5 = '\0';
                        src5--;
                    }

                    templ_flank_ptr = (FLANK *) smalloc( sizeof( FLANK ) );

                    templ_flank_ptr->left  = GetReverse( src1 );
                    templ_flank_ptr->right = strdup( src3 );

                    templ_flank_ptr->leftlen  = strlen( src1 );
                    templ_flank_ptr->rightlen = strlen( src3 );

                    templ_flank_ptr->leftlen =
                      min( templ_flank_ptr->leftlen, REFLEN );
                    templ_flank_ptr->rightlen =
                      min( templ_flank_ptr->rightlen, REFLEN );

                    templ_flank_ptr->leftcp =
                      GetComplement( templ_flank_ptr->left );
                    templ_flank_ptr->rightcp =
                      GetComplement( templ_flank_ptr->right );

                    templ_flank_ptr->sequence = strdup( src2 );
                    templ_flank_ptr->pattern  = strdup( src4 );

                    if ( !templ_flank_ptr->left || !templ_flank_ptr->right )
                        doCriticalErrorAndQuit(
                          "\n\nFlankAlign - memory error 1. Aborting!\n\n" );

                    templ_flank_ptr->patsize     = strlen( src4 );
                    templ_flank_ptr->id          = rid;
                    templ_flank_ptr->count       = 0;
                    templ_flank_ptr->scoresum    = 0;
                    templ_flank_ptr->comparisons = 0;
                    templ_flank_ptr->leftcon     = 0;
                    templ_flank_ptr->rightcon    = 0;

                    EasyListInsertTail( ref_list, templ_flank_ptr );

                } else if ( rid > 0 ) {

                    read_count++;

                    // printf("rid: %ld\n",rid);

                    /* SKIP

                    src1 = buffer;
                    while (*src1!=',' && *src1!='\0') src1++;
                    if (*src1==',') { *src1='\0'; src1++; }

                    src2 = src1;
                    while (*src2!=',' && *src2!='\0') src2++;
                    if (*src2==',') { *src2='\0'; src2++; }

                    src3 = src2;
                    while (*src3!=',' && *src3!='\0') { src3++; }
                    if (*src3==',') { *src3='\0'; src3++; }

                    src4 = src3;
                    while (*src4!=',' && *src4!='\0') { *src4=toupper(*src4); if
                    (NULL==strchr("ACGT",*src4)) { *src4='N'; } src4++; } if
                    (*src4==',') { *src4='\0'; src4++; }

                    src5 = src4 + strlen(src4) - 1;
                    while(*src5=='\n' || *src5=='\r') { *src5='\0'; src5--; }

                    templ_flank_ptr = (FLANK*)smalloc(sizeof(FLANK));

                    start = atoi(src1); end = atoi(src2); length = strlen(src3);

                    //printf("start: %d end: %d length: %d\n\n",start, end,
                    length);
                    //exit(1);

                    templ_flank_ptr->left = (char*)calloc(start+2,sizeof(char));
                    templ_flank_ptr->right =
                    (char*)calloc(length-end+2,sizeof(char));

                    if (!templ_flank_ptr->left || !templ_flank_ptr->right)
                       doCriticalErrorAndQuit("\n\nRefFlankAlign - memory
                    error 2. Aborting!\n\n");


                    dst = templ_flank_ptr->left;
                    for (i=start-2; i>=0; i--) { *dst=src3[i]; *dst++; }
                    *dst=0;

                    dst = templ_flank_ptr->right;
                    for (i=end; i<length; i++) { *dst=src3[i]; *dst++; }
                    *dst=0;

                    templ_flank_ptr->patsize = strlen(src4);
                    templ_flank_ptr->id = rid;

                    templ_flank_ptr->leftlen = strlen(templ_flank_ptr->left);
                    templ_flank_ptr->rightlen = strlen(templ_flank_ptr->right);

                    templ_flank_ptr->leftcp =
                    GetComplement(templ_flank_ptr->left);
                    templ_flank_ptr->rightcp =
                    GetComplement(templ_flank_ptr->right);


                    EasyListInsertTail(read_list,templ_flank_ptr);
                    */
                }

                if ( ( hasdata = fgets( buffer, MAXLINESIZE, rd_fp ) ) ==
                       NULL ||
                     buffer[0] == '@' )
                    break;
            }

            /* process cluster*/
            processed++;
            ThreadCounter++;

            // if (processed>=10) exit(0);

            pid_t pID = fork();

            if ( pID == 0 ) // child
            {
                // Code only executed by child process

                // BASICALIGNPAIR		*ba;
                FILE * fp;
                FLANK *read_ptr, *ref_ptr;
                int    comma;

                /* Redirect standard files to /dev/null */
                // freopen( "/dev/null", "r", stdin);
                // freopen( "/dev/null", "w", stdout);
                // freopen( "/dev/null", "w", stderr);

                EasyListQuickSort( ref_list, __patsizeCmp );

                sprintf( tempbuf, "%s/%d_%d.xml", xmldir, a1, a2 );
                fp = fopen( tempbuf, "w" );

                if ( !fp )
                    doCriticalErrorAndQuit(
                      "\n\nRefFlankAlign - can't open output XML file (%s). "
                      "Aborting!\n\n",
                      tempbuf );

                fprintf( fp, "<?xml-stylesheet type=\"text/xsl\" "
                             "href=\"../family.xsl\"?>\n" );
                fprintf( fp,
                  "<cluster trim=\"%d\" read_count=\"%lld\" ref_count=\"%zu\" "
                  "report=\"all\" maxedits=\"%d\" id=\"%d_%d\">\n",
                  REFLEN, read_count, ref_list->size, MAXERRORS, a1, a2 );

                if ( NULL == ref_list->head ) {
                    fprintf( fp,
                      "\n\nRefFlankAlign - empty reference list, this should "
                      "not happen(%d_%d). Aborting!\n\n",
                      a1, a2 );
                    _exit( 1 );
                }

                for ( nd1 = ref_list->head; nd1 != NULL; nd1 = nd1->next ) {

                    read_ptr = (FLANK *) EasyListItem( nd1 );

                    MinimumRepresentation( &( read_ptr->pattern ) );

                    // fprintf(fp,"%d - %d - %s - %s\n",read_ptr->leftlen,
                    // read_ptr->rightlen,read_ptr->left, read_ptr->right);
                    // fflush(fp);

                    comma = 0;
                    for ( nd2 = nd1->next; nd2 != NULL; nd2 = nd2->next ) {

                        int lerr, rerr;
                        int lerr1, rerr1;
                        int lerr2, rerr2;

                        int errors, errors1, errors2;
                        int sumerr1, sumerr2;

                        double sizeerror;

                        ref_ptr = (FLANK *) EasyListItem( nd2 );

                        // fprintf(fp,"\t%d - %d - %s - %s\n",ref_ptr->leftlen,
                        // ref_ptr->rightlen,ref_ptr->left, ref_ptr->right);
                        // fflush(fp);

                        /* if sizes are too different, stop */
                        if ( ref_ptr->patsize > read_ptr->patsize ) {
                            sizeerror =
                              ( ref_ptr->patsize / (double) read_ptr->patsize -
                                1.0 ) *
                              100;
                        } else {
                            sizeerror =
                              ( read_ptr->patsize / (double) ref_ptr->patsize -
                                1.0 ) *
                              100;
                        }

                        if ( sizeerror > PATLEN_SIZE_ERR ) {
                            break;
                        }

                        /* same direction */
                        lerr1 =
                          narrowbandUnitCostDistanceSameLengthBestScoreLastRowColumnLowMemUsesPointers(
                            read_ptr->left, ref_ptr->left, read_ptr->leftlen,
                            ( ref_ptr->leftlen > read_ptr->leftlen )
                              ? ref_ptr->leftlen
                              : read_ptr->leftlen,
                            MAXERRORS );
                        rerr1 =
                          narrowbandUnitCostDistanceSameLengthBestScoreLastRowColumnLowMemUsesPointers(
                            read_ptr->right, ref_ptr->right, read_ptr->rightlen,
                            ( ref_ptr->rightlen > read_ptr->rightlen )
                              ? ref_ptr->rightlen
                              : read_ptr->rightlen,
                            MAXERRORS );

                        // lerr1 = rerr1 = 100;
                        if ( lerr1 < 0 || rerr1 < 0 ) {
                            fprintf( fp,
                              "\n\nRefFlankAlign - wrong sizes for narrowband "
                              "alignment, seq2 must be larger. Aborting!\n\n" );
                            _exit( 1 );
                        }

                        /* try aligning to the compliment of the other flank
                         * instead */
                        lerr2 =
                          narrowbandUnitCostDistanceSameLengthBestScoreLastRowColumnLowMemUsesPointers(
                            read_ptr->rightcp, ref_ptr->left,
                            read_ptr->rightlen,
                            ( ref_ptr->leftlen > read_ptr->rightlen )
                              ? ref_ptr->leftlen
                              : read_ptr->rightlen,
                            MAXERRORS );
                        rerr2 =
                          narrowbandUnitCostDistanceSameLengthBestScoreLastRowColumnLowMemUsesPointers(
                            read_ptr->leftcp, ref_ptr->right, read_ptr->leftlen,
                            ( ref_ptr->rightlen > read_ptr->leftlen )
                              ? ref_ptr->rightlen
                              : read_ptr->leftlen,
                            MAXERRORS );

                        // lerr2 = rerr2 = 100;
                        if ( lerr2 < 0 || rerr2 < 0 ) {
                            fprintf( fp,
                              "\n\nRefFlankAlign - wrong sizes for narrowband "
                              "alignment, seq2 must be larger. Aborting!\n\n" );
                            _exit( 1 );
                        }

                        errors1 = max( lerr1, rerr1 );
                        errors2 = max( lerr2, rerr2 );
                        errors  = min( errors1, errors2 );

                        sumerr1 = lerr1 + rerr1;
                        sumerr2 = lerr2 + rerr2;

                        read_ptr->scoresum += min( sumerr1, sumerr2 );
                        read_ptr->comparisons++;

                        if ( lerr1 <= MAXERRORS || lerr2 <= MAXERRORS ) {
                            ref_ptr->leftcon  = 1;
                            read_ptr->leftcon = 1;
                        }

                        if ( rerr1 <= MAXERRORS || rerr2 <= MAXERRORS ) {
                            ref_ptr->rightcon  = 1;
                            read_ptr->rightcon = 1;
                        }

                        if ( errors <= MAXERRORS ) {

                            read_ptr->count++;
                        }
                    }
                }

                /* print components */
                distcount = undistcount = leftcount = rightcount = anycount = 0;
                fprintf( fp, "  <classes>\n" );

                // counts
                for ( nd1 = ref_list->head; nd1 != NULL; nd1 = nd1->next ) {
                    ref_ptr = (FLANK *) EasyListItem( nd1 );
                    if ( 1 == ref_ptr->leftcon )
                        leftcount++;
                    if ( 1 == ref_ptr->rightcon )
                        rightcount++;
                    if ( 0 == ref_ptr->leftcon && 0 == ref_ptr->rightcon )
                        distcount++;
                    if ( 1 == ref_ptr->leftcon && 1 == ref_ptr->rightcon )
                        undistcount++;
                    if ( 1 == ref_ptr->leftcon || 1 == ref_ptr->rightcon )
                        anycount++;
                }

                fprintf( fp,
                  "    <class component_count=\"%d\" name=\"distinguishable\" "
                  "ref_count=\"%d\" density=\"None\">\n",
                  distcount ? 1 : 0, distcount );
                if ( distcount ) {
                    fprintf( fp, "      <nodes count=\"%d\">", distcount );
                    for ( i = 0, nd1 = ref_list->head; nd1 != NULL;
                          nd1 = nd1->next ) {
                        ref_ptr = (FLANK *) EasyListItem( nd1 );

                        if ( 0 == ref_ptr->leftcon && 0 == ref_ptr->rightcon ) {
                            if ( i >= 1 ) {
                                fprintf( fp, " " );
                            }
                            fprintf( fp, "%lld", ref_ptr->id );
                            i++;
                        }
                    }
                    fprintf( fp, "</nodes>\n" );
                }
                fprintf( fp, "    </class>\n" );

                fprintf( fp,
                  "    <class name=\"left-connected\" ref_count=\"%d\" "
                  "density=\"None\" component_count=\"%d\">\n",
                  leftcount, leftcount ? 1 : 0 );
                if ( leftcount ) {
                    fprintf( fp, "      <component density=\"None\">\n" );
                    fprintf( fp, "        <nodes count=\"%d\">", leftcount );
                    for ( i = 0, nd1 = ref_list->head; nd1 != NULL;
                          nd1 = nd1->next ) {
                        ref_ptr = (FLANK *) EasyListItem( nd1 );

                        if ( 1 == ref_ptr->leftcon ) {
                            if ( i >= 1 ) {
                                fprintf( fp, " " );
                            }
                            fprintf( fp, "%lld", ref_ptr->id );
                            i++;
                        }
                    }
                    fprintf( fp, "</nodes>\n" );
                    fprintf( fp, "      </component>\n" );
                }
                fprintf( fp, "    </class>\n" );

                fprintf( fp,
                  "    <class name=\"right-connected\" ref_count=\"%d\" "
                  "density=\"None\" component_count=\"%d\">\n",
                  rightcount, rightcount ? 1 : 0 );
                if ( rightcount ) {
                    fprintf( fp, "      <component density=\"None\">\n" );
                    fprintf( fp, "        <nodes count=\"%d\">", rightcount );
                    for ( i = 0, nd1 = ref_list->head; nd1 != NULL;
                          nd1 = nd1->next ) {
                        ref_ptr = (FLANK *) EasyListItem( nd1 );

                        if ( 1 == ref_ptr->rightcon ) {
                            if ( i >= 1 ) {
                                fprintf( fp, " " );
                            }
                            fprintf( fp, "%lld", ref_ptr->id );
                            i++;
                        }
                    }
                    fprintf( fp, "</nodes>\n" );
                    fprintf( fp, "      </component>\n" );
                }
                fprintf( fp, "    </class>\n" );

                fprintf( fp,
                  "    <class name=\"left-or-right connected\" "
                  "ref_count=\"%d\" density=\"None\" component_count=\"%d\">\n",
                  anycount, anycount ? 1 : 0 );
                if ( anycount ) {
                    fprintf( fp, "      <component density=\"None\">\n" );
                    fprintf( fp, "        <nodes count=\"%d\">", anycount );
                    for ( i = 0, nd1 = ref_list->head; nd1 != NULL;
                          nd1 = nd1->next ) {
                        ref_ptr = (FLANK *) EasyListItem( nd1 );

                        if ( 1 == ref_ptr->leftcon || 1 == ref_ptr->rightcon ) {
                            if ( i >= 1 ) {
                                fprintf( fp, " " );
                            }
                            fprintf( fp, "%lld", ref_ptr->id );
                            i++;
                        }
                    }
                    fprintf( fp, "</nodes>\n" );
                    fprintf( fp, "      </component>\n" );
                }
                fprintf( fp, "    </class>\n" );

                fprintf( fp,
                  "    <class name=\"undistinguishable\" ref_count=\"%d\" "
                  "density=\"None\" component_count=\"%d\">\n",
                  undistcount, undistcount ? 1 : 0 );
                if ( undistcount ) {
                    fprintf( fp, "      <component density=\"None\">\n" );
                    fprintf( fp, "        <nodes count=\"%d\">", undistcount );
                    for ( i = 0, nd1 = ref_list->head; nd1 != NULL;
                          nd1 = nd1->next ) {
                        ref_ptr = (FLANK *) EasyListItem( nd1 );

                        if ( 1 == ref_ptr->leftcon && 1 == ref_ptr->rightcon ) {
                            if ( i >= 1 ) {
                                fprintf( fp, " " );
                            }
                            fprintf( fp, "%lld", ref_ptr->id );
                            i++;
                        }
                    }
                    fprintf( fp, "</nodes>\n" );
                    fprintf( fp, "      </component>\n" );
                }
                fprintf( fp, "    </class>\n" );

                fprintf( fp, "  </classes>\n" );

                /* print the representative with highest count to temp file */
                fprintf( fp, "  <references>\n" );
                largestCount = scoresum = 0;
                totlensum               = 0;
                largestPtr              = ref_list->head
                               ? (FLANK *) EasyListItem( ref_list->head )
                               : NULL;
                for ( nd1 = ref_list->head; nd1 != NULL; nd1 = nd1->next ) {

                    ref_ptr = (FLANK *) EasyListItem( nd1 );
                    if ( ref_ptr->count > largestCount ) {
                        largestCount = ref_ptr->count;
                        largestPtr   = ref_ptr;
                    }

                    scoresum += ref_ptr->scoresum;
                    totlensum += ref_ptr->comparisons *
                                 ( ref_ptr->leftlen + ref_ptr->rightlen );

                    fprintf( fp, "    <reference>\n" );
                    fprintf(
                      fp, "      <pattern>%s</pattern>\n", ref_ptr->pattern );
                    fprintf(
                      fp, "      <array>%s</array>\n", ref_ptr->sequence );
                    fprintf( fp, "    </reference>\n" );
                }

                if ( NULL == ref_list->head || NULL == largestPtr ) {
                    fprintf( fp,
                      "\n\nRefFlankAlign - empty reference list, this should "
                      "not happen(%d_%d). Aborting!\n\n",
                      a1, a2 );
                    _exit( 1 );
                }

                /* print references */
                if ( largestPtr ) {
                    char  scorebf[30];
                    FILE *fprtemp;
                    sprintf( tempbuf, "%s/%d_%d.rep", outdir, a1, a2 );
                    fprtemp = fopen( tempbuf, "w" );
                    if ( !fprtemp ) {
                        fprintf( fp, "\n\nRefFlankAlign - error opening temp "
                                     "representative file. Aborting!" );
                        _exit( 1 );
                    }

                    if ( 0 == totlensum )
                        strcpy( scorebf, "None" );
                    else
                        sprintf( scorebf, "%.2lf", 1.0 - scoresum / totlensum );

                    fprintf( fprtemp, "%d_%d\t%s\t%s", a1, a2, scorebf,
                      largestPtr->pattern );
                    fclose( fprtemp );
                }

                fprintf( fp, "  </references>\n" );
                fprintf( fp, "</cluster>\n" );
                fclose( fp );

                _exit( 0 );

            } else if ( pID < 0 ) // failed to fork
            {
                doCriticalErrorAndQuit( "\n\nRefFlankAlign - failed to fork "
                                        "code=%d (%s). Aborting!\n\n",
                  errno, strerror( errno ) );
                exit( 1 );
            }

            PARR[ThreadCounter - 1] = pID;

            // Code below only executed by parent process
            printf( "processing: %d_%d (pid: %d) threadcounter: %d\n", a1, a2,
              pID, ThreadCounter );

            /* child has its own copy now, so it's ok to destroy for parent */
            EasyListDestroy( ref_list );

            /* if full batch scheled, wait for one to finishe and schedule
             * another */
            if ( ThreadCounter >= MAXTHREADS ) {

                pid_t retid = waitpid( PARR[0], &status, 0 );

                if ( -1 == retid )
                    doCriticalErrorAndQuit( "\n\nRefFlankAlign - wait function "
                                            "returned -1. Aborting!" );
                else
                    printf( "\tterminated: %d\n", retid );

                for ( i = 1; i < ThreadCounter; i++ ) {
                    PARR[i - 1] = PARR[i];
                }

                ThreadCounter--;
            }

            /* end of file? */
            if ( !hasdata )
                break;

            /* read cluster info */
            a1 = a2 = 0;
            sscanf( buffer, "@(%d_%d)", &a1, &a2 );
            if ( a1 <= 0 || a2 <= 0 ) {
                printf( "Can't read cluster number. \n%s\n", buffer );
                exit( 1 );
            }
        }

    /* wait for all child processes to finish */
    for ( i = 0; i < ThreadCounter; i++ ) {

        pid_t retid = waitpid( PARR[i], &status, 0 );

        if ( -1 == retid )
            doCriticalErrorAndQuit(
              "\n\nRefFlankAlign - wait function returned -1. Aborting!" );
        else
            printf( "\tterminated: %d\n", retid );
    }
    ThreadCounter = 0;

    /* now print from temp files to representative file */
    rewind( rd_fp );
    while ( ( hasdata = fgets( buffer, MAXLINESIZE, rd_fp ) ) != NULL ) {

        if ( buffer[0] == '@' ) {

            a1 = a2 = 0;
            sscanf( buffer, "@(%d_%d)", &a1, &a2 );
            if ( a1 > 0 && a2 > 0 ) {
                FILE *fprtemp;
                sprintf( tempbuf, "%s/%d_%d.rep", outdir, a1, a2 );
                fprtemp = fopen( tempbuf, "r" );
                if ( !fprtemp )
                    doCriticalErrorAndQuit(
                      "\n\nRefFlankAlign - error opening temp representative "
                      "file (%s). Aborting!",
                      tempbuf );
                if ( fgets( tempbuf, 1000, fprtemp ) )
                    fprintf( RFP, "%s\n", tempbuf );
                fclose( fprtemp );
            }
        }
    }

    /* done */
    fclose( rd_fp );
    fclose( RFP );
    sprintf( tempbuf, "rm %s/*.rep", outdir );
    system( tempbuf );

    printf( "\ndone!!! processed: %lld, ThreadCounter: %d\n", processed,
      ThreadCounter );

    return 0;
}
