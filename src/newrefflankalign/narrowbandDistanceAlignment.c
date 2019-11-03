/*
 *  narrowbandDistanceAlignment.c
 *  NarrowbandDistanceAlignment
 *
 *  Created by Gary Benson on 1/20/12.
 *  Copyright 2012 Boston University. All rights reserved.
 *
 */

//#define debug

#include "narrowbandDistanceAlignment.h"
#include <stdio.h>
#include <stdlib.h>

#define min3( a, b, c ) \
    ( ( a <= b ) ? ( ( a <= c ) ? a : c ) : ( ( b <= c ) ? b : c ) )

int narrowbandUnitCostDistanceNoEndPenaltySeq2(
  char *seq1, char *seq2, int len1, int len2, int maxerr ) {
    // seq 1 is on left side
    // seq 2 is across top
    // no penalty at right end of seq 2
    // means look for best score in last row
    // allocates full array (len1 + 1) * (len2 + 1)
    // len(seq 2) must be >= len(seq 1), returns -1 otherwise

    int **S;
    int   score;
    int   row, col;
    int   match;
    int   left, up, diagonal;
    int   from, to;

    if ( len1 > len2 ) {
        return ( -1 );
    }

    S = (int **) calloc( ( len1 + 1 ), sizeof( int * ) );
    for ( row = 0; row <= len1; row++ ) {
        S[row] = (int *) calloc( len2 + 1, sizeof( int ) );
    }

    // standard unit cost alignment
    // zero row
    row  = 0;
    from = row - maxerr;
    to   = row + maxerr;
    if ( from < 0 )
        from = 0;
    if ( to > len2 )
        to = len2;

    for ( col = from; col <= to; col++ ) {
        S[0][col] = col;
    }
    // add boundary condition
    if ( to + 1 <= len2 ) {
        S[0][to + 1] = 1000;
    }

    // every other row
    for ( row = 1; row <= len1; row++ ) {
        from = row - maxerr;
        to   = row + maxerr;
        if ( from < 0 )
            from = 0;
        if ( to > len2 )
            to = len2;
        // add boundary condition
        if ( from > 0 )
            S[row][from - 1] = 1000;
        // columns
        for ( col = from; col <= to; col++ ) {
            if ( col == 0 )
                S[row][0] = row;
            else {
                match       = ( seq1[row - 1] == seq2[col - 1] ) ? 0 : 1;
                up          = S[row - 1][col] + 1;
                left        = S[row][col - 1] + 1;
                diagonal    = S[row - 1][col - 1] + match;
                S[row][col] = min3( diagonal, left, up );
                // if ((row==3)&&(col==3)) {
                //	printf("\nm:%d, d:%d, l:%d, u:%d
                //S[%d][%d]:%d,seq1[%d]:%c,seq2[%d]:%c",
                //		   match,diagonal,left,up,row,col,S[row][col],row-1,seq1[row-1],col-1,seq2[col-1]);
            }
        }
        // add boundary condition
        if ( to + 1 <= len2 ) {
            S[row][to + 1] = 1000;
        }
    }

#ifdef debug
    // print to check
    printf( "\n    -" );
    for ( col = 0; col < len2; col++ ) {
        printf( "  %c", seq2[col] );
    }

    for ( row = 0; row <= len1; row++ ) {
        if ( row == 0 )
            printf( "\n - " );
        else
            printf( "\n %c ", seq1[row - 1] );
        from = row - maxerr;
        to   = row + maxerr;
        for ( col = 0; col <= len2; col++ ) {
            if ( ( col < from - 1 ) || ( col > to + 1 ) )
                printf( " | " );
            else {
                if ( S[row][col] == 1000 )
                    printf( " X " );
                else
                    printf( "%2d ", S[row][col] );
            }
        }
    }
#endif

    // find minimum score
    row  = len1;
    from = row - maxerr;
    to   = row + maxerr;
    if ( from < 0 )
        from = 0;
    if ( to > len2 )
        to = len2;
    score = S[row][from];
    for ( col = from + 1; col <= to; col++ ) {
        if ( S[row][col] < score ) {
            score = S[row][col];
        }
    }

    // free memory
    free( S );
    return ( score );
}

int narrowbandUnitCostDistanceNoEndPenaltySeq2LowMem(
  char *seq1, char *seq2, int len1, int len2, int maxerr ) {
    // seq 1 is on left side
    // seq 2 is across top
    // no penalty at right end of seq 2
    // means look for best score in last row
    // allocates narrowband
    // len(seq 2) must be >= len(seq 1), returns -1 otherwise

    int **S;
    int   score;
    int   row, col;
    int   match;
    int   left, up, diagonal;
    int   from, to;
    int   rowlength;
    int   leftboundary, rightboundary;
    int   previousleftboundary;

    if ( len1 > len2 ) {
        return ( -1 );
    }

    S = (int **) calloc( ( len1 + 1 ), sizeof( int * ) );
    for ( row = 0; row <= len1; row++ ) {
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;
        S[row]    = (int *) calloc(
          leftboundary + rowlength + rightboundary, sizeof( int ) );
    }

    // standard unit cost alignment
    // zero row
    row          = 0;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;

    if ( leftboundary ) { // boundary condition left
        S[row][0] = 1000;
    }
    for ( col = leftboundary; col <= leftboundary + rowlength - 1; col++ ) {
        S[row][col] = col; // zero row fill
    }
    if ( rightboundary ) { // boundary condition right
        S[row][leftboundary + rowlength + rightboundary - 1] = 1000;
    }

    // save for next row
    previousleftboundary = leftboundary;

    // every other row
    for ( row = 1; row <= len1; row++ ) {
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;

        if ( leftboundary ) { // boundary condition left
            S[row][0] = 1000;
        }
        for ( col = leftboundary; col <= leftboundary + rowlength - 1; col++ ) {
            if ( col == 0 )
                S[row][0] = row; // column zero fill
            else {
                match = ( seq1[row - 1] == seq2[from + col - leftboundary - 1] )
                          ? 0
                          : 1;
                up   = S[row - 1][col + previousleftboundary] + 1;
                left = S[row][col - 1] + 1;
                diagonal =
                  S[row - 1][col - ( 1 - previousleftboundary )] + match;
                S[row][col] = min3( diagonal, left, up );
            }
            // if ((row==1)&&(col==1)) {
            //	printf("\nrow:%d, col:%d,seqleft:%c,
            //seqtop:%c\n",row,col,seq1[row-1],seq2[from+col-leftboundary-1]);
            //}
        }
        if ( rightboundary ) { // boundary condition right
            S[row][leftboundary + rowlength + rightboundary - 1] = 1000;
        }
        previousleftboundary = leftboundary;
    }

#ifdef debug
    // print to check
    printf( "\n    -" );
    for ( col = 0; col < len2; col++ ) {
        printf( "  %c", seq2[col] );
    }

    for ( row = 0; row <= len1; row++ ) {
        if ( row == 0 )
            printf( "\n - " );
        else
            printf( "\n %c ", seq1[row - 1] );
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;

        for ( col = 0; col < from - 1; col++ )
            printf( " | " );
        for ( col = 0; col <= leftboundary + rowlength + rightboundary - 1;
              col++ ) {
            if ( S[row][col] == 1000 )
                printf( " X " );
            else
                printf( "%2d ", S[row][col] );
        }
    }
#endif

    // find minimum score
    row          = len1;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;

    score = S[row][0];
    for ( col = 1; col <= leftboundary + rowlength + rightboundary - 1;
          col++ ) {
        if ( S[row][col] < score ) {
            score = S[row][col];
        }
    }

    // free memory
    for ( row = 0; row <= len1; row++ ) {
        free( S[row] );
    }
    free( S );

    return ( score );
}

int narrowbandUnitCostDistanceNoEndPenaltySeq2LowMemUsesPointers(
  char *seq1, char *seq2, int len1, int len2, int maxerr ) {
    // seq 1 is on left side
    // seq 2 is across top
    // no penalty at right end of seq 2
    // means look for best score in last row
    // allocates narrowband
    // uses pointers for score array rather than indexes
    // len(seq 2) must be >= len(seq 1), returns -1 otherwise

    int **S;
    int   score;
    int   row, col;
    int   match;
    int   left, up, diagonal;
    int   from, to;
    int   rowlength;
    int   leftboundary, rightboundary;
    int   previousleftboundary;
    int * Sp, *Spu, *Spl, *Spd;

    if ( len1 > len2 ) {
        return ( -1 );
    }

    S = (int **) calloc( ( len1 + 1 ), sizeof( int * ) );
    for ( row = 0; row <= len1; row++ ) {
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;
        S[row]    = (int *) calloc(
          leftboundary + rowlength + rightboundary, sizeof( int ) );
    }

    // standard unit cost alignment
    // zero row
    row          = 0;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;

    Sp = S[row];
    if ( leftboundary ) { // boundary condition left
        ( *Sp ) = 1000;
    }
    Sp += leftboundary;
    for ( col = leftboundary; col <= leftboundary + rowlength - 1; col++ ) {
        ( *Sp ) = col; // zero row fill
        Sp++;
    }
    if ( rightboundary ) { // boundary condition right
        ( *Sp ) = 1000;
    }

    // save for next row
    previousleftboundary = leftboundary;

    // every other row
    for ( row = 1; row <= len1; row++ ) {
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;

        Sp  = S[row];
        Spu = S[row - 1] + previousleftboundary;
        Spl = Sp - 1;
        Spd = Spu - 1;
        if ( leftboundary ) { // boundary condition left
            ( *Sp ) = 1000;
        }
        Sp += leftboundary;
        Spu += leftboundary;
        Spl += leftboundary;
        Spd += leftboundary;
        for ( col = leftboundary; col <= leftboundary + rowlength - 1; col++ ) {
            if ( col == 0 )
                ( *Sp ) = row; // column zero fill
            else {
                match = ( seq1[row - 1] == seq2[from + col - leftboundary - 1] )
                          ? 0
                          : 1;
                up       = ( *Spu ) + 1;
                left     = ( *Spl ) + 1;
                diagonal = ( *Spd ) + match;
                ( *Sp )  = min3( diagonal, left, up );
            }
            // if ((row==1)&&(col==1)) {
            //	printf("\nrow:%d, col:%d,seqleft:%c,
            //seqtop:%c\n",row,col,seq1[row-1],seq2[from+col-leftboundary-1]);
            //}
            Sp++;
            Spu++;
            Spl++;
            Spd++;
        }
        if ( rightboundary ) { // boundary condition right
            ( *Sp ) = 1000;
        }
        previousleftboundary = leftboundary;
    }

#ifdef debug
    // print to check
    printf( "\n    -" );
    for ( col = 0; col < len2; col++ ) {
        printf( "  %c", seq2[col] );
    }

    for ( row = 0; row <= len1; row++ ) {
        if ( row == 0 )
            printf( "\n - " );
        else
            printf( "\n %c ", seq1[row - 1] );
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;

        for ( col = 0; col < from - 1; col++ )
            printf( " | " );
        for ( col = 0; col <= leftboundary + rowlength + rightboundary - 1;
              col++ ) {
            if ( S[row][col] == 1000 )
                printf( " X " );
            else
                printf( "%2d ", S[row][col] );
        }
    }
#endif

    // find minimum score
    row          = len1;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;

    Sp    = S[row];
    score = ( *Sp );
    Sp++;
    for ( col = 1; col <= leftboundary + rowlength + rightboundary - 1;
          col++ ) {
        if ( ( *Sp ) < score ) {
            score = ( *Sp );
        }
        Sp++;
    }

    // free memory
    for ( row = 0; row <= len1; row++ ) {
        free( S[row] );
    }
    free( S );

    return ( score );
}

int narrowbandUnitCostDistanceSameLengthBestScoreLastRowColumnLowMemUsesPointers(
  char *seq1, char *seq2, int len1, int len2, int maxerr ) {
    // seq 1 is on left side
    // seq 2 is across top
    // assumes sequences are same length
    // len(seq 2) must be = len(seq 1), returns -1 otherwise
    // looks for best score in last row and last column
    // allocates narrowband
    // uses pointers for score array rather than indexes

    int **S;
    int   score;
    int   row, col;
    int   match;
    int   left, up, diagonal;
    int   from, to;
    int   rowlength;
    int   leftboundary, rightboundary;
    int   previousleftboundary;
    int * Sp, *Spu, *Spl, *Spd;

    if ( len1 != len2 ) {
        return ( -1 );
    }

    S = (int **) calloc( ( len1 + 1 ), sizeof( int * ) );
    for ( row = 0; row <= len1; row++ ) {
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;
        S[row]    = (int *) calloc(
          leftboundary + rowlength + rightboundary, sizeof( int ) );
    }

    // standard unit cost alignment
    // zero row
    row          = 0;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;

    Sp = S[row];
    if ( leftboundary ) { // boundary condition left
        ( *Sp ) = 1000;
    }
    Sp += leftboundary;
    for ( col = leftboundary; col <= leftboundary + rowlength - 1; col++ ) {
        ( *Sp ) = col; // zero row fill
        Sp++;
    }
    if ( rightboundary ) { // boundary condition right
        ( *Sp ) = 1000;
    }

    // save for next row
    previousleftboundary = leftboundary;

    // every other row
    for ( row = 1; row <= len1; row++ ) {
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;

        Sp  = S[row];
        Spu = S[row - 1] + previousleftboundary;
        Spl = Sp - 1;
        Spd = Spu - 1;
        if ( leftboundary ) { // boundary condition left
            ( *Sp ) = 1000;
        }
        Sp += leftboundary;
        Spu += leftboundary;
        Spl += leftboundary;
        Spd += leftboundary;
        for ( col = leftboundary; col <= leftboundary + rowlength - 1; col++ ) {
            if ( col == 0 )
                ( *Sp ) = row; // column zero fill
            else {
                match = ( seq1[row - 1] == seq2[from + col - leftboundary - 1] )
                          ? 0
                          : 1;
                up       = ( *Spu ) + 1;
                left     = ( *Spl ) + 1;
                diagonal = ( *Spd ) + match;
                ( *Sp )  = min3( diagonal, left, up );
            }
            // if ((row==1)&&(col==1)) {
            //	printf("\nrow:%d, col:%d,seqleft:%c,
            //seqtop:%c\n",row,col,seq1[row-1],seq2[from+col-leftboundary-1]);
            //}
            Sp++;
            Spu++;
            Spl++;
            Spd++;
        }
        if ( rightboundary ) { // boundary condition right
            ( *Sp ) = 1000;
        }
        previousleftboundary = leftboundary;
    }

#ifdef debug
    // print to check
    printf( "\n    -" );
    for ( col = 0; col < len2; col++ ) {
        printf( "  %c", seq2[col] );
    }

    for ( row = 0; row <= len1; row++ ) {
        if ( row == 0 )
            printf( "\n - " );
        else
            printf( "\n %c ", seq1[row - 1] );
        from         = row - maxerr;
        to           = row + maxerr;
        leftboundary = 0; // false
        if ( from <= 0 )
            from = 0;
        else
            leftboundary = 1; // true
        rightboundary = 0;    // false
        if ( to >= len2 )
            to = len2;
        else
            rightboundary = 1; // true
        rowlength = to - from + 1;

        for ( col = 0; col < from - 1; col++ )
            printf( " | " );
        for ( col = 0; col <= leftboundary + rowlength + rightboundary - 1;
              col++ ) {
            if ( S[row][col] == 1000 )
                printf( " X " );
            else
                printf( "%2d ", S[row][col] );
        }
    }
#endif

    // find minimum score in last row
    row          = len1;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;

    Sp    = S[row];
    score = ( *Sp );
    Sp++;
    for ( col = 1; col <= leftboundary + rowlength + rightboundary - 1;
          col++ ) {
        if ( ( *Sp ) < score ) {
            score = ( *Sp );
        }
        Sp++;
    }

    // find minimum score in last column
    row          = len1;
    from         = row - maxerr;
    to           = row + maxerr;
    leftboundary = 0; // false
    if ( from <= 0 )
        from = 0;
    else
        leftboundary = 1; // true
    rightboundary = 0;    // false
    if ( to >= len2 )
        to = len2;
    else
        rightboundary = 1; // true
    rowlength = to - from + 1;
    // printf("\n");
    while ( ( row >= 0 ) && ( to == len2 ) ) {
        Sp = S[row] + leftboundary + rowlength + rightboundary - 1;
        // printf("%d ",(*Sp));
        if ( ( *Sp ) < score ) {
            score = ( *Sp );
        }
        row--;
        if ( row >= 0 ) {
            from         = row - maxerr;
            to           = row + maxerr;
            leftboundary = 0; // false
            if ( from <= 0 )
                from = 0;
            else
                leftboundary = 1; // true
            rightboundary = 0;    // false
            if ( to >= len2 )
                to = len2;
            else
                rightboundary = 1; // true
            rowlength = to - from + 1;
        }
    }
    // free memory
    for ( row = 0; row <= len1; row++ ) {
        free( S[row] );
    }
    free( S );

    // minimum score from last row and column
    return ( score );
}
