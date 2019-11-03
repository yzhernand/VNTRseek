/*
 *  bitwise LCS.c
 *  linear time bitwise LCS
 *
 *  Created by Gary Benson on 8/15/11.
 *  Revised by Gary Benson on 7/17/12
 *  Copyright 2011 Boston University. All rights reserved.
 *
 */

#include "bitwise LCS single word.h"
#include "convert to bitstring64.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define wordSize 64

int LCS_single_word( char *string1, char *string2, int n, int m ) {
    // the procedure computes fast bitwise LCS
    // n is the length of string1
    // m is the length of string2
    // these must both be <= wordSize
    // n and m are not determined by strlen because they
    // may be part of longer strings
    // returns =1 if error

    int                    i;
    unsigned long long int matchA;
    unsigned long long int matchC;
    unsigned long long int matchG;
    unsigned long long int matchT;
    unsigned long long int matchN;
    unsigned long long int bitmask;
    unsigned long long int complement;
    unsigned long long int matchString;
    unsigned long long int onlyOnesNotInOriginal;
    unsigned long long int addResult;
    unsigned long long int xorResult;
    unsigned long long int complementErase;
    int                    countOneBits;
    int                    junkBits;
    unsigned long long int junkBitsMask;

    // UINT_MAX is             4294967295 which is 32 bits long
    // ULLONG_MAX is 18446744073709551615 which is 64 bits long
    // printf("\n%u %llu",UINT_MAX,ULLONG_MAX);
    // printf("\n");

    // length of strings for LCS is n and m
    if ( n > wordSize - 1 ) {
        printf( "\nError, string1 is longer than %d characters. n=%d",
          wordSize - 1, n );
        return ( -1 );
    }

    junkBits     = wordSize - 1 - n;
    junkBitsMask = 0xFFFFFFFFFFFFFFFF >> junkBits;
    // printf("\njunkBitsMask >>> %s",convertToBitString64(junkBitsMask));

    //*************************encode match strings A C G T N for string1
    // loop through string1 and store bits in matchA, matchC, etc.
    // position zero corresponds to column zero in the score matrix, i.e., no
    // character so we start with i = 0 and bitmask = 2
    bitmask = 0x0000000000000002; // 2
    matchA  = 0x0000000000000000;
    matchC  = 0x0000000000000000;
    matchG  = 0x0000000000000000;
    matchT  = 0x0000000000000000;
    matchN  = 0x0000000000000000;
    for ( i = 0; i < n; i++ ) {
        // printf("\nbitmask %6llu",bitmask);
        switch ( string1[i] ) {
        case 'A':
            matchA |= bitmask;
            break;
        case 'C':
            matchC |= bitmask;
            break;
        case 'G':
            matchG |= bitmask;
            break;
        case 'T':
            matchT |= bitmask;
            break;
        case 'N':
            matchN |= bitmask;
            break;
        default:
            printf(
              "\nError, non-ACGTN character read in string:%c", string1[i] );
            return ( -1 );
            break;
        }
        bitmask <<= 1; // bitmask = bitmask<<1
    }

    // debug
    // printf("\nmatchA: %s",convertToBitString64(matchA));
    // printf("\nmatchC: %s",convertToBitString64(matchC));
    // printf("\nmatchG: %s",convertToBitString64(matchG));
    // printf("\nmatchT: %s",convertToBitString64(matchT));
    // printf("\nmatchN: %s",convertToBitString64(matchN));
    // printf("\nhere");

    //**********************load complement of row zero in the score matrix
    //which is all 1s
    complement = ~0x0000000000000000;

    // debug
    // printf("\n\nRow 0");
    // printf("\n\ncomplement:            %s",convertToBitString64(complement));

    // loop for each letter in string2
    // row zero in the score matrix corresponds to the initial value of
    // complement
    for ( i = 0; i < m; i++ ) {
        switch ( string2[i] ) {
        case 'A':
            matchString = matchA;
            break;
        case 'C':
            matchString = matchC;
            break;
        case 'G':
            matchString = matchG;
            break;
        case 'T':
            matchString = matchT;
            break;
        case 'N':
            matchString = matchN;
            break;
        default:
            printf(
              "\nError, non-ACGTN character read in string:%c", string2[i] );
            return ( -1 );
            break;
        }

        // AND complement and matchString
        onlyOnesNotInOriginal = complement & matchString;

        // ADD complement and onlyOnesNotInOriginal
        addResult = complement + onlyOnesNotInOriginal;

        // XOR complement and onlyOnesNotInOriginal
        xorResult = complement ^ onlyOnesNotInOriginal;

        // OR
        complement = addResult | xorResult;

        // debug
        // printf("\n\n                     %s",string1);
        // printf("\n                    %2d %c",i+1,string2[i]);
        // printf("\nmatchString %s",convertToBitString64(matchString));
        // printf("\nonlyOnesNotInOriginal:
        // %s",convertToBitString64(onlyOnesNotInOriginal)); printf("\naddResult:
        // %s",convertToBitString64(addResult)); printf("\nxorResult:
        // %s",convertToBitString64(xorResult)); printf("\ncomplement:
        // %s",convertToBitString64(complement));
    }

    // debug
    // printf("\n\nFinal complement:      %s",convertToBitString64(complement));
    // printf("\nFinal ~complement:     %s",convertToBitString64(~complement));

    complement = ~complement;

    // find length of LCS
    // time proportional to number of ones in bit string
    complementErase = complement;

    // mask bits past N in final word using JunkBits
    complementErase &= junkBitsMask;

    countOneBits = 0;
    while ( complementErase ) // stop when zero
    {
        countOneBits++;
        complementErase &= ( complementErase - 1 ); // removes last one bit
    }
    // debug
    // printf("\n\nNumber of one bits = LCS = %d\n",countOneBits);

    return ( countOneBits );
}