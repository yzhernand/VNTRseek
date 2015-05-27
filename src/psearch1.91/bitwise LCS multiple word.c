/*
 *  bitwise LCS multiple word.c
 *  linear time bitwise LCS
 *
 *  Created by Gary Benson on 8/15/11.
 *  Revised by Gary Benson on 7/17/12
 *  Copyright 2011 Boston University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "bitwise LCS multiple word.h"
#include "bitwise LCS single word.h"
#include "convert to bitstring64.h"

#define wordSize 64

#define freeMultipleWords \
	free(matchA);\
	free(matchC);\
	free(matchG);\
	free(matchT);\
	free(matchN);\
	free(Complement);

int LCS_multiple_word(char *string1, char *string2, int n, int m){
	//the procedure computes fast bitwise LCS
	//n is the length of string1
	//m is the length of string2
	//n and m are *not* determined by strlen because they 
	//may be part of longer strings
	//returns -1 if error
	
	int i,j;
	unsigned long long int *matchA;
	unsigned long long int *matchC;
	unsigned long long int *matchG;
	unsigned long long int *matchT;
	unsigned long long int *matchN;
	unsigned long long int *matchVector; //To hold one of matchA, matchC, etc.  no memory allocated
	unsigned long long int *Complement;
	unsigned long long int bitmask;
	unsigned long long int complement;
	unsigned long long int matchString;
	unsigned long long int onlyOnesNotInOriginal;
	unsigned long long int addResult;
	unsigned long long int xorResult;
	unsigned long long int complementErase;
	unsigned long long int carryBit;
	int countOneBits;
	int integerPart;
	int N,M;
	char *stringN, *stringM;
	int NWords;
	
	//UINT_MAX is             4294967295 which is 32 bits long
	//ULLONG_MAX is 18446744073709551615 which is 64 bits long
	//printf("\n%u %llu",UINT_MAX,ULLONG_MAX);
	//printf("\n");
	
	
	//length of strings for LCS is n and m
	//number of wordSize-bit words for longer string is nWords
	//printf("\n n = %d, m = %d",n,m);
	//printf("\nnWords = %d, mWords = %d",nWords,mWords);
	
	//choose the longer string for the horizontal calculation: complement, matchString, etc.
	if(n>=m){
		N=n;
		M=m;
		stringN = string1;
		stringM = string2;
	}
	else{
		N=m;
		M=n;
		stringN = string2;
		stringM = string1;
	}
	
	//compute number of wordSize-bit words required.  All computation is done in wordSize-bit words.
	if(N>wordSize-1){
		integerPart = (N-wordSize+1)/wordSize;
		NWords = (N-wordSize+1-wordSize*integerPart)>0?integerPart+2:integerPart+1;
	}
	else return(LCS_single_word(stringN, stringM, N, M));
	
	//debug
	/*
	printf("\nN  M  NWords");
	printf("\n%d  %d  %d",N,M,NWords);
	*/
	
	//*************************encode match strings A C G T N for string1 
	matchA = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchC = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchG = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchT = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchN = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	Complement = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));	
	
	//loop through stringN and store bits in matchA, matchC, etc.
	//column zero in the score matrix is ignored
	//so we start with i = 0 and bitmask = 1
	
	for(j=0;j<NWords;j++)
	{
		bitmask=0x0000000000000001;
		
		matchA[j]=0x0000000000000000;
		matchC[j]=0x0000000000000000;
		matchG[j]=0x0000000000000000;
		matchT[j]=0x0000000000000000;
		matchN[j]=0x0000000000000000;	
		for(i=j*wordSize;i<(j+1)*wordSize;i++)
		{
			if (i<=N) 
			{
				if(i||j)
				{
					//if both i and j are zero, it means we are at the zero position in the row
					//don't process because it doesn't correspond to a character in the string
					switch (stringN[i-1]) {
						case 'A':
							matchA[j] |= bitmask;
							break;
						case 'C':
							matchC[j] |= bitmask;
							break;
						case 'G':
							matchG[j] |= bitmask;
							break;
						case 'T':
							matchT[j] |= bitmask;
							break;
						case 'N':
							matchN[j] |= bitmask;
							break;
						default:
							printf("\nError, non-ACGTN character read at position %d in string:%c",i,stringN[i]);
							freeMultipleWords;
							return(-1);
							break;
					}
				}
				bitmask <<=1; //bitmask = bitmask<<1
			}
		}
	}		

	//debug
	/*
	printf("\nmatchA:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchA[j]));}
	printf("\nmatchC:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchC[j]));}
	printf("\nmatchG:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchG[j]));}
	printf("\nmatchT:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchT[j]));}
	printf("\nmatchN:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(matchN[j]));}
	*/
	//**********************load complement of row zero in the score matrix which is all 1s
	
	for(j=0;j<NWords;j++)
	{
		Complement[j] = ~0x0000000000000000;
	}	
	//debug
	/*
	printf("\n\nRow 0");
	printf("\ncomplement:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(Complement[j]));}
	*/
	
	//loop for each letter in string2
	//row zero in the score matrix corresponds to the initial value of complement
	for(i=0;i<M;i++)
	{
		carryBit = 0x0000000000000000;
		
		switch (stringM[i]) {
			case 'A':
				matchVector = matchA;
				break;
			case 'C':
				matchVector = matchC;
				break;
			case 'G':
				matchVector = matchG;
				break;
			case 'T':
				matchVector = matchT;
				break;
			case 'N':
				matchVector = matchN;
				break;
			default:
				printf("\nError, non-ACGTN character read in string:%c",stringM[i]);
				freeMultipleWords;
				return(-1);
				break;
		}
		
		for(j=0;j<NWords;j++)
		{
			matchString = matchVector[j];
			
			complement = Complement[j];
			
			//AND complement and matchString
			onlyOnesNotInOriginal = complement&matchString;
			
			//ADD complement and onlyOnesNotInOriginal 
			//(and carryBit for strings longer than wordSize bits)
			//test if add sets carryBit
			//saw this online. if a+b causes overflow, a+b<a, a+b<b.  Is it correct?
			if(carryBit)
			{
				addResult = complement+carryBit;
				carryBit = (addResult==0)?1:0;  //if the carryBit is one, then only carry if whole number rolls over to zero
				addResult += onlyOnesNotInOriginal;
				if(addResult<onlyOnesNotInOriginal) carryBit = 1;  //must test again in case carryBit not reset to one on first add				
			}
			else //carryBit not set from last NWord
			{
				addResult = complement+onlyOnesNotInOriginal;
				carryBit = (addResult<onlyOnesNotInOriginal)?1:0;
			}
							
			//XOR complement and onlyOnesNotInOriginal
			xorResult = complement^onlyOnesNotInOriginal;
			
			//OR 
			Complement[j] = addResult|xorResult;
			
			//debug
			/*
			printf("\n\nRow %d, Column %d",i+1,j);
			printf("\n\n                     %s",string1);
			printf("\n                    %2d %c",i+1,string2[i]);
			printf("\nmatchString:           %s",convertToBitString64(matchString));
			printf("\nonlyOnesNotInOriginal: %s",convertToBitString64(onlyOnesNotInOriginal));
			printf("\naddResult:             %s",convertToBitString64(addResult));
			printf("\nxorResult:             %s",convertToBitString64(xorResult));
			printf("\ncomplement:            %s",convertToBitString64(Complement[j]));
			*/
			
		}
	}
		
	//debug
	/*
	printf("\n\nFinal:");
	printf("\ncomplement:         ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(Complement[j]));}
	printf("\n~complement:        ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(~Complement[j]));}
	*/
	 
	//find length of LCS
	countOneBits = 0 ;
	for(j=0;j<NWords;j++)
	{
		Complement[j] = ~Complement[j];
	
		//time proportional to number of ones in bit string
		complementErase = Complement[j];
		while (complementErase)//stop when zero
		{
			countOneBits++;
			complementErase &= (complementErase - 1); //removes last one bit
		}
	}
	
	//debug
	/*
	printf("\nNumber of one bits = LCS = %d",countOneBits);
   */
	
	
	freeMultipleWords;
	return(countOneBits);
	
	
}
