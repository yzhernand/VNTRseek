/*
 *  bitwise edit distance alignment multiple word no end penalty.c
 *  bitwise edit distance alignment
 *
 *  Created by Gary Benson on 7/20/12.
 *  Copyright 2012 Boston University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "bitwise edit distance alignment multiple word no end penalty.h"
#include "convert to bitstring64.h"

#define wordSize 64
#define wordSizeMinusOne 63

#define freeMultipleWords \
free(matchA);\
free(matchC);\
free(matchG);\
free(matchT);\
free(matchN);\
free(P_D);\
free(P_S_or_P_D);

int Edit_Distance_multiple_word_NoEndPenaltySeq1(char *stringN, char *stringM, int N, int M)
{
	
	//the procedure computes fast bitwise edit distance
	//***longer string is stringN
	//***no end penalty on the right end of stringN
	//n is the length of stringN
	//m is the length of stringM
	//n and m are *not* determined by strlen because they 
	//may be part of longer strings
	//does not call Edit_Distance_single_word even when N<wordsize
	
	//returns zero if M is zero
	
	//returns -1 if error
	
	
	int i,j;
	unsigned long long int *matchA;
	unsigned long long int *matchC;
	unsigned long long int *matchG;
	unsigned long long int *matchT;
	unsigned long long int *matchN;
	unsigned long long int *matchVector; //to hold one of matchA, matchC, etc.  no memory allocated
	unsigned long long int bitmask;
	unsigned long long int matchString;
	unsigned long long int *P_D;
	unsigned long long int *P_S_or_P_D;
	unsigned long long int M_m;
	unsigned long long int R_IS;
	unsigned long long int not_M_I;
	unsigned long long int not_M_I_xor_P_D;
	unsigned long long int VC_0;
	unsigned long long int VC_plus_1;
	unsigned long long int VC_0_shift;
	unsigned long long int VC_plus_1_shift;
	unsigned long long int sum;
	unsigned long long int sum_and_R_IS;
	unsigned long long int highBitMask64 = 0x8000000000000000;
	unsigned long long int carryBitSum;
	unsigned long long int carryBitVC_plus_1_shift;
	unsigned long long int carryBitVC_0_shift;
	unsigned long long int oldCarryBitVC_plus_1_shift;
	unsigned long long int oldCarryBitVC_0_shift;
	
	
	unsigned long long int P_DErase;
	unsigned long long int P_IErase;	
	unsigned long long int mask;
	int Score;
	int bestFinalRowScore;
	
	int nWords;
	int integerPart;
	int NWords;
	
	//***beginning of identical code to Edit_Distance_multiple_word
	
	//UINT_MAX is             4294967295 which is 32 bits long
	//ULLONG_MAX is 18446744073709551615 which is 64 bits long
	//printf("\n%u %llu",UINT_MAX,ULLONG_MAX);
	//printf("\n");
	
	//compute number of wordSize-bit words required.  
	//All computation is done in wordSize-bit words.
	//First bit word can only hold wordSize-1 string positions
	//so there is space for the zero column
	
	if(M==0) return(0);
	
	if(N>wordSize-1){
		integerPart = (N-wordSize+1)/wordSize;
		nWords = (N-wordSize+1-wordSize*integerPart)>0?integerPart+2:integerPart+1;
	}
	else nWords=1;  //does not call Edit_Distance_single_word
	
	
	//length of strings for edit distance is n and m
	//number of wordSize-bit words is nWords
	//debug
	NWords = nWords;
	
	//debug
	/*
	 printf("\nN  M  NWords");
	 printf("\n%d  %d  %d",N,M,NWords);
	 */
	
	//storage allocation
	matchA = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchC = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchG = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchT = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	matchN = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	P_D =  (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	P_S_or_P_D = (unsigned long long int *)calloc(NWords,sizeof(unsigned long long int));
	
	//*************************encode match strings A C G T N for string1 
	//loop through stringN and store bits in matchA, matchC, etc.
	//column zero in the score matrix is ignored
	//so we start with i = 0 and bitmask = 1 in the first nWord only
	
	
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
				bitmask <<=1; //bitmask = bitmask<<1, moves set bit one position higher
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
	
	//initialize PD and PS||PD for row zero
	//all of row zero is increase, except zero position
	//which must be S or D so that an initial first R_IS run is computed correctly
	//here it is set to D for global alignment so that the vertical class value
	//below it will be set to +1 (first column increasing)
	//otherwise, set to S for semi-global alignment so no penalty 
	//for starting row (first column all zeros)
	
	P_D[0] = 1;        //puts 1 in zero position, as required above for global alignment
	P_S_or_P_D[0] = 1; //required for consistency
	
	//debug
	/*
	 printf("\n\nRow 0");
	 printf("\nP_D:         ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_D[j]));}
	 printf("\nP_S_or_P_D:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_S_or_P_D[j]));}
	 */
	
	
	//loop for each letter in string2
	//row zero in the score matrix corresponds to the initial value of complement
	for(i=0;i<M;i++)
	{
		carryBitSum = 0x0000000000000000;
		carryBitVC_plus_1_shift = 0x0000000000000001; //set up for adding 1 to first position after shift
		carryBitVC_0_shift = 0x0000000000000000;
		
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
			oldCarryBitVC_plus_1_shift = carryBitVC_plus_1_shift;
			oldCarryBitVC_0_shift = carryBitVC_0_shift;
			
			matchString = matchVector[j];
			
			//compute bitstrings
			//match subsets
			M_m = matchString|P_D[j];
			R_IS = ~M_m;
			not_M_I = R_IS|P_S_or_P_D[j];
			not_M_I_xor_P_D = not_M_I^P_D[j];
			
			//the add (with carry)
			sum = not_M_I+P_S_or_P_D[j]+carryBitSum;
			//get carryBitSum for next round
			//saw this online. if a+b causes overflow (carry), a+b<a, a+b<b.  
			if(sum<not_M_I) carryBitSum = 1;  
			sum_and_R_IS = sum&R_IS;
			
			//the new vertical classes
			VC_0 = sum_and_R_IS^not_M_I_xor_P_D;
			VC_plus_1 = P_D[j]|(sum_and_R_IS&P_S_or_P_D[j]);
			
			//get shift carry bits
			//mask and shift 
			carryBitVC_0_shift = (VC_0&highBitMask64)>>wordSizeMinusOne;
			carryBitVC_plus_1_shift = (VC_plus_1&highBitMask64)>>wordSizeMinusOne;
			
			//the vertical class shifts (with carry)
			VC_0_shift=(VC_0<<1)+oldCarryBitVC_0_shift;
			VC_plus_1_shift=(VC_plus_1<<1)+oldCarryBitVC_plus_1_shift;
			
			//reset first position in VC_plus_1_shift for boundary condition
			//this is done with oldCarryBitVC__plus_1_shift above
			//VC_plus_1_shift+=1;  
			
			//debug
			/*
			 printf("\n\nRow %d, Column %d",i+1,j);
			 printf("\nP_D:                    %s",convertToBitString64(P_D[j]));
			 printf("\nP_S_or_P_D:             %s",convertToBitString64(P_S_or_P_D[j]));
			 printf("\nmatchString:            %s",convertToBitString64(matchString));
			 printf("\nM_m:                    %s",convertToBitString64(M_m));
			 printf("\nR_IS                    %s",convertToBitString64(R_IS));
			 printf("\nnot_M_I                 %s",convertToBitString64(not_M_I));
			 printf("\nnot_M_I_xor_P_D         %s",convertToBitString64(not_M_I_xor_P_D));
			 printf("\nsum                     %s",convertToBitString64(sum));
			 printf("\nsum_and_R_IS            %s",convertToBitString64(sum_and_R_IS));
			 printf("\nVC_0                    %s",convertToBitString64(VC_0));
			 printf("\nVC_0_shift              %s",convertToBitString64(VC_0_shift));
			 printf("\nVC_plus_1               %s",convertToBitString64(VC_plus_1));
			 printf("\nVC_plus_1_shift         %s",convertToBitString64(VC_plus_1_shift));
			 
			 printf("\n\ncarryBitSum             %s",convertToBitString64(carryBitSum));
			 printf("\ncarryBitVC_0_shift      %s",convertToBitString64(carryBitVC_0_shift));
			 printf("\ncarryBitVC_plus_1_shift %s",convertToBitString64(carryBitVC_plus_1_shift));
			 */
			
			//new pair classes
			P_D[j] = M_m&VC_plus_1_shift;
			P_S_or_P_D[j] = VC_plus_1_shift|(M_m&VC_0_shift);
			
			//debug
			/*
			 printf("\n\nNew P_D:                %s",convertToBitString64(P_D[j]));
			 printf("\nNew P_S_or_P_D:         %s",convertToBitString64(P_S_or_P_D[j]));
			 */
		}
	}
	
	//debug
	/*
	 printf("\n\nFinal:");
	 printf("\nP_D:         ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_D[j]));}
	 printf("\nP_S_or_P_D:  ");for(j=0;j<NWords;j++){printf("%s ",convertToBitString64(P_S_or_P_D[j]));}
	 */
	
	//find best edit distance score in final row
	mask=0x0000000000000001;
	bestFinalRowScore=M;
	Score=M;
	for(j=0;j<NWords;j++)
	{	
		P_DErase=P_D[j];
		P_IErase=~(P_S_or_P_D[j]);
		for (i=0;i<wordSize;i++) {
			//do nothing if first NWord, first bit.  
			//This was set to 1 artifically in both PD and PSorPD
			if(i|j) 
			{
				Score-=(int)P_DErase&mask;
				Score+=(int)P_IErase&mask;
				if (Score<bestFinalRowScore) bestFinalRowScore=Score;
			}
			P_DErase>>=1;
			P_IErase>>=1;
			
			//debug
			/*
			 printf("\n\nP_DErase: %s",convertToBitString64(P_DErase));
			 printf("\nP_IErase: %s",convertToBitString64(P_IErase));
			 */
		}
		
	}
	
	//***end of identical code to Edit_Distance_multiple_word
	
	//debug
	
	/*
	printf("\nBest Final Row Score is %d",bestFinalRowScore);
	*/ 
	
	freeMultipleWords;
	return(bestFinalRowScore);
	
	
}

