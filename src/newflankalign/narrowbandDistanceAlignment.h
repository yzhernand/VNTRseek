/*
 *  narrowbandDistanceAlignment.h
 *  NarrowbandDistanceAlignment
 *
 *  Created by Gary Benson on 1/20/12.
 *  Copyright 2012 Boston University. All rights reserved.
 *
 */

int narrowbandUnitCostDistanceNoEndPenaltySeq2(char *seq1, char *seq2, int len1, int len2, int maxerr);
int narrowbandUnitCostDistanceNoEndPenaltySeq2LowMem(char *seq1, char *seq2, int len1, int len2, int maxerr);
int narrowbandUnitCostDistanceNoEndPenaltySeq2LowMemUsesPointers(char *seq1, char *seq2, int len1, int len2, int maxerr);


