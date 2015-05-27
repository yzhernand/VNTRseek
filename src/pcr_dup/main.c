#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LinkedList.h"
#include "bitwise LCS single word.h"
#include "bitwise LCS multiple word.h"

#define MAXREADLENGTH 1000
#define FLANKLENGTH 20
#define MAXDIFFFORPCRDUPLICATE 2
//may want to use a combination of length and mismatch/indel differences

struct read {
	int readid;
	int length;
	int TRstart;
	int	TRend;
	char *sequence;
	char *fastaHeader1;
	char *fastaHeader2;
};

char *linkedListInScanRemoveCurrentItemData(struct linkedList *linkedList)
{
	//removes the item pointed to by whileLoopPointer
	//resets whileLoopPointer to PREVIOUS item or head if none
	//does not delete the dataItem
	//returns a pointer to the dataItem (must be cast after return)
	//so that caller may delete if required
	//returns NULL if the whileLoopPointer points to NULL or to the head 
	
	struct linkedListItem *currentItem;
	struct linkedListItem *previousItem;
	struct linkedListItem *nextItem;
	char *currentItemData;
	
	currentItem = linkedList->whileLoopPointer;
	if ((currentItem == NULL)||(currentItem == linkedList->head)) return(NULL); //at head or null item
	
	previousItem = currentItem->previous;
	nextItem = currentItem->next;
	
	
	if (previousItem!=NULL)	{
		previousItem->next = currentItem->next;
		linkedList->whileLoopPointer = previousItem;
	}
	else {//first item in list
		linkedList->head->next = currentItem->next;
		linkedList->whileLoopPointer = linkedList->head;
	}
	if (nextItem!=NULL){
		nextItem->previous = currentItem->previous;
	}
	else {//last item in list
		linkedList->tail->next = currentItem->previous;
	}
	currentItemData = currentItem->itemData;
	free(currentItem);
	return(currentItemData);
	
}
void linkedListMergeSorted(struct linkedList *linkedList1, struct linkedList *linkedList2)
{
	struct read *readItemData1;
	struct read *readItemData2;
	//merges two linked lists, 
	//assumes both lists are sorted on the TRend field
	//puts items in the second list into the first list, maintaining the 
	//sort on the TRend field
	
	
	linkedListInitializeForScan(linkedList1);
	linkedListInitializeForScan(linkedList2);
	readItemData1 = (struct read *)linkedListNextItemDataInScan(linkedList1);
	readItemData2 = (struct read *)linkedListNextItemDataInScan(linkedList2);
	while ((readItemData1!=NULL)||(readItemData2!=NULL))
	{
		if (readItemData1 == NULL)//empty list 1
		{
			linkedListInsertItemAtTail(linkedList1,(char *)readItemData2);
			readItemData2  = (struct read *)linkedListNextItemDataInScan(linkedList2);
		}
		else if (readItemData2 == NULL)//empty list 2 
		{
			//done
			return;
		}
		else //both lists not empty, compare TRends
		{
			if (readItemData2->TRend < readItemData1->TRend) {
				linkedListInScanInsertItemDataBeforeCurrent(linkedList1,(char *)readItemData2);
				//at this point,readItemData2 is on two lists
				
				//get next readItemData2
				readItemData2 = (struct read *)linkedListNextItemDataInScan(linkedList2);
			}
			else {
				readItemData1 = (struct read *)linkedListNextItemDataInScan(linkedList1);
			}

		}
	}
}



void linkedListPrintItemData(struct linkedList *linkedList)
{
	struct read *readItemData;

	linkedListInitializeForScan(linkedList);
	while ((readItemData = (struct read *)linkedListNextItemDataInScan(linkedList))!=NULL) {
		
		printf("\n%s %d %d",readItemData->fastaHeader1,readItemData->TRstart,readItemData->TRend);
	}
}
void linkedListBubbleSort(struct linkedList *linkedList) //need to sort on field in ItemData
{
	int numItems;
	int counter,second,outerloopcounter,change;
	struct read *readItemDataA, *readItemDataB;
	struct linkedListItem *pointsToB, *pointsToA;
	
	linkedListInitializeForScan(linkedList);
	//count number of items in linked list
	numItems=0;
	while (linkedListNextItemDataInScan(linkedList)!=NULL) {
		numItems++;
	}
	printf(" numItems: %d",numItems);
	
	if (counter>=2) {
		change = 1;
		outerloopcounter = 0;
		while(change){
			change = 0;
			//read first two items
			linkedListInitializeForScan(linkedList);
			readItemDataA = (struct read *)linkedListNextItemDataInScan(linkedList);
			readItemDataB = (struct read *)linkedListNextItemDataInScan(linkedList);
			second=1; //counts position in list of readItemData2 (zero based)
			while(second < numItems - outerloopcounter){
				if(readItemDataA->TRend > readItemDataB->TRend){
					change = 1;
					//swap
					pointsToB = linkedList->whileLoopPointer;
					pointsToA = pointsToB->previous;
					//next two should be mutually exclusive because if A is the first item
					//in the list, then its previous pointer is NULL
					if (linkedList->head->next == pointsToA) {
						linkedList->head->next = pointsToB;
					}
					if (pointsToA->previous!=NULL) {
						pointsToA->previous->next = pointsToB;
					}
					//next two should be mutually exclusive because if B is the last item
					//in the list, then its next pointer is NULL
					if (linkedList->tail->next == pointsToB) {
						linkedList->tail->next = pointsToA;
					}
					if (pointsToB->next!=NULL) {
						pointsToB->next->previous = pointsToA;
					}
					pointsToA->next = pointsToB->next;
					pointsToB->next = pointsToA;
					pointsToB->previous = pointsToA->previous;
					pointsToA->previous = pointsToB;
					//what ends up in linkedlist while loop counter?
					linkedList->whileLoopPointer = pointsToA;
					
				}
				else {
					readItemDataA = readItemDataB;
				}

				readItemDataB = (struct read *)linkedListNextItemDataInScan(linkedList);
				second++;
			}
			outerloopcounter++;
		}
		
	}
}	


int main (int argc, const char * argv[]) {

	int h,i,j;
	char line[100000];
	char fastaHeader1[1000], fastaHeader2[1000], pattern[10000], sequence[100000];
	int TRstart, TRend;
	long long int readid;
	int lineCounter;
	int sscanfreturnvalue;
	struct read *readItemData;
	int stringlength;
	struct linkedList **trStartList;
	struct read *readItemDataA, *readItemDataB;
	int endVariance, numDifferences;
	struct linkedListItem *currentItemPtr;
	int withinEndVariance;
	int trEndDelta, trStartDelta, readEndDelta;
	int *coverageInRead;
	int start, end;
	struct linkedList *windowList;
	int trailing;
	int cantSpan;
	int readLength;
	float copies;
	int patternsize;
	int bothEndsCantSpan;
	int comparisonCount;
	int lcsScore;
	int highLCSScoreCount;
	int shortestReadLengthInPair;
	
	FILE *readFilep, *outFile;
	
	printf("\n");
	for(i=0;i<argc;i++)
		printf("%s ",argv[i]);
	
	printf("\n%d %d",sizeof(char *), sizeof(long long int));
	
	if(argc!=5)
	{
		printf("\n\nPlease use: %s ReadFile OutFile EndVariance NumDifferences", argv[0]);
		printf("\n\nWhere:");
		printf("\n  ReadFile is a file containing the reads and associated information regarding contained TRs");
		printf("\n     Note the expected format is .index.seq");
		printf("\n     Note that the same read may appear several times in the file because of different TRs");
		printf("\n  OutFile is a text file which contains the PCR duplicate cluster lists composed offasta headers");
		printf("\n  EndVariance is the maximum allowed difference in the TRstart and TRend positions (individually, not combined)"); 
		printf("\n  NumDifferences is the number of differences allowed between two read sequences to declare them duplicates");
		printf("\n     Note, NumDifferences includes the end variance");
		printf("\n     Note.  Right now, the program uses MAXDIFFFORPCRDUPLICATE instead of the input parameter.");
		printf("\n");
		printf("\nThis program finds PCR duplicates in reads that have already been processed and found to contain TRs");
		printf("\nThe reads are grouped by starting and ending location of the TRs so all-pairs comparison is not required");
		printf("\nComparison is done using fast LCS");
		printf("\nFormat of output is: ???");
		printf("\n ");
		exit(-1);
	}
	
	endVariance = atoi(argv[3]);
	numDifferences = atoi(argv[4]);
	
	//allocate structure for "read" data
	trStartList = (struct linkedList **) calloc(MAXREADLENGTH+1,sizeof(struct linkedList *));
	for (i=0; i<=MAXREADLENGTH; i++) {
		trStartList[i]=linkedListAlloc();
	}
	
	//open out file
	outFile = fopen(argv[2], "w");
	if (outFile==NULL) {
		printf("\nCan't open for writing file %s",argv[2]);
		exit(0);
	}

	fprintf(outFile,"\n");
	for(i=0;i<argc;i++)
		fprintf(outFile,"%s ",argv[i]);
	
	//read "read" data and store each read on a linked list with the same TRstart
	if(strcmp(argv[1],"stdin")!=0)
	{
		readFilep = fopen(argv[1], "r");
		if (readFilep==NULL) {
			printf("\nCan't open file %s",argv[1]);
			exit(0);
		}
	}
	else {
		readFilep = stdin;
	}

	lineCounter = 0;
	while (!feof(readFilep)) {
		lineCounter++;
		if(fgets(line, 10000, readFilep)!=NULL)	//NULL if reading error or blank end line
		{
			//the following is the format of a .index.seq file
			sscanfreturnvalue = sscanf(line,"%lld %s %s %d %d %f %d %s %s",&readid,fastaHeader1,fastaHeader2,&TRstart,&TRend,&copies,&patternsize,pattern,sequence); //%lld is 64 bit signed integer
			if (sscanfreturnvalue == 9) //not  if blank line or misformatted line
			{
				//printf("\n%lld %s %s %d %d %f %d %s %s",readid,fastaHeader1,fastaHeader2,TRstart,TRend,copies,patternsize,pattern,sequence); //%lld is 64 bit signed integer
					   
				//store read data in readItemData
				readItemData = (struct read *)calloc(1, sizeof(struct read));
				stringlength = strlen(fastaHeader1);
				readItemData->fastaHeader1 = (char *)calloc(stringlength+1, sizeof(char));
				stringlength = strlen(fastaHeader2);
				readItemData->fastaHeader2 = (char *)calloc(stringlength+1, sizeof(char));
				stringlength = strlen(sequence);
				readItemData->sequence = (char *)calloc(stringlength+1, sizeof(char));
				strcpy(readItemData->fastaHeader1,fastaHeader1);
				strcpy(readItemData->fastaHeader2,fastaHeader2);
				strcpy(readItemData->sequence,sequence);
				readItemData->readid=readid;
				readItemData->TRstart=TRstart;
				readItemData->TRend=TRend;
				readItemData->length=strlen(sequence);
				
				//store readItemData on linked list by TRstart
				linkedListInsertItemAtHead(trStartList[TRstart], (char *)readItemData);
			}
			
			else printf("\nerror: %s",line);
		}
	}
	printf("\nNumber of reads = %d",lineCounter);
	fprintf(outFile,"\nNumber of reads = %d",lineCounter);
	

	//reorder within the lists so they are sorted by TRend
	printf("\nSorting");
	for (i=0; i<=MAXREADLENGTH; i++) {
		if (trStartList[i]->head->next!=NULL) {
			//printf("\nTRstart position = %d",i);
			printf("\n%d",i);
			//linkedListPrintItemData(trStartList[i]);
			//printf("\n***");
			linkedListBubbleSort(trStartList[i]);
			//printf("\n***");
			//linkedListPrintItemData(trStartList[i]);
		}
	}
	
	//calculate TR coverage within all reads to look for end bias
	coverageInRead = (int *)calloc(MAXREADLENGTH+1,sizeof(int));
	for (i=0; i<=MAXREADLENGTH; i++) {
		linkedListInitializeForScan(trStartList[i]);
		start = i; //start of TR coverage in read
		while((readItemData = (struct read *)linkedListNextItemDataInScan(trStartList[i]))!=NULL) 
		{
			end = readItemData->TRend; //end of TR coverage in read
			for (j=start; j<=end; j++) {
				coverageInRead[j]++;
			}
		}
	}
	printf("\n\nCoverage:");
	for (i=0; i<=MAXREADLENGTH; i++) printf("\n%3d %7d",i,coverageInRead[i]);
	

	//calculate which reads cannot be spanning because either or both ends are too short
	cantSpan=0;
	bothEndsCantSpan=0;
	for (i=0; i<=MAXREADLENGTH; i++) {
		linkedListInitializeForScan(trStartList[i]);
		start = i; //start of TR coverage in read
		while((readItemData = (struct read *)linkedListNextItemDataInScan(trStartList[i]))!=NULL) 
		{
			end = readItemData->TRend; //end of TR coverage in read
			readLength = readItemData->length;
			if ((start<FLANKLENGTH)||(readLength-end<FLANKLENGTH)) {
				cantSpan++;
			}
			if ((start<FLANKLENGTH)&&(readLength-end<FLANKLENGTH)) {
				bothEndsCantSpan++;
			}
		}
	}
	printf("\nNumber of reads that can't span at one end or the other for flank length %d = %d",FLANKLENGTH,cantSpan);
	printf("\nNumber of reads that can't span at either end for flank length %d = %d",FLANKLENGTH,bothEndsCantSpan);
	
	//exit(-1);
	
	
//	//create windowList to hold merged lists for TRend comparison
//	windowList = linkedListAlloc();
//	
//	//testing
//	linkedListMergeSorted(windowList,trStartList[110]);
//	linkedListPrintItemData(windowList);
//	linkedListMergeSorted(windowList,trStartList[112]);
//	linkedListPrintItemData(windowList);
//	linkedListInitializeForScan(windowList);
//	while((readItemData = (struct read *)linkedListNextItemDataInScan(windowList))!=NULL)
//	{
//		if ((readItemData->TRstart == 110)||(readItemData->TRstart == 112)) {
//			linkedListInScanRemoveCurrentItemData(windowList);
//		}
//	}
//	printf("\n\n");
//	linkedListPrintItemData(windowList);
//	printf("\n\n");
//	linkedListPrintItemData(trStartList[110]);
//	linkedListPrintItemData(trStartList[112]);
		
	windowList = linkedListAlloc();
	comparisonCount=0;
	highLCSScoreCount=0;
	for (i=1; i<=MAXREADLENGTH; i++) {
		//load leading trStartList into windowlinkedListInitializeForScan(trStartList[i]);
		printf("\nAdding %d",i);
		linkedListMergeSorted(windowList,trStartList[i]);
		//remove trailing trStartList items if i-endVariance-1 not less than one
		trailing = i - endVariance - 1;
		if (trailing>=1) {
			printf(" removing %d",trailing);
			linkedListInitializeForScan(windowList);
			while ((readItemData = (struct read *)linkedListNextItemDataInScan(windowList))!=NULL) {
				if (readItemData->TRstart==trailing) {
					linkedListInScanRemoveCurrentItemData(windowList);
				}
			}			
		}
		//for debugging
		/*
		 if ((i==110)||(i==112))
		{
			
			linkedListPrintItemData(windowList);
		}
		*/ 
		//now scan through windowList
		linkedListInitializeForScan(windowList);
		while((readItemDataA = (struct read *)linkedListNextItemDataInScan(windowList))!=NULL) 
		{
			if (readItemDataA->TRstart==i)
			{
				//readItemDataA is from trStartList[i]
				//it has to be compared with all - endVariance from its position
				//up and down the list.  Tests for + endVariance are done in 
				//later iterations.
				
				for (h=1; h<=2; h++) 
				{
					//h==1 means up towards head
					//h==2 means down towards tail
					currentItemPtr = windowList->whileLoopPointer;
					withinEndVariance = 1;
					while (((h==1)&&(currentItemPtr->previous!=NULL)&&(withinEndVariance))
						   ||((h==2)&&(currentItemPtr->next!=NULL)&&(withinEndVariance)))
					{
						if (h==1) currentItemPtr = currentItemPtr->previous;
						else currentItemPtr = currentItemPtr->next;
						
						readItemDataB = (struct read *)currentItemPtr->itemData;
						trEndDelta = readItemDataA->TRend - readItemDataB->TRend;
						if (trEndDelta<0) trEndDelta = -trEndDelta; //abs value of TRend difference
						if (trEndDelta<=endVariance) 
						{
							trStartDelta = readItemDataA->TRstart - readItemDataB->TRstart;
							if (trStartDelta<0) trStartDelta = -trStartDelta; //abs value of TRstart difference
							//is the next comparison unnecessary?
							if (trStartDelta<=endVariance)
							{
								readEndDelta = readItemDataA->length - readItemDataB->length;
								if (readEndDelta<0) readEndDelta = - readEndDelta;
								if (readEndDelta<=endVariance) {
								//is this all comparisons?
									if ((h==1)||(trStartDelta>0)) //only do comparison on same TRstart going up
									{
										//compare sequences
										comparisonCount++;
										//printf("\ncompare: %s %d %d|%s %d %d",
										//readItemDataA->fastaHeader1,readItemDataA->TRstart,readItemDataA->TRend,
										//readItemDataB->fastaHeader1,readItemDataB->TRstart,readItemDataB->TRend);
										lcsScore=LCS_multiple_word(readItemDataA->sequence, readItemDataB->sequence, readItemDataA->length,readItemDataB->length);
										if (lcsScore==-1)//error during LCS
										{
											printf("\nError during comparison: %s %d %d %d|%s %d %d %d",
												   readItemDataA->fastaHeader1,readItemDataA->TRstart,readItemDataA->TRend,readItemDataA->length,
												   readItemDataB->fastaHeader1,readItemDataB->TRstart,readItemDataB->TRend,readItemDataB->length);
										}
										shortestReadLengthInPair = (readItemDataA->length <= readItemDataB->length)?readItemDataA->length:readItemDataB->length;
										if(lcsScore >= shortestReadLengthInPair-MAXDIFFFORPCRDUPLICATE)
										{
											highLCSScoreCount++;
											/* gelfand, aug 11, changed to be like before to have readit and not fastaheader */
											fprintf(outFile,"\ncompare: %d %d %d %d|%d %d %d %d|LCS: %d",
												   readItemDataA->readid,readItemDataA->TRstart,readItemDataA->TRend,readItemDataA->length,
												   readItemDataB->readid,readItemDataB->TRstart,readItemDataB->TRend,readItemDataB->length,
												   lcsScore);
											//fprintf(outFile,"\n%s\n%s",readItemDataA->sequence,readItemDataB->sequence);
											
											/*printf("\ncompare: %s %d %d %d|%s %d %d %d|LCS: %d",
												   readItemDataA->fastaHeader1,readItemDataA->TRstart,readItemDataA->TRend,readItemDataA->length,
												   readItemDataB->fastaHeader1,readItemDataB->TRstart,readItemDataB->TRend,readItemDataB->length,
												   lcsScore);
											 */
										}
										
									}
								}
							}
						}
						else withinEndVariance = 0;
					}
				}
			}			
		}
	}
	printf("\nNumber of Comparisons is %d",comparisonCount);
	printf("\nNumber of High LCS Scores is %d",highLCSScoreCount);
	fprintf(outFile,"\nNumber of Comparisons is %d",comparisonCount);
	fprintf(outFile,"\nNumber of High LCS Scores is %d",highLCSScoreCount);
	
	//store reads determined to be duplicates in clusters
	//report clusters
	
    printf("\nDone.");
	fclose(readFilep);
	fclose(outFile);
	return 0;
}
