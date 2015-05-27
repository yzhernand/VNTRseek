/*
 *  LinkedList.c
 *  Large_Indels
 *
 *  Created by Gary Benson on 5/23/10.
 *  Copyright 2010 Boston University. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "LinkedList.h"

struct linkedList *linkedListAlloc(void){
	//The returned linkedList points to a head and tail which are pointers to linkedListItems with null pointers. 
	//This makes the programming easier. 
	struct linkedList *newLinkedList = (struct linkedList *)calloc(1, sizeof(struct linkedList));
	struct linkedListItem *newHead = (struct linkedListItem *)calloc(1, sizeof(struct linkedListItem));
	struct linkedListItem *newTail = (struct linkedListItem *)calloc(1, sizeof(struct linkedListItem));
	newHead->next = NULL;
	newTail->next = NULL;
	newLinkedList->head = newHead;
	newLinkedList->tail = newTail;
	newLinkedList->whileLoopPointer = newHead->next;
	return(newLinkedList);
}

void linkedListInsertItemDataAtIndex(struct linkedList *linkedList, char *itemData, int index){
	//The inserted itemData pointer is wrapped in a linkedListItem which is given an itemIndex equal to index 
	//even if there are missing indexes between zero and index.  
	//An itemData pointer already at index is replaced with the new itemData pointer
	struct linkedListItem *current;
	struct linkedListItem *previous;
	struct linkedListItem *newItem;
  unsigned int notFound;	
	
  if((linkedList->head->next == NULL)&&(linkedList->tail->next == NULL))
	{	//linkedList is empty
		newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
		newItem->itemData = itemData;
		newItem->itemIndex = index;
		newItem->previous = NULL;
		newItem->next = NULL;
		linkedList->head->next = newItem;
		linkedList->tail->next = newItem;
	}
	else if((linkedList->head->next!=NULL)&&(linkedList->tail->next!=NULL))
	{
		//linkedList is not empty
		previous = linkedList->head;
		current = linkedList->head->next;
		notFound=1;
		while((current!=NULL)&&(current->itemIndex<=index)&&(notFound))
		{
			if(current->itemIndex == index)
			{
				//link in new itemData in place of old itemData
				current->itemData = itemData;
				notFound = 0;
			}
			else
			{
				previous = current;
				current=current->next;
			}
		}
		if(current==NULL)
		{
			//link in new item at end of linkedList
			newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
			newItem->itemData = itemData;
			newItem->itemIndex = index;
			newItem->next = NULL;
			newItem->previous = previous;
			previous->next = newItem;
			linkedList->tail->next = newItem;			
		}
		else if(current->itemIndex>index)
		{
			//link in new item between previous and current
			newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
			newItem->itemData = itemData;
			newItem->itemIndex = index;
			newItem->next = current;
			newItem->previous = previous;
			previous->next = newItem;
			
		}
	}
	else{
		//not found
		printf("\nError in adding item at index to linkedList, only one of head or tail is pointing to NULL");
	}
}

char *linkedListReturnItemDataAtIndex(struct linkedList *linkedList, int index){
	
	struct linkedListItem *current = linkedList->head->next;
	while((current!=NULL)&&(current->itemIndex<=index))
	{
		if(current->itemIndex == index) return((char *)(current->itemData));
		else current=current->next;
	}
	return(NULL); //not found
	
}

	
void linkedListInitializeForWhileLoop(struct linkedList *linkedList){
	//stores pointer to item after current item in whileLoopPointer
	linkedList->whileLoopPointer = linkedList->head->next;
}

char *linkedListNextItemDataInWhileLoop(struct linkedList *linkedList){
	//stores pointer to item after current item in whileLoopPointer	
	char *itemData;
	
	if(linkedList->whileLoopPointer!=NULL)
	{
		itemData = linkedList->whileLoopPointer->itemData;
		linkedList->whileLoopPointer = linkedList->whileLoopPointer->next;
		return(itemData);
	}
	else return(NULL);
}

void linkedListInitializeForScan(struct linkedList *linkedList){
	//retains current item pointer in whileLoopPointer
	linkedList->whileLoopPointer = linkedList->head;
}

char *linkedListNextItemDataInScan(struct linkedList *linkedList){
	//retains current item pointer in whileLoopPointer
	char *itemData;
	
	if(linkedList->whileLoopPointer->next!=NULL)
	{
		linkedList->whileLoopPointer = linkedList->whileLoopPointer->next;
		itemData = linkedList->whileLoopPointer->itemData;
		return(itemData);		
	}
	else return (NULL);
}

void linkedListInScanInsertItemDataBeforeCurrent(struct linkedList *linkedList, char *itemData){
	//retains current item pointer in whileLoopPointer

	//The inserted itemData pointer is wrapped in a linkedListItem which is then inserted before the item pointed to by 
	//whileLoopPointer

	struct linkedListItem *current;
	struct linkedListItem *newItem;
	struct linkedListItem *previous;
	if((linkedList->head->next == NULL)&&(linkedList->tail->next == NULL))
	{	
		//linkedList is empty
		//do nothing
		//but notify 
		printf("Error in linkedListInScanInsertItemDataBeforeCurrent.  Linked list is empty.  There is no current element."); 
	}
	else if((linkedList->head->next!=NULL)&&(linkedList->tail->next!=NULL))
	{
		//linked list is not empty
		if(linkedList->head->next == linkedList->whileLoopPointer)
		{
			//current points to first list element so insert at head
			//and leave whileLoopPointer alone
			linkedListInsertItemAtHead(linkedList, itemData);
		}
		else
		{
			//current is not first list element
			current = linkedList->whileLoopPointer;
			previous = current->previous;
			//link in new item between previous-> and current->
			newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
			newItem->itemData = itemData;
			newItem->next = current;
			newItem->previous = previous;
			current->previous = newItem;
			previous->next = newItem;
			
		}
	}
	else{
		//not found
		printf("\nError in linkedListInScanInsertItemDataBeforeCurrent, only one of head or tail is pointing to NULL");
	}
		
	
}

char *linkedListPreviousItemDataFromHereInScan(struct linkedList *linkedList)
//does not change whileLoopPointer
{
	return(linkedList->whileLoopPointer->previous->itemData);
}


void linkedListInsertItemAtHead(struct linkedList *linkedList, char *itemData){
	
	//The inserted itemData pointer is wrapped in a linkedListItem which is then inserted at the head of the list
	//index of each item is not set and is ignored
	
	struct linkedListItem *current;
	struct linkedListItem *newItem;
	
	if((linkedList->head->next == NULL)&&(linkedList->tail->next == NULL))
	{	//linkedList is empty
		newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
		newItem->itemData = itemData;
		newItem->next = NULL;
		newItem->previous = NULL;
		linkedList->head->next = newItem;
		linkedList->tail->next = newItem;
	}
	else if((linkedList->head->next!=NULL)&&(linkedList->tail->next!=NULL))
	{
		//linkedList is not empty
		current = linkedList->head->next;
		//link in new item between head-> and current->
		newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
		newItem->itemData = itemData;
		newItem->next = current;
		newItem->previous = NULL;
		current->previous = newItem;
		linkedList->head->next = newItem;
	}
	else{
		//not found
		printf("\nError in adding item at head to linkedList, only one of head or tail is pointing to NULL");
	}
}

void linkedListInsertItemAtTail(struct linkedList *linkedList, char *itemData){

	//The inserted itemData pointer is wrapped in a linkedListItem which is then inserted at the tail of the list
	//index of each item is not set and is ignored
	
	struct linkedListItem *current;
	struct linkedListItem *newItem;

	if((linkedList->head->next == NULL)&&(linkedList->tail->next == NULL))
	{	//linkedList is empty
		newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
		newItem->itemData = itemData;
		newItem->next = NULL;
		newItem->previous = NULL;
		linkedList->head->next = newItem;
		linkedList->tail->next = newItem;
	}
	else if((linkedList->head->next!=NULL)&&(linkedList->tail->next!=NULL))
	{
		//linkedList is not empty
		current = linkedList->tail->next;
		//link in new item between head-> and current->
		newItem = (struct linkedListItem *)calloc(1,sizeof(struct linkedListItem));
		newItem->itemData = itemData;
		newItem->next = NULL;
		newItem->previous = current;
		current->next = newItem;
		linkedList->tail->next = newItem;
	}
	else{
		//not found
		printf("\nError in adding item at tail to linkedList, only one of head or tail is pointing to NULL");
	}
	
	
}


char *linkedListReturnItemDataAtTail(struct linkedList *linkedList){

  return (linkedList->tail->next->itemData);
}