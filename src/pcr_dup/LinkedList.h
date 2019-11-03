/*
 *  LinkedList.h
 *  Large_Indels
 *
 *  Created by Gary Benson on 5/23/10.
 *  Copyright 2010 Boston University. All rights reserved.
 *
 */

#ifndef LinkedListHeader
#define LinkedListHeader

// This is the linked list
// It has an itemIndex field so it can mimic an array.  It is like a mutable
// array (size mutable) It supports allocation, inserting data, returning data,
// and stepping through the list one item at a time The data is not part of the
// linked list.  It is a separate data structure which is pointed to by each
// element of the linked list Look closely at the procedures and use (in
// ReadGroup) of the two functions for While Loop The initializer sets the
// whileLoopPointer to the first element in the list The NextItem retrieves the
// data and updates the whileLoopPointer to the next element in the list

struct linkedListItem {
    char *                 itemData;
    int                    itemIndex;
    struct linkedListItem *next;
    struct linkedListItem *previous;
};

struct linkedList {
    struct linkedListItem *head;
    struct linkedListItem *tail;
    struct linkedListItem *whileLoopPointer;
};

struct linkedList *linkedListAlloc( void );
void               linkedListInsertItemDataAtIndex(
                struct linkedList *linkedList, char *itemData, int index );
char *linkedListReturnItemDataAtIndex(
  struct linkedList *linkedList, int index );
void  linkedListInitializeForWhileLoop( struct linkedList *linkedList );
char *linkedListNextItemDataInWhileLoop( struct linkedList *linkedList );
void  linkedListInsertItemAtHead(
   struct linkedList *linkedList, char *itemData );
void linkedListInsertItemAtTail(
  struct linkedList *linkedList, char *itemData );
void  linkedListInitializeForScan( struct linkedList *linkedList );
char *linkedListNextItemDataInScan( struct linkedList *linkedList );
void  linkedListInScanInsertItemDataBeforeCurrent(
   struct linkedList *linkedList, char *itemData );
char *linkedListPreviousItemDataFromHereInScan( struct linkedList *linkedList );
char *linkedListReturnItemDataAtTail( struct linkedList *linkedList );
#endif