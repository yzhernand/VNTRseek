/****************************************************************
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
 *   Please direct all questions related to this program or bug reports or
 *updates to Yevgeniy Gelfand (ygelfand@bu.edu)
 *
 ****************************************************************/

/*************************************************
 *   This module implements 1 key hash and  2 key
 *   hash
 **************************************************/

#ifndef GHASH_H
#define GHASH_H

#include "../libs/easylife/easylife.h"
#include <stdio.h>
#include <stdlib.h>

/* 1key function prototypes */
typedef struct tagGSHITEM {
    unsigned int       key;
    void *             data;
    struct tagGSHITEM *next;
} GSHITEM;

typedef struct {
    unsigned int size;
    GSHITEM **   rack;
} GSHASH;

GSHASH *CreateSingleHash( unsigned int size );
void    DestroySingleHash( GSHASH *hash );
void *  GetSingleHashItem( GSHASH *hash, unsigned int key );
int     SetSingleHashItem( GSHASH *hash, unsigned int key, void *data );
int     ClearSingleHashItem( GSHASH *hash, unsigned int key );
void    PrintSingleHash( GSHASH *hash );

/* 2key function prototypes */
typedef struct tagGDHITEM {
    unsigned int       key1;
    unsigned int       key2;
    void *             data;
    struct tagGDHITEM *next;
} GDHITEM;

typedef struct {
    unsigned int size;
    GDHITEM **   rack;
} GDHASH;

GDHASH *CreateDoubleHash( unsigned int size );
void    DestroyDoubleHash( GDHASH *hash, void ( *destroy )( void *item ) );
void *  GetDoubleHashItem( GDHASH *hash, unsigned int key1, unsigned int key2 );
int     SetDoubleHashItem(
      GDHASH *hash, unsigned int key1, unsigned int key2, void *data );
int  ClearDoubleHashItem( GDHASH *hash, unsigned int key1, unsigned int key2 );
void PrintDoubleHash( GDHASH *hash );

/* 3key function prototypes */
typedef struct tagGTHITEM {
    unsigned int       key1;
    unsigned int       key2;
    unsigned int       key3;
    void *             data;
    struct tagGTHITEM *next;
} GTHITEM;

typedef struct {
    unsigned int size;
    GTHITEM **   rack;
    EASY_LIST *  walker;
} GTHASH;

GTHASH *CreateTripleHash( unsigned int size );
void    DestroyTripleHash( GTHASH *hash, void ( *destroy )( void *item ) );
void *  GetTripleHashItem(
    GTHASH *hash, unsigned int key1, unsigned int key2, unsigned int key3 );
int SetTripleHashItem( GTHASH *hash, unsigned int key1, unsigned int key2,
  unsigned int key3, void *data );

/****************************************************
 *               function definition
 *****************************************************/

int EasyLifeGetPrime( int size );

/* Robert Jenkins' 32 bit integer hash function */
unsigned int rjhash( unsigned int a ) {
    a = ( a + 0x7ed55d16 ) + ( a << 12 );
    a = ( a ^ 0xc761c23c ) ^ ( a >> 19 );
    a = ( a + 0x165667b1 ) + ( a << 5 );
    a = ( a + 0xd3a2646c ) ^ ( a << 9 );
    a = ( a + 0xfd7046c5 ) + ( a << 3 );
    a = ( a ^ 0xb55a4f09 ) ^ ( a >> 16 );
    return a;
}

/****************************************************
 *               1 key hash
 *****************************************************/

GSHASH *CreateSingleHash( unsigned int size ) {
    GSHASH *hash;

    /* allocate the hash structure */
    hash = (GSHASH *) malloc( sizeof( GSHASH ) );
    if ( hash == NULL )
        return NULL;

    /* if zero size use default size */
    if ( size <= 0 )
        size = 14591;
    else {
        size = EasyLifeGetPrime( size );
    }

    /* allocate the item rack with zeroes */
    hash->rack = (GSHITEM **) calloc( size, sizeof( GSHITEM * ) );
    if ( hash->rack == NULL ) {
        free( hash );
        return NULL;
    }

    /* set value and return */
    hash->size = size;

    return hash;
}

void DestroySingleHash( GSHASH *hash ) {
    unsigned int i;
    GSHITEM *    curr;

    if ( NULL == hash )
        return;

    /* step trough and free all items in rack */
    for ( i = 0; i < hash->size; i++ ) {
        /* free all items in that slot of the rack */
        while ( hash->rack[i] ) {
            curr          = hash->rack[i];
            hash->rack[i] = curr->next;
            free( curr );
        }
    }

    /* free rack and structure */
    free( hash->rack );
    free( hash );

    return;
}

void *GetSingleHashItem( GSHASH *hash, unsigned int key ) {
    GSHITEM *    hitem;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key ) ) % hash->size;

    for ( hitem = hash->rack[hkey]; hitem != NULL; hitem = hitem->next ) {
        if ( hitem->key == key )
            break;
    }

    return ( hitem == NULL ? NULL : hitem->data );
}

int SetSingleHashItem( GSHASH *hash, unsigned int key, void *data ) {
    GSHITEM *    hitem;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key ) ) % hash->size;

    for ( hitem = hash->rack[hkey]; hitem != NULL; hitem = hitem->next ) {
        if ( hitem->key == key )
            break;
    }

    /* if item is not in list insert */
    if ( hitem == NULL ) {
        hitem = (GSHITEM *) malloc( sizeof( GSHITEM ) );
        if ( hitem == NULL )
            return -1;
        hitem->key       = key;
        hitem->next      = hash->rack[hkey];
        hash->rack[hkey] = hitem;
    }
    hitem->data = data;

    return 0;
}

int ClearSingleHashItem( GSHASH *hash, unsigned int key ) {
    GSHITEM *    hitem, *prev;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key ) ) % hash->size;

    prev = NULL;
    for ( hitem = hash->rack[hkey]; hitem != NULL; ) {
        if ( hitem->key == key )
            break;
        prev  = hitem;
        hitem = hitem->next;
    }

    /* if item was found remove */
    if ( hitem != NULL ) {
        if ( prev == NULL ) {
            hash->rack[hkey] = hitem->next;
            free( hitem );
        } else {
            prev->next = hitem->next;
            free( hitem );
        }
    }

    return 0;
}

void PrintSingleHash( GSHASH *hash ) {
    GSHITEM *    hitem;
    unsigned int i;

    /* compute the hash key */
    for ( i = 0; i < hash->size; i++ ) {
        /* print rack id */
        printf( "\n%d :", i );

        for ( hitem = hash->rack[i]; hitem != NULL; hitem = hitem->next ) {
            printf( " %d", hitem->key );
        }
    }

    return;
}

/****************************************************
 *               2 key hash
 *****************************************************/

GDHASH *CreateDoubleHash( unsigned int size ) {
    GDHASH *hash;

    /* allocate the hash structure */
    hash = (GDHASH *) malloc( sizeof( GDHASH ) );
    if ( hash == NULL )
        return NULL;

    /* if zero size use default size */
    if ( size <= 0 )
        size = 14591;
    else {
        size = EasyLifeGetPrime( size );
    }

    /* allocate the item rack with zeroes */
    hash->rack = (GDHITEM **) calloc( size, sizeof( GDHITEM * ) );
    if ( hash->rack == NULL ) {
        free( hash );
        return NULL;
    }

    /* set value and return */
    hash->size = size;

    return hash;
}

void DestroyDoubleHash( GDHASH *hash, void ( *destroy )( void *item ) ) {
    unsigned int i;
    GDHITEM *    curr;

    if ( NULL == hash )
        return;

    /* step trough and free all items in rack */
    for ( i = 0; i < hash->size; i++ ) {
        /* free all items in that slot of the rack */
        while ( hash->rack[i] ) {
            curr          = hash->rack[i];
            hash->rack[i] = curr->next;
            if ( destroy )
                destroy( curr->data );
            free( curr );
        }
    }

    /* free rack and structure */
    free( hash->rack );
    free( hash );

    return;
}

void *GetDoubleHashItem( GDHASH *hash, unsigned int key1, unsigned int key2 ) {
    GDHITEM *    hitem;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key1 ) + rjhash( key2 ) ) % hash->size;

    for ( hitem = hash->rack[hkey]; hitem != NULL; hitem = hitem->next ) {
        if ( hitem->key1 == key1 && hitem->key2 == key2 )
            break;
    }

    return ( hitem == NULL ? 0 : hitem->data );
}

int SetDoubleHashItem(
  GDHASH *hash, unsigned int key1, unsigned int key2, void *data ) {
    GDHITEM *    hitem;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key1 ) + rjhash( key2 ) ) % hash->size;

    for ( hitem = hash->rack[hkey]; hitem != NULL; hitem = hitem->next ) {
        if ( hitem->key1 == key1 && hitem->key2 == key2 )
            break;
    }

    /* if item is not in list insert */
    if ( hitem == NULL ) {
        hitem = (GDHITEM *) malloc( sizeof( GDHITEM ) );
        if ( hitem == NULL )
            return -1;
        hitem->key1      = key1;
        hitem->key2      = key2;
        hitem->next      = hash->rack[hkey];
        hash->rack[hkey] = hitem;
    }
    hitem->data = data;

    return 0;
}

int ClearDoubleHashItem( GDHASH *hash, unsigned int key1, unsigned int key2 ) {
    GDHITEM *    hitem, *prev;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key1 ) + rjhash( key2 ) ) % hash->size;

    prev = NULL;
    for ( hitem = hash->rack[hkey]; hitem != NULL; ) {
        if ( hitem->key1 == key1 && hitem->key2 == key2 )
            break;
        prev  = hitem;
        hitem = hitem->next;
    }

    /* if item was found remove */
    if ( hitem != NULL ) {
        if ( prev == NULL ) {
            hash->rack[hkey] = hitem->next;
            free( hitem->data );
            free( hitem );
        } else {
            prev->next = hitem->next;
            free( hitem->data );
            free( hitem );
        }
    }

    return 0;
}

void PrintDoubleHash( GDHASH *hash ) {
    GDHITEM *    hitem;
    unsigned int i;

    /* compute the hash key */
    for ( i = 0; i < hash->size; i++ ) {
        /* print rack id */
        printf( "\n%d :", i );

        for ( hitem = hash->rack[i]; hitem != NULL; hitem = hitem->next ) {
            printf( " %d-%d", hitem->key1, hitem->key2 );
        }
    }

    return;
}

/****************************************************
 *               3 key hash
 *****************************************************/

GTHASH *CreateTripleHash( unsigned int size ) {
    GTHASH *hash;

    /* allocate the hash structure */
    hash = (GTHASH *) malloc( sizeof( GTHASH ) );
    if ( hash == NULL )
        return NULL;

    /* if zero size use default size */
    if ( size <= 0 )
        size = 14591;
    else {
        size = EasyLifeGetPrime( size );
    }

    /* allocate the item rack with zeroes */
    hash->rack = (GTHITEM **) calloc( size, sizeof( GTHITEM * ) );
    if ( hash->rack == NULL ) {
        free( hash );
        return NULL;
    }

    /* set value and return */
    hash->size = size;

    /* create a list to keep track of rack items (for faster clearing and
     * walkthough) */
    hash->walker = EasyListCreate( NULL, NULL );

    return hash;
}

void DestroyTripleHash( GTHASH *hash, void ( *destroy )( void *item ) ) {
    unsigned int *i;
    GTHITEM *     curr;
    EASY_NODE *   iter;

    if ( NULL == hash )
        return;

    /* step trough and free all items in rack */
    for ( iter = hash->walker->head; iter != NULL; iter = iter->next ) {
        i = (unsigned int *) EasyListItem( iter );

        /* free all items in that slot of the rack */
        while ( hash->rack[*i] ) {
            curr           = hash->rack[*i];
            hash->rack[*i] = curr->next;
            if ( destroy )
                destroy( curr->data );
            free( curr );
        }
    }

    /* free rack and structure */
    EasyListDestroy( hash->walker );
    free( hash->rack );
    free( hash );

    return;
}

void ClearTripleHash( GTHASH *hash, void ( *destroy )( void *item ) ) {
    unsigned int *i;
    GTHITEM *     curr;
    EASY_NODE *   iter;

    if ( NULL == hash )
        return;

    /* step trough and free all items in rack */
    // printf("\n\nwalker size: %d hash->size: %d: ",hash->walker->size,
    // hash->size ); fflush(stdout);
    for ( iter = hash->walker->head; iter != NULL; iter = iter->next )
    // for(i=0;i<hash->size;i++)
    {
        i = (unsigned int *) EasyListItem( iter );
        // printf(" %d(",i); fflush(stdout);
        // printf("size: %d, item %d\n",hash->size,i); fflush(stdout);
        // free all items in that slot of the rack
        while ( hash->rack[*i] ) {
            curr           = hash->rack[*i];
            hash->rack[*i] = curr->next;
            // printf(" zz"); fflush(stdout);
            if ( destroy )
                destroy( curr->data );
            free( curr );
        }
        // printf(")"); fflush(stdout);
    }
    EasyListDestroy( hash->walker );
    hash->walker = EasyListCreate( NULL, NULL );
}

void *GetTripleHashItem(
  GTHASH *hash, unsigned int key1, unsigned int key2, unsigned int key3 ) {
    GTHITEM *    hitem;
    unsigned int hkey;

    /* compute the hash key */
    hkey = ( rjhash( key1 ) + rjhash( key2 ) + rjhash( key3 ) ) % hash->size;

    for ( hitem = hash->rack[hkey]; hitem != NULL; hitem = hitem->next ) {
        if ( hitem->key1 == key1 && hitem->key2 == key2 && hitem->key3 == key3 )
            break;
    }

    return ( hitem == NULL ? 0 : hitem->data );
}

int SetTripleHashItem( GTHASH *hash, unsigned int key1, unsigned int key2,
  unsigned int key3, void *data ) {
    GTHITEM *     hitem;
    unsigned int *hkey;

    /* compute the hash key */
    *hkey = ( rjhash( key1 ) + rjhash( key2 ) + rjhash( key3 ) ) % hash->size;

    /* set walker key on first rack item */
    if ( NULL == hash->rack[*hkey] )
        EasyListInsertHead( hash->walker, hkey );

    for ( hitem = hash->rack[*hkey]; hitem != NULL; hitem = hitem->next ) {
        if ( hitem->key1 == key1 && hitem->key2 == key2 && hitem->key3 == key3 )
            break;
    }

    /* if item is not in list insert */
    if ( hitem == NULL ) {
        hitem = (GTHITEM *) malloc( sizeof( GTHITEM ) );
        if ( hitem == NULL )
            return -1;
        hitem->key1       = key1;
        hitem->key2       = key2;
        hitem->key3       = key3;
        hitem->next       = hash->rack[*hkey];
        hash->rack[*hkey] = hitem;
    }
    hitem->data = data;

    return 0;
}

#endif
