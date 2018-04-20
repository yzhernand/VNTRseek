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
*   Please direct all questions related to this program or bug reports or updates to
*   Yevgeniy Gelfand (ygelfand@bu.edu)
*
****************************************************************/

/***************************************************************
    cluster.h : Implementation of connected components algorithm
                used for clustering.
                Items inserted must have a unique integer key.
****************************************************************/



#ifndef CLUSTER_H
#define CLUSTER_H

#include <stdio.h>
#include <float.h>
#include "profile.h"

                            /* INTERFACE */

typedef struct tagCLUSTERBASE CLUSTERBASE;

CLUSTERBASE*    ClusterBaseCreate(void);
int             ClusterBaseAddLink(CLUSTERBASE* base, int key1, int key2, int direction);
int             ClusterBaseAddLinkProfile(CLUSTERBASE* base, int key1 , int key2, PROFILE *prof1, PROFILE *prof1rc,PROFILE *prof2, PROFILE *prof2rc, int dir);
int				ClusterBaseAddUnique(CLUSTERBASE* base, int key1); // added to support 1 repeat clusters
int             ClusterBasePrint(CLUSTERBASE* base, FILE* fp);
int             ClusterBaseDestroy(CLUSTERBASE* base);
int             ClusterBaseAskLink(CLUSTERBASE* base, int key1, int key2);
int             ClusterBaseRestoreDirection(CLUSTERBASE* base, int key, int direction);



                        /* IMPLEMENTATION */


#define CLUSTERBASESIZE     1000000
#define CLUSTERMAGICNUM     15423

typedef struct tagCLUSTERITEM
{
    struct tagCLUSTERKEY* key;
    struct tagCLUSTERITEM* next;
}   CLUSTERITEM;


typedef struct tagCLUSTERHEAD
{
    int count;
    struct tagCLUSTERHEAD* prevcluster;
    struct tagCLUSTERHEAD* nextcluster;
    CLUSTERITEM* first;
    CLUSTERITEM* last;
    char *comment;
} CLUSTERHEAD;


typedef struct tagCLUSTERKEY
{
    int key;
    int direction;
    char leftcon,rightcon;
    PROFILE *prof,*profrc;
    struct tagCLUSTERKEY* next;
    CLUSTERHEAD* cluster;

} CLUSTERKEY;


struct tagCLUSTERBASE
{
    int magic; // set to CLUSTERMAGICNUM on valid structures
    CLUSTERKEY** keys;
    CLUSTERHEAD* firstcluster;
    CLUSTERHEAD* lastcluster;
    int size;
    unsigned long long int trcount;
};



CLUSTERKEY*     cl_findkey(CLUSTERBASE* base, int key);
CLUSTERHEAD*    cl_newcluster(CLUSTERBASE* base);
int     cl_insertkey(CLUSTERBASE* base, CLUSTERHEAD* cluster, int key, int dir);
int     cl_mergeclusters(CLUSTERBASE* base, CLUSTERHEAD* c1, CLUSTERHEAD* c2, int flip);





CLUSTERBASE*    ClusterBaseCreate(void)
{
    CLUSTERBASE* base;
    CLUSTERKEY** keys;

    // allocate a clusterbase structure
    base = (CLUSTERBASE*) malloc(sizeof(CLUSTERBASE));
    if(base==NULL) return NULL;

    // allocate the rack with 'zero-out' memory
    keys = (CLUSTERKEY**) calloc(CLUSTERBASESIZE, sizeof(CLUSTERKEY*));
    if(keys==NULL)
    {
        free(base);
        return NULL;
    }

    // assemble the parts
    base->magic = CLUSTERMAGICNUM;
    base->keys = keys;
    base->firstcluster = NULL;
    base->lastcluster = NULL;
    base->size = CLUSTERBASESIZE;
    base->trcount = 0;

    return base;
}

int     ClusterBaseAddComment(CLUSTERBASE* base, int key1, char *comment) {
      
    CLUSTERKEY *pk1;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;  
    
    pk1 = cl_findkey(base,key1);
    if(pk1==NULL) return -1;

    pk1->cluster->comment = strdup(comment);

    if (pk1->cluster->comment)
        return 0;
    else
        return -1;
}

// this function added to support 1 repeat clusters, Gelfand Jan 5, 2007
int		ClusterBaseAddUnique(CLUSTERBASE* base, int key1) {
   
    CLUSTERHEAD *cluster;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;

    cluster = cl_newcluster(base);
    if(cluster==NULL) return -1;
    if(cl_insertkey(base,cluster,key1,0)) return -1;

	return 0;
}

/* this function is called after PAM/CLARA to assign original directions inside clusterr */
int ClusterBaseRestoreDirection(CLUSTERBASE* base, int key, int direction) 
{

    CLUSTERKEY *pk1;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;

    pk1 = cl_findkey(base,key);
    if(pk1==NULL) 
      return -1;
    else
      pk1->direction = direction;
  
    return 0;
}

int     ClusterBaseAddLinkProfile(CLUSTERBASE* base, int key1 , int key2, PROFILE *prof1, PROFILE *prof1rc,PROFILE *prof2, PROFILE *prof2rc, int dir)
{
    CLUSTERKEY *pk1,*pk2;
    CLUSTERHEAD *cluster;
    int flip;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;

    // make sure keys are different
    if(key1==key2) return -1;

    pk1 = cl_findkey(base,key1);
    pk2 = cl_findkey(base,key2);

    // if both keys are new
    if(pk1==NULL && pk2==NULL)
    {
        cluster = cl_newcluster(base);
        if(cluster==NULL) return -1;
        if(cl_insertkey(base,cluster,key1,0)) return -1;
        if(cl_insertkey(base,cluster,key2,dir)) return -1;
    }
    // if both keys have been added before
    else if(pk1!=NULL && pk2!=NULL)
    {
        // make sure the clusters are different
        if(pk1->cluster!=pk2->cluster)
        {
            flip = pk1->direction ^ pk2->direction ^ dir;
            if(cl_mergeclusters(base,pk1->cluster,pk2->cluster,flip)) return -1;
        }
    }
    // if one has been added but not the other
    else
    {
        if(pk1!=NULL)
        {
            if(pk1->direction==0)
            {
                if(cl_insertkey(base,pk1->cluster,key2,dir)) return -1;
            }
            else
            {
                if(cl_insertkey(base,pk1->cluster,key2,!dir)) return -1;
            }
        }
        else
        {
            if(pk2->direction==0)
            {
                if(cl_insertkey(base,pk2->cluster,key1,dir)) return -1;
            }
            else
            {
                if(cl_insertkey(base,pk2->cluster,key1,!dir)) return -1;
            }
        }
    }

	// this is needed for pam
    pk1 = cl_findkey(base,key1);
    pk2 = cl_findkey(base,key2);
	pk1->prof = prof1;
	pk1->profrc = prof1rc;
	pk2->prof = prof2;
	pk2->profrc = prof2rc;

    return 0;
}

int     ClusterBaseAddLink(CLUSTERBASE* base, int key1 , int key2, int dir)
{
    CLUSTERKEY *pk1,*pk2;
    CLUSTERHEAD *cluster;
    int flip;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;

    // make sure keys are different
    if(key1==key2) return -1;

    pk1 = cl_findkey(base,key1);
    pk2 = cl_findkey(base,key2);

    // if both keys are new
    if(pk1==NULL && pk2==NULL)
    {
        cluster = cl_newcluster(base);
        if(cluster==NULL) return -1;
        if(cl_insertkey(base,cluster,key1,0)) return -1;
        if(cl_insertkey(base,cluster,key2,dir)) return -1;
    }
    // if both keys have been added before
    else if(pk1!=NULL && pk2!=NULL)
    {
        // make sure the clusters are different
        if(pk1->cluster!=pk2->cluster)
        {
            flip = pk1->direction ^ pk2->direction ^ dir;
            if(cl_mergeclusters(base,pk1->cluster,pk2->cluster,flip)) return -1;
        }
    }
    // if one has been added but not the other
    else
    {
        if(pk1!=NULL)
        {
            if(pk1->direction==0)
            {
                if(cl_insertkey(base,pk1->cluster,key2,dir)) return -1;
            }
            else
            {
                if(cl_insertkey(base,pk1->cluster,key2,!dir)) return -1;
            }
        }
        else
        {
            if(pk2->direction==0)
            {
                if(cl_insertkey(base,pk2->cluster,key1,dir)) return -1;
            }
            else
            {
                if(cl_insertkey(base,pk2->cluster,key1,!dir)) return -1;
            }
        }
    }

    return 0;
}

int     ClusterBaseAskLink(CLUSTERBASE* base, int key1, int key2)
{
    CLUSTERKEY *pk1,*pk2;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;

    pk1 = cl_findkey(base,key1);
    pk2 = cl_findkey(base,key2);

    // if both keys have been added before
    if(pk1!=NULL && pk2!=NULL)
    {
        // if both keys point to same cluster return 1
        if(pk1->cluster==pk2->cluster) return 1;
    }
    return 0;
}



CLUSTERKEY*     cl_findkey(CLUSTERBASE* base, int key)
{
    int pos;
    CLUSTERKEY* pk;

    // locate key in hash rack
    pos = abs(key)%base->size;
    for(pk=base->keys[pos];pk!=NULL;pk=pk->next) if(pk->key==key) break;

    return pk;
}


CLUSTERHEAD*    cl_newcluster(CLUSTERBASE* base)
{
    CLUSTERHEAD* cluster;

    // allocate cluster structure
    cluster = (CLUSTERHEAD*) malloc(sizeof(CLUSTERHEAD));
    if(cluster==NULL) return NULL;

    // set cluster members
    cluster->count = 0;
    cluster->prevcluster = NULL;
    cluster->nextcluster = NULL;
    cluster->first = NULL;
    cluster->last = NULL;
    cluster->comment = NULL;

    // place at end of cluster list
    if(base->firstcluster==NULL)
    {
        base->firstcluster = cluster;
        base->lastcluster = cluster;
    }
    else
    {
        base->lastcluster->nextcluster = cluster;
        cluster->prevcluster = base->lastcluster;
        base->lastcluster = cluster;
    }

    return cluster;
}

int     cl_insertkey(CLUSTERBASE* base, CLUSTERHEAD* cluster, int key, int dir)
{
    int pos;
    CLUSTERKEY* pk;
    CLUSTERITEM* pi;

    // compute position of key
    pos = abs(key)%base->size;

    // allocate a clusterkey structure
    pk = (CLUSTERKEY*) malloc(sizeof(CLUSTERKEY));
    if(pk==NULL) return -1;

    // allocate a clusteritem structure
    pi = (CLUSTERITEM*) malloc(sizeof(CLUSTERITEM));
    if(pi==NULL){free(pk);return -1;}

    // set clusterkey members and put in list
    pk->key = key;
    pk->direction = dir;
    pk->leftcon = 0;
    pk->rightcon = 0;
    pk->cluster = cluster;
    pk->next = base->keys[pos];
	pk->prof = NULL;
	pk->profrc = NULL;
    base->keys[pos] = pk;

    // set clusteritem member and put in list
    pi->key = pk;
    pi->next = NULL;
    if(cluster->first==NULL)
    {
        cluster->first = pi;
        cluster->last = pi;
    }
    else
    {
        cluster->last->next = pi;
        cluster->last = pi;
    }
    cluster->count++;
    base->trcount++;

    return 0;
}

int     cl_mergeclusters(CLUSTERBASE* base, CLUSTERHEAD* c1, CLUSTERHEAD* c2, int flip)
{
    CLUSTERHEAD *hold;
    CLUSTERITEM *item;

    // make c1 the one with more items
    if(c1->count<c2->count)
    {
        hold = c1;
        c1 = c2;
        c2 = hold;
    }

    // make all keys in c2 point to c1
    for(item = c2->first;item!=NULL;item=item->next)
    {
        item->key->cluster = c1;

        // also flip direction if flip is 1
        item->key->direction ^= flip;
    }

    // append c2's item list to c1
    c1->last->next = c2->first;
    c1->last = c2->last;
    c1->count += c2->count;

    // remove c2 from list and free it.
    if(c2==base->firstcluster)
    {
        base->firstcluster = c2->nextcluster;
        c2->nextcluster->prevcluster = NULL;
    }
    else if(c2==base->lastcluster)
    {
        base->lastcluster = c2->prevcluster;
        c2->prevcluster->nextcluster = NULL;
    }
    else //cluster is in the middle
    {
        c2->nextcluster->prevcluster = c2->prevcluster;
        c2->prevcluster->nextcluster = c2->nextcluster;
    }

    free(c2);

    return 0;
}


int     ClusterBasePrint(CLUSTERBASE* base, FILE* fp)
{
    CLUSTERHEAD* cluster;
    CLUSTERITEM* item;

    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;


    for(cluster=base->firstcluster;cluster!=NULL;cluster=cluster->nextcluster)
    {

       {

            fprintf(fp,"@");
            for(item=cluster->first;item!=NULL;item=item->next)
            {

		/* for ref to ref */
		if (OPTION == 'R') {

                  if(item->key->leftcon || item->key->rightcon)
                      fprintf(fp," %d+",item->key->key);	// UNDISTINGUISHABLE
                  else
                      fprintf(fp," %d-",item->key->key);	// DISTINGUISHABLE

		} else {

                  if(item->key->direction)
                      fprintf(fp," %d\"",item->key->key);
                  else
                      fprintf(fp," %d'",item->key->key);
		}

            }
            if (cluster->comment) 
                fprintf(fp," #%s\n", cluster->comment );
            else
                fprintf(fp,"\n");
        }
    }

    return 0;
}

int             ClusterBaseDestroy(CLUSTERBASE* base)
{
    CLUSTERHEAD *clus;
    CLUSTERITEM *item;


    // verify if valid cluster base
    if(base==NULL) return -1;
    if(base->magic!=CLUSTERMAGICNUM) return -1;


    // free clusterheads, clusteritems, and clusterkeys
    while(base->firstcluster!=NULL)
    {
        clus = base->firstcluster;
        base->firstcluster = clus->nextcluster;

        while(clus->first!=NULL)
        {
            item = clus->first;
            clus->first = item->next;

            free(item->key);
            free(item);
        }

        free(clus);
    }

    // free hash rack and base structure
    free(base->keys);
    base->magic = 0;
    free(base);

    return 0;
}

#endif



