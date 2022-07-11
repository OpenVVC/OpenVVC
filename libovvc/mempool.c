/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
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
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include "ovmem.h"

#include "mempool_internal.h"
#include "mempool.h"


static void ovmempool_free(MemPool *mpool);


MemPool *
ovmempool_init(size_t elem_size)
{
    struct MemPool *mpool = ov_mallocz(sizeof(*mpool));
    if (!mpool) {
        goto failalloc;
    }

    mpool->elem_size = elem_size;

    /* The pool keeps a ref to itself so we avoid freeing it
       while some of its elements can still point to it
       the pool will be freed only when all allocated elements
       have returned to it and ovmem_pool_uninit() has been called*/
    mpool->nb_ref = 1;

    pthread_mutex_init(&mpool->pool_mtx, NULL);
failalloc:
    return mpool;
}

MemPoolElem *
ovmempool_popelem(MemPool *mpool)
{
    pthread_mutex_lock(&mpool->pool_mtx);

    MemPoolElem *elem = mpool->stack_elem;
    if (elem) {
        mpool->stack_elem = elem->next_elem;
        elem->next_elem = NULL;
    } else {
        elem = ov_mallocz(sizeof(*elem));
        if (!elem) {
            goto ret_and_unlock;
        }
        /* Keep track of parent pool so elem can be released
           without knowledge of responsible mempool */
        elem->mempool = mpool;

        elem->data = ov_mallocz(mpool->elem_size);
        if (!elem->data) {
            ov_freep(&elem);
            goto ret_and_unlock;
        }
    }

    /* Keep track of ref in use so we avoid freeing the
       mempool if some of its elements did not return */
    mpool->nb_ref++;

ret_and_unlock:
    pthread_mutex_unlock(&mpool->pool_mtx);
    return elem;
}

static void
ovmempool_free(MemPool *mpool)
{
    MemPoolElem *stack_elem = mpool->stack_elem;

    while (stack_elem) {
        /* Free all elements directly available in the pool */
        MemPoolElem *next_to_free = stack_elem->next_elem;

        ov_free(stack_elem->data);
        ov_free(stack_elem);

        stack_elem = next_to_free;
    }

    ov_freep(&mpool);
}

void
ovmempool_pushelem(MemPoolElem *released_elem)
{  
    if (released_elem) {
        MemPool *mpool = released_elem->mempool;
        pthread_mutex_lock(&mpool->pool_mtx);

        released_elem->next_elem = mpool->stack_elem;
        mpool->stack_elem = released_elem;
        mpool->nb_ref--;
        if (!mpool->nb_ref) {
            pthread_mutex_unlock(&mpool->pool_mtx);
            ovmempool_free(mpool);
            return;
        }
        pthread_mutex_unlock(&mpool->pool_mtx);
    }
}

void
ovmempool_uninit(MemPool **mpool_p)
{
    MemPool *mpool = *mpool_p;

    mpool->nb_ref--;

    pthread_mutex_destroy(&mpool->pool_mtx);

    if (!mpool->nb_ref) {
        ovmempool_free(mpool);
    }

    *mpool_p = NULL;
}
