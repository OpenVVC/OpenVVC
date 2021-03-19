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

failalloc:
    return mpool;
}

MemPoolElem *
ovmempool_popelem(MemPool *mpool)
{
    MemPoolElem *elem = mpool->stack_elem;

    if (elem) {
        mpool->stack_elem = elem->next_elem;
        elem->next_elem = NULL;
    } else {
        elem = ov_mallocz(sizeof(*elem));
        if (!elem) {
            return elem;
        }
        /* Keep track of parent pool so elem can be released
           without knowledge of responsible mempool */
        elem->mempool = mpool;

        elem->data = ov_mallocz(mpool->elem_size);
        if (!elem->data) {
            ov_freep(&elem);
            return elem;
        }
    }

    /* Keep track of ref in use so we avoid freeing the
       mempool if some of its elements did not return */
    mpool->nb_ref++;

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

       released_elem->next_elem = mpool->stack_elem;
       mpool->stack_elem = released_elem;
       mpool->nb_ref--;
       if (!mpool->nb_ref) {
           ovmempool_free(mpool);
       }
   }
}

void
ovmempool_uninit(MemPool **mpool_p)
{
    MemPool *mpool = *mpool_p;

    mpool->nb_ref--;

    if (!mpool->nb_ref) {
        ovmempool_free(mpool);
    }

    *mpool_p = NULL;
}
