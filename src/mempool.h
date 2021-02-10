#ifndef OV_MEMPOOL_H
#define OV_MEMPOOL_H

#include <stddef.h>

typedef struct MemPoolElem MemPoolElem;
typedef struct MemPool MemPool;


MemPool *ovmempool_init(size_t elem_size);

MemPoolElem *ovmempool_popelem(MemPool *mpool);

void ovmempool_pushelem(MemPoolElem *released_elem);

void ovmempool_uninit(MemPool **mpool_p);

#endif
