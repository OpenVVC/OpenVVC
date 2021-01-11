#ifndef OV_MEMPOOL_H
#define OV_MEMPOOL_H
typedef struct MemPoolElem MemPoolElem;
typedef struct MemPool MemPool;


MemPool *ovmempool_init(size_t elem_size);


void ovmempool_uninit(MemPool **mpool_p);

#endif
