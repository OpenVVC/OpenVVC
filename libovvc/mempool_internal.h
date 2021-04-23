#ifndef OV_MEMPOOL_INTERNAL_H
#define OV_MEMPOOL_INTERNAL_H
#include <stddef.h>
#include <pthread.h>

struct MemPoolElem
{
    struct MemPool *mempool;
    void *data;
    struct MemPoolElem *next_elem;
};

struct MemPool
{
    struct MemPoolElem *stack_elem;
    size_t elem_size;
    int nb_ref;
    pthread_mutex_t pool_mtx;
};

#endif
