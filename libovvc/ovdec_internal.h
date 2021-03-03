#ifndef OVDEC_INTERNAL_H
#define OVDEC_INTERNAL_H

/* Private decoder functions and structures declarations and generic
 */

#include "ovdefs.h"
#include "dec_structures.h"
#include "mempool.h"

struct MVPlane
{

    OVMV *mvs;
    uint64_t *dirs;

    /* Pool elems */
    void *dir_elem;
    void *mv_elem;
};

struct MVPool
{
    MemPool *dir_pool;
    MemPool *mv_pool;
    /* FIXME dimension info ?*/
};

struct PicPartInfo;

int mvpool_init(struct MVPool **mv_pool_p, const struct PicPartInfo *const pinfo);

void mvpool_uninit(struct MVPool **mv_pool_p);

int mvpool_request_mv_plane(struct MVPool *mv_pool, struct MVPlane *mv_plane);

void mvpool_release_mv_plane(struct MVPlane *mv_plane);

#endif
