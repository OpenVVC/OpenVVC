#include "mempool.h"
#include "mempool_internal.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"

#include "ovframe.h"

static void
framepool_release_planes(OVFrame *const frame)
{
    const int nb_comp = 3;
    int i;

    for (i = 0; i < nb_comp; ++i) {
        struct MemPoolElem *pool_elem;

        pool_elem = frame->internal.pool_elem[i];

        if (pool_elem) {
            ovmempool_pushelem(pool_elem);
        }

        frame->internal.pool_elem[i] = NULL;
        frame->data[i] = NULL;
    }
}

int
framepool_request_planes(OVFrame *const frame, struct FramePool *const fp)
{
    const int nb_comp = 3;
    int i;

    for (i = 0; i < nb_comp; ++i) {
        MemPool *pool = fp->plane_pool[i];
        struct MemPoolElem *pool_elem;
        struct PlaneProp *prop = &fp->plane_prop[i];
        pool_elem = ovmempool_popelem(pool);
        if (!pool_elem) {
             goto failpop;
        }

        frame->internal.pool_elem[i] = pool_elem;

        frame->data[i]     = pool_elem->data;

        frame->width[i]    = prop->width;
        frame->height[i]   = prop->height;
        frame->linesize[i] = prop->stride;
    }

    frame->internal.ref_count = 1;

    return 0;

failpop:
    framepool_release_planes(frame);
    return OVVC_ENOMEM;
}

int
ovframe_new_ref(OVFrame **dst, OVFrame *src)
{
     if (!src) {
         return -1;
     }

     src->internal.ref_count++;

     *dst = src;

     return 0;
}

void
ovframe_unref(OVFrame **frame)
{
    if (!*frame){
        ov_log(NULL, OVLOG_ERROR, "Trying to unref NULL frame\n");
        return;
    }

    (*frame)->internal.ref_count--;

    if (!(*frame)->internal.ref_count) {
        framepool_release_planes(*frame);
        ov_freep(frame);
    }

    *frame = NULL;
}
