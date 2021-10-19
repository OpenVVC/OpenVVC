#include "mempool.h"
#include "mempool_internal.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"

#include "ovframe.h"
#include "ovframepool.h"

int
ovframe_new_ref(OVFrame **dst, OVFrame *src)
{
     if (!src) {
         return -1;
     }

    unsigned ref_count = atomic_fetch_add_explicit(&src->internal.ref_count, 1, memory_order_acq_rel);
    ov_log(NULL, OVLOG_TRACE, "NewRef Frame %ld ref_count: %d\n", src, ref_count);

    *dst = src;

    return 0;
}

void
ovframe_unref(OVFrame **frame_p)
{
    if (!frame_p)
        return;

    if (!*frame_p){
        ov_log(NULL, OVLOG_ERROR, "Trying to unref NULL frame\n");
        return;
    }

    unsigned ref_count = atomic_fetch_add_explicit(&(*frame_p)->internal.ref_count, -1, memory_order_acq_rel);
    ov_log(NULL, OVLOG_TRACE, "Unref Frame %ld ref_count: %d\n", *frame_p, ref_count);

    if (!ref_count) {
        ovframepool_release_frame(frame_p);
    }
}
