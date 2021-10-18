#include <string.h>

#include "mempool.h"
#include "mempool_internal.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"

#include "ovframepool.h"
#include "nvcl.h"
#include "nvcl_structures.h"
#include "ovframe.h"
#include "ovdpb.h"

void
dpbpriv_uninit_framepool(struct DPBInternal *dpb_priv)
{
    ovframepool_uninit(&dpb_priv->frame_pool);
}

int
dpbpriv_init_framepool(struct DPBInternal *dpb_priv, const OVSPS *const sps)
{
    int ret;
    ret = ovframepool_init(&dpb_priv->frame_pool, sps->sps_chroma_format_idc,
                           sps->sps_pic_width_max_in_luma_samples,
                           sps->sps_pic_height_max_in_luma_samples);
    if (ret < 0) {
        goto fail_init;
    }

    return 0;
fail_init:
    return ret;
}

int
dpbpriv_request_frame(struct DPBInternal *dpb_priv, OVFrame **frame)
{
    int ret;

    *frame = ov_mallocz(sizeof(**frame));
    if (!(*frame)) {
        return OVVC_ENOMEM;
    }

    ret = ovframepool_request_planes(*frame, dpb_priv->frame_pool);
    if (ret < 0) {
        goto failrequest;
    }

    return 0;

failrequest:
    ov_freep(frame);
    return OVVC_ENOMEM;
}

