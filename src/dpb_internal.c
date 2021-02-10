#include <string.h>

#include "mempool.h"
#include "mempool_internal.h"
#include "overror.h"
#include "ovmem.h"

#include "nvcl.h"
#include "nvcl_structures.h"
#include "ovframe.h"
#include "ovdpb.h"

struct DPBInternal
{
    struct FramePool frame_pool;
};

void
dpbpriv_uninit_framepool(struct DPBInternal *dpb_priv)
{
    const int nb_comp = 3;
    struct FramePool *fp = &dpb_priv->frame_pool;
    int i;
    for (i = 0; i < nb_comp; ++i) {
        if (fp->plane_pool[i]) {
            ovmempool_uninit(&fp->plane_pool[i]);
        }
    }

    memset(&fp->plane_prop, 0, sizeof(fp->plane_prop));
}

int
dpbpriv_init_framepool(struct DPBInternal *dpb_priv, const OVSPS *const sps)
{
    static const uint8_t comp_shift[3] = {0, 2, 2};
    size_t pic_w = (size_t) sps->sps_pic_width_max_in_luma_samples;
    size_t pic_h = (size_t) sps->sps_pic_height_max_in_luma_samples;
    uint8_t bd_shift = !!sps->sps_bitdepth_minus8;
    /* TODO non 420 chromat_formats */
    uint8_t nb_comp = 3;
    struct FramePool *fp = &dpb_priv->frame_pool;

    /* FIXME allocation size overflow */
    size_t pic_size = (pic_w * pic_h) << bd_shift;

    int i;

    for (i = 0; i < nb_comp; ++i) {
        size_t elem_size = pic_size >> comp_shift[i];

        fp->plane_pool[i] = ovmempool_init(elem_size);

        if (!fp->plane_pool[i]) {
            goto fail_poolinit;
        }

        fp->plane_prop[i].stride = (pic_w << bd_shift) >> comp_shift[i];
        fp->plane_prop[i].height = pic_h >> comp_shift[i];
        fp->plane_prop[i].width  = pic_w >> comp_shift[i];
        fp->plane_prop[i].depth  = bd_shift;
    }

    return 0;

fail_poolinit :
    dpbpriv_uninit_framepool(dpb_priv);
    return OVVC_ENOMEM;
}

int
dpbpriv_request_frame(struct DPBInternal *dpb_priv, OVFrame **frame)
{
    int ret;

    *frame = ov_mallocz(sizeof(**frame));
    if (!(*frame)) {
        return OVVC_ENOMEM;
    }

    ret = framepool_request_planes(*frame, &dpb_priv->frame_pool);
    if (ret < 0) {
        goto failrequest;
    }

    /* Copy frame_info */


    return 0;

failrequest:
    ov_freep(frame);
    return OVVC_ENOMEM;
}
