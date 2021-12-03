#include "mempool_internal.h"
#include "mempool.h"
#include "ovmem.h"
#include "ovframe.h"
#include "overror.h"
#include "ovutils.h"
#include "ovframepool.h"

struct PlaneProp
{
    uint16_t stride;
    uint16_t width;
    uint16_t height;
    uint16_t depth;
};

struct FramePool
{
    struct MemPool *frame_pool;
    struct MemPool *plane_pool[4];
    struct PlaneProp plane_prop[4];
    const struct ChromaFmtInfo *fmt_info;
};

struct FrameProperties
{
   uint8_t chromat_format; 
};

struct ChromaFmtInfo
{
   uint8_t nb_comp;
   uint8_t bd_shift;
   uint8_t shift_h[3];
   uint8_t shift_v[3];
};

static const struct ChromaFmtInfo yuv420_8 = {
  .nb_comp = 3,
  .bd_shift = 0,
  .shift_h = {0, 1, 1},
  .shift_v = {0, 1, 1},
};

static const struct ChromaFmtInfo yuv420_10 = {
  .nb_comp = 3,
  .bd_shift = 1,
  .shift_h = {0, 1, 1},
  .shift_v = {0, 1, 1},
};

static const struct ChromaFmtInfo *const select_frame_format(enum ChromaFmt fmt, uint8_t bitdepth_min8)
{
    if (bitdepth_min8) {
        return &yuv420_10;
    } else {
        return &yuv420_8;
    }
}

static void set_plane_properties(struct PlaneProp *const pln, const struct ChromaFmtInfo *const fmt_info,
                                 uint8_t comp_idx, uint16_t pic_w, uint16_t pic_h)
{
    pln->stride = (pic_w << fmt_info->bd_shift) >> fmt_info->shift_v[comp_idx];
    pln->height = pic_h >> fmt_info->shift_v[comp_idx];
    pln->width  = pic_w >> fmt_info->shift_h[comp_idx];
    pln->depth  = fmt_info->bd_shift;
}

void
ovframepool_uninit(struct FramePool **fpool_p)
{
    struct FramePool *fpool = *fpool_p;
    const struct  ChromaFmtInfo *const fmt_info = fpool->fmt_info;
    int i;
    for (i = 0; i < fmt_info->nb_comp; ++i) {
        if (fpool->plane_pool[i]) {
            ovmempool_uninit(&fpool->plane_pool[i]);
        }
    }

    if (fpool->frame_pool) {
        ovmempool_uninit(&fpool->frame_pool);
    }

    ov_freep(fpool_p);
}

int
ovframepool_init(struct FramePool **fpool_p, uint8_t fmt, uint8_t bitdepth_min8, uint16_t pic_w, uint16_t pic_h)
{
    const struct ChromaFmtInfo *const fmt_info = select_frame_format(fmt, bitdepth_min8);

    /* FIXME allocation size overflow */
    size_t pic_size = (pic_w * pic_h) << fmt_info->bd_shift;
    struct FramePool *fpool;

    int i;

    *fpool_p = ov_mallocz(sizeof(struct FramePool));
    if (!*fpool_p) {
        goto fail_alloc;
    }

    fpool = *fpool_p;

    fpool->fmt_info = fmt_info;

    fpool->frame_pool = ovmempool_init(sizeof(OVFrame));
    if (!fpool->frame_pool) {
        goto fail_poolinit;
    }

    for (i = 0; i < fmt_info->nb_comp; ++i) {
        uint8_t comp_shift = fmt_info->shift_h[i] + fmt_info->shift_v[i];
        size_t elem_size = pic_size >> comp_shift;

        fpool->plane_pool[i] = ovmempool_init(elem_size);

        if (!fpool->plane_pool[i]) {
            goto fail_poolinit;
        }

        set_plane_properties(&fpool->plane_prop[i], fmt_info, i, pic_w, pic_h);
    }

    return 0;

fail_poolinit :
    ov_log(NULL, OVLOG_ERROR, "Failed frame pool alloc\n");
    ovframepool_uninit(fpool_p);
    return OVVC_ENOMEM;
fail_alloc:
     *fpool_p = NULL;
    return OVVC_ENOMEM;
}

static void
ovframepool_release_planes(OVFrame *const frame)
{
    const int nb_comp = frame->internal.frame_pool->fmt_info->nb_comp;
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

static int
ovframepool_request_planes(OVFrame *const frame, struct FramePool *const fpool)
{
    const int nb_comp = fpool->fmt_info->nb_comp;
    int i;

    frame->internal.frame_pool = fpool;

    for (i = 0; i < nb_comp; ++i) {
        MemPool *pool = fpool->plane_pool[i];
        struct MemPoolElem *pool_elem;
        struct PlaneProp *prop = &fpool->plane_prop[i];
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

    atomic_init(&frame->internal.ref_count, 0);

    return 0;

failpop:
    ovframepool_release_planes(frame);
    return OVVC_ENOMEM;
}

static OVFrame *
create_new_frame(struct FramePool *fpool)
{
    OVFrame *frame;
    int ret;
    struct MemPoolElem *felem = ovmempool_popelem(fpool->frame_pool);
    if (!felem) {
        return NULL;
    }

    frame = felem->data;

    frame->internal.felem = felem;

    ret = ovframepool_request_planes(frame, fpool);
    if (ret < 0) {
        goto failrequest;
    }

    return frame;

failrequest:
    ovmempool_pushelem(felem);
    return NULL;
}

OVFrame *
ovframepool_request_frame(struct FramePool *fpool)
{
    OVFrame *frame = create_new_frame(fpool);
    return frame;
}

void
ovframepool_release_frame(OVFrame **frame_p)
{
    ovframepool_release_planes(*frame_p);
    ovmempool_pushelem((*frame_p)->internal.felem);
    *frame_p = NULL;
}

