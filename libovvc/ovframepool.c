#include "mempool_internal.h"
#include "mempool.h"
#include "ovmem.h"
#include "ovframe.h"
#include "overror.h"
#include "ovutils.h"
#include "ovframepool.h"

struct FramePool
{
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

static const struct ChromaFmtInfo yuv420_10 = {
  .nb_comp = 3,
  .bd_shift = 1,
  .shift_h = {0, 1, 1},
  .shift_v = {0, 1, 1},
};

static const struct ChromaFmtInfo *const select_frame_format(enum ChromaFmt fmt)
{
    return &yuv420_10;
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
    struct FramePool *fp = *fpool_p;
    const struct  ChromaFmtInfo *const fmt_info = fp->fmt_info;
    int i;
    for (i = 0; i < fmt_info->nb_comp; ++i) {
        if (fp->plane_pool[i]) {
            ovmempool_uninit(&fp->plane_pool[i]);
        }
    }

    ov_freep(fpool_p);
}

int
ovframepool_init(struct FramePool **fpool_p, uint8_t fmt, uint16_t pic_w, uint16_t pic_h)
{
    const struct ChromaFmtInfo *const fmt_info = select_frame_format(fmt);

    /* FIXME allocation size overflow */
    size_t pic_size = (pic_w * pic_h) << fmt_info->bd_shift;
    struct FramePool *fp;

    int i;

    *fpool_p = ov_mallocz(sizeof(struct FramePool));
    if (!*fpool_p) {
        goto fail_alloc;
    }

    fp = *fpool_p;

    fp->fmt_info = fmt_info;

    for (i = 0; i < fmt_info->nb_comp; ++i) {
        uint8_t comp_shift = fmt_info->shift_h[i] + fmt_info->shift_v[i];
        size_t elem_size = pic_size >> comp_shift;

        fp->plane_pool[i] = ovmempool_init(elem_size);

        if (!fp->plane_pool[i]) {
            goto fail_poolinit;
        }

        set_plane_properties(&fp->plane_prop[i], fmt_info, i, pic_w, pic_h);
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

void
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

int
ovframepool_request_planes(OVFrame *const frame, struct FramePool *const fp)
{
    const int nb_comp = fp->fmt_info->nb_comp;
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

    atomic_init(&frame->internal.ref_count, 0);

    return 0;

failpop:
    ovframepool_release_planes(frame);
    return OVVC_ENOMEM;
}

