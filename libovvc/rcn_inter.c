#include <string.h>
#include <stdlib.h>
#include "ovdpb.h"
#include "dec_structures.h"
#include "ovdefs.h"
#include "ovutils.h"
#include "ctudec.h"
#include "rcn_structures.h"
#include "rcn.h"
#include "rcn_mc.h"
#include "drv.h"
#include "rcn.h"
#include "drv_utils.h"

#define MAX_PB_SIZE 128

#define REF_PADDING_C 1
#define EPEL_EXTRA_AFTER  2
#define EPEL_EXTRA REF_PADDING_C + EPEL_EXTRA_AFTER

#define REF_PADDING_L 3
#define QPEL_EXTRA_BEFORE 3
#define QPEL_EXTRA_AFTER  4
#define QPEL_EXTRA REF_PADDING_L + QPEL_EXTRA_AFTER

#define GRAD_SHIFT 6

/* Link to GRAD SHIFT */
#define PROF_DMV_MAX ((1 << 5) - 1)

#define PROF_MV_SHIFT 8
#define PROF_MV_RND (1 << (PROF_MV_SHIFT - 1))

#define BITDEPTH 10
#define PROF_SMP_SHIFT (14 - BITDEPTH)
#define PROF_SMP_RND (1 << (14 - 1))
#define PROF_SMP_OFFSET (1 << (PROF_SMP_SHIFT + 1) - 1)

#define PROF_BUFF_PADD_H 1
#define PROF_BUFF_PADD_W 1
#define SB_H 4
#define SB_W 4

#define BDOF_WGT_LIMIT ((1 << 4) - 1)
#define BDOF_SHIFT   (14 + 1 - BITDEPTH)
#define BDOF_OFFSET  ((1 << (BDOF_SHIFT - 1)))

static const int8_t dmvr_mv_x[25 + 25] = {
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    -2, -1, 0, 1, 2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
    2, 1, 0, -1, -2,
};

static const int8_t dmvr_mv_y[25 + 25] = {
    -2, -2, -2, -2, -2,
    -1, -1, -1, -1, -1,
     0,  0,  0,  0,  0,
     1,  1,  1,  1,  1,
     2,  2,  2,  2,  2,
     2,  2,  2,  2,  2,
     1,  1,  1,  1,  1,
     0,  0,  0,  0,  0,
    -1, -1, -1, -1, -1,
    -2, -2, -2, -2, -2,
};
struct OVDMV {
    int32_t x;
    int32_t y;
};

#define BIT_DEPTH 10
#define ov_clip_pixel(a) ov_clip_uintp2(a, BIT_DEPTH)

enum CUMode {
    OV_NA = 0xFF,
    OV_INTER = 1,
    OV_INTRA = 2,
    OV_INTER_SKIP = 3,
    OV_MIP = 4,
};
static void
rcn_bdof(int16_t *dst, int dst_stride,
         const int16_t *ref_bdof0, const int16_t *ref_bdof1, int ref_stride,
         const int16_t *grad_x0, const int16_t *grad_y0,
         const int16_t *grad_x1, const int16_t *grad_y1,
         int grad_stride, uint8_t pb_w, uint8_t pb_h);

static OVMV
clip_mv(int pos_x, int pos_y, int pic_w, int pic_h, int pb_w, int pb_h, OVMV mv)
{
    int x_max  = (pic_w + 2 - pos_x) << 4;
    int y_max  = (pic_h + 2 - pos_y) << 4;
    int x_min  = (-pb_w - 3 - pos_x) << 4;
    int y_min  = (-pb_h - 3 - pos_y) << 4;

    mv.x = ov_clip(mv.x, x_min, x_max);
    mv.y = ov_clip(mv.y, y_min, y_max);

    return mv;
}

void
compute_prof_grad(const uint16_t* src, int src_stride, int sb_w, int sb_h,
                  int grad_stride, int16_t* grad_x, int16_t* grad_y)
{
    int y, x;
    const int nb_smp_h = sb_h;
    const int nb_smp_w = sb_w;

    src += src_stride + 1;

    for (y = 0; y < nb_smp_h; ++y) {
        for (x = 0; x < nb_smp_w; ++x) {
            grad_y[x]  = (((int16_t)src[x + src_stride] - (1 << 13)) >> GRAD_SHIFT);
            grad_y[x] -= (((int16_t)src[x - src_stride] - (1 << 13)) >> GRAD_SHIFT);
            grad_x[x]  = (((int16_t)src[x + 1] - (1 << 13)) >> GRAD_SHIFT);
            grad_x[x] -= (((int16_t)src[x - 1] - (1 << 13)) >> GRAD_SHIFT);
        }
        grad_x += grad_stride;
        grad_y += grad_stride;
        src += src_stride;
    }
}

#define PROF_DELTA_LIMIT (1 << (BITDEPTH + 3))
void
rcn_prof(uint16_t* dst, int dst_stride, const uint16_t* src, int src_stride,
         const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
         const int32_t* dmv_scale_h, const int32_t* dmv_scale_v,
         uint8_t bidir)
{
    int idx = 0;
    int x, y;

    for (y = 0; y < SB_H; ++y) {
        for (x = 0; x < SB_W; ++x) {
            int32_t add = dmv_scale_h[idx] * grad_x[x] + dmv_scale_v[idx] * grad_y[x];
            int16_t val;

            add = ov_clip(add, -PROF_DELTA_LIMIT, PROF_DELTA_LIMIT - 1);

            val = (int16_t)src[x] -(1 << 13);
            val += add;

            /* Clipping if not bi directional */
            if (!bidir) {
                val = (val + 8200 /*+ PROF_SMP_OFFSET*/) >> PROF_SMP_SHIFT;
                dst[x] = ov_clip(val, 0, 1023);
            } else {
                dst[x] = val + (1 << 13);
            }

            idx++;
        }

        grad_x += grad_stride;
        grad_y += grad_stride;

        dst += dst_stride;
        src += src_stride;
    }
}

/* FIXME check edge_emulation OK */

static void
extend_bdof_buff(const uint16_t *const src, uint16_t *dst_prof,
                 int16_t ref_stride, int16_t pb_w, int16_t pb_h,
                 uint8_t ext_x, uint8_t ext_y)
{
    int16_t tmp_prof_stride = 128;
    const uint16_t *ref = src  - ref_stride  - 1;

    uint16_t     *dst = dst_prof;
    uint16_t *dst_lst = dst_prof + (pb_h + 1) * tmp_prof_stride;
    int i, j;

    /* Position ref according to precision */
    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    const uint16_t *ref_lst = ref + (pb_h + 1) * ref_stride;

    /* Copy or extend upper and lower ref_line */
    for (i = 0; i < pb_w + 2; ++i) {
        dst[i]     = (ref[i]     << PROF_SMP_SHIFT);
        dst_lst[i] = (ref_lst[i] << PROF_SMP_SHIFT);
    }

    dst += tmp_prof_stride;
    dst_lst = dst + pb_w + 1;

    ref = src - 1;

    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    ref_lst = ref + pb_w + 1;

    /* Copy or extend left and right column*/
    for (j = 0; j < pb_h; ++j) {
        dst[0]     = (ref[0]     << PROF_SMP_SHIFT);
        dst_lst[0] = (ref_lst[0] << PROF_SMP_SHIFT);

        ref     += ref_stride;
        ref_lst += ref_stride;
        dst     += tmp_prof_stride;
        dst_lst += tmp_prof_stride;
    }
}

static void
extend_prof_buff(const uint16_t *const src, uint16_t *dst_prof, int16_t ref_stride, uint8_t ext_x, uint8_t ext_y)
{
    int16_t tmp_prof_stride = (128);
    const uint16_t *ref = src  - ref_stride  - 1;
    uint16_t       *dst = dst_prof;
    uint16_t *dst_lst = dst_prof + (SB_H + 1) * tmp_prof_stride;
    int i, j;

    /* Position ref according to precision */
    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    const uint16_t *ref_lst = ref + (SB_H + 1) * ref_stride;

    /* Copy or extend upper and lower ref_line */
    for (i = 0; i < SB_W + 2; ++i) {
        dst[i]     = (ref[i]     << PROF_SMP_SHIFT);
        dst_lst[i] = (ref_lst[i] << PROF_SMP_SHIFT);
    }

    dst += tmp_prof_stride;
    dst_lst = dst + SB_W + 1;

    ref = src - 1;

    if (ext_x) {
        ref += 1;
    }

    if (ext_y) {
        ref += ref_stride;
    }

    ref_lst = ref + SB_W + 1;

    /* Copy or extend left and right column*/
    for (j = 0; j < SB_H; ++j) {
        dst[0]     = (ref[0]     << PROF_SMP_SHIFT);
        dst_lst[0] = (ref_lst[0] << PROF_SMP_SHIFT);
        #if 0
        for (i = 1; i< 5; ++i) {
            dst[i] = dst[i];
        }
        #endif

        ref += ref_stride;
        ref_lst += ref_stride;
        dst += tmp_prof_stride;
        dst_lst += tmp_prof_stride;
    }
}

static void
extend_bdof_grad(uint16_t *dst_grad, int16_t grad_stride, int16_t pb_w, int16_t pb_h)
{
    uint16_t       *dst = dst_grad + grad_stride;
    const uint16_t *ref = dst + 1;

    uint16_t       *dst_lst = (uint16_t*)ref + pb_w;
    const uint16_t *ref_lst = dst_lst - 1;

    int j;

    /* Copy or extend left and right column*/
    for (j = 0; j < pb_h; ++j) {
        #if 1
        dst[0]     = ref[0];
        dst_lst[0] = ref_lst[0];
        #else
        //dst[0] = 0;
        dst[0]     = ref[0];
        dst_lst[0] = 0;
        #endif

        ref += grad_stride;
        dst += grad_stride;
        ref_lst += grad_stride;
        dst_lst += grad_stride;
    }

    /* Copy or extend upper and lower ref_line */
    dst = dst_grad;
    ref = dst + grad_stride;
    ref_lst = dst + (pb_h) * grad_stride;
    dst_lst = (uint16_t*)ref_lst + grad_stride;

    memcpy(dst,     ref    , sizeof(*ref) * (pb_w + 2));
    memcpy(dst_lst, ref_lst, sizeof(*ref) * (pb_w + 2));
}

static void
emulate_block_border(uint16_t *buf, const uint16_t *src,
                     ptrdiff_t buf_linesize,
                     ptrdiff_t src_linesize,
                     int block_w, int block_h,
                     int src_x, int src_y, int w, int h)
{
    int x, y;
    int start_y, start_x, end_y, end_x;

    if (!w || !h)
        return;

    if (src_y >= h) {
        src -= src_y * src_linesize;
        src += (h - 1) * src_linesize;
        src_y = h - 1;
    } else if (src_y <= -block_h) {
        src -= src_y * src_linesize;
        src += (1 - block_h) * src_linesize;
        src_y = 1 - block_h;
    }

    if (src_x >= w) {
        src  += (w - 1 - src_x);
        src_x = w - 1;
    } else if (src_x <= -block_w) {
        src  += (1 - block_w - src_x);
        src_x = 1 - block_w;
    }

    start_y = OVMAX(0, -src_y);
    start_x = OVMAX(0, -src_x);

    end_y = OVMIN(block_h, h-src_y);
    end_x = OVMIN(block_w, w-src_x);

    w    = end_x - start_x;
    src += start_y * src_linesize + start_x;
    buf += start_x;

    // top
    for (y = 0; y < start_y; y++) {
        memcpy(buf, src, w * sizeof(uint16_t));
        buf += buf_linesize;
    }

    // copy existing part
    for (; y < end_y; y++) {
        memcpy(buf, src, w * sizeof(uint16_t));
        src += src_linesize;
        buf += buf_linesize;
    }

    // bottom
    src -= src_linesize;
    for (; y < block_h; y++) {
        memcpy(buf, src, w * sizeof(uint16_t));
        buf += buf_linesize;
    }

    buf -= block_h * buf_linesize + start_x;

    while (block_h--) {
        uint16_t *bufp = (uint16_t *) buf;

        // left
        for(x = 0; x < start_x; x++) {
            bufp[x] = bufp[start_x];
        }

        // right
        for (x = end_x; x < block_w; x++) {
            bufp[x] = bufp[end_x - 1];
        }
        buf += buf_linesize;
    }
}

static uint8_t
test_for_edge_emulation_c(int pb_x, int pb_y, int pic_w, int pic_h,
                          int pb_w, int pb_h)
{
    uint8_t emulate_edge = 0;
    emulate_edge  =      pb_x - REF_PADDING_C < 0;
    emulate_edge |= 2 * (pb_y - REF_PADDING_C < 0);
    emulate_edge |= 4 * (pb_x >= pic_w);
    emulate_edge |= 8 * (pb_y >= pic_h);
    emulate_edge |= 4 * ((pb_x + pb_w + EPEL_EXTRA_AFTER) >= pic_w);
    emulate_edge |= 8 * ((pb_y + pb_h + EPEL_EXTRA_AFTER) >= pic_h);
    return emulate_edge;
}

static uint8_t
test_for_edge_emulation(int pb_x, int pb_y, int pic_w, int pic_h,
                        int pu_w, int pu_h)
{
    /* FIXME thi could be simplified */
    uint8_t emulate_edge = 0;
    emulate_edge =       pb_x - REF_PADDING_L < 0;
    emulate_edge |= 2 * (pb_y - REF_PADDING_L < 0);

    emulate_edge |= 4 * (pb_x >= pic_w);
    emulate_edge |= 8 * (pb_y >= pic_h);

    emulate_edge |= 4 * ((pb_x + pu_w + QPEL_EXTRA_AFTER) >= pic_w);
    emulate_edge |= 8 * ((pb_y + pu_h + QPEL_EXTRA_AFTER) >= pic_h);
    return emulate_edge;
}

static struct OVBuffInfo
derive_ref_buf_c(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                 uint16_t *edge_buff0, uint16_t *edge_buff1,
                 int log2_pu_w, int log2_pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    uint16_t *const ref_cb  = (uint16_t *) ref_pic->frame->data[1];
    uint16_t *const ref_cr  = (uint16_t *) ref_pic->frame->data[2];

    int src_stride = ref_pic->frame->linesize[1] >> 1;

    /*FIXME check buff side derivation */
    int ref_pos_x = pos_x + (mv.x >> 5);
    int ref_pos_y = pos_y + (mv.y >> 5);

    const int pu_w = (1 << log2_pu_w) >> 1;
    const int pu_h = (1 << log2_pu_h) >> 1;

    const int pic_w = ref_pic->frame->width[0]  >> 1;
    const int pic_h = ref_pic->frame->height[0] >> 1;

    uint16_t *src_cb  = &ref_cb[ref_pos_x + ref_pos_y * src_stride];
    uint16_t *src_cr  = &ref_cr[ref_pos_x + ref_pos_y * src_stride];

    uint8_t emulate_edge = test_for_edge_emulation_c(ref_pos_x, ref_pos_y, pic_w, pic_h,
                                                     pu_w, pu_h);;

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (src_stride) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        int cpy_w = pu_w + EPEL_EXTRA;
        int cpy_h = pu_h + EPEL_EXTRA;
        /*FIXME clip to frame?*/
        int start_pos_x = ref_pos_x - REF_PADDING_C;
        int start_pos_y = ref_pos_y - REF_PADDING_C;

        emulate_block_border(edge_buff0, (src_cb - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             cpy_w, cpy_h, start_pos_x, start_pos_y,
                             pic_w, pic_h);
        emulate_block_border(edge_buff1, (src_cr - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             cpy_w, cpy_h, start_pos_x, start_pos_y,
                             pic_w, pic_h);

        ref_buff.cb = edge_buff0 + buff_off;
        ref_buff.cr = edge_buff1 + buff_off;
        ref_buff.stride_c = RCN_CTB_STRIDE;

    } else {
        ref_buff.cb = src_cb;
        ref_buff.cr = src_cr;
        ref_buff.stride_c = src_stride;
    }
    return ref_buff;
}

static void
padd_dmvr_c(int16_t *const _ref, int16_t stride, uint8_t pu_w, uint8_t pu_h)
{
    int i;

    int16_t *ref = _ref;
    ref -= REF_PADDING_C + REF_PADDING_C * stride;

    for (i = 0; i < pu_h + EPEL_EXTRA; ++i) {
        ref[-1] = ref[0];
        ref[-2] = ref[0];
        ref[pu_w + EPEL_EXTRA    ] = ref[pu_w + EPEL_EXTRA - 1];
        ref[pu_w + EPEL_EXTRA + 1] = ref[pu_w + EPEL_EXTRA - 1];
        ref += stride;
    }

    ref  = _ref;
    ref -= REF_PADDING_C + REF_PADDING_C * stride + 2;
    memcpy(ref -     stride, ref, sizeof(*ref) * (pu_w + EPEL_EXTRA + 4));
    memcpy(ref - 2 * stride, ref, sizeof(*ref) * (pu_w + EPEL_EXTRA + 4));

    ref += stride * (pu_h + EPEL_EXTRA - 1);

    memcpy(ref + 1 * stride, ref, sizeof(*ref) * (pu_w + EPEL_EXTRA + 4));
    memcpy(ref + 2 * stride, ref, sizeof(*ref) * (pu_w + EPEL_EXTRA + 4));
}

static void
padd_dmvr(int16_t *const _ref, int16_t stride, uint8_t pu_w, uint8_t pu_h)
{
    int i;

    int16_t *ref = _ref;
    ref -= QPEL_EXTRA_BEFORE + QPEL_EXTRA_BEFORE * stride;

    for (i = 0; i < pu_h + QPEL_EXTRA; ++i) {
        ref[-1] = ref[0];
        ref[-2] = ref[0];
        ref[pu_w + QPEL_EXTRA    ] = ref[pu_w + QPEL_EXTRA - 1];
        ref[pu_w + QPEL_EXTRA + 1] = ref[pu_w + QPEL_EXTRA - 1];
        ref += stride;
    }

    ref  = _ref;
    ref -= QPEL_EXTRA_BEFORE + QPEL_EXTRA_BEFORE * stride + 2;
    memcpy(ref -     stride, ref, sizeof(*ref) * (pu_w + QPEL_EXTRA + 4));
    memcpy(ref - 2 * stride, ref, sizeof(*ref) * (pu_w + QPEL_EXTRA + 4));

    ref += stride * (pu_h + QPEL_EXTRA - 1);

    memcpy(ref + 1 * stride, ref, sizeof(*ref) * (pu_w + QPEL_EXTRA + 4));
    memcpy(ref + 2 * stride, ref, sizeof(*ref) * (pu_w + QPEL_EXTRA + 4));
}

static struct OVBuffInfo
derive_ref_buf_y(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                uint16_t *edge_buff, int log2_pu_w, int log2_pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    uint16_t *const ref_y  = (uint16_t *) ref_pic->frame->data[0];

    int src_stride = ref_pic->frame->linesize[0] >> 1;

    int ref_pos_x = pos_x + (mv.x >> 4);
    int ref_pos_y = pos_y + (mv.y >> 4);

    int pu_w = 1 << log2_pu_w;
    int pu_h = 1 << log2_pu_h;

    const int pic_w = ref_pic->frame->width[0];
    const int pic_h = ref_pic->frame->height[0];

    uint8_t emulate_edge = test_for_edge_emulation(ref_pos_x, ref_pos_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    /* FIXME Frame thread synchronization here to ensure data is available
     */

    if (emulate_edge){
        const uint16_t *src_y  = &ref_y[ref_pos_x + ref_pos_y * src_stride];
        int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
        int buff_off = REF_PADDING_L * (RCN_CTB_STRIDE) + (REF_PADDING_L);
        int cpy_w = pu_w + QPEL_EXTRA;
        int cpy_h = pu_h + QPEL_EXTRA;

        /*FIXME clip to frame?*/
        int start_pos_x = ref_pos_x - REF_PADDING_L;
        int start_pos_y = ref_pos_y - REF_PADDING_L;

        emulate_block_border(edge_buff + 2 * RCN_CTB_STRIDE + 2, (src_y - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             cpy_w, cpy_h, start_pos_x, start_pos_y,
                             pic_w, pic_h);
        ref_buff.y  = edge_buff + buff_off + 2 * RCN_CTB_STRIDE + 2;
        ref_buff.stride = RCN_CTB_STRIDE;

    } else {

        ref_buff.y  = &ref_y[ref_pos_x + ref_pos_y * src_stride];
        ref_buff.stride = src_stride;
    }
    return ref_buff;
}

static struct OVBuffInfo
derive_dmvr_ref_buf_y(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                      uint16_t *edge_buff, int pu_w, int pu_h)
{
    struct OVBuffInfo ref_buff;
    uint16_t *const ref_y  = (uint16_t *) ref_pic->frame->data[0];

    const int pic_w = ref_pic->frame->width[0];
    const int pic_h = ref_pic->frame->height[0];

    int src_stride = ref_pic->frame->linesize[0] >> 1;

    OVMV mv_clipped = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_pos_x = pos_x + (mv_clipped.x >> 4);
    int ref_pos_y = pos_y + (mv_clipped.y >> 4);

    const uint16_t *src_y = &ref_y[ref_pos_x + ref_pos_y * src_stride];

    int src_off  = (REF_PADDING_L * src_stride) + (REF_PADDING_L);
    int buff_off = (REF_PADDING_L * RCN_CTB_STRIDE) + (REF_PADDING_L);

    int cpy_w = pu_w + QPEL_EXTRA;
    int cpy_h = pu_h + QPEL_EXTRA;

    int start_pos_x = ref_pos_x - REF_PADDING_L;
    int start_pos_y = ref_pos_y - REF_PADDING_L;

    emulate_block_border(edge_buff + 2 * RCN_CTB_STRIDE + 2, (src_y - src_off),
                         RCN_CTB_STRIDE, src_stride,
                         cpy_w, cpy_h, start_pos_x, start_pos_y,
                         pic_w, pic_h);

    ref_buff.y  = edge_buff + buff_off + 2 * RCN_CTB_STRIDE + 2;
    ref_buff.stride = RCN_CTB_STRIDE;

    return ref_buff;
}

static struct OVBuffInfo
derive_dmvr_ref_buf_c(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                      uint16_t *edge_buff0, uint16_t *edge_buff1, int pu_w, int pu_h)
{
    struct OVBuffInfo ref_buff;
    uint16_t *const ref_cb  = (uint16_t *) ref_pic->frame->data[1];
    uint16_t *const ref_cr  = (uint16_t *) ref_pic->frame->data[2];

    const int pic_w = ref_pic->frame->width[1];
    const int pic_h = ref_pic->frame->height[1];

    int src_stride = ref_pic->frame->linesize[1] >> 1;

    OVMV mv_clipped = clip_mv(pos_x << 1, pos_y << 1, pic_w << 1, pic_h << 1, pu_w << 1, pu_h << 1, mv);

    int ref_pos_x = pos_x + (mv_clipped.x >> 5);
    int ref_pos_y = pos_y + (mv_clipped.y >> 5);

    const uint16_t *src_cb = &ref_cb[ref_pos_x + ref_pos_y * src_stride];
    const uint16_t *src_cr = &ref_cr[ref_pos_x + ref_pos_y * src_stride];

    int src_off  = (REF_PADDING_C * src_stride) + (REF_PADDING_C);
    int buff_off = (REF_PADDING_C * RCN_CTB_STRIDE) + (REF_PADDING_C);

    int cpy_w = pu_w + EPEL_EXTRA;
    int cpy_h = pu_h + EPEL_EXTRA;

    int start_pos_x = ref_pos_x - REF_PADDING_C;
    int start_pos_y = ref_pos_y - REF_PADDING_C;

    emulate_block_border(edge_buff0 + 2 * RCN_CTB_STRIDE + 2, (src_cb - src_off),
                         RCN_CTB_STRIDE, src_stride,
                         cpy_w, cpy_h, start_pos_x, start_pos_y,
                         pic_w, pic_h);

    emulate_block_border(edge_buff1 + 2 * RCN_CTB_STRIDE + 2, (src_cr - src_off),
                         RCN_CTB_STRIDE, src_stride,
                         cpy_w, cpy_h, start_pos_x, start_pos_y,
                         pic_w, pic_h);

    ref_buff.cb  = edge_buff0 + buff_off + 2 * RCN_CTB_STRIDE + 2;
    ref_buff.cr  = edge_buff1 + buff_off + 2 * RCN_CTB_STRIDE + 2;
    ref_buff.stride_c = RCN_CTB_STRIDE;
    //padd_dmvr(edge_buff + buff_off + 2 * RCN_CTB_STRIDE + 2, ref_buff.stride, pu_w, pu_h);

    return ref_buff;
}

static void
rcn_motion_compensation_b(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                          uint8_t x0, uint8_t y0,
                          uint8_t log2_pu_w, uint8_t log2_pu_h,
                          OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = ref_idx0;
    uint8_t ref_idx_1 = ref_idx1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx_0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx_1];

    int16_t bcw_weights[5] = { -2, 3, 4, 5, 10 };
    int16_t wt0, wt1;

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];
    uint16_t edge_buff0_1[RCN_CTB_SIZE];
    uint16_t edge_buff1_1[RCN_CTB_SIZE];
    int16_t tmp_buff[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);

    struct OVBuffInfo ref0_b = derive_ref_buf_y(ref0, mv0, pos_x, pos_y, edge_buff0,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    struct OVBuffInfo ref1_b = derive_ref_buf_y(ref1, mv1, pos_x, pos_y, edge_buff1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    uint8_t prec_x0 = (mv0.x) & 0xF;
    uint8_t prec_y0 = (mv0.y) & 0xF;
    uint8_t prec_x1 = (mv1.x) & 0xF;
    uint8_t prec_y1 = (mv1.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    dst.y  += x0 + y0 * dst.stride;
    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;


    mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride, pu_h, prec_x0, prec_y0, pu_w);
    
    if( mv0.bcw_idx_plus1 == 0 || mv0.bcw_idx_plus1 == 3){
        mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride, tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
    }
    else{
        wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        wt0 = 8 - wt1;
        // put_vvc_bi_w_qpel_h4_10_sse(uint8_t* _dst,
        //                              ptrdiff_t _dststride,
        //                              uint8_t* _src,
        //                              ptrdiff_t _srcstride,
        //                              int16_t* src2,
        //                              ptrdiff_t src2stride,
        //                              int height,
        //                              int denom,
        //                              int _wx0,
        //                              int _wx1,
        //                              int _ox0,
        //                              int _ox1,
        //                              intptr_t mx,
        //                              intptr_t my,
        //                              int width)
        // int denom = 3;
        // int ox0 = 0; int ox1 = 0;
        // mc_l->bidir_w[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride, tmp_buff,  MAX_PB_SIZE,
        //                                 pu_h, denom, wt0, wt1, ox0, ox1, prec_x1, prec_y1, pu_w);

        mc_l->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride, tmp_buff, 
                                        pu_h, prec_x1, prec_y1, pu_w, wt0, wt1);
    }

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);

    const struct OVBuffInfo ref0_c = derive_ref_buf_c(ref0, mv0,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff0, edge_buff0_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_c = derive_ref_buf_c(ref1, mv1,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff1, edge_buff1_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);
    prec_x0 = (mv0.x) & 0x1F;
    prec_y0 = (mv0.y) & 0x1F;

    prec_x1 = (mv1.x) & 0x1F;
    prec_y1 = (mv1.y) & 0x1F;

    prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int16_t* ref_data0 = tmp_buff;
    int16_t* ref_data1 = tmp_buff + MAX_PB_SIZE / 2;

    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);

    if( mv0.bcw_idx_plus1 == 0 || mv0.bcw_idx_plus1 == 3 ){
        mc_c->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE, ref1_c.cb, ref1_c.stride_c, ref_data0, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
        mc_c->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE, ref1_c.cr, ref1_c.stride_c, ref_data1, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
    }
    else{
        mc_c->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE, ref1_c.cb, ref1_c.stride_c, ref_data0, 
                                                    pu_h >> 1, prec_x1, prec_y1, pu_w >> 1, wt0, wt1);
        mc_c->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE, ref1_c.cr, ref1_c.stride_c, ref_data1, 
                                                    pu_h >> 1, prec_x1, prec_y1, pu_w >> 1, wt0, wt1);
    }

}

static void
rcn_motion_compensation_b_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                            uint8_t x0, uint8_t y0,
                            uint8_t log2_pu_w, uint8_t log2_pu_h,
                            OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = ref_idx0;
    uint8_t ref_idx_1 = ref_idx1;
    
    int16_t bcw_weights[5] = { -2, 3, 4, 5, 10 };
    int16_t wt0, wt1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx_0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx_1];

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];
    int16_t tmp_buff[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    struct OVBuffInfo ref0_b = derive_ref_buf_y(ref0, mv0, pos_x, pos_y, edge_buff0,
                                                log2_pu_w, log2_pu_h, log2_ctb_s);

    struct OVBuffInfo ref1_b = derive_ref_buf_y(ref1, mv1, pos_x, pos_y, edge_buff1,
                                                log2_pu_w, log2_pu_h, log2_ctb_s);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    uint8_t prec_x0 = (mv0.x) & 0xF;
    uint8_t prec_y0 = (mv0.y) & 0xF;
    uint8_t prec_x1 = (mv1.x) & 0xF;
    uint8_t prec_y1 = (mv1.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    dst.y  += x0 + y0 * dst.stride;

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride,
                                                pu_h, prec_x0, prec_y0, pu_w);

    if( mv0.bcw_idx_plus1 == 0 || mv0.bcw_idx_plus1 == 3){
        mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride, 
                                                    tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
    }
    else{
        wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        wt0 = 8 - wt1;
        mc_l->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride, tmp_buff, 
                                        pu_h, prec_x1, prec_y1, pu_w, wt0, wt1);
    }

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE,
                                  ctudec->lmcs_info.lmcs_lut_fwd_luma,
                                  pu_w, pu_h);

}

#define DMVR_NUM_ITERATION 2
#define DMVR_SAD_STRIDE  ((2 * DMVR_NUM_ITERATION) + 1)
#define DMVR_NB_IDX (DMVR_SAD_STRIDE * DMVR_SAD_STRIDE)
#define DMVR_SB_PXL_LVL 4


struct DMVRDelta
{
    int32_t h;
    int32_t v;
};

static uint64_t
rcn_dmvr_sad(const int16_t *ref0, const int16_t *ref1,
             int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
  uint64_t sum = 0;
  int i;

  if (pb_w == 16) {
      for (i = 0; i < (pb_h >> 1); ++i) {
          sum += abs(ref0[0]  - ref1[0]);
          sum += abs(ref0[1]  - ref1[1]);
          sum += abs(ref0[2]  - ref1[2]);
          sum += abs(ref0[3]  - ref1[3]);
          sum += abs(ref0[4]  - ref1[4]);
          sum += abs(ref0[5]  - ref1[5]);
          sum += abs(ref0[6]  - ref1[6]);
          sum += abs(ref0[7]  - ref1[7]);
          sum += abs(ref0[8]  - ref1[8]);
          sum += abs(ref0[9]  - ref1[9]);
          sum += abs(ref0[10] - ref1[10]);
          sum += abs(ref0[11] - ref1[11]);
          sum += abs(ref0[12] - ref1[12]);
          sum += abs(ref0[13] - ref1[13]);
          sum += abs(ref0[14] - ref1[14]);
          sum += abs(ref0[15] - ref1[15]);

          ref0 += dmvr_stride << 1;
          ref1 += dmvr_stride << 1;
      }
  } else {
      for (i = 0; i < (pb_h >> 1); ++i) {
          sum += abs(ref0[0]  - ref1[0]);
          sum += abs(ref0[1]  - ref1[1]);
          sum += abs(ref0[2]  - ref1[2]);
          sum += abs(ref0[3]  - ref1[3]);
          sum += abs(ref0[4]  - ref1[4]);
          sum += abs(ref0[5]  - ref1[5]);
          sum += abs(ref0[6]  - ref1[6]);
          sum += abs(ref0[7]  - ref1[7]);

          ref0 += dmvr_stride << 1;
          ref1 += dmvr_stride << 1;
      }
  }
  return sum;
}

/*FIXME return min_dmvr_idx; */
static uint8_t
dmvr_compute_sads(const int16_t *ref0, const int16_t *ref1,
                  uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const int16_t *const ref0_start = ref0;
    const int16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx;

    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        /* TODO skip 12 already filled by activation first check */
        uint64_t cost = rcn_dmvr_sad(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
        if (idx == 12) {
            cost -= cost >> 2;
        }

        sad_array[idx] = cost;

        if (cost < min_cost || (idx == 12 && cost <= min_cost)) {
            min_cost = cost;
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

/* FIXME understand this */
static inline int32_t
div_for_maxq7(int64_t num, int64_t den)
{
  int32_t sign, q;
  sign = 0;

  if (num < 0) {
    sign = 1;
    num = -num;
  }

  q = 0;
  den = (den << 3);

  if (num >= den) {
    num -= den;
    q++;
  }

  q = (q << 1);

  den = (den >> 1);

  if (num >= den) {
    num -= den;
    q++;
  }

  q = (q << 1);

  if (num >= (den >> 1)) {
    q++;
  }

  if (sign) {
    return -q;
  }

  return q;
}


static struct DMVRDelta
refine_mv(uint64_t *sad_buff)
{
    struct DMVRDelta delta = {0};
    uint64_t sad_lc[5];
    sad_lc[2] = sad_buff[-DMVR_SAD_STRIDE];
    sad_lc[1] = sad_buff[              -1];
    sad_lc[0] = sad_buff[               0];
    sad_lc[3] = sad_buff[               1];
    sad_lc[4] = sad_buff[ DMVR_SAD_STRIDE];

    /* delta mv horizontal */
    int64_t den_h = (int64_t)(sad_lc[1] + sad_lc[3] - (sad_lc[0] << 1));
    int64_t den_v = (int64_t)(sad_lc[2] + sad_lc[4] - (sad_lc[0] << 1));

    if (0 != den_h) {
        if ((sad_lc[1] != sad_lc[0]) && (sad_lc[3] != sad_lc[0])) {
            int64_t num = (int64_t)(sad_lc[1] - sad_lc[3]) << DMVR_SB_PXL_LVL;
            delta.h = div_for_maxq7(num, den_h);
        } else {
            delta.h = (sad_lc[1] == sad_lc[0]) ? -8 : 8;
        }
    }

    /* delta mv vertical */
    if (0 != den_v) {
        if ((sad_lc[2] != sad_lc[0]) && (sad_lc[4] != sad_lc[0])) {
            int64_t num = (int64_t)((sad_lc[2] - sad_lc[4]) << DMVR_SB_PXL_LVL);
            delta.v = div_for_maxq7(num, den_v);
        } else {
            delta.v = (sad_lc[2] == sad_lc[0]) ? -8 : 8;
        }
    }
    return delta;
}

#define MV_FRACTIONAL_BITS_INTERNAL 4
#define DMVR_REF_PADD 2
#define MV_MIN -(1 << 17)
#define MV_MAX ((1 << 17) - 1)

uint8_t
rcn_dmvr_mv_refine(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_pu_w, uint8_t log2_pu_h,
                   OVMV *mv0, OVMV *mv1, uint8_t ref_idx0, uint8_t ref_idx1, uint8_t 
                   apply_bdof)
{
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx1];

    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];

    /*FIXME permit smaller stride to reduce tables */
    int16_t ref_dmvr0[(16 + 2 * DMVR_REF_PADD) * (128 + 2 * DMVR_REF_PADD)] = {0};
    int16_t ref_dmvr1[(16 + 2 * DMVR_REF_PADD) * (128 + 2 * DMVR_REF_PADD)] = {0};

    int16_t tmp_buff[RCN_CTB_SIZE];

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;
    OVMV tmp0 = *mv0;
    OVMV tmp1 = *mv1;

    struct OVBuffInfo ref0_b = derive_dmvr_ref_buf_y(ref0, *mv0, pos_x, pos_y, edge_buff0,
                                                     pu_w, pu_h);

    struct OVBuffInfo ref1_b = derive_dmvr_ref_buf_y(ref1, *mv1, pos_x, pos_y, edge_buff1,
                                                     pu_w, pu_h);

    uint8_t prec_x0 = (mv0->x) & 0xF;
    uint8_t prec_y0 = (mv0->y) & 0xF;

    uint8_t prec_x1 = (mv1->x) & 0xF;
    uint8_t prec_y1 = (mv1->y) & 0xF;

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int16_t dmvr_stride = 128 + 4;
    uint64_t dmvr_sad;
    uint64_t min_cost; 


    /* Interpolate ref 0/1 with 2 additional samples before and after */
    mc_l->bilinear[prec_0_mc_type]((uint16_t *)ref_dmvr0, dmvr_stride,
                                   ref0_b.y -2 -2 * ref0_b.stride, ref0_b.stride, pu_h + 4,
                                   prec_x0, prec_y0, pu_w + 4);

    mc_l->bilinear[prec_1_mc_type]((uint16_t *)ref_dmvr1, dmvr_stride,
                                   ref1_b.y -2 -2 * ref1_b.stride, ref1_b.stride, pu_h + 4,
                                   prec_x1, prec_y1, pu_w + 4);

    /* Compute SAD on center part */
    dmvr_sad = rcn_dmvr_sad(ref_dmvr0 + 2 + 2 * dmvr_stride,
                            ref_dmvr1 + 2 + 2 * dmvr_stride,
                            dmvr_stride, pu_w, pu_h);

    min_cost = (dmvr_sad - (dmvr_sad >> 2));

    /* skip MV refinement if cost is small or zero */
    if (min_cost >= (pu_w * pu_h)) {
        uint64_t sad[25];
        uint8_t dmvr_idx = dmvr_compute_sads(ref_dmvr0 + 2 + 2 * dmvr_stride,
                                             ref_dmvr1 + 2 + 2 * dmvr_stride,
                                             sad, pu_w, pu_h);

        int32_t delta_h = dmvr_mv_x[dmvr_idx] << 4;
        int32_t delta_v = dmvr_mv_y[dmvr_idx] << 4;

        min_cost = sad[dmvr_idx];

        if ((abs(delta_h) != (2 << MV_FRACTIONAL_BITS_INTERNAL)) &&
            (abs(delta_v) != (2 << MV_FRACTIONAL_BITS_INTERNAL))) {
            struct DMVRDelta add = refine_mv(&sad[dmvr_idx]);

            delta_h += add.h;
            delta_v += add.v;
        }

        /* FIXME if dmvr_idx != 12 force DBF ?*/

        mv0->x += delta_h;
        mv0->y += delta_v;

        mv1->x -= delta_h;
        mv1->y -= delta_v;

        /* FIXME clipping should be done on copy 
         * since we use actual MVs for TMVP ? */
        mv0->x = ov_clip(mv0->x, MV_MIN, MV_MAX);
        mv0->y = ov_clip(mv0->y, MV_MIN, MV_MAX);

        mv1->x = ov_clip(mv1->x, MV_MIN, MV_MAX);
        mv1->y = ov_clip(mv1->y, MV_MIN, MV_MAX);
    }

    /* TODO Call classic reconstruction */

    padd_dmvr((int16_t *)ref0_b.y, ref0_b.stride, pu_w, pu_h);
    padd_dmvr((int16_t *)ref1_b.y, ref1_b.stride, pu_w, pu_h);

    prec_x0 = (mv0->x) & 0xF;
    prec_y0 = (mv0->y) & 0xF;

    prec_x1 = (mv1->x) & 0xF;
    prec_y1 = (mv1->y) & 0xF;

    if(inter_ctx->prec_amvr == 3){
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    dst.y  += x0 + y0 * dst.stride;

    prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int delta_h2 = (mv0->x >> 4) - (tmp0.x >> 4);
    int delta_v2 = (mv0->y >> 4) - (tmp0.y >> 4);

    int delta_h3 = (mv1->x >> 4) - (tmp1.x >> 4);
    int delta_v3 = (mv1->y >> 4) - (tmp1.y >> 4);

    ref0_b.y += delta_h2 + (int16_t)ref0_b.stride * delta_v2;
    ref1_b.y += delta_h3 + (int16_t)ref1_b.stride * delta_v3;

    uint8_t disable_bdof = apply_bdof ? min_cost < 2 * (pu_w * pu_h) : 1;

    if (disable_bdof) {
        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride,
                                                    pu_h, prec_x0, prec_y0, pu_w);
        mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride,
                                                    tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
    } else {
        int16_t grad_x0[(16 + 2) * (16 + 2)];
        int16_t grad_y0[(16 + 2) * (16 + 2)];
        int16_t grad_x1[(16 + 2) * (16 + 2)];
        int16_t grad_y1[(16 + 2) * (16 + 2)];

        int16_t tmp_buff1[RCN_CTB_SIZE];

        int16_t ref_stride = 128;
        int16_t grad_stride = pu_w + 2;

        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff + 128 + 1,
                                                ref0_b.y, ref0_b.stride, pu_h,
                                                prec_x0, prec_y0, pu_w);

        mc_l->bidir0[prec_1_mc_type][log2_pu_w - 1](tmp_buff1 + 128 + 1,
                                                ref1_b.y, ref1_b.stride, pu_h,
                                                prec_x1, prec_y1, pu_w);

        /* Padding for grad derivation */
        extend_bdof_buff(ref0_b.y, (uint16_t*)tmp_buff, ref0_b.stride, pu_w, pu_h, prec_x0 >> 3, prec_y0 >> 3);
        extend_bdof_buff(ref1_b.y, (uint16_t*)tmp_buff1, ref1_b.stride, pu_w, pu_h, prec_x1 >> 3, prec_y1 >> 3);

        compute_prof_grad((uint16_t *)tmp_buff, ref_stride, pu_w, pu_h, grad_stride,
                          grad_x0 + grad_stride + 1, grad_y0 + grad_stride + 1);

        compute_prof_grad((uint16_t *)tmp_buff1, ref_stride, pu_w, pu_h, grad_stride,
                          grad_x1 + grad_stride + 1, grad_y1 + grad_stride + 1);

        /* Grad padding */
        extend_bdof_grad((uint16_t *)grad_x0, grad_stride, pu_w, pu_h);
        extend_bdof_grad((uint16_t *)grad_y0, grad_stride, pu_w, pu_h);
        extend_bdof_grad((uint16_t *)grad_x1, grad_stride, pu_w, pu_h);
        extend_bdof_grad((uint16_t *)grad_y1, grad_stride, pu_w, pu_h);

        /* Reference padding overwrite for weights derivation */
        extend_bdof_grad((uint16_t *)tmp_buff, ref_stride, pu_w, pu_h);
        extend_bdof_grad((uint16_t *)tmp_buff1, ref_stride, pu_w, pu_h);

        /* Split into 4x4 subblocks for BDOF computation */
        rcn_bdof((int16_t *)dst.y, dst.stride, tmp_buff + 128 + 1, tmp_buff1 + 128 + 1,
                 ref_stride, grad_x0, grad_y0, grad_x1, grad_y1,
                 grad_stride, pu_w, pu_h);

    }

    if (ctudec->lmcs_info.lmcs_enabled_flag){
        rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE,
                                        ctudec->lmcs_info.lmcs_lut_fwd_luma,
                                        pu_w, pu_h);
    }

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;
    //uint16_t *edge_buff0_1 = edge_buff0 + 24 * RCN_CTB_STRIDE;
    //uint16_t *edge_buff1_1 = edge_buff1 + 24 * RCN_CTB_STRIDE;
    uint16_t edge_buff0_1[RCN_CTB_SIZE];
    uint16_t edge_buff1_1[RCN_CTB_SIZE];

    struct OVBuffInfo ref0_c = derive_dmvr_ref_buf_c(ref0, tmp0,
                                                     pos_x >> 1, pos_y >> 1,
                                                     edge_buff0, edge_buff0_1,
                                                     pu_w >> 1, pu_h >> 1);

    struct OVBuffInfo ref1_c = derive_dmvr_ref_buf_c(ref1, tmp1,
                                                     pos_x >> 1, pos_y >> 1,
                                                     edge_buff1, edge_buff1_1,
                                                     pu_w >> 1, pu_h >> 1);

    padd_dmvr_c((int16_t *)ref0_c.cb, ref0_c.stride_c, pu_w >> 1, pu_h >> 1);
    padd_dmvr_c((int16_t *)ref0_c.cr, ref0_c.stride_c, pu_w >> 1, pu_h >> 1);

    padd_dmvr_c((int16_t *)ref1_c.cb, ref1_c.stride_c, pu_w >> 1, pu_h >> 1);
    padd_dmvr_c((int16_t *)ref1_c.cr, ref1_c.stride_c, pu_w >> 1, pu_h >> 1);

    prec_x0 = (mv0->x) & 0x1F;
    prec_y0 = (mv0->y) & 0x1F;

    prec_x1 = (mv1->x) & 0x1F;
    prec_y1 = (mv1->y) & 0x1F;

    prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int16_t* ref_data0 = tmp_buff;
    int16_t* ref_data1 = tmp_buff + MAX_PB_SIZE / 2;

    delta_h2 = (mv0->x >> 5) - (tmp0.x >> 5);
    delta_v2 = (mv0->y >> 5) - (tmp0.y >> 5);

    delta_h3 = (mv1->x >> 5) - (tmp1.x >> 5);
    delta_v3 = (mv1->y >> 5) - (tmp1.y >> 5);

    ref0_c.cb += (delta_h2) + ((int16_t)ref0_c.stride_c * (delta_v2));
    ref1_c.cb += (delta_h3) + ((int16_t)ref1_c.stride_c * (delta_v3));

    ref0_c.cr += (delta_h2) + ((int16_t)ref0_c.stride_c * (delta_v2));
    ref1_c.cr += (delta_h3) + ((int16_t)ref1_c.stride_c * (delta_v3));

    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 2](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 2](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);

    mc_c->bidir1[prec_1_mc_type][log2_pu_w - 2](dst.cb, RCN_CTB_STRIDE, ref1_c.cb, ref1_c.stride_c, ref_data0, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
    mc_c->bidir1[prec_1_mc_type][log2_pu_w - 2](dst.cr, RCN_CTB_STRIDE, ref1_c.cr, ref1_c.stride_c, ref_data1, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);

    return disable_bdof;
}

struct PROFInfo
{
    int32_t dmv_scale_h_0[16];
    int32_t dmv_scale_v_0[16];
    int32_t dmv_scale_h_1[16];
    int32_t dmv_scale_v_1[16];
};

#define ov_clip_pixel2(a) ov_clip_uintp2(a, BITDEPTH)
static void
tmp_mrg(uint16_t* _dst, ptrdiff_t _dststride,
        const uint16_t* _src0, ptrdiff_t _srcstride,
        const int16_t* _src1, int height, intptr_t mx,
        intptr_t my, int width)
{
    int x, y;
    const int16_t* src0 = (int16_t *)_src0;
    const int16_t* src1 = (int16_t *)_src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 1;
    int offset = 16400;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_clip_pixel2((src0[x]-(1<<13) + src1[x]-(1<<13) + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
tmp_mrg_w(uint16_t* _dst, ptrdiff_t _dststride,
        const uint16_t* _src0, ptrdiff_t _srcstride,
        const int16_t* _src1, int height, intptr_t mx,
        intptr_t my, int width, int wt0, int wt1)
{
    int x, y;
    const int16_t* src0 = (int16_t *)_src0;
    const int16_t* src1 = (int16_t *)_src1;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int log_weights = floor_log2(wt0 + wt1);
    int shift = 14 - BITDEPTH + log_weights;
    int offset = 16400 << (log_weights - 1) ;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_clip_pixel2(((src0[x] - (1<<13)) * wt0 + (src1[x]-(1<<13)) * wt1 + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
rcn_apply_bdof_subblock(const int16_t* src0, int src0_stride,
                        const int16_t* src1, int src1_stride,
                        int16_t *dst, int dst_stride,
                        const int16_t *gradX0, const int16_t *gradX1,
                        const int16_t *gradY0, const int16_t *gradY1, int grad_stride,
                        int wgt_x, int wgt_y)
{
    int i/*, j*/;

    for (i = 0; i < SB_H; i++) {
        #if 0
        for (j = 0; j < SB_W; j += 4) {
        #endif
            int32_t b0, b1, b2, b3;
            int16_t val0, val1, val2, val3;

            b0 = wgt_x * (gradX0[0] - gradX1[0]) + wgt_y * (gradY0[0] - gradY1[0]);
            b1 = wgt_x * (gradX0[1] - gradX1[1]) + wgt_y * (gradY0[1] - gradY1[1]);
            b2 = wgt_x * (gradX0[2] - gradX1[2]) + wgt_y * (gradY0[2] - gradY1[2]);
            b3 = wgt_x * (gradX0[3] - gradX1[3]) + wgt_y * (gradY0[3] - gradY1[3]);

            val0 = (int16_t)((src0[0] + src1[0] + b0 + BDOF_OFFSET) >> BDOF_SHIFT);
            val1 = (int16_t)((src0[1] + src1[1] + b1 + BDOF_OFFSET) >> BDOF_SHIFT);
            val2 = (int16_t)((src0[2] + src1[2] + b2 + BDOF_OFFSET) >> BDOF_SHIFT);
            val3 = (int16_t)((src0[3] + src1[3] + b3 + BDOF_OFFSET) >> BDOF_SHIFT);

            dst[0] = ov_clip(val0, 0, 1023);
            dst[1] = ov_clip(val1, 0, 1023);
            dst[2] = ov_clip(val2, 0, 1023);
            dst[3] = ov_clip(val3, 0, 1023);
        #if 0
        }
        #endif

        dst += dst_stride;

        src0 += src0_stride;
        src1 += src1_stride;

        gradX0 += grad_stride;
        gradX1 += grad_stride;

        gradY0 += grad_stride;
        gradY1 += grad_stride;
    }
}

static void
derive_bdof_weights(const int16_t* ref0, const int16_t* ref1,
                    const int16_t* grad_x0, const int16_t* grad_x1,
                    const int16_t* grad_y0, const int16_t* grad_y1,
                    const int src0_stride, const int src1_stride,
                    const int grad_stride,
                    int *weight_x, int *weight_y)
{
    int sum_avg_x = 0;
    int sum_avg_y = 0;

    int sum_delta_x = 0;
    int sum_delta_y = 0;
    int wgt_x = 0;
    int wgt_y = 0;

    /* FIXME understand this part */
    /* sum / substract avg_grad_x based on avg_grad_y signs*/
    int sum_avg_x_y_signs = 0;

    int i, j;

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            int32_t avg_grad_x = (grad_x0[j] + grad_x1[j]) >> 1;
            int32_t avg_grad_y = (grad_y0[j] + grad_y1[j]) >> 1;

            int32_t delta_ref = ((ref1[j] -(1<<13)) >> 4) - ((ref0[j] -(1<<13)) >> 4);

            sum_avg_x += abs(avg_grad_x);
            sum_avg_y += abs(avg_grad_y);

            sum_avg_x_y_signs += (avg_grad_y < 0 ? -avg_grad_x : (avg_grad_y == 0 ? 0 : avg_grad_x));

            sum_delta_x += (avg_grad_x < 0 ? -delta_ref : (avg_grad_x == 0 ? 0 : delta_ref));
            sum_delta_y += (avg_grad_y < 0 ? -delta_ref : (avg_grad_y == 0 ? 0 : delta_ref));
        }

        ref1 += src1_stride;
        ref0 += src0_stride;

        grad_x0 += grad_stride;
        grad_x1 += grad_stride;

        grad_y0 += grad_stride;
        grad_y1 += grad_stride;
    }

    if (sum_avg_x) {
        int log2_renorm_x = floor_log2(sum_avg_x);

        wgt_x = (sum_delta_x << 2) >> log2_renorm_x;
        wgt_x = ov_clip(wgt_x, -BDOF_WGT_LIMIT, BDOF_WGT_LIMIT);
        *weight_x = wgt_x;
    }

    if (sum_avg_y) {
        int log2_renorm_y = floor_log2(sum_avg_y);
        int x_offset = 0;

        if (wgt_x) {
            /* FIXME understand this part */
            int high = sum_avg_x_y_signs >> 12;
            int low  = sum_avg_x_y_signs & ((1 << 12) - 1);
            x_offset = (((wgt_x * high) << 12) + (wgt_x * low)) >> 1;
        }

        wgt_y = ((sum_delta_y << 2) - x_offset) >> log2_renorm_y;
        wgt_y = ov_clip(wgt_y, -BDOF_WGT_LIMIT, BDOF_WGT_LIMIT);
        *weight_y = wgt_y;
    }
}

static void
rcn_bdof(int16_t *dst, int dst_stride,
         const int16_t *ref_bdof0, const int16_t *ref_bdof1, int ref_stride,
         const int16_t *grad_x0, const int16_t *grad_y0,
         const int16_t *grad_x1, const int16_t *grad_y1,
         int grad_stride, uint8_t pb_w, uint8_t pb_h)
{
    int nb_sb_w = (pb_w >> 2);
    int nb_sb_h = (pb_h >> 2);

    const int16_t *grad_x0_ln = grad_x0;
    const int16_t *grad_y0_ln = grad_y0;
    const int16_t *grad_x1_ln = grad_x1;
    const int16_t *grad_y1_ln = grad_y1;

    const int16_t *ref0_ln = ref_bdof0 - 128 - 1;
    const int16_t *ref1_ln = ref_bdof1 - 128 - 1;
    int16_t *dst_ln = dst;

    int i, j;

    for (i = 0; i < nb_sb_h; i++) {
        const int16_t *ref0_tmp = ref0_ln;
        const int16_t *ref1_tmp = ref1_ln;

        grad_x0 = grad_x0_ln;
        grad_y0 = grad_y0_ln;
        grad_x1 = grad_x1_ln;
        grad_y1 = grad_y1_ln;

        dst = dst_ln;

        for (j = 0; j < nb_sb_w; j++) {
            int wgt_x = 0;
            int wgt_y = 0;

            derive_bdof_weights(ref0_tmp, ref1_tmp,
                                grad_x0, grad_x1, grad_y0, grad_y1,
                                ref_stride, ref_stride,
                                grad_stride,
                                &wgt_x, &wgt_y);

            rcn_apply_bdof_subblock(ref0_tmp + ref_stride + 1, ref_stride,
                                    ref1_tmp + ref_stride + 1, ref_stride,
                                    dst, dst_stride,
                                    grad_x0 + grad_stride + 1, grad_x1 + grad_stride + 1,
                                    grad_y0 + grad_stride + 1, grad_y1 + grad_stride + 1,
                                    grad_stride,
                                    wgt_x, wgt_y);

            grad_x0 += 1 << 2;
            grad_x1 += 1 << 2;
            grad_y0 += 1 << 2;
            grad_y1 += 1 << 2;
            ref0_tmp += 1 << 2;
            ref1_tmp += 1 << 2;
            dst    += 1 << 2;
        }

        grad_x0_ln += grad_stride << 2;
        grad_y0_ln += grad_stride << 2;
        grad_x1_ln += grad_stride << 2;
        grad_y1_ln += grad_stride << 2;

        ref0_ln += ref_stride << 2;
        ref1_ln += ref_stride << 2;

        dst_ln += dst_stride << 2;
    }
}

void
rcn_bdof_mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
               uint8_t x0, uint8_t y0,
               uint8_t log2_pu_w, uint8_t log2_pu_h,
               OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = ref_idx0;
    uint8_t ref_idx_1 = ref_idx1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx_0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx_1];

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    struct OVBuffInfo ref0_b = derive_ref_buf_y(ref0, mv0, pos_x, pos_y, edge_buff0,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    struct OVBuffInfo ref1_b = derive_ref_buf_y(ref1, mv1, pos_x, pos_y, edge_buff1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    uint8_t prec_x0 = (mv0.x) & 0xF;
    uint8_t prec_y0 = (mv0.y) & 0xF;
    uint8_t prec_x1 = (mv1.x) & 0xF;
    uint8_t prec_y1 = (mv1.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }
    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    /*FIXME permit smaller stride to reduce tables */
    int16_t ref_bdof0[(16 + 2 * PROF_BUFF_PADD_H) * (128 + 2 * PROF_BUFF_PADD_W)] = {0};
    int16_t ref_bdof1[(16 + 2 * PROF_BUFF_PADD_H) * (128 + 2 * PROF_BUFF_PADD_W)] = {0};

    /* FIXME tab dim */
    int16_t grad_x0[(16 + 2) * (16 + 2)];
    int16_t grad_y0[(16 + 2) * (16 + 2)];
    int16_t grad_x1[(16 + 2) * (16 + 2)];
    int16_t grad_y1[(16 + 2) * (16 + 2)];

    int16_t ref_stride = 128;
    int16_t grad_stride = pu_w + 2;

    int pb_w = 1 << log2_pu_w;
    int pb_h = 1 << log2_pu_h;

    mc_l->bidir0[prec_0_mc_type][log2_pu_w](ref_bdof0 + 128 + 1,
                                            ref0_b.y, ref0_b.stride, pu_h,
                                            prec_x0, prec_y0, pu_w);

    mc_l->bidir0[prec_1_mc_type][log2_pu_w](ref_bdof1 + 128 + 1,
                                            ref1_b.y, ref1_b.stride, pu_h,
                                            prec_x1, prec_y1, pu_w);

    /* Padding for grad derivation */
    extend_bdof_buff(ref0_b.y, (uint16_t *)ref_bdof0, ref0_b.stride, pb_w, pb_h, prec_x0 >> 3, prec_y0 >> 3);
    extend_bdof_buff(ref1_b.y, (uint16_t *)ref_bdof1, ref1_b.stride, pb_w, pb_h, prec_x1 >> 3, prec_y1 >> 3);

    compute_prof_grad((uint16_t *)ref_bdof0, ref_stride, pb_w, pb_h, grad_stride,
                      grad_x0 + grad_stride + 1, grad_y0 + grad_stride + 1);

    compute_prof_grad((uint16_t *)ref_bdof1, ref_stride, pb_w, pb_h, grad_stride,
                      grad_x1 + grad_stride + 1, grad_y1 + grad_stride + 1);

    /* Grad padding */
    extend_bdof_grad((uint16_t *)grad_x0, grad_stride, pb_w, pb_h);
    extend_bdof_grad((uint16_t *)grad_y0, grad_stride, pb_w, pb_h);
    extend_bdof_grad((uint16_t *)grad_x1, grad_stride, pb_w, pb_h);
    extend_bdof_grad((uint16_t *)grad_y1, grad_stride, pb_w, pb_h);

    /* Reference padding overwrite for weights derivation */
    extend_bdof_grad((uint16_t *)ref_bdof0, ref_stride, pb_w, pb_h);
    extend_bdof_grad((uint16_t *)ref_bdof1, ref_stride, pb_w, pb_h);
    dst.y += x0;
    dst.y += y0 * RCN_CTB_STRIDE;

    /* Split into 4x4 subblocks for BDOF computation */
    rcn_bdof((int16_t *)dst.y, dst.stride, ref_bdof0 + 128 + 1, ref_bdof1 + 128 + 1,
             ref_stride, grad_x0, grad_y0, grad_x1, grad_y1,
             grad_stride, pb_w, pb_h);

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE,
                                  ctudec->lmcs_info.lmcs_lut_fwd_luma,
                                  pu_w, pu_h);
}


static void
rcn_prof_motion_compensation_b_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                                 uint8_t x0, uint8_t y0,
                                 uint8_t log2_pu_w, uint8_t log2_pu_h,
                                 OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1,
                                 uint8_t prof_dir, const struct PROFInfo *const prof_info)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = ref_idx0;
    uint8_t ref_idx_1 = ref_idx1;
    
    int16_t bcw_weights[5] = { -2, 3, 4, 5, 10 };
    int16_t wt0, wt1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx_0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx_1];

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];
    int16_t tmp_buff[RCN_CTB_SIZE];
    int16_t tmp_buff1[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    const struct OVBuffInfo ref0_b = derive_ref_buf_y(ref0, mv0, pos_x, pos_y, edge_buff0,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_b = derive_ref_buf_y(ref1, mv1, pos_x, pos_y, edge_buff1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    uint8_t prec_x0 = (mv0.x) & 0xF;
    uint8_t prec_y0 = (mv0.y) & 0xF;
    uint8_t prec_x1 = (mv1.x) & 0xF;
    uint8_t prec_y1 = (mv1.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    dst.y  += x0 + y0 * dst.stride;

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    if (prof_dir & 0x1) {
        int16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (128 + 2 * PROF_BUFF_PADD_W)];
        int16_t tmp_prof_stride = (128);

        int16_t tmp_grad_x[16];
        int16_t tmp_grad_y[16];
        mc_l->bidir0[prec_0_mc_type][log2_pu_w](tmp_prof + 128 + 1,
                                              ref0_b.y, ref0_b.stride, pu_h,
                                              prec_x0, prec_y0, pu_w);

        extend_prof_buff(ref0_b.y, (uint16_t *)tmp_prof, ref0_b.stride, prec_x0 >> 3, prec_y0 >> 3);

        compute_prof_grad((uint16_t *)tmp_prof, tmp_prof_stride, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

        rcn_prof((uint16_t *)tmp_buff, MAX_PB_SIZE, (uint16_t *)tmp_prof + 128 + 1, tmp_prof_stride, tmp_grad_x, tmp_grad_y,
                 4, prof_info->dmv_scale_h_0, prof_info->dmv_scale_v_0, 1);
    } else {
        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride,
                                                    pu_h, prec_x0, prec_y0, pu_w);
    }

    if (prof_dir & 0x2) {
        int16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (128 + 2 * PROF_BUFF_PADD_W)];
        int16_t tmp_prof_stride = (128);

        int16_t tmp_grad_x[16];
        int16_t tmp_grad_y[16];
        mc_l->bidir0[prec_1_mc_type][log2_pu_w](tmp_prof + 128 + 1,
                                                ref1_b.y, ref1_b.stride, pu_h,
                                                prec_x1, prec_y1, pu_w);

        extend_prof_buff(ref1_b.y, (uint16_t *)tmp_prof, ref1_b.stride, prec_x1 >> 3, prec_y1 >> 3);

        compute_prof_grad((uint16_t *)tmp_prof, tmp_prof_stride, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

        rcn_prof((uint16_t *)tmp_buff1, MAX_PB_SIZE, (uint16_t *)tmp_prof + 128 + 1, tmp_prof_stride,
                 tmp_grad_x, tmp_grad_y,
                 4, prof_info->dmv_scale_h_1, prof_info->dmv_scale_v_1, 1);
                 /*FIXME merge */
        if( mv0.bcw_idx_plus1 == 0 || mv0.bcw_idx_plus1 == 3){
            tmp_mrg(dst.y, RCN_CTB_STRIDE, (uint16_t *)tmp_buff1, MAX_PB_SIZE,
                tmp_buff, pu_h, 0, 0, pu_w);
        }
        else{
            wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
            wt0 = 8 - wt1;
            tmp_mrg_w(dst.y, RCN_CTB_STRIDE, (uint16_t *)tmp_buff1, MAX_PB_SIZE,
                tmp_buff, pu_h, 0, 0, pu_w, wt1, wt0);
        }


    } else {
        if( mv0.bcw_idx_plus1 == 0 || mv0.bcw_idx_plus1 == 3){
            mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride,
                                                    tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
        }
        else{
            wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
            wt0 = 8 - wt1;
            mc_l->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.y, RCN_CTB_STRIDE, ref1_b.y, ref1_b.stride,
                                                    tmp_buff, pu_h, prec_x1, prec_y1, pu_w, wt0, wt1);
        }
    }

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE,
                                  ctudec->lmcs_info.lmcs_lut_fwd_luma,
                                  pu_w, pu_h);
}

static void
rcn_motion_compensation_b_c(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                            uint8_t x0, uint8_t y0,
                            uint8_t log2_pu_w, uint8_t log2_pu_h,
                            OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = ref_idx0;
    uint8_t ref_idx_1 = ref_idx1;
    
    int16_t bcw_weights[5] = { -2, 3, 4, 5, 10 };
    int16_t wt0, wt1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx_0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx_1];

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];
    uint16_t edge_buff0_1[RCN_CTB_SIZE];
    uint16_t edge_buff1_1[RCN_CTB_SIZE];
    int16_t tmp_buff[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    const struct OVBuffInfo ref0_c = derive_ref_buf_c(ref0, mv0,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff0, edge_buff0_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_c = derive_ref_buf_c(ref1, mv1,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff1, edge_buff1_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);
    uint8_t prec_x0 = (mv0.x) & 0x1F;
    uint8_t prec_y0 = (mv0.y) & 0x1F;

    uint8_t prec_x1 = (mv1.x) & 0x1F;
    uint8_t prec_y1 = (mv1.y) & 0x1F;

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int16_t* ref_data0 = tmp_buff;
    int16_t* ref_data1 = tmp_buff + MAX_PB_SIZE / 2;

    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);

    if( mv0.bcw_idx_plus1 == 0 || mv0.bcw_idx_plus1 == 3 ){
        mc_c->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE, ref1_c.cb, ref1_c.stride_c, ref_data0, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
        mc_c->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE, ref1_c.cr, ref1_c.stride_c, ref_data1, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
    }
    else{
        wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        wt0 = 8 - wt1;
        mc_c->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE, ref1_c.cb, ref1_c.stride_c, ref_data0, 
                                                    pu_h >> 1, prec_x1, prec_y1, pu_w >> 1, wt0, wt1);
        mc_c->bidir_w2[prec_1_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE, ref1_c.cr, ref1_c.stride_c, ref_data1, 
                                                    pu_h >> 1, prec_x1, prec_y1, pu_w >> 1, wt0, wt1);
    }
}

void
rcn_mcp(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
        OVMV mv, uint8_t type, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.y  += x0 + y0 * dst.stride;
    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    uint16_t tmp_buff [RCN_CTB_SIZE];

    const OVFrame *const frame0 =  type ? ref1->frame : ref0->frame;

    const uint16_t *const ref0_y  = (uint16_t *) frame0->data[0];
    const uint16_t *const ref0_cb = (uint16_t *) frame0->data[1];
    const uint16_t *const ref0_cr = (uint16_t *) frame0->data[2];

    int src_stride   = frame0->linesize[0] >> 1;
    int src_stride_c = frame0->linesize[1] >> 1;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width[0];
    const int pic_h = frame0->height[0];

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }
    uint8_t prec_x_c = (mv.x) & 0x1F;
    uint8_t prec_y_c = (mv.y) & 0x1F;

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);
    int prec_c_mc_type = (prec_x_c > 0) + ((prec_y_c > 0) << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const uint16_t *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];
    const uint16_t *src_cb = &ref0_cb[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];
    const uint16_t *src_cr = &ref0_cr[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];

    /* FIXME
     * Thread synchronization to ensure data is available before usage
     */

    if (emulate_edge){
        int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
        int buff_off = REF_PADDING_L * (RCN_CTB_STRIDE) + (REF_PADDING_L);

        emulate_block_border(tmp_buff, (src_y - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             pu_w + QPEL_EXTRA, pu_h + QPEL_EXTRA,
                             ref_x - REF_PADDING_L, ref_y - REF_PADDING_L,
                             pic_w, pic_h);

        src_y = tmp_buff + buff_off;
        src_stride = RCN_CTB_STRIDE;
    }

    mc_l->unidir[prec_mc_type][log2_pu_w](dst.y, RCN_CTB_STRIDE,
                                          src_y, src_stride, pu_h,
                                          prec_x, prec_y, pu_w);

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);


    emulate_edge = test_for_edge_emulation_c(ref_x >> 1, ref_y >> 1, pic_w >> 1, pic_h >> 1,
                                             pu_w >> 1, pu_h >> 1);;

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cb - src_off),
                             RCN_CTB_STRIDE, src_stride_c,
                             (pu_w >> 1)  + EPEL_EXTRA, (pu_h >> 1) + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cb = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    mc_c->unidir[prec_c_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE,
                                                src_cb, src_stride_c,
                                                pu_h >> 1, prec_x_c, prec_y_c, pu_w >> 1);

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (frame0->linesize[1] >> 1) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cr - src_off),
                             RCN_CTB_STRIDE, frame0->linesize[1] >> 1,
                             (pu_w >> 1) + EPEL_EXTRA, (pu_h >> 1) + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cr = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }
    mc_c->unidir[prec_c_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE,
                                                src_cr, src_stride_c,
                                                pu_h >> 1, prec_x_c, prec_y_c, pu_w >> 1);
}

void
rcn_mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
          OVMV mv, uint8_t type, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.y  += x0 + y0 * dst.stride;

    uint16_t tmp_buff [RCN_CTB_SIZE];

    const OVFrame *const frame0 =  type ? ref1->frame : ref0->frame;

    const uint16_t *const ref0_y  = (uint16_t *) frame0->data[0];

    int src_stride   = frame0->linesize[0] >> 1;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width[0];
    const int pic_h = frame0->height[0];

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const uint16_t *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];

    /* FIXME
     * Thread synchronization to ensure data is available before usage
     */

    if (emulate_edge){
        int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
        int buff_off = REF_PADDING_L * (RCN_CTB_STRIDE) + (REF_PADDING_L);

        emulate_block_border(tmp_buff, (src_y - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             pu_w + QPEL_EXTRA, pu_h + QPEL_EXTRA,
                             ref_x - REF_PADDING_L, ref_y - REF_PADDING_L,
                             pic_w, pic_h);

        src_y = tmp_buff + buff_off;
        src_stride = RCN_CTB_STRIDE;
    }

    mc_l->unidir[prec_mc_type][log2_pu_w](dst.y, RCN_CTB_STRIDE,
                                          src_y, src_stride, pu_h,
                                          prec_x, prec_y, pu_w);

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);
}

void
rcn_prof_mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0,
               int log2_pu_w, int log2_pu_h,
               OVMV mv, uint8_t type,
               uint8_t ref_idx,
               const int32_t *dmv_scale_h, const int32_t *dmv_scale_v)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.y  += x0 + y0 * dst.stride;

    uint16_t tmp_buff [RCN_CTB_SIZE];

    const OVFrame *const frame0 =  type ? ref1->frame : ref0->frame;

    const uint16_t *const ref0_y  = (uint16_t *) frame0->data[0];

    int src_stride   = frame0->linesize[0] >> 1;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width[0];
    const int pic_h = frame0->height[0];

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if(inter_ctx->prec_amvr == 3){
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const uint16_t *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];

    /* FIXME
     * Thread synchronization to ensure data is available before usage
     */
    int16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (128 + 2 * PROF_BUFF_PADD_W)];
    int16_t tmp_prof_stride = (128);

    int16_t tmp_grad_x[16];
    int16_t tmp_grad_y[16];

    if (emulate_edge){
        int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
        int buff_off = REF_PADDING_L * (RCN_CTB_STRIDE) + (REF_PADDING_L);

        emulate_block_border(tmp_buff, (src_y - src_off),
                             RCN_CTB_STRIDE, src_stride,
                             pu_w + QPEL_EXTRA, pu_h + QPEL_EXTRA,
                             ref_x - REF_PADDING_L, ref_y - REF_PADDING_L,
                             pic_w, pic_h);

        src_y = tmp_buff + buff_off;
        src_stride = RCN_CTB_STRIDE;
    }

    /* FIXME specialize and reduce buff sizes */
    mc_l->bidir0[prec_mc_type][log2_pu_w](tmp_prof + 128 + 1,
                                          src_y, src_stride, pu_h,
                                          prec_x, prec_y, pu_w);

    extend_prof_buff(src_y, (uint16_t *)tmp_prof, src_stride, prec_x >> 3, prec_y >> 3);

    compute_prof_grad((uint16_t *)tmp_prof, tmp_prof_stride, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

    rcn_prof(dst.y, dst.stride, (uint16_t *)tmp_prof + 128 + 1, tmp_prof_stride, tmp_grad_x, tmp_grad_y,
             4, dmv_scale_h, dmv_scale_v, 0);

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);
}

void
rcn_mcp_c(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
          OVMV mv, uint8_t type, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    uint16_t tmp_buff [RCN_CTB_SIZE];

    const OVFrame *const frame0 =  type ? ref1->frame : ref0->frame;

    const uint16_t *const ref0_cb = (uint16_t *) frame0->data[1];
    const uint16_t *const ref0_cr = (uint16_t *) frame0->data[2];

    int src_stride_c = frame0->linesize[1] >> 1;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width[0];
    const int pic_h = frame0->height[0];

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x_c = (mv.x) & 0x1F;
    uint8_t prec_y_c = (mv.y) & 0x1F;

    int prec_c_mc_type = (prec_x_c > 0) + ((prec_y_c > 0) << 1);

    const uint16_t *src_cb = &ref0_cb[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];
    const uint16_t *src_cr = &ref0_cr[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];

    /* FIXME
     * Thread synchronization to ensure data is available before usage
     */

    uint8_t emulate_edge = test_for_edge_emulation_c(ref_x >> 1, ref_y >> 1, pic_w >> 1, pic_h >> 1,
                                                     pu_w >> 1, pu_h >> 1);;

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cb - src_off),
                             RCN_CTB_STRIDE, src_stride_c,
                             (pu_w >> 1)  + EPEL_EXTRA, (pu_h >> 1) + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cb = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    mc_c->unidir[prec_c_mc_type][log2_pu_w - 1](dst.cb, RCN_CTB_STRIDE,
                                                src_cb, src_stride_c,
                                                pu_h >> 1, prec_x_c, prec_y_c, pu_w >> 1);

    if (emulate_edge){
        int src_off  = REF_PADDING_C * (frame0->linesize[1] >> 1) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cr - src_off),
                             RCN_CTB_STRIDE, frame0->linesize[1] >> 1,
                             (pu_w >> 1) + EPEL_EXTRA, (pu_h >> 1) + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cr = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }
    mc_c->unidir[prec_c_mc_type][log2_pu_w - 1](dst.cr, RCN_CTB_STRIDE,
                                                src_cr, src_stride_c,
                                                pu_h >> 1, prec_x_c, prec_y_c, pu_w >> 1);
}

void
rcn_mcp_b(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
          const OVPartInfo *const part_ctx,
          const OVMV mv0, const OVMV mv1,
          unsigned int x0, unsigned int y0,
          unsigned int log2_pb_w, unsigned int log2_pb_h,
          uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    if (inter_dir == 3) {

        rcn_motion_compensation_b(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);

    } else if (inter_dir & 0x2) {

        rcn_mcp(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);

    } else if (inter_dir & 0x1) {

        rcn_mcp(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);

    }
}

void
rcn_mcp_b_l(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
            const OVPartInfo *const part_ctx,
            const OVMV mv0, const OVMV mv1,
            unsigned int x0, unsigned int y0,
            unsigned int log2_pb_w, unsigned int log2_pb_h,
            uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    if (inter_dir == 3) {

        rcn_motion_compensation_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);

    } else if (inter_dir & 0x2) {

        rcn_mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);

    } else if (inter_dir & 0x1) {

        rcn_mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);

    }
}


void
rcn_prof_mcp_b_l(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
            const OVPartInfo *const part_ctx,
            const OVMV mv0, const OVMV mv1,
            unsigned int x0, unsigned int y0,
            unsigned int log2_pb_w, unsigned int log2_pb_h,
            uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
            uint8_t prof_dir, const struct PROFInfo *const prof_info)
{
    if (inter_dir == 3) {

        rcn_prof_motion_compensation_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1, prof_dir, prof_info);

    } else if (inter_dir & 0x2) {

        rcn_prof_mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1,
                       prof_info->dmv_scale_h_1,
                       prof_info->dmv_scale_v_1);

    } else if (inter_dir & 0x1) {

        rcn_prof_mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0,
                       prof_info->dmv_scale_h_0,
                       prof_info->dmv_scale_v_0);

    }
}

void
rcn_mcp_b_c(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
            const OVPartInfo *const part_ctx,
            const OVMV mv0, const OVMV mv1,
            unsigned int x0, unsigned int y0,
            unsigned int log2_pb_w, unsigned int log2_pb_h,
            uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    if (inter_dir == 3) {

        rcn_motion_compensation_b_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);

    } else if (inter_dir & 0x2) {

        rcn_mcp_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);

    } else if (inter_dir & 0x1) {

        rcn_mcp_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);

    }
}

void rcn_ciip_weighted_sum(OVCTUDec*const ctudec, struct OVBuffInfo* tmp_intra, struct OVBuffInfo* tmp_inter,
                            unsigned int x0, unsigned int y0,
                            unsigned int log2_pb_w, unsigned int log2_pb_h)
{
    //Compute weight in function of neighboring coding modes
    int x_right  = x0 + (1 << log2_pb_w) - 1;
    int y_bottom = y0 + (1 << log2_pb_h) - 1;
    int mode_abv = ctudec->part_map.cu_mode_x[x_right  >> ctudec->part_ctx->log2_min_cb_s];
    int mode_lft = ctudec->part_map.cu_mode_y[y_bottom >> ctudec->part_ctx->log2_min_cb_s];
    uint8_t cu_intra_abv = mode_abv == OV_INTRA || mode_abv == OV_MIP;
    uint8_t cu_intra_lft = mode_lft == OV_INTRA || mode_lft == OV_MIP;
    int wt = 1;
    wt += (cu_intra_abv + cu_intra_lft) ;

    //Apply weighted sum to the final CIIP predicted block
    struct OVBuffInfo dst = ctudec->rcn_ctx.ctu_buff;
    dst.y  += x0 + y0 * dst.stride;
    tmp_intra->y  += x0 + y0 * tmp_intra->stride;
    tmp_inter->y  += x0 + y0 * tmp_inter->stride;
    put_weighted_ciip_pixels(dst.y, dst.stride, tmp_intra->y, tmp_inter->y, tmp_inter->stride,
                       1 << log2_pb_w, 1 << log2_pb_h, wt);

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    tmp_intra->cb += (x0 >> 1) + (y0 >> 1) * tmp_intra->stride_c;
    tmp_inter->cb += (x0 >> 1) + (y0 >> 1) * tmp_inter->stride_c;

    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    tmp_intra->cr += (x0 >> 1) + (y0 >> 1) * tmp_intra->stride_c;
    tmp_inter->cr += (x0 >> 1) + (y0 >> 1) * tmp_inter->stride_c;

    struct MCFunctions *mc_c = &ctudec->rcn_ctx.rcn_funcs.mc_c;
    if (log2_pb_w <= 2){
        mc_c->unidir[0][0](dst.cb, dst.stride_c, tmp_inter->cb, tmp_inter->stride_c, 1 << (log2_pb_h - 1), 0, 0, 1 << (log2_pb_w - 1));
        mc_c->unidir[0][0](dst.cr, dst.stride_c, tmp_inter->cr, tmp_inter->stride_c, 1 << (log2_pb_h - 1), 0, 0, 1 << (log2_pb_w - 1));
    }
    else{
        put_weighted_ciip_pixels(dst.cb, dst.stride_c, tmp_intra->cb, tmp_inter->cb, tmp_inter->stride_c,
                           1 << (log2_pb_w - 1), 1 << (log2_pb_h - 1), wt);
        put_weighted_ciip_pixels(dst.cr, dst.stride_c, tmp_intra->cr, tmp_inter->cr, tmp_inter->stride_c,
                           1 << (log2_pb_w - 1), 1 << (log2_pb_h - 1), wt);
    }
}

void
rcn_ciip_b(OVCTUDec*const ctudec,
           const OVMV mv0, const OVMV mv1,
           unsigned int x0, unsigned int y0,
           unsigned int log2_pb_w, unsigned int log2_pb_h,
           uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    //Inter merge mode
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    const OVPartInfo *const part_ctx = ctudec->part_ctx;
    struct OVBuffInfo tmp_inter;
    uint16_t tmp_inter_l [RCN_CTB_SIZE], tmp_inter_cb[RCN_CTB_SIZE], tmp_inter_cr[RCN_CTB_SIZE] ;
    tmp_inter.y  = &tmp_inter_l [RCN_CTB_PADDING];
    tmp_inter.cb = &tmp_inter_cb[RCN_CTB_PADDING];
    tmp_inter.cr = &tmp_inter_cr[RCN_CTB_PADDING];
    tmp_inter.stride   = RCN_CTB_STRIDE;
    tmp_inter.stride_c = RCN_CTB_STRIDE;
    rcn_mcp_b(ctudec, tmp_inter, inter_ctx, part_ctx, mv0, mv1, x0, y0, log2_pb_w, log2_pb_h,
        inter_dir, ref_idx0, ref_idx1);

    //Intra Planar mode
    struct OVBuffInfo tmp_intra;
    uint16_t tmp_intra_l [RCN_CTB_SIZE], tmp_intra_cb[RCN_CTB_SIZE], tmp_intra_cr[RCN_CTB_SIZE] ;
    tmp_intra.y  = &tmp_intra_l [RCN_CTB_PADDING];
    tmp_intra.cb = &tmp_intra_cb[RCN_CTB_PADDING];
    tmp_intra.cr = &tmp_intra_cr[RCN_CTB_PADDING];
    tmp_intra.stride   = RCN_CTB_STRIDE;
    tmp_intra.stride_c = RCN_CTB_STRIDE;
    vvc_intra_pred(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0, y0, log2_pb_w, log2_pb_h);
    vvc_intra_pred_chroma(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0 >> 1, y0 >> 1, log2_pb_w - 1, log2_pb_h - 1);

    rcn_ciip_weighted_sum(ctudec, &tmp_intra, &tmp_inter, x0, y0, log2_pb_w, log2_pb_h);
}

void
rcn_ciip(OVCTUDec *const ctudec,
         int x0, int y0, int log2_pb_w, int log2_pb_h,
         OVMV mv, uint8_t ref_idx)
{
    //Inter merge mode
    struct OVBuffInfo tmp_inter;
    uint16_t tmp_inter_l [RCN_CTB_SIZE], tmp_inter_cb[RCN_CTB_SIZE], tmp_inter_cr[RCN_CTB_SIZE] ;
    tmp_inter.y  = &tmp_inter_l [RCN_CTB_PADDING];
    tmp_inter.cb = &tmp_inter_cb[RCN_CTB_PADDING];
    tmp_inter.cr = &tmp_inter_cr[RCN_CTB_PADDING];
    tmp_inter.stride   = RCN_CTB_STRIDE;
    tmp_inter.stride_c = RCN_CTB_STRIDE;
    rcn_mcp(ctudec, tmp_inter, x0, y0, log2_pb_w, log2_pb_h, mv, 0, ref_idx);

    //Intra Planar mode
    struct OVBuffInfo tmp_intra;
    uint16_t tmp_intra_l [RCN_CTB_SIZE], tmp_intra_cb[RCN_CTB_SIZE], tmp_intra_cr[RCN_CTB_SIZE] ;
    tmp_intra.y  = &tmp_intra_l [RCN_CTB_PADDING];
    tmp_intra.cb = &tmp_intra_cb[RCN_CTB_PADDING];
    tmp_intra.cr = &tmp_intra_cr[RCN_CTB_PADDING];
    tmp_intra.stride   = RCN_CTB_STRIDE;
    tmp_intra.stride_c = RCN_CTB_STRIDE;
    vvc_intra_pred(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0, y0, log2_pb_w, log2_pb_h);
    vvc_intra_pred_chroma(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0 >> 1, y0 >> 1, log2_pb_w - 1, log2_pb_h - 1);

    rcn_ciip_weighted_sum(ctudec, &tmp_intra, &tmp_inter, x0, y0, log2_pb_w, log2_pb_h);

}

static int16_t g_globalGeoWeights [GEO_NUM_PRESTORED_MASK][GEO_WEIGHT_MASK_SIZE * GEO_WEIGHT_MASK_SIZE];
static int16_t g_weightOffset     [GEO_NUM_PARTITION_MODE][GEO_NUM_CU_SIZE][GEO_NUM_CU_SIZE][2];

const int16_t g_GeoParams[GEO_NUM_PARTITION_MODE][2] =
{
    {0, 1},
    {0, 3},

    {2, 0},
    {2, 1},
    {2, 2},
    {2, 3},

    {3, 0},
    {3, 1},
    {3, 2},
    {3, 3},

    {4, 0},
    {4, 1},
    {4, 2},
    {4, 3},

    {5, 0},
    {5, 1},
    {5, 2},
    {5, 3},

    {8, 1},
    {8, 3},

    {11, 0},
    {11, 1},
    {11, 2},
    {11, 3},

    {12, 0},
    {12, 1},
    {12, 2},
    {12, 3},

    {13, 0},
    {13, 1},
    {13, 2},
    {13, 3},

    {14, 0},
    {14, 1},
    {14, 2},
    {14, 3},

    {16, 1},
    {16, 3},

    {18, 1},
    {18, 2},
    {18, 3},

    {19, 1},
    {19, 2},
    {19, 3},

    {20, 1},
    {20, 2},
    {20, 3},

    {21, 1},
    {21, 2},
    {21, 3},

    {24, 1},
    {24, 3},

    {27, 1},
    {27, 2},
    {27, 3},

    {28, 1},
    {28, 2},
    {28, 3},

    {29, 1},
    {29, 2},
    {29, 3},

    {30, 1},
    {30, 2},
    {30, 3}
};

const int8_t g_Dis[GEO_NUM_ANGLES] = { 
    8,  8,  8,  8,  4,  4,  2,  1,
    0, -1, -2, -4, -4, -8, -8, -8,

   -8, -8, -8, -8, -4, -4, -2, -1,
    0,  1,  2,  4,  4,  8,  8,  8
};

static const int8_t g_angle2mask[GEO_NUM_ANGLES] = {
    0, -1,  1,  2,  3,  4, -1, -1,
    5, -1, -1,  4,  3,  2,  1, -1,

    0, -1,  1,  2,  3,  4, -1, -1,
    5, -1, -1,  4,  3,  2,  1, -1
};

void rcn_init_gpm_params(){
  int angle_idx;
  for (angle_idx = 0; angle_idx < (GEO_NUM_ANGLES >> 2) + 1; angle_idx++){
    int x, y;

    if (g_angle2mask[angle_idx] == -1){
      continue;
    }

    int dist_x = angle_idx;
    int dist_y = (dist_x + (GEO_NUM_ANGLES >> 2)) % GEO_NUM_ANGLES;
    int16_t rho = ((int)g_Dis[dist_x] << (GEO_MAX_CU_LOG2 + 1)) +
                  ((int)g_Dis[dist_y] << (GEO_MAX_CU_LOG2 + 1));

    static const int16_t offset_msk = (2 * GEO_MAX_CU_SIZE - GEO_WEIGHT_MASK_SIZE) >> 1;
    int idx = 0;
    for (y = 0; y < GEO_WEIGHT_MASK_SIZE; y++) {

      int16_t lookUpY = (((y + offset_msk) << 1) + 1) * g_Dis[dist_y];

      for (x = 0; x < GEO_WEIGHT_MASK_SIZE; x++, idx++) {
        int16_t sx_i = ((x + offset_msk) << 1) + 1;
        int16_t weightIdx = sx_i * g_Dis[dist_x] + lookUpY - rho;
        int weightLinearIdx = 32 + weightIdx;

        g_globalGeoWeights[g_angle2mask[angle_idx]][idx] = ov_clip((weightLinearIdx + 4) >> 3, 0, 8);
      }
    }
  }

  for( int hIdx = 0; hIdx < GEO_NUM_CU_SIZE; hIdx++ )
  {
    int16_t height = 1 << ( hIdx + GEO_MIN_CU_LOG2);
    for( int wIdx = 0; wIdx < GEO_NUM_CU_SIZE; wIdx++ ){
      int16_t width = 1 << (wIdx + GEO_MIN_CU_LOG2);
      for( int splitDir = 0; splitDir < GEO_NUM_PARTITION_MODE; splitDir++ ){
        int16_t angle         = g_GeoParams[splitDir][0];
        int16_t distance      = g_GeoParams[splitDir][1];
        int16_t offsetX       = (GEO_WEIGHT_MASK_SIZE - width) >> 1;
        int16_t offsetY       = (GEO_WEIGHT_MASK_SIZE - height) >> 1;
        if( distance > 0 ){
          if( angle % 16 == 8 || (angle % 16 != 0 && height >= width) ){
            offsetY += angle < 16 ? ((distance * (int32_t)height) >> 3) : -((distance * (int32_t)height) >> 3);
          }
          else{
            offsetX += angle < 16 ? ((distance * (int32_t)width) >> 3) : -((distance * (int32_t)width) >> 3);
          }
        }
        g_weightOffset[splitDir][hIdx][wIdx][0] = offsetX;
        g_weightOffset[splitDir][hIdx][wIdx][1] = offsetY;
      }
    }
  }
}


static void
rcn_gpm_weights_and_steps(int split_dir, int log2_pb_w_l, int log2_pb_h_l, int* step_x, int* step_y, 
                            int16_t** weight, int cr_scale)
{
    static const int8_t g_angle2mirror[GEO_NUM_ANGLES] = {
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 2, 2, 2,

        0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 2, 2, 2, 2
    };

  int16_t angle = g_GeoParams[split_dir][0];
  int16_t wIdx = log2_pb_w_l - GEO_MIN_CU_LOG2;
  int16_t hIdx = log2_pb_h_l - GEO_MIN_CU_LOG2;
  *step_x = 1 << cr_scale;
  *step_y = 0;

  if (g_angle2mirror[angle] == 2){
      int weight_offset = (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[split_dir][hIdx][wIdx][1])
                        * GEO_WEIGHT_MASK_SIZE
                        + g_weightOffset[split_dir][hIdx][wIdx][0];

      *step_y = -(int)((GEO_WEIGHT_MASK_SIZE << cr_scale) + (1 << log2_pb_w_l));
      *weight = &g_globalGeoWeights[g_angle2mask[angle]][weight_offset];
  } else if (g_angle2mirror[angle] == 1) {
      int weight_offset = g_weightOffset[split_dir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[split_dir][hIdx][wIdx][0]);
      *step_x = -1 << cr_scale;
      *step_y = (GEO_WEIGHT_MASK_SIZE << cr_scale) + (1 << log2_pb_w_l);
      *weight = &g_globalGeoWeights[g_angle2mask[angle]][weight_offset];
  } else {
      int weight_offset = g_weightOffset[split_dir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[split_dir][hIdx][wIdx][0];
      *step_y = (GEO_WEIGHT_MASK_SIZE << cr_scale) - (1 << log2_pb_w_l);
      *weight = &g_globalGeoWeights[g_angle2mask[angle]][weight_offset];
  }
}


static void
rcn_gpm_mc(OVCTUDec *const ctudec, struct OVBuffInfo dst, int split_dir,
              uint8_t x0, uint8_t y0, uint8_t log2_pu_w, uint8_t log2_pu_h, 
              int type0, OVMV mv0, int type1, OVMV mv1)

{

    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &rcn_ctx->rcn_funcs.mc_l;
    struct MCFunctions *mc_c = &rcn_ctx->rcn_funcs.mc_c;
    /* FIXME derive ref_idx */
    uint8_t ref_idx_0 = mv0.ref_idx;
    uint8_t ref_idx_1 = mv1.ref_idx;

    OVPicture *ref0 = type0 == 1 ? inter_ctx->rpl0[ref_idx_0]: inter_ctx->rpl1[ref_idx_0];
    OVPicture *ref1 = type1 == 1 ? inter_ctx->rpl0[ref_idx_1]: inter_ctx->rpl1[ref_idx_1];

    /* TMP buffers for edge emulation
     * FIXME use tmp buffers in local contexts
     */
    uint16_t edge_buff0[RCN_CTB_SIZE];
    uint16_t edge_buff1[RCN_CTB_SIZE];
    uint16_t edge_buff0_1[RCN_CTB_SIZE];
    uint16_t edge_buff1_1[RCN_CTB_SIZE];
    int16_t tmp_buff0[RCN_CTB_SIZE];
    int16_t tmp_buff1[RCN_CTB_SIZE];

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    /* FIXME we should not need ctb_x/y
     * it could be retrieved from position in frame buff
     */
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width[0],
                  ref0->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width[0],
                  ref1->frame->height[0], 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    const struct OVBuffInfo ref0_b = derive_ref_buf_y(ref0, mv0, pos_x, pos_y, edge_buff0,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_b = derive_ref_buf_y(ref1, mv1, pos_x, pos_y, edge_buff1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    uint8_t prec_x0 = (mv0.x) & 0xF;
    uint8_t prec_y0 = (mv0.y) & 0xF;
    uint8_t prec_x1 = (mv1.x) & 0xF;
    uint8_t prec_y1 = (mv1.y) & 0xF;

    if(inter_ctx->prec_amvr == 3){
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    dst.y  += x0 + y0 * dst.stride;
    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff0, ref0_b.y, ref0_b.stride, pu_h, prec_x0, prec_y0, pu_w);
    mc_l->bidir0[prec_1_mc_type][log2_pu_w - 1](tmp_buff1, ref1_b.y, ref1_b.stride, pu_h, prec_x1, prec_y1, pu_w);

    int16_t* weight;
    int step_x, step_y;
    rcn_gpm_weights_and_steps(split_dir, log2_pu_w, log2_pu_h, &step_x, &step_y, &weight, 0);
    put_gpm_pel_bi_pixels(dst.y, RCN_CTB_STRIDE, tmp_buff1, MAX_PB_SIZE, tmp_buff0, pu_h, prec_x1, prec_y1, pu_w,
                            step_x, step_y, weight);

    rcn_ctx->rcn_funcs.lmcs_reshape(dst.y, RCN_CTB_STRIDE, ctudec->lmcs_info.lmcs_lut_fwd_luma, pu_w, pu_h);


    const struct OVBuffInfo ref0_c = derive_ref_buf_c(ref0, mv0,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff0, edge_buff0_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);

    const struct OVBuffInfo ref1_c = derive_ref_buf_c(ref1, mv1,
                                                      pos_x >> 1, pos_y >> 1,
                                                      edge_buff1, edge_buff1_1,
                                                      log2_pu_w, log2_pu_h, log2_ctb_s);
    prec_x0 = (mv0.x) & 0x1F;
    prec_y0 = (mv0.y) & 0x1F;

    prec_x1 = (mv1.x) & 0x1F;
    prec_y1 = (mv1.y) & 0x1F;

    prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    int16_t* ref_data0 = tmp_buff0;
    int16_t* ref_data1 = tmp_buff0 + MAX_PB_SIZE / 2;
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h >> 1, prec_x0, prec_y0, pu_w >> 1);

    int16_t* ref_data01 = tmp_buff1;
    int16_t* ref_data11 = tmp_buff1 + MAX_PB_SIZE / 2;
    mc_c->bidir0[prec_1_mc_type][log2_pu_w - 1](ref_data01, ref1_c.cb, ref1_c.stride_c, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);
    mc_c->bidir0[prec_1_mc_type][log2_pu_w - 1](ref_data11, ref1_c.cr, ref1_c.stride_c, pu_h >> 1, prec_x1, prec_y1, pu_w >> 1);

    rcn_gpm_weights_and_steps(split_dir, log2_pu_w, log2_pu_h, &step_x, &step_y, &weight, 1);
    put_gpm_pel_bi_pixels(dst.cb, RCN_CTB_STRIDE, ref_data01, MAX_PB_SIZE, ref_data0, 
                            pu_h >> 1, prec_x1, prec_y1, pu_w >> 1, step_x, step_y, weight);
    put_gpm_pel_bi_pixels(dst.cr, RCN_CTB_STRIDE, ref_data11, MAX_PB_SIZE, ref_data1,
                          pu_h >> 1, prec_x1, prec_y1, pu_w >> 1, step_x, step_y, weight);
}


void
rcn_gpm_b(OVCTUDec *const ctudec, struct VVCGPM* gpm_ctx,
            int x0, int y0, int log2_pb_w, int log2_pb_h)
{   

    int type0 = gpm_ctx->inter_dir0 ;
    int type1 = gpm_ctx->inter_dir1 ;
    struct OVBuffInfo dst = ctudec->rcn_ctx.ctu_buff;

    rcn_gpm_mc(ctudec, dst, gpm_ctx->split_dir, x0, y0, log2_pb_w, log2_pb_h, 
                          type0, gpm_ctx->mv0, type1, gpm_ctx->mv1);

}
