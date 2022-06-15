/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include <string.h>
#include <stdlib.h>
#include "ovdpb.h"
#include "dec_structures.h"
#include "ovdefs.h"
#include "ovutils.h"
#include "ctudec.h"
#include "rcn_structures.h"
#include "rcn.h"
#include "drv.h"
#include "drv_utils.h"

#define MAX_PB_SIZE 128

#define REF_PADDING_C 1
#define EPEL_EXTRA_AFTER  2
#define EPEL_EXTRA REF_PADDING_C + EPEL_EXTRA_AFTER

#define REF_PADDING_L 3
#define QPEL_EXTRA_BEFORE 3
#define QPEL_EXTRA_AFTER  4
#define QPEL_EXTRA REF_PADDING_L + QPEL_EXTRA_AFTER

#define PROF_BUFF_STRIDE 128
#define PROF_BUFF_PADD_H 1
#define PROF_BUFF_PADD_W 1
#define SB_H 4
#define SB_W 4

#define BCW_DENOM 2
#define RPR_NO_SCALE (1 << RPR_SCALE_BITS)

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

static int16_t bcw_weights[5] = { -2, 3, 4, 5, 10 };

struct OVDMV {
    int32_t x;
    int32_t y;
};

static OVMV
clip_mv(int pos_x, int pos_y, int pic_w, int pic_h, int pb_w, int pb_h, OVMV mv)
{
    int x_max  = (pic_w + 2 - pos_x) << 4;
    int y_max  = (pic_h + 2 - pos_y) << 4;
    int x_min  = -((pb_w + 3 + pos_x) << 4);
    int y_min  = -((pb_h + 3 + pos_y) << 4);

    mv.x = ov_clip(mv.x, x_min, x_max);
    mv.y = ov_clip(mv.y, y_min, y_max);

    return mv;
}

static OVMV
clip_mv_rpr(int pos_x, int pos_y, int pic_w, int pic_h, int pb_w, int pb_h, OVMV mv)
{
    int prec_x = mv.x & 0xF;
    int prec_y = mv.y & 0xF;
    int x_max  = (pic_w + 3 - pos_x) << 4;
    int y_max  = (pic_h + 3 - pos_y) << 4;
    int x_min  = -((pb_w + 3 + pos_x) << 4);
    int y_min  = -((pb_h + 3 + pos_y) << 4);

    mv.x = ov_clip(mv.x, x_min + prec_x, x_max + prec_x);
    mv.y = ov_clip(mv.y, y_min + prec_y, y_max + prec_y);

    return mv;
}

static void
gpm_weights_and_steps(int split_dir, int log2_pb_w_l, int log2_pb_h_l, int* step_x, int* step_y,
                      int16_t** weight, int cr_scale);

static void
rcn_inter_synchronization(const OVPicture *ref_pic, int ref_pos_x, int ref_pos_y, int pu_w, int pu_h, int log2_ctu_s)
{
    const int pic_w = ref_pic->frame->width;
    const int pic_h = ref_pic->frame->height;

    /*Frame thread synchronization to ensure data is available
    */
    int nb_ctb_pic_w = (pic_w + ((1 << log2_ctu_s) - 1)) >> log2_ctu_s;
    int nb_ctb_pic_h = (pic_h + ((1 << log2_ctu_s) - 1)) >> log2_ctu_s;
    int tl_ctu_y = ov_clip((ref_pos_y - QPEL_EXTRA_BEFORE) >> log2_ctu_s, 0, nb_ctb_pic_h - 1);
    int tl_ctu_x = ov_clip((ref_pos_x - QPEL_EXTRA_BEFORE) >> log2_ctu_s, 0, nb_ctb_pic_w - 1);
    int br_ctu_y = ov_clip((ref_pos_y + QPEL_EXTRA_AFTER + pu_h) >> log2_ctu_s, 0, nb_ctb_pic_h - 1);
    int br_ctu_x = ov_clip((ref_pos_x + QPEL_EXTRA_AFTER + pu_w) >> log2_ctu_s, 0, nb_ctb_pic_w - 1);
    FrameSynchroFunction sync_func = (FrameSynchroFunction)atomic_load(ref_pic->sync.func);

    sync_func(ref_pic, tl_ctu_x, tl_ctu_y, br_ctu_x, br_ctu_y);
}

static void
emulate_block_border(OVSample *buf, const OVSample *src,
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
        memcpy(buf, src, w * sizeof(OVSample));
        buf += buf_linesize;
    }

    // copy existing part
    for (; y < end_y; y++) {
        memcpy(buf, src, w * sizeof(OVSample));
        src += src_linesize;
        buf += buf_linesize;
    }

    // bottom
    src -= src_linesize;
    for (; y < block_h; y++) {
        memcpy(buf, src, w * sizeof(OVSample));
        buf += buf_linesize;
    }

    buf -= block_h * buf_linesize + start_x;

    while (block_h--) {
        OVSample *bufp = (OVSample *) buf;

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
    emulate_edge |= 4 * ((pb_x + pb_w + EPEL_EXTRA_AFTER) > pic_w);
    emulate_edge |= 8 * ((pb_y + pb_h + EPEL_EXTRA_AFTER) > pic_h);
    return emulate_edge;
}

static uint8_t
test_for_edge_emulation(int pb_x, int pb_y, int pic_w, int pic_h,
                        int pu_w, int pu_h)
{
    uint8_t emulate_edge = 0;
    emulate_edge  =      pb_x - REF_PADDING_L < 0;
    emulate_edge |= 2 * (pb_y - REF_PADDING_L < 0);
    emulate_edge |= 4 * (pb_x >= pic_w);
    emulate_edge |= 8 * (pb_y >= pic_h);
    emulate_edge |= 4 * ((pb_x + pu_w + QPEL_EXTRA_AFTER) > pic_w);
    emulate_edge |= 8 * ((pb_y + pu_h + QPEL_EXTRA_AFTER) > pic_h);
    return emulate_edge;
}


static uint8_t
check_identical_motion(struct InterDRVCtx *const inter_ctx, int inter_dir, const OVMV mv0, const OVMV mv1,
                       uint8_t ref_idx0, uint8_t ref_idx1)
{
    if (inter_dir != 3)
        return 0;

    uint16_t poc0 = inter_ctx->rpl0[ref_idx0] ? inter_ctx->rpl0[ref_idx0]->poc : -1;
    uint16_t poc1 = inter_ctx->rpl1[ref_idx1] ? inter_ctx->rpl1[ref_idx1]->poc : -1;

    return (poc0 == poc1 && mv0.x == mv1.x && mv0.y == mv1.y);
}


static struct OVBuffInfo
derive_ref_buf_c(const OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                 OVSample *edge_buff0, OVSample *edge_buff1,
                 int log2_pu_w, int log2_pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    OVSample *const ref_cb  = (OVSample *) ref_pic->frame->data[1];
    OVSample *const ref_cr  = (OVSample *) ref_pic->frame->data[2];

    int src_stride = ref_pic->frame->linesize[1]/sizeof(OVSample);
    const int pic_w = ref_pic->frame->width >> 1;
    const int pic_h = ref_pic->frame->height >> 1;

    /*FIXME check buff side derivation */
    int ref_pos_x = pos_x + (mv.x >> 5);
    int ref_pos_y = pos_y + (mv.y >> 5);

    const int pu_w = (1 << log2_pu_w) >> 1;
    const int pu_h = (1 << log2_pu_h) >> 1;

    OVSample *src_cb  = &ref_cb[ref_pos_x + ref_pos_y * src_stride];
    OVSample *src_cr  = &ref_cr[ref_pos_x + ref_pos_y * src_stride];

    uint8_t emulate_edge = test_for_edge_emulation_c(ref_pos_x, ref_pos_y, pic_w, pic_h,
                                                     pu_w, pu_h);;

    if (emulate_edge) {
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
padd_dmvr_c(OVSample *const _ref, int16_t stride, uint8_t pu_w, uint8_t pu_h)
{
    int i;

    OVSample *ref = _ref;
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
padd_dmvr(OVSample *const _ref, int16_t stride, uint8_t pu_w, uint8_t pu_h)
{
    int i;

    OVSample *ref = _ref;
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
derive_ref_buf_y(OVPicture *const ref_pic, OVMV mv, int pos_x, int pos_y,
                 OVSample *edge_buff, int log2_pu_w, int log2_pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    OVSample *const ref_y  = (OVSample *) ref_pic->frame->data[0];

    int src_stride = ref_pic->frame->linesize[0]/sizeof(OVSample);

    int ref_pos_x = pos_x + (mv.x >> 4);
    int ref_pos_y = pos_y + (mv.y >> 4);

    int pu_w = 1 << log2_pu_w;
    int pu_h = 1 << log2_pu_h;

    const int pic_w = ref_pic->frame->width;
    const int pic_h = ref_pic->frame->height;

    uint8_t emulate_edge = test_for_edge_emulation(ref_pos_x, ref_pos_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    /*Frame thread synchronization to ensure data is available
    */
    rcn_inter_synchronization(ref_pic, ref_pos_x, ref_pos_y, pu_w, pu_h, log2_ctu_s);

    if (emulate_edge) {
        const OVSample *src_y  = &ref_y[ref_pos_x + ref_pos_y * src_stride];
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
                      OVSample *edge_buff, int pu_w, int pu_h, int log2_ctu_s)
{
    struct OVBuffInfo ref_buff;
    OVSample *const ref_y  = (OVSample *) ref_pic->frame->data[0];

    const int pic_w = ref_pic->frame->width;
    const int pic_h = ref_pic->frame->height;

    int src_stride = ref_pic->frame->linesize[0]/sizeof(OVSample);

    OVMV mv_clipped = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_pos_x = pos_x + (mv_clipped.x >> 4);
    int ref_pos_y = pos_y + (mv_clipped.y >> 4);

    const OVSample *src_y = &ref_y[ref_pos_x + ref_pos_y * src_stride];

    int src_off  = (REF_PADDING_L * src_stride) + (REF_PADDING_L);
    int buff_off = (REF_PADDING_L * RCN_CTB_STRIDE) + (REF_PADDING_L);

    int cpy_w = pu_w + QPEL_EXTRA;
    int cpy_h = pu_h + QPEL_EXTRA;

    int start_pos_x = ref_pos_x - REF_PADDING_L;
    int start_pos_y = ref_pos_y - REF_PADDING_L;

    /* Frame thread synchronization to ensure data is available */
    rcn_inter_synchronization(ref_pic, ref_pos_x, ref_pos_y, pu_w, pu_h, log2_ctu_s);

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
                      OVSample *edge_buff0, OVSample *edge_buff1, int pu_w, int pu_h)
{
    struct OVBuffInfo ref_buff;
    OVSample *const ref_cb  = (OVSample *) ref_pic->frame->data[1];
    OVSample *const ref_cr  = (OVSample *) ref_pic->frame->data[2];

    const int pic_w = ref_pic->frame->width >> 1;
    const int pic_h = ref_pic->frame->height >> 1;

    int src_stride = ref_pic->frame->linesize[1]/sizeof(OVSample);

    OVMV mv_clipped = clip_mv(pos_x << 1, pos_y << 1, pic_w << 1, pic_h << 1, pu_w << 1, pu_h << 1, mv);

    int ref_pos_x = pos_x + (mv_clipped.x >> 5);
    int ref_pos_y = pos_y + (mv_clipped.y >> 5);

    const OVSample *src_cb = &ref_cb[ref_pos_x + ref_pos_y * src_stride];
    const OVSample *src_cr = &ref_cr[ref_pos_x + ref_pos_y * src_stride];

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

    return ref_buff;
}

static void
motion_compensation_b_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pu_w, uint8_t log2_pu_h,
                        OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    uint8_t use_alt_filter = mv0.prec_amvr == MV_PRECISION_HALF;
    uint8_t use_bcw = mv0.bcw_idx_plus1 != 0 && mv0.bcw_idx_plus1 != 3;

    int16_t wt0, wt1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx1];
    if (!ref0 || !ref1) return;

    OVSample *edge_buff0 = (OVSample *)rcn_ctx->data.edge_buff0;
    OVSample *edge_buff1 = (OVSample *)rcn_ctx->data.edge_buff1;

    int16_t *tmp_buff = (int16_t *) rcn_ctx->data.tmp_buff;

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width,
                  ref0->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width,
                  ref1->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv1);


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

    if (use_alt_filter) {
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    dst.y  += x0 + y0 * dst.stride;

    mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride,
                                                pu_h, prec_x0, prec_y0, pu_w);

    if (!use_bcw && !inter_ctx->weighted_pred_status) {
        mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                    tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
    } else if (use_bcw) {
        wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        wt0 = 8 - wt1;
        mc_l->bidir_w[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                     tmp_buff,  MAX_PB_SIZE, pu_h, BCW_DENOM, wt0, wt1, 0, 0,
                                                     prec_x1, prec_y1, pu_w);
    } else {
        int16_t offset0 = inter_ctx->wp_info0[ref_idx0].offset_y;
        int16_t offset1 = inter_ctx->wp_info1[ref_idx1].offset_y;
        wt1 = inter_ctx->wp_info1[ref_idx1].weight_y;
        wt0 = inter_ctx->wp_info0[ref_idx0].weight_y;
        mc_l->bidir_w[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                     tmp_buff,  MAX_PB_SIZE, pu_h, inter_ctx->weighted_denom, wt0, wt1, offset0, offset1,
                                                     prec_x1, prec_y1, pu_w);
    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride, ctudec->lmcs_info.luts,
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
rcn_dmvr_sad_16(const uint16_t *ref0, const uint16_t *ref1,
                int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
    uint64_t sum = 0;
    int i;
    for (i = 0; i < (pb_h >> 1); ++i) {
        sum += abs((int16_t)ref0[0]  - (int16_t)ref1[0]);
        sum += abs((int16_t)ref0[1]  - (int16_t)ref1[1]);
        sum += abs((int16_t)ref0[2]  - (int16_t)ref1[2]);
        sum += abs((int16_t)ref0[3]  - (int16_t)ref1[3]);
        sum += abs((int16_t)ref0[4]  - (int16_t)ref1[4]);
        sum += abs((int16_t)ref0[5]  - (int16_t)ref1[5]);
        sum += abs((int16_t)ref0[6]  - (int16_t)ref1[6]);
        sum += abs((int16_t)ref0[7]  - (int16_t)ref1[7]);
        sum += abs((int16_t)ref0[8]  - (int16_t)ref1[8]);
        sum += abs((int16_t)ref0[9]  - (int16_t)ref1[9]);
        sum += abs((int16_t)ref0[10] - (int16_t)ref1[10]);
        sum += abs((int16_t)ref0[11] - (int16_t)ref1[11]);
        sum += abs((int16_t)ref0[12] - (int16_t)ref1[12]);
        sum += abs((int16_t)ref0[13] - (int16_t)ref1[13]);
        sum += abs((int16_t)ref0[14] - (int16_t)ref1[14]);
        sum += abs((int16_t)ref0[15] - (int16_t)ref1[15]);

        ref0 += dmvr_stride << 1;
        ref1 += dmvr_stride << 1;
    }
    return sum;
}

static uint64_t
rcn_dmvr_sad_8(const uint16_t *ref0, const uint16_t *ref1,
               int16_t dmvr_stride, int16_t pb_w, int16_t pb_h)
{
    uint64_t sum = 0;
    int i;
    for (i = 0; i < (pb_h >> 1); ++i) {
        sum += abs((int16_t)ref0[0] - (int16_t)ref1[0]);
        sum += abs((int16_t)ref0[1] - (int16_t)ref1[1]);
        sum += abs((int16_t)ref0[2] - (int16_t)ref1[2]);
        sum += abs((int16_t)ref0[3] - (int16_t)ref1[3]);
        sum += abs((int16_t)ref0[4] - (int16_t)ref1[4]);
        sum += abs((int16_t)ref0[5] - (int16_t)ref1[5]);
        sum += abs((int16_t)ref0[6] - (int16_t)ref1[6]);
        sum += abs((int16_t)ref0[7] - (int16_t)ref1[7]);

        ref0 += dmvr_stride << 1;
        ref1 += dmvr_stride << 1;
    }
    return sum;
}

static uint8_t
dmvr_compute_sads_16(const uint16_t *ref0, const uint16_t *ref1,
                     uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const uint16_t *const ref0_start = ref0;
    const uint16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx = 12;

    for (idx = 0; idx < 12; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_16(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }

    for (idx = 13; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_16(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }
    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        if (sad_array[idx] < min_cost || (idx == 12 && sad_array[idx] <= min_cost)) {
            min_cost = sad_array[idx];
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

static uint8_t
dmvr_compute_sads_8(const uint16_t *ref0, const uint16_t *ref1,
                    uint64_t *sad_array, int sb_w, int sb_h)
{
    const int32_t stride_l0 = 128 + 4;
    const int32_t stride_l1 = 128 + 4;

    const uint16_t *const ref0_start = ref0;
    const uint16_t *const ref1_start = ref1;
    uint64_t min_cost = (uint64_t) -1;

    uint8_t idx;
    uint8_t dmvr_idx = 12;

    for (idx = 0; idx < 12; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_8(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }

    for (idx = 13; idx < DMVR_NB_IDX; ++idx) {
        ref0 = ref0_start + (int16_t)dmvr_mv_x[idx]
                          + (int16_t)dmvr_mv_y[idx] * stride_l0;

        ref1 = ref1_start - (int16_t)dmvr_mv_x[idx]
                          - (int16_t)dmvr_mv_y[idx] * stride_l1;

        sad_array[idx] = rcn_dmvr_sad_8(ref0, ref1, stride_l1,
                                     sb_w, sb_h);
    }
    for (idx = 0; idx < DMVR_NB_IDX; ++idx) {
        if (sad_array[idx] < min_cost || (idx == 12 && sad_array[idx] <= min_cost)) {
            min_cost = sad_array[idx];
            dmvr_idx = idx;
        }
    }

    return dmvr_idx;
}

/* FIXME understand this */
static inline int32_t
div_for_maxq7(int64_t num, int64_t den)
{
    uint64_t sign = -(num < 0);
    int32_t q = 0;
    uint64_t msk;

    num = (num ^ sign) - sign;
    den <<= 3;

    msk = -(num >= den);
    num -= msk & den;
    q |= msk & 0x1;

    q   <<= 1;
    den >>= 1;

    msk = -(num >= den);
    num -= msk & den;
    q |= msk & 0x1;

    q   <<= 1;
    den >>= 1;

    msk = -(num >= den);
    q |= msk & 0x1;

    return (q ^ sign) - sign;
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
            int64_t num = (int64_t)(sad_lc[1]  << DMVR_SB_PXL_LVL) - (int64_t)(sad_lc[3] << DMVR_SB_PXL_LVL);
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

static void
extend_bdof_grad(int16_t *dst_grad, int16_t grad_stride, int16_t pb_w, int16_t pb_h)
{
    int16_t       *dst = dst_grad + grad_stride;
    const int16_t *ref = dst + 1;

    int16_t       *dst_lst = (int16_t*)ref + pb_w;
    const int16_t *ref_lst = dst_lst - 1;

    int j;

    /* Copy or extend left and right column*/
    for (j = 0; j < pb_h; ++j) {
        dst[0]     = ref[0];
        dst_lst[0] = ref_lst[0];

        ref     += grad_stride;
        dst     += grad_stride;
        ref_lst += grad_stride;
        dst_lst += grad_stride;
    }

    /* Copy or extend upper and lower ref_line */
    dst = dst_grad;
    ref = dst + grad_stride;
    ref_lst = dst + (pb_h) * grad_stride;
    dst_lst = (int16_t*)ref_lst + grad_stride;

    memcpy(dst,     ref    , sizeof(*ref) * (pb_w + 2));
    memcpy(dst_lst, ref_lst, sizeof(*ref) * (pb_w + 2));
}

static uint8_t
rcn_dmvr_mv_refine(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_pu_w, uint8_t log2_pu_h,
                   OVMV *mv0, OVMV *mv1, uint8_t ref_idx0, uint8_t ref_idx1, uint8_t
                   apply_bdof)
{
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    struct DMVRFunctions *dmvr = &ctudec->rcn_funcs.dmvr;
    struct BDOFFunctions *bdof = &ctudec->rcn_funcs.bdof;
    uint8_t use_alt_filter = mv0->prec_amvr == MV_PRECISION_HALF;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx1];
    if (!ref0 | !ref1) return 0;

    OVSample edge_buff0[RCN_CTB_SIZE];
    OVSample edge_buff1[RCN_CTB_SIZE];

    /*FIXME permit smaller stride to reduce tables */
    uint16_t ref_dmvr0[(16 + 2 * DMVR_REF_PADD) * (128 + 2 * DMVR_REF_PADD)] = {0};
    uint16_t ref_dmvr1[(16 + 2 * DMVR_REF_PADD) * (128 + 2 * DMVR_REF_PADD)] = {0};

    int16_t *tmp_buff = (int16_t *) rcn_ctx->data.tmp_buff;

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pu_w = 1 << log2_pu_w;
    int pu_h = 1 << log2_pu_h;

    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0->x = ov_clip(mv0->x, MV_MIN, MV_MAX);
    mv0->y = ov_clip(mv0->y, MV_MIN, MV_MAX);

    mv1->x = ov_clip(mv1->x, MV_MIN, MV_MAX);
    mv1->y = ov_clip(mv1->y, MV_MIN, MV_MAX);

    OVMV tmp0 = *mv0;
    OVMV tmp1 = *mv1;

    struct OVBuffInfo ref0_b = derive_dmvr_ref_buf_y(ref0, *mv0, pos_x, pos_y, edge_buff0,
                                                     pu_w, pu_h, ctudec->part_ctx->log2_ctu_s);

    struct OVBuffInfo ref1_b = derive_dmvr_ref_buf_y(ref1, *mv1, pos_x, pos_y, edge_buff1,
                                                     pu_w, pu_h, ctudec->part_ctx->log2_ctu_s);

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
    mc_l->bilinear[prec_0_mc_type][log2_pu_w-1](ref_dmvr0, dmvr_stride,
                                                ref0_b.y -2 -2 * ref0_b.stride,
                                                ref0_b.stride, pu_h + 4,
                                                prec_x0, prec_y0, pu_w + 4);

    mc_l->bilinear[prec_1_mc_type][log2_pu_w-1](ref_dmvr1, dmvr_stride,
                                                ref1_b.y -2 -2 * ref1_b.stride,
                                                ref1_b.stride, pu_h + 4,
                                                prec_x1, prec_y1, pu_w + 4);

    /* Compute SAD on center part */
    dmvr_sad = dmvr->sad[pu_w==16](ref_dmvr0 + 2 + 2 * dmvr_stride,
                                   ref_dmvr1 + 2 + 2 * dmvr_stride,
                                   dmvr_stride, pu_w, pu_h);

    min_cost = (dmvr_sad - (dmvr_sad >> 2));

    /* skip MV refinement if cost is small or zero */
    if (min_cost >= (pu_w * pu_h)) {
        uint64_t sad[25];
        sad[12] = min_cost;
        uint8_t dmvr_idx = dmvr->computeSB[pu_w==16](ref_dmvr0 + 2 + 2 * dmvr_stride,
                                                     ref_dmvr1 + 2 + 2 * dmvr_stride,
                                                     sad, pu_w, pu_h);

        int32_t delta_h = dmvr_mv_x[dmvr_idx] * 16;
        int32_t delta_v = dmvr_mv_y[dmvr_idx] * 16;

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

    padd_dmvr(ref0_b.y, ref0_b.stride, pu_w, pu_h);
    padd_dmvr(ref1_b.y, ref1_b.stride, pu_w, pu_h);

    prec_x0 = (mv0->x) & 0xF;
    prec_y0 = (mv0->y) & 0xF;

    prec_x1 = (mv1->x) & 0xF;
    prec_y1 = (mv1->y) & 0xF;

    if (use_alt_filter) {
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

    ref0_b.y += delta_h2 + ref0_b.stride * delta_v2;
    ref1_b.y += delta_h3 + ref1_b.stride * delta_v3;

    uint8_t disable_bdof = apply_bdof ? min_cost < 2 * (pu_w * pu_h) : 1;

    if (disable_bdof) {
        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride,
                                                    pu_h, prec_x0, prec_y0, pu_w);
        mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                    tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
    } else {
        int16_t grad_x0[(16 + 2) * (16 + 3)];
        int16_t grad_y0[(16 + 2) * (16 + 3)];
        int16_t grad_x1[(16 + 2) * (16 + 3)];
        int16_t grad_y1[(16 + 2) * (16 + 3)];

        int16_t *tmp_buff1 = (int16_t*) rcn_ctx->data.tmp_buff1;

        int16_t ref_stride = 128;
        int16_t grad_stride = pu_w + 2;

        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff + 128 + 1,
                                                    ref0_b.y, ref0_b.stride, pu_h,
                                                    prec_x0, prec_y0, pu_w);

        mc_l->bidir0[prec_1_mc_type][log2_pu_w - 1](tmp_buff1 + 128 + 1,
                                                    ref1_b.y, ref1_b.stride, pu_h,
                                                    prec_x1, prec_y1, pu_w);

        /* Padding for grad derivation */
        bdof->extend_bdof_buff(ref0_b.y, tmp_buff, ref0_b.stride, pu_w, pu_h, prec_x0 >> 3, prec_y0 >> 3);
        bdof->extend_bdof_buff(ref1_b.y, tmp_buff1, ref1_b.stride, pu_w, pu_h, prec_x1 >> 3, prec_y1 >> 3);

        bdof->grad(tmp_buff, ref_stride, pu_w, pu_h, grad_stride,
                   grad_x0 + grad_stride + 1, grad_y0 + grad_stride + 1);

        bdof->grad(tmp_buff1, ref_stride, pu_w, pu_h, grad_stride,
                   grad_x1 + grad_stride + 1, grad_y1 + grad_stride + 1);

        /* Grad padding */
        extend_bdof_grad(grad_x0, grad_stride, pu_w, pu_h);
        extend_bdof_grad(grad_y0, grad_stride, pu_w, pu_h);
        extend_bdof_grad(grad_x1, grad_stride, pu_w, pu_h);
        extend_bdof_grad(grad_y1, grad_stride, pu_w, pu_h);

        /* Reference padding overwrite for weights derivation */
        extend_bdof_grad(tmp_buff, ref_stride, pu_w, pu_h);
        extend_bdof_grad(tmp_buff1, ref_stride, pu_w, pu_h);

        /* Split into 4x4 subblocks for BDOF computation */
        bdof->rcn_bdof(bdof, dst.y, dst.stride, tmp_buff + 128 + 1, tmp_buff1 + 128 + 1,
                       ref_stride, grad_x0, grad_y0, grad_x1, grad_y1,
                       grad_stride, pu_w, pu_h);

    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride,
                                           ctudec->lmcs_info.luts,
                                           pu_w, pu_h);

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    uint8_t use_wp_c = inter_ctx->weighted_pred_status && (inter_ctx->wp_info1[ref_idx1].flag_c |
                                                           inter_ctx->wp_info0[ref_idx0].flag_c);
    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;

    OVSample *edge_buff0_1 = (OVSample *)rcn_ctx->data.edge_buff0;
    OVSample *edge_buff1_1 = (OVSample *)rcn_ctx->data.edge_buff1;

    pu_w = pu_w >> 1;
    pu_h = pu_h >> 1;

    struct OVBuffInfo ref0_c = derive_dmvr_ref_buf_c(ref0, tmp0,
                                                     pos_x >> 1, pos_y >> 1,
                                                     edge_buff0, edge_buff0_1,
                                                     pu_w, pu_h);

    struct OVBuffInfo ref1_c = derive_dmvr_ref_buf_c(ref1, tmp1,
                                                     pos_x >> 1, pos_y >> 1,
                                                     edge_buff1, edge_buff1_1,
                                                     pu_w, pu_h);

    padd_dmvr_c(ref0_c.cb, ref0_c.stride_c, pu_w, pu_h);
    padd_dmvr_c(ref0_c.cr, ref0_c.stride_c, pu_w, pu_h);

    padd_dmvr_c(ref1_c.cb, ref1_c.stride_c, pu_w, pu_h);
    padd_dmvr_c(ref1_c.cr, ref1_c.stride_c, pu_w, pu_h);

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

    ref0_c.cb += delta_h2 + (ref0_c.stride_c * delta_v2);
    ref1_c.cb += delta_h3 + (ref1_c.stride_c * delta_v3);

    ref0_c.cr += delta_h2 + (ref0_c.stride_c * delta_v2);
    ref1_c.cr += delta_h3 + (ref1_c.stride_c * delta_v3);

    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 2](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h, prec_x0, prec_y0, pu_w);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 2](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h, prec_x0, prec_y0, pu_w);

    if (!use_wp_c) {
    mc_c->bidir1[prec_1_mc_type][log2_pu_w - 2](dst.cb, dst.stride_c, ref1_c.cb, ref1_c.stride_c, ref_data0, pu_h, prec_x1, prec_y1, pu_w);
    mc_c->bidir1[prec_1_mc_type][log2_pu_w - 2](dst.cr, dst.stride_c, ref1_c.cr, ref1_c.stride_c, ref_data1, pu_h, prec_x1, prec_y1, pu_w);
    } else {
    }

    return disable_bdof;
}

struct PROFInfo
{
    int16_t dmv_scale_h_0[16];
    int16_t dmv_scale_v_0[16];
    int16_t dmv_scale_h_1[16];
    int16_t dmv_scale_v_1[16];
};

static void
rcn_bdof_mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
               uint8_t x0, uint8_t y0, uint8_t log2_pu_w, uint8_t log2_pu_h,
               OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1)
{
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    struct BDOFFunctions *bdof = &ctudec->rcn_funcs.bdof;
    uint8_t use_alt_filter = mv0.prec_amvr == MV_PRECISION_HALF;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx1];
    if (!ref0 | !ref1) return;

    OVSample edge_buff0[RCN_CTB_SIZE];
    OVSample edge_buff1[RCN_CTB_SIZE];

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width,
                  ref0->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width,
                  ref1->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv1);


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

    if (use_alt_filter) {
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
    int16_t grad_x0[(16 + 2) * (16 + 3)];
    int16_t grad_y0[(16 + 2) * (16 + 3)];
    int16_t grad_x1[(16 + 2) * (16 + 3)];
    int16_t grad_y1[(16 + 2) * (16 + 3)];

    int16_t ref_stride = 128;
    int16_t grad_stride = pu_w + 2;

    int pb_w = 1 << log2_pu_w;
    int pb_h = 1 << log2_pu_h;

    mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](ref_bdof0 + 128 + 1,
                                                ref0_b.y, ref0_b.stride, pu_h,
                                                prec_x0, prec_y0, pu_w);

    mc_l->bidir0[prec_1_mc_type][log2_pu_w - 1](ref_bdof1 + 128 + 1,
                                                ref1_b.y, ref1_b.stride, pu_h,
                                                prec_x1, prec_y1, pu_w);

    /* Padding for grad derivation */
    bdof->extend_bdof_buff(ref0_b.y, ref_bdof0, ref0_b.stride, pb_w, pb_h, prec_x0 >> 3, prec_y0 >> 3);
    bdof->extend_bdof_buff(ref1_b.y, ref_bdof1, ref1_b.stride, pb_w, pb_h, prec_x1 >> 3, prec_y1 >> 3);

    bdof->grad(ref_bdof0, ref_stride, pb_w, pb_h, grad_stride,
               grad_x0 + grad_stride + 1, grad_y0 + grad_stride + 1);

    bdof->grad(ref_bdof1, ref_stride, pb_w, pb_h, grad_stride,
               grad_x1 + grad_stride + 1, grad_y1 + grad_stride + 1);

    /* Grad padding */
    extend_bdof_grad(grad_x0, grad_stride, pb_w, pb_h);
    extend_bdof_grad(grad_y0, grad_stride, pb_w, pb_h);
    extend_bdof_grad(grad_x1, grad_stride, pb_w, pb_h);
    extend_bdof_grad(grad_y1, grad_stride, pb_w, pb_h);

    /* Reference padding overwrite for weights derivation */
    extend_bdof_grad(ref_bdof0, ref_stride, pb_w, pb_h);
    extend_bdof_grad(ref_bdof1, ref_stride, pb_w, pb_h);

    dst.y += x0 + y0 * dst.stride;

    /* Split into 4x4 subblocks for BDOF computation */
    bdof->rcn_bdof(bdof, dst.y, dst.stride, ref_bdof0 + 128 + 1, ref_bdof1 + 128 + 1,
                   ref_stride, grad_x0, grad_y0, grad_x1, grad_y1,
                   grad_stride, pb_w, pb_h);

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride,
                                           ctudec->lmcs_info.luts,
                                           pu_w, pu_h);
}

static void
prof_motion_compensation_b_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                             uint8_t x0, uint8_t y0,
                             uint8_t log2_pu_w, uint8_t log2_pu_h,
                             OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1,
                             uint8_t prof_dir, const struct PROFInfo *const prof_info)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    struct PROFFunctions *prof = &ctudec->rcn_funcs.prof;
    uint8_t use_alt_filter = mv0.prec_amvr == MV_PRECISION_HALF;
    uint8_t use_bcw = mv0.bcw_idx_plus1 != 0 && mv0.bcw_idx_plus1 != 3;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx1];
    if (!ref0 || !ref1) return;


    OVSample *edge_buff0 = (OVSample *)rcn_ctx->data.edge_buff0;
    OVSample *edge_buff1 = (OVSample *)rcn_ctx->data.edge_buff1;

    int16_t *tmp_buff1 = (int16_t*) rcn_ctx->data.tmp_buff1;

    int16_t *tmp_buff = (int16_t *) rcn_ctx->data.tmp_buff;


    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width,
                  ref0->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width,
                  ref1->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv1);


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

    if (use_alt_filter) {
        prec_x0 += (prec_x0 == 8) ? 8 : 0;
        prec_y0 += (prec_y0 == 8) ? 8 : 0;
        prec_x1 += (prec_x1 == 8) ? 8 : 0;
        prec_y1 += (prec_y1 == 8) ? 8 : 0;
    }

    dst.y  += x0 + y0 * dst.stride;

    uint8_t prec_0_mc_type = (prec_x0 > 0) + ((prec_y0 > 0) << 1);
    uint8_t prec_1_mc_type = (prec_x1 > 0) + ((prec_y1 > 0) << 1);

    if (prof_dir & 0x1) {
        int16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (PROF_BUFF_STRIDE + 2 * PROF_BUFF_PADD_W)];

        int16_t tmp_grad_x[16];
        int16_t tmp_grad_y[16];
        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_prof + PROF_BUFF_STRIDE + 1,
                                                    ref0_b.y, ref0_b.stride, pu_h,
                                                    prec_x0, prec_y0, pu_w);

        prof->extend_prof_buff(ref0_b.y, tmp_prof, ref0_b.stride, prec_x0 >> 3, prec_y0 >> 3);

        prof->grad(tmp_prof, PROF_BUFF_STRIDE, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

        prof->rcn1((OVSample *)tmp_buff, MAX_PB_SIZE, tmp_prof + PROF_BUFF_STRIDE + 1, PROF_BUFF_STRIDE,
                  tmp_grad_x, tmp_grad_y,
                  4, prof_info->dmv_scale_h_0, prof_info->dmv_scale_v_0);
    } else {
        mc_l->bidir0[prec_0_mc_type][log2_pu_w - 1](tmp_buff, ref0_b.y, ref0_b.stride,
                                                    pu_h, prec_x0, prec_y0, pu_w);
    }

    if (prof_dir & 0x2) {
        int16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (PROF_BUFF_STRIDE + 2 * PROF_BUFF_PADD_W)];

        int16_t tmp_grad_x[16];
        int16_t tmp_grad_y[16];

        mc_l->bidir0[prec_1_mc_type][log2_pu_w - 1](tmp_prof + PROF_BUFF_STRIDE + 1,
                                                    ref1_b.y, ref1_b.stride, pu_h,
                                                    prec_x1, prec_y1, pu_w);

        prof->extend_prof_buff(ref1_b.y, tmp_prof, ref1_b.stride, prec_x1 >> 3, prec_y1 >> 3);

        prof->grad(tmp_prof, PROF_BUFF_STRIDE, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

        prof->rcn1((OVSample *)tmp_buff1, MAX_PB_SIZE, tmp_prof + PROF_BUFF_STRIDE + 1, PROF_BUFF_STRIDE,
                  tmp_grad_x, tmp_grad_y, 4,
                  prof_info->dmv_scale_h_1, prof_info->dmv_scale_v_1);

        if (!use_bcw && !inter_ctx->weighted_pred_status) {
            prof->tmp_prof_mrg(dst.y, dst.stride, tmp_buff1, MAX_PB_SIZE,
                               tmp_buff, pu_h, 0, 0, pu_w);
        } else if (use_bcw) {
            int16_t wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
            int16_t wt0 = 8 - wt1;
            prof->tmp_prof_mrg_w(dst.y, dst.stride, tmp_buff1, MAX_PB_SIZE,
                                 tmp_buff, pu_h, 0, 0, pu_w, BCW_DENOM, wt1, wt0, 0, 0);
        } else {
            int16_t offset0 = inter_ctx->wp_info0[ref_idx0].offset_y;
            int16_t offset1 = inter_ctx->wp_info1[ref_idx1].offset_y;
            int16_t wt1 = inter_ctx->wp_info1[ref_idx1].weight_y;
            int16_t wt0 = inter_ctx->wp_info0[ref_idx0].weight_y;
            prof->tmp_prof_mrg_w(dst.y, dst.stride, tmp_buff1, MAX_PB_SIZE,
                                 tmp_buff, pu_h, 0, 0, pu_w, inter_ctx->weighted_denom, wt1, wt0, offset1, offset0);
        }

    } else {
        int wt0, wt1;
        if (!use_bcw && !inter_ctx->weighted_pred_status) {
            mc_l->bidir1[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                        tmp_buff, pu_h, prec_x1, prec_y1, pu_w);
        } else if (use_bcw) {
            wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
            wt0 = 8 - wt1;
            mc_l->bidir_w[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                         tmp_buff,  MAX_PB_SIZE, pu_h, BCW_DENOM, wt0, wt1, 0, 0,
                                                         prec_x1, prec_y1, pu_w);
        } else {
            int16_t offset0 = inter_ctx->wp_info0[ref_idx0].offset_y;
            int16_t offset1 = inter_ctx->wp_info1[ref_idx1].offset_y;
            wt1 = inter_ctx->wp_info1[ref_idx1].weight_y;
            wt0 = inter_ctx->wp_info0[ref_idx0].weight_y;
            mc_l->bidir_w[prec_1_mc_type][log2_pu_w - 1](dst.y, dst.stride, ref1_b.y, ref1_b.stride,
                                                         tmp_buff,  MAX_PB_SIZE, pu_h, inter_ctx->weighted_denom, wt0, wt1, offset0, offset1,
                                                         prec_x1, prec_y1, pu_w);
        }
    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride,
                                           ctudec->lmcs_info.luts,
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
    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;
    uint8_t use_wp_c = inter_ctx->weighted_pred_status && (inter_ctx->wp_info1[ref_idx1].flag_c |
                                                           inter_ctx->wp_info0[ref_idx0].flag_c);
    uint8_t use_bcw = mv0.bcw_idx_plus1 != 0 && mv0.bcw_idx_plus1 != 3;
    int16_t wt0, wt1;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx0];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx1];
    if (!ref0 || !ref1) return;

    OVSample *edge_buff0 = (OVSample *)rcn_ctx->data.edge_buff0;
    OVSample *edge_buff1 = (OVSample *)rcn_ctx->data.edge_buff1;

    OVSample *edge_buff0_1 = (OVSample *)rcn_ctx->data.edge_buff0_1;
    OVSample *edge_buff1_1 = (OVSample *)rcn_ctx->data.edge_buff1_1;

    int16_t *tmp_buff = (int16_t *) rcn_ctx->data.tmp_buff;

    /*FIXME we suppose here both refs possess the same size*/

    const int log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    mv0 = clip_mv(pos_x, pos_y, ref0->frame->width,
                  ref0->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv0);

    mv1 = clip_mv(pos_x, pos_y, ref1->frame->width,
                  ref1->frame->height, 1 << log2_pu_w, 1 << log2_pu_h, mv1);


    const int pu_w = 1 << (log2_pu_w - 1);
    const int pu_h = 1 << (log2_pu_h - 1);

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

    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 2](ref_data0, ref0_c.cb, ref0_c.stride_c, pu_h, prec_x0, prec_y0, pu_w);
    mc_c->bidir0[prec_0_mc_type][log2_pu_w - 2](ref_data1, ref0_c.cr, ref0_c.stride_c, pu_h, prec_x0, prec_y0, pu_w);

    if (!use_bcw  && !use_wp_c) {
        mc_c->bidir1[prec_1_mc_type][log2_pu_w - 2](dst.cb, dst.stride_c, ref1_c.cb, ref1_c.stride_c, ref_data0, pu_h, prec_x1, prec_y1, pu_w);
        mc_c->bidir1[prec_1_mc_type][log2_pu_w - 2](dst.cr, dst.stride_c, ref1_c.cr, ref1_c.stride_c, ref_data1, pu_h, prec_x1, prec_y1, pu_w);
    } else if (use_bcw) {
        wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        wt0 = 8 - wt1;
        mc_c->bidir_w[prec_1_mc_type][log2_pu_w - 2](dst.cb, dst.stride_c, ref1_c.cb, ref1_c.stride_c,
                                                     ref_data0,  MAX_PB_SIZE, pu_h, BCW_DENOM, wt0, wt1, 0, 0,
                                                     prec_x1, prec_y1, pu_w);
        mc_c->bidir_w[prec_1_mc_type][log2_pu_w - 2](dst.cr, dst.stride_c, ref1_c.cr, ref1_c.stride_c,
                                                     ref_data1,  MAX_PB_SIZE, pu_h, BCW_DENOM, wt0, wt1, 0, 0,
                                                     prec_x1, prec_y1, pu_w);
    } else {
        int16_t denom = inter_ctx->weighted_denom_c;
        int16_t wt1 = inter_ctx->wp_info1[ref_idx1].weight_cb;
        int16_t wt0 = inter_ctx->wp_info0[ref_idx0].weight_cb;
        int16_t o1 = inter_ctx->wp_info1[ref_idx1].offset_cb;
        int16_t o0 = inter_ctx->wp_info0[ref_idx0].offset_cb;
        mc_c->bidir_w[prec_1_mc_type][log2_pu_w - 2](dst.cb, dst.stride_c, ref1_c.cb, ref1_c.stride_c,
                                                     ref_data0,  MAX_PB_SIZE, pu_h, denom, wt0, wt1, o0, o1,
                                                     prec_x1, prec_y1, pu_w);
        wt1 = inter_ctx->wp_info1[ref_idx1].weight_cr;
        wt0 = inter_ctx->wp_info0[ref_idx0].weight_cr;
        o1 = inter_ctx->wp_info1[ref_idx1].offset_cr;
        o0 = inter_ctx->wp_info0[ref_idx0].offset_cr;
        mc_c->bidir_w[prec_1_mc_type][log2_pu_w - 2](dst.cr, dst.stride_c, ref1_c.cr, ref1_c.stride_c,
                                                     ref_data1,  MAX_PB_SIZE, pu_h, denom, wt0, wt1, o0, o1,
                                                     prec_x1, prec_y1, pu_w);
    }
}

static void
mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
      OVMV mv, uint8_t rpl_idx, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    uint8_t use_alt_filter = mv.prec_amvr == MV_PRECISION_HALF;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.y  += x0 + y0 * dst.stride;

    OVSample *tmp_buff = (OVSample *)rcn_ctx->data.tmp_buff;
    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;

    const OVSample *const ref0_y  = (OVSample *) frame0->data[0];

    int src_stride   = frame0->linesize[0]/sizeof(OVSample);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width;
    const int pic_h = frame0->height;

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if (use_alt_filter) {
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const OVSample *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];

    /*
     * Thread synchronization to ensure data is available before usage
     */
    rcn_inter_synchronization(ref_pic, ref_x, ref_y, pu_w, pu_h, log2_ctb_s);

    if (emulate_edge) {
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

    if (!inter_ctx->weighted_pred_status) {
        mc_l->unidir[prec_mc_type][log2_pu_w - 1](dst.y, dst.stride,
                                                  src_y, src_stride, pu_h,
                                                  prec_x, prec_y, pu_w);
    } else {
        int16_t offset0 = rpl_idx ? inter_ctx->wp_info1[ref_idx].offset_y : inter_ctx->wp_info0[ref_idx].offset_y;
        int wx = rpl_idx ? inter_ctx->wp_info1[ref_idx].weight_y : inter_ctx->wp_info0[ref_idx].weight_y;
        mc_l->unidir_w[prec_mc_type][log2_pu_w - 1](dst.y, dst.stride,
                                                    src_y, src_stride, pu_h,
                                                    inter_ctx->weighted_denom, wx, offset0,
                                                    prec_x, prec_y, pu_w);
    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride, ctudec->lmcs_info.luts,
                                           pu_w, pu_h);
}

static void
rcn_mcp_bidir0_l(OVCTUDec *const ctudec, uint16_t* dst, int dst_stride, int x0, int y0,
                 int log2_pu_w, int log2_pu_h, OVMV mv, uint8_t rpl_idx, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    uint8_t use_alt_filter = mv.prec_amvr == MV_PRECISION_HALF;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];
    if (!ref0 | !ref1) return;

    dst  += x0 + y0 * dst_stride;

    OVSample* tmp_buff = (OVSample*)rcn_ctx->data.tmp_buff;

    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;

    const OVSample *const ref0_y  = (OVSample *) frame0->data[0];

    int src_stride   = frame0->linesize[0] /sizeof(OVSample);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width;
    const int pic_h = frame0->height;

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if (use_alt_filter) {
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const OVSample *src_y  = &ref0_y [ ref_x + ref_y * src_stride];

    /*
     * Thread synchronization to ensure data is available before usage
     */
    rcn_inter_synchronization(ref_pic, ref_x, ref_y, pu_w, pu_h, log2_ctb_s);

    if (emulate_edge) {
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
    mc_l->bidir0[prec_mc_type][log2_pu_w - 1]((int16_t*)dst,
                                          src_y, src_stride, pu_h,
                                          prec_x, prec_y, pu_w);
}

static void
rcn_prof_mcp_l(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0,
               int log2_pu_w, int log2_pu_h,
               OVMV mv, uint8_t rpl_idx, uint8_t ref_idx,
               const int16_t *dmv_scale_h, const int16_t *dmv_scale_v)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    struct PROFFunctions *prof = &ctudec->rcn_funcs.prof;
    uint8_t use_alt_filter = mv.prec_amvr == MV_PRECISION_HALF;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.y  += x0 + y0 * dst.stride;

    OVSample *tmp_buff = (OVSample *)rcn_ctx->data.tmp_buff;

    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;

    const OVSample *const ref0_y  = (OVSample *) frame0->data[0];

    int src_stride   = frame0->linesize[0]/sizeof(OVSample);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width;
    const int pic_h = frame0->height;

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if (use_alt_filter) {
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const OVSample *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];

    /*
     * Thread synchronization to ensure data is available before usage
     */
    rcn_inter_synchronization(ref_pic, ref_x, ref_y, pu_w, pu_h, log2_ctb_s);

    uint16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (PROF_BUFF_STRIDE + 2 * PROF_BUFF_PADD_W)];

    int16_t tmp_grad_x[16];
    int16_t tmp_grad_y[16];

    if (emulate_edge) {
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
    mc_l->bidir0[prec_mc_type][log2_pu_w - 1](tmp_prof + PROF_BUFF_STRIDE + 1,
                                              src_y, src_stride, pu_h,
                                              prec_x, prec_y, pu_w);

    prof->extend_prof_buff(src_y, tmp_prof, src_stride, prec_x >> 3, prec_y >> 3);

    prof->grad(tmp_prof, PROF_BUFF_STRIDE, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

    if (!inter_ctx->weighted_pred_status) {
        prof->rcn0(dst.y, dst.stride, tmp_prof + PROF_BUFF_STRIDE + 1, PROF_BUFF_STRIDE, tmp_grad_x, tmp_grad_y,
                   4, dmv_scale_h, dmv_scale_v);
    } else {
        int16_t offset = rpl_idx ? inter_ctx->wp_info1[ref_idx].offset_y : inter_ctx->wp_info0[ref_idx].offset_y;
        int16_t wx = rpl_idx ? inter_ctx->wp_info1[ref_idx].weight_y : inter_ctx->wp_info0[ref_idx].weight_y;

        prof->rcn0_w(dst.y, dst.stride, tmp_prof + PROF_BUFF_STRIDE + 1, PROF_BUFF_STRIDE, tmp_grad_x, tmp_grad_y,
                     4, dmv_scale_h, dmv_scale_v, inter_ctx->weighted_denom, wx, offset);

    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride, ctudec->lmcs_info.luts,
                                           pu_w, pu_h);

}

static void
rcn_prof_mcp_bi_l(OVCTUDec *const ctudec, uint16_t* dst, uint16_t dst_stride, int x0, int y0,
                  int log2_pu_w, int log2_pu_h,
                  OVMV mv, uint8_t rpl_idx, uint8_t ref_idx,
                  const int16_t *dmv_scale_h, const int16_t *dmv_scale_v)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    struct PROFFunctions *prof = &ctudec->rcn_funcs.prof;
    uint8_t use_alt_filter = mv.prec_amvr == MV_PRECISION_HALF;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst  += x0 + y0 * dst_stride;

    OVSample *tmp_buff = (OVSample *)rcn_ctx->data.tmp_buff;

    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;

    const OVSample *const ref0_y  = (OVSample *) frame0->data[0];

    int src_stride   = frame0->linesize[0]/sizeof(OVSample);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    const int pic_w = frame0->width;
    const int pic_h = frame0->height;

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w, pu_h, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x   = (mv.x) & 0xF;
    uint8_t prec_y   = (mv.y) & 0xF;
    if (use_alt_filter) {
        prec_x += (prec_x == 8) ? 8 : 0;
        prec_y += (prec_y == 8) ? 8 : 0;
    }

    int prec_mc_type   = (prec_x  > 0) + ((prec_y > 0)   << 1);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, pic_w, pic_h,
                                                   pu_w, pu_h);;

    const OVSample *src_y  = &ref0_y [ ref_x       + ref_y        * src_stride];

    /*
     * Thread synchronization to ensure data is available before usage
     */
    rcn_inter_synchronization(ref_pic, ref_x, ref_y, pu_w, pu_h, log2_ctb_s);

    int16_t tmp_prof[(SB_H + 2 * PROF_BUFF_PADD_H) * (PROF_BUFF_STRIDE + 2 * PROF_BUFF_PADD_W)];

    int16_t tmp_grad_x[16];
    int16_t tmp_grad_y[16];

    if (emulate_edge) {
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
    mc_l->bidir0[prec_mc_type][log2_pu_w - 1](tmp_prof + PROF_BUFF_STRIDE + 1,
                                              src_y, src_stride, pu_h,
                                              prec_x, prec_y, pu_w);

    prof->extend_prof_buff(src_y, tmp_prof, src_stride, prec_x >> 3, prec_y >> 3);

    prof->grad(tmp_prof, PROF_BUFF_STRIDE, SB_W, SB_H, 4, tmp_grad_x, tmp_grad_y);

    prof->rcn1((OVSample*) dst, dst_stride, (uint16_t *)tmp_prof + PROF_BUFF_STRIDE + 1, PROF_BUFF_STRIDE, tmp_grad_x, tmp_grad_y,
              4, dmv_scale_h, dmv_scale_v);
}

static void
mcp_c(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
      OVMV mv, uint8_t rpl_idx, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;
    uint8_t use_wp_c = inter_ctx->weighted_pred_status && (rpl_idx ? inter_ctx->wp_info1[ref_idx].flag_c
                                                                   : inter_ctx->wp_info0[ref_idx].flag_c);

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    OVSample *tmp_buff = (OVSample *)rcn_ctx->data.tmp_buff;

    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;

    const OVSample *const ref0_cb = (OVSample *) frame0->data[1];
    const OVSample *const ref0_cr = (OVSample *) frame0->data[2];

    int src_stride_c = frame0->linesize[1]/sizeof(OVSample);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << (log2_pu_w - 1);
    const int pu_h = 1 << (log2_pu_h - 1);

    const int pic_w = frame0->width;
    const int pic_h = frame0->height;

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w << 1, pu_h << 1, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x_c = mv.x & 0x1F;
    uint8_t prec_y_c = mv.y & 0x1F;

    int prec_c_mc_type = (prec_x_c > 0) + ((prec_y_c > 0) << 1);

    const OVSample *src_cb = &ref0_cb[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];
    const OVSample *src_cr = &ref0_cr[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];

    uint8_t emulate_edge = test_for_edge_emulation_c(ref_x >> 1, ref_y >> 1, pic_w >> 1, pic_h >> 1,
                                                     pu_w, pu_h);

    if (emulate_edge) {
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cb - src_off), RCN_CTB_STRIDE, src_stride_c,
                             pu_w  + EPEL_EXTRA, pu_h + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cb = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    if (!use_wp_c) {
    mc_c->unidir[prec_c_mc_type][log2_pu_w - 2](dst.cb, dst.stride_c, src_cb, src_stride_c,
                                                pu_h, prec_x_c, prec_y_c, pu_w);
    } else {
        int16_t denom = inter_ctx->weighted_denom_c;
        int16_t wx = rpl_idx ? inter_ctx->wp_info1[ref_idx].weight_cb
                             : inter_ctx->wp_info0[ref_idx].weight_cb;
        int16_t ox = rpl_idx ? inter_ctx->wp_info1[ref_idx].offset_cb
                             : inter_ctx->wp_info0[ref_idx].offset_cb;
        mc_c->unidir_w[prec_c_mc_type][log2_pu_w - 2](dst.cb, dst.stride_c, src_cb, src_stride_c,
                                                      pu_h, denom, wx, ox, prec_x_c, prec_y_c, pu_w);
    }

    if (emulate_edge) {
        int src_off  = REF_PADDING_C * (frame0->linesize[1]/sizeof(OVSample)) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cr - src_off),
                             RCN_CTB_STRIDE, frame0->linesize[1]/sizeof(OVSample),
                             pu_w + EPEL_EXTRA, pu_h + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cr = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    if (!use_wp_c) {
    mc_c->unidir[prec_c_mc_type][log2_pu_w - 2](dst.cr, dst.stride_c, src_cr, src_stride_c,
                                                pu_h, prec_x_c, prec_y_c, pu_w);
    } else {
        int16_t denom = inter_ctx->weighted_denom_c;
        int16_t wx = rpl_idx ? inter_ctx->wp_info1[ref_idx].weight_cr
                             : inter_ctx->wp_info0[ref_idx].weight_cr;
        int16_t ox = rpl_idx ? inter_ctx->wp_info1[ref_idx].offset_cr
                             : inter_ctx->wp_info0[ref_idx].offset_cr;
        mc_c->unidir_w[prec_c_mc_type][log2_pu_w - 2](dst.cr, dst.stride_c, src_cr, src_stride_c,
                                                      pu_h, denom, wx, ox, prec_x_c, prec_y_c, pu_w);
    }
}

static void
rcn_mcp_bidir0_c(OVCTUDec *const ctudec, uint16_t* dst_cb, uint16_t* dst_cr, int dst_stride, int x0, int y0,
                 int log2_pu_w, int log2_pu_h, OVMV mv, uint8_t rpl_idx, uint8_t ref_idx)
{
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];
    if (!ref0 | !ref1) return;

    dst_cb += (x0 >> 1) + (y0 >> 1) * dst_stride;
    dst_cr += (x0 >> 1) + (y0 >> 1) * dst_stride;

    OVSample *tmp_buff = (OVSample *)rcn_ctx->data.tmp_buff;

    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;

    const OVSample *const ref0_cb = (OVSample *) frame0->data[1];
    const OVSample *const ref0_cr = (OVSample *) frame0->data[2];

    int src_stride_c = frame0->linesize[1] /sizeof(OVSample);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;

    const int pu_w = 1 << (log2_pu_w - 1);
    const int pu_h = 1 << (log2_pu_h - 1);

    const int pic_w = frame0->width;
    const int pic_h = frame0->height;

    mv = clip_mv(pos_x, pos_y, pic_w, pic_h, pu_w << 1, pu_h << 1, mv);

    int ref_x = pos_x + (mv.x >> 4);
    int ref_y = pos_y + (mv.y >> 4);

    uint8_t prec_x_c = (mv.x) & 0x1F;
    uint8_t prec_y_c = (mv.y) & 0x1F;

    int prec_c_mc_type = (prec_x_c > 0) + ((prec_y_c > 0) << 1);

    const OVSample *src_cb = &ref0_cb[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];
    const OVSample *src_cr = &ref0_cr[(ref_x >> 1) + (ref_y >> 1) * src_stride_c];

    uint8_t emulate_edge = test_for_edge_emulation_c(ref_x >> 1, ref_y >> 1, pic_w >> 1, pic_h >> 1,
                                                     pu_w, pu_h);;

    if (emulate_edge) {
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cb - src_off), RCN_CTB_STRIDE, src_stride_c,
                             pu_w  + EPEL_EXTRA, pu_h + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cb = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    mc_c->bidir0[prec_c_mc_type][log2_pu_w - 2]((int16_t*) dst_cb,
                                                src_cb, src_stride_c,
                                                pu_h, prec_x_c, prec_y_c, pu_w);

    src_stride_c = frame0->linesize[1] /sizeof(OVSample);
    if (emulate_edge) {
        int src_off  = REF_PADDING_C * src_stride_c + (REF_PADDING_C);
        int buff_off = REF_PADDING_C * (RCN_CTB_STRIDE) + (REF_PADDING_C);
        emulate_block_border(tmp_buff, (src_cr - src_off), RCN_CTB_STRIDE, src_stride_c,
                             pu_w + EPEL_EXTRA, pu_h + EPEL_EXTRA,
                             (pos_x >> 1) + (mv.x >> 5) - REF_PADDING_C, (pos_y >> 1) + (mv.y >> 5) - REF_PADDING_C,
                             (pic_w >> 1), (pic_h >> 1));
        src_cr = tmp_buff + buff_off;
        src_stride_c = RCN_CTB_STRIDE;
    }

    mc_c->bidir0[prec_c_mc_type][log2_pu_w - 2]((int16_t*) dst_cr, src_cr, src_stride_c,
                                                pu_h, prec_x_c, prec_y_c, pu_w);
}

static uint8_t
compute_rpr_filter_idx(int scale_factor, int flag_4x4)
{
    const int rpr_thres_1 = (5 << RPR_SCALE_BITS) >> 2;
    const int rpr_thres_2 = (7 << RPR_SCALE_BITS) >> 2;

    //TODOrpr: use MC filtre 4x4 for affine (and use RPR affine filters)
    int filter_idx = flag_4x4 ? 3 : 0 ;

    if (scale_factor > rpr_thres_2) {
        filter_idx += 2;
    } else if (scale_factor > rpr_thres_1) {
        filter_idx += 1;
    }

    return filter_idx;
}


static void
clip_rpr_position(int* pos_x, int* pos_y, int pic_w, int pic_h, int pb_w, int pb_h, int shift_pos)
{
    int prec_x = *pos_x & ((1 << shift_pos) - 1);
    int prec_y = *pos_y & ((1 << shift_pos) - 1);

    int x_max  = (pic_w + 3) << shift_pos;
    int y_max  = (pic_h + 3) << shift_pos;
    int x_min  = -((pb_w + 4) << shift_pos);
    int y_min  = -((pb_h + 4) << shift_pos);

    *pos_x = ov_clip(*pos_x, x_min + prec_x, x_max + prec_x);
    *pos_y = ov_clip(*pos_y, y_min + prec_y, y_max + prec_y);
}


static void
mcp_rpr_l(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
          OVMV mv, uint8_t rpl_idx, uint8_t ref_idx, int scaling_hor, int scaling_ver)
{
    struct OVRCNCtx *const rcn_ctx   = &ctudec->rcn_ctx;
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;
    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    OVSample* tmp_emul = (OVSample *) rcn_ctx->data.tmp_buff;
    uint16_t* tmp_rpr  = (uint16_t *) rcn_ctx->data.tmp_rpr;

    uint16_t tmp_emul_str = 2 * RCN_CTB_STRIDE;
    uint16_t tmp_rpr_str  = 2 * RCN_CTB_STRIDE;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];
    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;

    const OVFrame *const frame0 = ref_pic->frame;
    const OVSample *const ref0_y  = (OVSample *)frame0->data[0];
    dst.y  += x0 + y0 * dst.stride;
    int src_stride   = frame0->linesize[0] / sizeof(OVSample);

    const int ref_pic_w = frame0->width;
    const int ref_pic_h = frame0->height;

    //MV precision is 4 bits for luma
    int shift_mv  = 4;
    int shift_pos = RPR_SCALE_BITS + shift_mv;
    int offset    = 1 << (RPR_SCALE_BITS - 1);
    uint8_t flag_4x4     = (log2_pu_w == 2 && log2_pu_h == 2);
    uint8_t filter_idx_h = compute_rpr_filter_idx(scaling_hor, flag_4x4);
    uint8_t filter_idx_v = compute_rpr_filter_idx(scaling_ver, flag_4x4);
    int stepX = ((scaling_hor + 8) >> 4) << 4;
    int stepY = ((scaling_ver + 8) >> 4) << 4;

    int32_t ref_pos_x = (((pos_x << shift_mv) + mv.x) * (int32_t)scaling_hor) + (1 << 7);
    int32_t ref_pos_y = (((pos_y << shift_mv) + mv.y) * (int32_t)scaling_ver) + (1 << 7);

    int     ref_x     = (ref_pos_x + offset) >> shift_pos;
    int     ref_y     = (ref_pos_y + offset) >> shift_pos;

    int ref_pu_w = ((ref_pos_x + (((pu_w - 1) * stepX) << shift_mv) + offset) >> shift_pos) - ref_x + 1;
    int ref_pu_h = ((ref_pos_y + (((pu_h - 1) * stepY) << shift_mv) + offset) >> shift_pos) - ref_y + 1;
    ref_pu_h = OVMAX(1, ref_pu_h);

    //Clip ref position now that ref_pu_w and ref_pu_h are computed
    clip_rpr_position(&ref_pos_x, &ref_pos_y, ref_pic_w, ref_pic_h, ref_pu_w, ref_pu_h, shift_pos);
    ref_x = (ref_pos_x + offset) >> shift_pos;
    ref_y = (ref_pos_y + offset) >> shift_pos;

    /*
     * Thread synchronization to ensure data is available before usage
     */
    rcn_inter_synchronization(ref_pic, ref_x, ref_y, ref_pu_w, ref_pu_h, log2_ctb_s);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, ref_pic_w, ref_pic_h,
                                                   ref_pu_w + 1, ref_pu_h + 1);;

    const OVSample *src  = &ref0_y [ ref_x + ref_y * src_stride];
    int buff_off = REF_PADDING_L * (tmp_emul_str) + (REF_PADDING_L);
    int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
    if (emulate_edge) {
        emulate_block_border(tmp_emul, (src - src_off),
                             tmp_emul_str, src_stride,
                             ref_pu_w + QPEL_EXTRA, ref_pu_h + QPEL_EXTRA,
                             ref_x - REF_PADDING_L, ref_y - REF_PADDING_L,
                             ref_pic_w, ref_pic_h);

        src = tmp_emul + REF_PADDING_L;
        src_stride = tmp_emul_str;
    } else {
        src = src - src_off + REF_PADDING_L;
    }

    int32_t   pos_mv_x, pos_mv_y;
    const OVSample* p_src = src ;
    uint16_t* p_tmp_rpr = tmp_rpr + REF_PADDING_L;
    uint8_t   prec_x, prec_y, prec_type;
    for(int col = 0; col < pu_w; col++)
    {
        pos_mv_x  = (ref_pos_x + ((col * stepX) << shift_mv) + offset) >>  RPR_SCALE_BITS ;
        prec_x    = pos_mv_x & 0xF;
        prec_y    = 0;
        prec_type = filter_idx_h || prec_x;
        p_src     = src + (pos_mv_x >> shift_mv) - ref_x;

        mc_l->rpr_h[prec_type][1](p_tmp_rpr, tmp_rpr_str, p_src, src_stride, ref_pu_h + QPEL_EXTRA + 1,
                                  prec_x, prec_y, 1, filter_idx_h);
        p_tmp_rpr   +=  1;
    }

    OVSample* p_dst     = dst.y;
    for(int row = 0; row < pu_h; row++ ) {
        pos_mv_y  = ( ref_pos_y + ((row * stepY) << shift_mv) + offset ) >> RPR_SCALE_BITS;
        prec_x    = 0;
        prec_y    = pos_mv_y & 0xF;
        prec_type = filter_idx_v || prec_y ;
        p_tmp_rpr = tmp_rpr + buff_off + ((pos_mv_y >> shift_mv) - ref_y) * tmp_rpr_str;

        mc_l->rpr_v_uni[prec_type][log2_pu_w-1](p_dst, dst.stride, p_tmp_rpr, tmp_rpr_str, 1,
                                                prec_x, prec_y, pu_w, filter_idx_v);
        p_dst    +=  dst.stride;
    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride,
                                           ctudec->lmcs_info.luts,
                                           pu_w, pu_h);
}


static void
rcn_mcp_rpr_bi_l(OVCTUDec *const ctudec, uint16_t* dst, uint16_t dst_stride, int x0, int y0,
                 int log2_pu_w, int log2_pu_h, OVMV mv, uint8_t rpl_idx, uint8_t ref_idx,
                 int scaling_hor, int scaling_ver)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;
    const int pu_w = 1 << log2_pu_w;
    const int pu_h = 1 << log2_pu_h;

    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    OVSample* tmp_emul = (OVSample *) rcn_ctx->data.tmp_buff;
    uint16_t* tmp_rpr  = (uint16_t *) rcn_ctx->data.tmp_rpr;
    uint16_t tmp_emul_str = 2*RCN_CTB_STRIDE;
    uint16_t tmp_rpr_str  = 2*RCN_CTB_STRIDE;

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];
    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;

    const OVFrame *const frame0 = ref_pic->frame;
    const OVSample *const ref0_y  = (OVSample *) frame0->data[0];
    dst  += x0 + y0 * dst_stride;
    int src_stride   = frame0->linesize[0] /sizeof(OVSample);

    const int ref_pic_w = frame0->width;
    const int ref_pic_h = frame0->height;

    //MV precision is 4 bits for luma
    int shift_mv  = 4;
    int offset    = 1 << (RPR_SCALE_BITS - 1);
    int shift_pos  = RPR_SCALE_BITS + shift_mv;
    uint8_t flag_4x4     = (log2_pu_w == 2 && log2_pu_h == 2);
    uint8_t filter_idx_h = compute_rpr_filter_idx(scaling_hor, flag_4x4);
    uint8_t filter_idx_v = compute_rpr_filter_idx(scaling_ver, flag_4x4);
    int stepX = (( scaling_hor + 8 ) >> 4) << 4;
    int stepY = (( scaling_ver + 8 ) >> 4) << 4;

    int32_t ref_pos_x = ((( pos_x << shift_mv)  + mv.x ) * (int32_t)scaling_hor) + (1<<7);
    int     ref_x     = (ref_pos_x + offset)  >> shift_pos;
    int32_t ref_pos_y = ((( pos_y << shift_mv ) + mv.y ) * (int32_t)scaling_ver) + (1<<7);
    int     ref_y     = (ref_pos_y + offset) >> shift_pos;
    int ref_pu_w = ((ref_pos_x + (((pu_w-1) * stepX) << shift_mv) + offset) >> shift_pos) - ref_x + 1 ;
    int ref_pu_h = ((ref_pos_y + (((pu_h-1) * stepY) << shift_mv) + offset) >> shift_pos) - ref_y + 1;
    ref_pu_h = OVMAX(1, ref_pu_h);

    //Clip ref position now that ref_pu_w and ref_pu_h are computed
    clip_rpr_position(&ref_pos_x, &ref_pos_y, ref_pic_w, ref_pic_h, ref_pu_w, ref_pu_h, shift_pos);
    ref_x = (ref_pos_x + offset) >> shift_pos;
    ref_y = (ref_pos_y + offset) >> shift_pos;

    /*
     * Thread synchronization to ensure data is available before usage
     */
    rcn_inter_synchronization(ref_pic, ref_x, ref_y, ref_pu_w, ref_pu_h, log2_ctb_s);

    uint8_t emulate_edge = test_for_edge_emulation(ref_x, ref_y, ref_pic_w, ref_pic_h,
                                                   ref_pu_w + 1, ref_pu_h + 1);;

    const OVSample *src  = &ref0_y [ ref_x + ref_y * src_stride];
    int buff_off = REF_PADDING_L * (tmp_emul_str) + (REF_PADDING_L);
    int src_off  = REF_PADDING_L * (src_stride) + (REF_PADDING_L);
    if (emulate_edge) {
        emulate_block_border(tmp_emul, (src - src_off),
                             tmp_emul_str, src_stride,
                             ref_pu_w + QPEL_EXTRA, ref_pu_h + QPEL_EXTRA,
                             ref_x - REF_PADDING_L, ref_y - REF_PADDING_L,
                             ref_pic_w, ref_pic_h);

        src = tmp_emul + REF_PADDING_L;
        src_stride = tmp_emul_str;
    } else {
        src = src - src_off + REF_PADDING_L ;
    }

    int32_t   pos_mv_x, pos_mv_y;
    const OVSample* p_src = src ;
    uint16_t* p_tmp_rpr = tmp_rpr + REF_PADDING_L;
    uint8_t   prec_x, prec_y, prec_type;
    for(int col = 0; col < pu_w; col++)
    {
        pos_mv_x  = (ref_pos_x + ((col * stepX) << shift_mv) + offset) >>  RPR_SCALE_BITS ;
        prec_x    = pos_mv_x & 0xF;
        prec_y    = 0;
        prec_type = filter_idx_h || prec_x;
        p_src     = src + (pos_mv_x >> shift_mv) - ref_x;

        mc_l->rpr_h[prec_type][1](p_tmp_rpr, tmp_rpr_str, p_src, src_stride, ref_pu_h + QPEL_EXTRA + 1,
                                  prec_x, prec_y, 1, filter_idx_h);
        p_tmp_rpr   +=  1;
    }

    uint16_t* p_dst     = dst;
    for(int row = 0; row < pu_h; row++ ) {
        pos_mv_y  = ( ref_pos_y + ((row * stepY) << shift_mv) + offset ) >> RPR_SCALE_BITS;
        prec_x    = 0;
        prec_y    = pos_mv_y & 0xF;
        prec_type = filter_idx_v || prec_y;
        p_tmp_rpr = tmp_rpr + buff_off + ((pos_mv_y >> shift_mv) - ref_y) * tmp_rpr_str;
        mc_l->rpr_v_bi[prec_type][log2_pu_w-1](p_dst, dst_stride, p_tmp_rpr, tmp_rpr_str, 1,
                                               prec_x, prec_y, pu_w, filter_idx_v);
        p_dst    +=  dst_stride;
    }

}

static void
mcp_rpr_c(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
          OVMV mv, uint8_t rpl_idx, uint8_t ref_idx, int scaling_hor, int scaling_ver)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;
    int pu_w = 1 << log2_pu_w;
    int pu_h = 1 << log2_pu_h;

    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    OVSample* tmp_emul = (OVSample *) rcn_ctx->data.tmp_buff;
    uint16_t* tmp_rpr  = (uint16_t *) rcn_ctx->data.tmp_rpr;
    uint16_t tmp_emul_str = 2*RCN_CTB_STRIDE;
    uint16_t tmp_rpr_str  = 2*RCN_CTB_STRIDE;

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];
    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    if (!ref_pic) return;
    const OVFrame *const frame0 = ref_pic->frame;
    const OVSample *const ref0_cb = (OVSample *) frame0->data[1];
    const OVSample *const ref0_cr = (OVSample *) frame0->data[2];
    int src_stride_c = frame0->linesize[1] /sizeof(OVSample);

    int ref_pic_w = frame0->width >> 1;
    int ref_pic_h = frame0->height >> 1;
    pos_x = pos_x >> 1;
    pos_y = pos_y >> 1;
    pu_w = pu_w >> 1;
    pu_h = pu_h >> 1;

    //MV precision in 5 bits for chroma
    int shift_mv  = 5;
    int offset    = 1 << (RPR_SCALE_BITS - 1);
    int shift_pos  = RPR_SCALE_BITS + shift_mv;
    uint8_t flag_4x4     = (log2_pu_w == 2 && log2_pu_h == 2);
    uint8_t filter_idx_h = compute_rpr_filter_idx(scaling_hor, flag_4x4);
    uint8_t filter_idx_v = compute_rpr_filter_idx(scaling_ver, flag_4x4);
    int stepX = (( scaling_hor + 8 ) >> 4) << 4;
    int stepY = (( scaling_ver + 8 ) >> 4) << 4;

    int32_t add_x = (1 - ref0->scale_info.chroma_hor_col_flag) * 8 * ( scaling_hor - (1<<RPR_SCALE_BITS));
    int32_t add_y = (1 - ref0->scale_info.chroma_ver_col_flag) * 8 * ( scaling_ver - (1<<RPR_SCALE_BITS));
    int32_t ref_pos_x = (( pos_x << shift_mv)  + mv.x ) * (int32_t)scaling_hor + add_x + (1 << 8);
    int     ref_x     = (ref_pos_x + offset)  >> shift_pos;
    int32_t ref_pos_y = (( pos_y << shift_mv ) + mv.y ) * (int32_t)scaling_ver + add_y + (1 << 8);
    int     ref_y     = (ref_pos_y + offset) >> shift_pos;
    int ref_pu_w = ((ref_pos_x + (((pu_w-1) * stepX) << shift_mv) + offset) >> shift_pos) - ref_x + 1 ;
    int ref_pu_h = ((ref_pos_y + (((pu_h-1) * stepY) << shift_mv) + offset) >> shift_pos) - ref_y + 1;
    ref_pu_h = OVMAX(1, ref_pu_h);

    //Clip ref position now that ref_pu_w and ref_pu_h are computed
    clip_rpr_position(&ref_pos_x, &ref_pos_y, ref_pic_w, ref_pic_h, ref_pu_w, ref_pu_h, shift_pos);
    ref_x = (ref_pos_x + offset)  >> shift_pos;
    ref_y = (ref_pos_y + offset) >> shift_pos;

    const OVSample *src_cb = &ref0_cb[ref_x + ref_y * src_stride_c];
    const OVSample *src_cr = &ref0_cr[ref_x + ref_y * src_stride_c];

    uint8_t prec_x, prec_y, prec_type;
    const OVSample* p_src;
    OVSample* p_dst;
    for(int comp_c = 0;  comp_c < 2; comp_c ++)
    {
        src_stride_c = frame0->linesize[1] /sizeof(OVSample);
        const OVSample * src_c = src_cb;
        OVSample* dst_c = dst.cb;
        if (comp_c == 1) {
            src_c = src_cr;
            dst_c = dst.cr;
        }

        uint8_t emulate_edge = test_for_edge_emulation_c(ref_x, ref_y, ref_pic_w, ref_pic_h,
                                                         ref_pu_w + 2, ref_pu_h + 2);;

        int buff_off = REF_PADDING_C * (tmp_emul_str) + (REF_PADDING_C);
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        if (emulate_edge) {
            emulate_block_border(tmp_emul, (src_c - src_off),
                                 tmp_emul_str, src_stride_c,
                                 (ref_pu_w)  + EPEL_EXTRA, (ref_pu_h) + EPEL_EXTRA,
                                 (ref_x) - REF_PADDING_C, (ref_y) - REF_PADDING_C,
                                 (ref_pic_w), (ref_pic_h));
            src_c = tmp_emul + REF_PADDING_C;
            src_stride_c = tmp_emul_str ;
        } else {
            src_c = src_c - src_off + REF_PADDING_C ;
        }

        int pos_mv_x, pos_mv_y;
        uint16_t* p_tmp_rpr = tmp_rpr + REF_PADDING_C;
        for(int col = 0; col < pu_w; col++)
        {
            pos_mv_x  = ( ref_pos_x + ((col * stepX) << shift_mv) + offset ) >> RPR_SCALE_BITS;
            prec_x    = pos_mv_x & 0x1F;
            prec_y    = 0;
            prec_type = filter_idx_h || prec_x;
            p_src     = src_c + (pos_mv_x >> shift_mv) - ref_x;

            mc_c->rpr_h[prec_type][1](p_tmp_rpr, tmp_rpr_str, p_src, src_stride_c,
                                      ref_pu_h + EPEL_EXTRA + 2, prec_x, prec_y, 1, filter_idx_h);
            p_tmp_rpr   +=  1;
        }

        p_dst        = dst_c;
        for(int row = 0; row < pu_h; row++ )
        {
            pos_mv_y  = ( ref_pos_y + ((row * stepY) << shift_mv) + offset ) >> RPR_SCALE_BITS;
            prec_x    = 0;
            prec_y    = pos_mv_y & 0x1F;
            prec_type = filter_idx_v || prec_y;
            p_tmp_rpr = tmp_rpr + buff_off + ((pos_mv_y >> shift_mv) - ref_y) * tmp_rpr_str;

            mc_c->rpr_v_uni[prec_type][log2_pu_w-2](p_dst, dst.stride_c, p_tmp_rpr, tmp_rpr_str,
                                                    1, prec_x, prec_y, pu_w, filter_idx_v);
            p_dst    +=  dst.stride_c;
        }
    }
}

static void
rcn_mcp_rpr_bi_c(OVCTUDec *const ctudec, uint16_t* dst_cb, uint16_t* dst_cr, uint16_t dst_stride, int x0, int y0,
                 int log2_pu_w, int log2_pu_h, OVMV mv, uint8_t rpl_idx, uint8_t ref_idx,
                 int scaling_hor, int scaling_ver)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int pos_x = (ctudec->ctb_x << log2_ctb_s) + x0;
    int pos_y = (ctudec->ctb_y << log2_ctb_s) + y0;
    int pu_w = 1 << log2_pu_w;
    int pu_h = 1 << log2_pu_h;

    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    OVSample* tmp_emul = (OVSample *) rcn_ctx->data.tmp_buff;
    uint16_t* tmp_rpr  = (uint16_t *) rcn_ctx->data.tmp_rpr;
    uint16_t tmp_emul_str = 2*RCN_CTB_STRIDE;
    uint16_t tmp_rpr_str  = 2*RCN_CTB_STRIDE;

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;

    OVPicture *ref0 = inter_ctx->rpl0[ref_idx];
    OVPicture *ref1 = inter_ctx->rpl1[ref_idx];
    OVPicture *ref_pic =  rpl_idx ? ref1 : ref0;
    if (!ref_pic) return;

    dst_cb += (x0 >> 1) + (y0 >> 1) * dst_stride;
    dst_cr += (x0 >> 1) + (y0 >> 1) * dst_stride;

    const OVFrame *const frame0 = ref_pic->frame;
    const OVSample *const ref0_cb = (OVSample *) frame0->data[1];
    const OVSample *const ref0_cr = (OVSample *) frame0->data[2];
    int src_stride_c = frame0->linesize[1] /sizeof(OVSample);

    int ref_pic_w = frame0->width >> 1;
    int ref_pic_h = frame0->height >> 1;
    pos_x = pos_x >> 1;
    pos_y = pos_y >> 1;
    pu_w = pu_w >> 1;
    pu_h = pu_h >> 1;

    //MV precision in 5 bits for chroma
    int shift_mv  = 5;
    int offset    = 1 << (RPR_SCALE_BITS - 1);
    int shift_pos  = RPR_SCALE_BITS + shift_mv;
    uint8_t flag_4x4     = (log2_pu_w == 2 && log2_pu_h == 2);
    uint8_t filter_idx_h = compute_rpr_filter_idx(scaling_hor, flag_4x4);
    uint8_t filter_idx_v = compute_rpr_filter_idx(scaling_ver, flag_4x4);
    int stepX = (( scaling_hor + 8 ) >> 4) << 4;
    int stepY = (( scaling_ver + 8 ) >> 4) << 4;

    int32_t add_x = (1 - ref0->scale_info.chroma_hor_col_flag) * 8 * ( scaling_hor - (1<<RPR_SCALE_BITS));
    int32_t add_y = (1 - ref0->scale_info.chroma_ver_col_flag) * 8 * ( scaling_ver - (1<<RPR_SCALE_BITS));
    int32_t ref_pos_x = (( pos_x << shift_mv)  + mv.x ) * (int32_t)scaling_hor + add_x + (1 << 8);
    int     ref_x     = (ref_pos_x + offset)  >> shift_pos;
    int32_t ref_pos_y = (( pos_y << shift_mv ) + mv.y ) * (int32_t)scaling_ver + add_y + (1 << 8);
    int     ref_y     = (ref_pos_y + offset) >> shift_pos;
    int ref_pu_w = ((ref_pos_x + (((pu_w-1) * stepX) << shift_mv) + offset) >> shift_pos) - ref_x + 1 ;
    int ref_pu_h = ((ref_pos_y + (((pu_h-1) * stepY) << shift_mv) + offset) >> shift_pos) - ref_y + 1;
    ref_pu_h = OVMAX(1, ref_pu_h);

    //Clip ref position now that ref_pu_w and ref_pu_h are computed
    clip_rpr_position(&ref_pos_x, &ref_pos_y, ref_pic_w, ref_pic_h, ref_pu_w, ref_pu_h, shift_pos);
    ref_x = (ref_pos_x + offset)  >> shift_pos;
    ref_y = (ref_pos_y + offset) >> shift_pos;

    const OVSample *src_cb = &ref0_cb[ref_x + ref_y * src_stride_c];
    const OVSample *src_cr = &ref0_cr[ref_x + ref_y * src_stride_c];

    uint8_t prec_x, prec_y, prec_type;
    const OVSample* p_src;
    uint16_t* p_dst;
    for(int comp_c = 0;  comp_c < 2; comp_c ++)
    {
        src_stride_c = frame0->linesize[1]/sizeof(OVSample);
        const OVSample * src_c = src_cb;
        uint16_t * dst_c = dst_cb;
        if (comp_c == 1) {
            src_c = src_cr;
            dst_c = dst_cr;
        }

        uint8_t emulate_edge = test_for_edge_emulation_c(ref_x, ref_y, ref_pic_w, ref_pic_h,
                                                         ref_pu_w + 2, ref_pu_h + 2);;

        int buff_off = REF_PADDING_C * (tmp_emul_str) + (REF_PADDING_C);
        int src_off  = REF_PADDING_C * (src_stride_c) + (REF_PADDING_C);
        if (emulate_edge) {
            emulate_block_border(tmp_emul, (src_c - src_off),
                                 tmp_emul_str, src_stride_c,
                                 (ref_pu_w)  + EPEL_EXTRA, (ref_pu_h) + EPEL_EXTRA,
                                 (ref_x) - REF_PADDING_C, (ref_y) - REF_PADDING_C,
                                 (ref_pic_w), (ref_pic_h));
            src_c = tmp_emul + REF_PADDING_C;
            src_stride_c = tmp_emul_str ;
        } else {
            src_c = src_c - src_off + REF_PADDING_C ;
        }

        int pos_mv_x, pos_mv_y;
        uint16_t* p_tmp_rpr = tmp_rpr + REF_PADDING_C;
        for(int col = 0; col < pu_w; col++)
        {
            pos_mv_x  = ( ref_pos_x + ((col * stepX) << shift_mv) + offset ) >> RPR_SCALE_BITS;
            prec_x    = pos_mv_x & 0x1F;
            prec_y    = 0;
            prec_type = filter_idx_h || prec_x;
            p_src     = src_c + (pos_mv_x >> shift_mv) - ref_x;

            mc_c->rpr_h[prec_type][1](p_tmp_rpr, tmp_rpr_str, p_src, src_stride_c,
                                      ref_pu_h + EPEL_EXTRA + 2, prec_x, prec_y, 1, filter_idx_h);
            p_tmp_rpr   +=  1;
        }

        p_dst        = dst_c;
        for(int row = 0; row < pu_h; row++ )
        {
            pos_mv_y  = ( ref_pos_y + ((row * stepY) << shift_mv) + offset ) >> RPR_SCALE_BITS;
            prec_x    = 0;
            prec_y    = pos_mv_y & 0x1F;
            prec_type = filter_idx_v || prec_y ;
            p_tmp_rpr = tmp_rpr + buff_off + ((pos_mv_y >> shift_mv) - ref_y) * tmp_rpr_str;

            mc_c->rpr_v_bi[prec_type][log2_pu_w-2](p_dst, dst_stride, p_tmp_rpr, tmp_rpr_str,
                                                   1, prec_x, prec_y, pu_w, filter_idx_v);
            p_dst    +=  dst_stride;
        }
    }
}


static void
mcp_rpr(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
        OVMV mv, uint8_t rpl_idx, uint8_t ref_idx, int scaling_hor, int scaling_ver)
{
    mcp_rpr_l(ctudec, dst, x0, y0, log2_pu_w, log2_pu_h, mv, rpl_idx, ref_idx, scaling_hor, scaling_ver);
    mcp_rpr_c(ctudec, dst, x0, y0, log2_pu_w, log2_pu_h, mv, rpl_idx, ref_idx, scaling_hor, scaling_ver);
}


static void
mc_rpr_b_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
           uint8_t x0, uint8_t y0, uint8_t log2_pb_w, uint8_t log2_pb_h,
           OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1,
           int rpr_scale_h0, int rpr_scale_v0,
           int rpr_scale_h1, int rpr_scale_v1, struct VVCGPM* gpm_ctx )
{
    uint8_t no_rpr_l0 = rpr_scale_h0 == RPR_NO_SCALE && rpr_scale_v0 == RPR_NO_SCALE;
    uint8_t no_rpr_l1 = rpr_scale_h1 == RPR_NO_SCALE && rpr_scale_v1 == RPR_NO_SCALE;
    int type0 = 0;
    int type1 = 1;
    uint8_t use_bcw = mv0.bcw_idx_plus1 != 0 && mv0.bcw_idx_plus1 != 3;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    if (gpm_ctx) {
        type0 = gpm_ctx->inter_dir0 - 1;
        type1 = gpm_ctx->inter_dir1 - 1;
    }

    struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    uint16_t* tmp_rpl0 = rcn_ctx->data.tmp_bi_mrg0;
    int tmp_rpl0_stride = RCN_CTB_STRIDE;
    if (no_rpr_l0) {
        tmp_rpl0_stride = MAX_PB_SIZE;
        rcn_mcp_bidir0_l(ctudec, tmp_rpl0, tmp_rpl0_stride, x0, y0, log2_pb_w, log2_pb_h, mv0, type0, ref_idx0);
    } else {
        rcn_mcp_rpr_bi_l(ctudec, tmp_rpl0, tmp_rpl0_stride, x0, y0, log2_pb_w, log2_pb_h, mv0, type0, ref_idx0,
                         rpr_scale_h0, rpr_scale_v0);
    }

    uint16_t* tmp_rpl1 = rcn_ctx->data.tmp_bi_mrg1;
    int tmp_rpl1_stride = RCN_CTB_STRIDE;
    if (no_rpr_l1) {
        tmp_rpl1_stride = MAX_PB_SIZE;
        rcn_mcp_bidir0_l(ctudec, tmp_rpl1, tmp_rpl1_stride, x0, y0, log2_pb_w, log2_pb_h, mv1, type1, ref_idx1);
    } else {
        rcn_mcp_rpr_bi_l(ctudec, tmp_rpl1, tmp_rpl1_stride, x0, y0, log2_pb_w, log2_pb_h, mv1, type1, ref_idx1,
                         rpr_scale_h1, rpr_scale_v1);
    }

    dst.y       += x0 + y0 * dst.stride;
    tmp_rpl0    += x0 + y0 * tmp_rpl0_stride;
    tmp_rpl1    += x0 + y0 * tmp_rpl1_stride;
    int pu_w = 1 <<log2_pb_w;
    int pu_h = 1 <<log2_pb_h;

    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    if (!use_bcw) {
        if (gpm_ctx) {
            int16_t* weight;
            int step_x, step_y;
            uint8_t chroma_flag = 0;
            gpm_weights_and_steps(gpm_ctx->split_dir, log2_pb_w, log2_pb_h, &step_x, &step_y, &weight, chroma_flag);
            mc_l->gpm_weighted(dst.y, dst.stride, (int16_t*) tmp_rpl0, tmp_rpl0_stride,
                               (int16_t*) tmp_rpl1, tmp_rpl1_stride, pu_h, pu_w, step_x, step_y, weight);
        } else {
            mc_l->rpr_sum(dst.y, dst.stride, tmp_rpl0, tmp_rpl0_stride,
                          tmp_rpl1, tmp_rpl1_stride, pu_h, 0, 0, pu_w);
        }
    } else {
        int wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        int wt0 = 8 - wt1;
        mc_l->rpr_w[log2_pb_w-1](dst.y, dst.stride, tmp_rpl0, tmp_rpl0_stride,
                                 tmp_rpl1, tmp_rpl1_stride, pu_h, BCW_DENOM, wt0, wt1, 0, 0, pu_w);
    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride, ctudec->lmcs_info.luts, pu_w, pu_h);
}

static void
mc_rpr_prof_b_l(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                uint8_t x0, uint8_t y0, uint8_t log2_pb_w, uint8_t log2_pb_h,
                OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1,
                int rpr_scale_h0, int rpr_scale_v0,
                int rpr_scale_h1, int rpr_scale_v1,
                uint8_t prof_dir, const struct PROFInfo *const prof_info )
{
    uint8_t no_rpr_l0 = rpr_scale_h0 == RPR_NO_SCALE && rpr_scale_v0 == RPR_NO_SCALE;
    uint8_t no_rpr_l1 = rpr_scale_h1 == RPR_NO_SCALE && rpr_scale_v1 == RPR_NO_SCALE;
    uint8_t use_bcw = mv0.bcw_idx_plus1 != 0 && mv0.bcw_idx_plus1 != 3;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;

    struct OVRCNCtx *const rcn_ctx = &ctudec->rcn_ctx;
    uint16_t *tmp_rpl0 = &rcn_ctx->data.tmp_bi_mrg0[0];
    uint16_t tmp_rpl0_stride = RCN_CTB_STRIDE;
    if (no_rpr_l0) {
        rcn_prof_mcp_bi_l(ctudec, tmp_rpl0, tmp_rpl0_stride, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0,
                          prof_info->dmv_scale_h_0, prof_info->dmv_scale_v_0);
    } else {
        rcn_mcp_rpr_bi_l(ctudec, tmp_rpl0, tmp_rpl0_stride, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0,
                         rpr_scale_h0, rpr_scale_v0);
    }

    uint16_t *tmp_rpl1 = &rcn_ctx->data.tmp_bi_mrg1[0];
    uint16_t tmp_rpl1_stride = RCN_CTB_STRIDE;
    if (no_rpr_l1) {
        rcn_prof_mcp_bi_l(ctudec, tmp_rpl1, tmp_rpl1_stride, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1,
                          prof_info->dmv_scale_h_1, prof_info->dmv_scale_v_1);
    } else {
        rcn_mcp_rpr_bi_l(ctudec, tmp_rpl1, tmp_rpl1_stride, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1,
                         rpr_scale_h1, rpr_scale_v1);
    }

    dst.y     += x0 + y0 * dst.stride;
    tmp_rpl0  += x0 + y0 * tmp_rpl0_stride;
    tmp_rpl1  += x0 + y0 * tmp_rpl1_stride;
    int pu_w = 1 <<log2_pb_w;
    int pu_h = 1 <<log2_pb_h;

    struct MCFunctions *mc_l = &ctudec->rcn_funcs.mc_l;
    if (!use_bcw) {
        mc_l->rpr_sum(dst.y, dst.stride, tmp_rpl0, tmp_rpl0_stride,
                      tmp_rpl1, tmp_rpl1_stride, pu_h, 0, 0, pu_w);
    } else {
        int wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        int wt0 = 8 - wt1;
        mc_l->rpr_w[log2_pb_w-1](dst.y, dst.stride, tmp_rpl0, tmp_rpl0_stride,
                                 tmp_rpl1, tmp_rpl1_stride, pu_h, BCW_DENOM, wt0, wt1, 0, 0, pu_w);
    }

    ctudec->rcn_funcs.lmcs_reshape_forward(dst.y, dst.stride, ctudec->lmcs_info.luts, pu_w, pu_h);
}


static void
mc_rpr_b_c(OVCTUDec *const ctudec, struct OVBuffInfo dst,
           uint8_t x0, uint8_t y0, uint8_t log2_pb_w, uint8_t log2_pb_h,
           OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1,
           int rpr_scale_h0, int rpr_scale_v0,
           int rpr_scale_h1, int rpr_scale_v1, struct VVCGPM* gpm_ctx )
{
    uint8_t no_rpr_l0 = rpr_scale_h0 == RPR_NO_SCALE && rpr_scale_v0 == RPR_NO_SCALE;
    uint8_t no_rpr_l1 = rpr_scale_h1 == RPR_NO_SCALE && rpr_scale_v1 == RPR_NO_SCALE;
    uint8_t use_bcw = mv0.bcw_idx_plus1 != 0 && mv0.bcw_idx_plus1 != 3;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    int type0 = 0;
    int type1 = 1;
    if (gpm_ctx) {
        type0 = gpm_ctx->inter_dir0 - 1;
        type1 = gpm_ctx->inter_dir1 - 1;
    }

    struct OVRCNCtx *const rcn_ctx   = &ctudec->rcn_ctx;
    uint16_t* tmp_rpl0_cb = rcn_ctx->data.tmp_bi_mrg0;
    uint16_t* tmp_rpl0_cr = rcn_ctx->data.tmp_bi_mrg1;
    uint16_t tmp_rpl0_stride = dst.stride_c;
    if (no_rpr_l0) {
        tmp_rpl0_stride = MAX_PB_SIZE;
        rcn_mcp_bidir0_c(ctudec, tmp_rpl0_cb, tmp_rpl0_cr, tmp_rpl0_stride, x0, y0,
                         log2_pb_w, log2_pb_h, mv0, type0, ref_idx0);
    } else {
        rcn_mcp_rpr_bi_c(ctudec, tmp_rpl0_cb, tmp_rpl0_cr, tmp_rpl0_stride, x0, y0,
                         log2_pb_w, log2_pb_h, mv0, type0, ref_idx0, rpr_scale_h0, rpr_scale_v0);
    }

    uint16_t* tmp_rpl1_cb = rcn_ctx->data.tmp_bi_mrg2;
    uint16_t* tmp_rpl1_cr = rcn_ctx->data.tmp_buff0;
    uint16_t tmp_rpl1_stride = dst.stride_c;
    if (no_rpr_l1) {
        tmp_rpl1_stride = MAX_PB_SIZE;
        rcn_mcp_bidir0_c(ctudec, tmp_rpl1_cb, tmp_rpl1_cr, tmp_rpl1_stride, x0, y0,
                         log2_pb_w, log2_pb_h, mv1, type1, ref_idx1);
    } else {
        rcn_mcp_rpr_bi_c(ctudec, tmp_rpl1_cb, tmp_rpl1_cr, tmp_rpl1_stride, x0, y0,
                         log2_pb_w, log2_pb_h, mv1, type1, ref_idx1, rpr_scale_h1, rpr_scale_v1);
    }

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;

    tmp_rpl0_cb += (x0 >> 1) + (y0 >> 1) * tmp_rpl0_stride;
    tmp_rpl1_cb += (x0 >> 1) + (y0 >> 1) * tmp_rpl1_stride;
    tmp_rpl0_cr += (x0 >> 1) + (y0 >> 1) * tmp_rpl0_stride;
    tmp_rpl1_cr += (x0 >> 1) + (y0 >> 1) * tmp_rpl1_stride;

    int pu_w = 1 << (log2_pb_w - 1);
    int pu_h = 1 << (log2_pb_h - 1);

    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;

    if (!use_bcw) {
        if (gpm_ctx) {
            int16_t* weight;
            int step_x, step_y;
            uint8_t chroma_flag = 1;
            gpm_weights_and_steps(gpm_ctx->split_dir, log2_pb_w, log2_pb_h, &step_x, &step_y, &weight, chroma_flag);
            mc_c->gpm_weighted(dst.cb, dst.stride_c, (int16_t*) tmp_rpl0_cb, tmp_rpl0_stride,
                               (int16_t*) tmp_rpl1_cb, tmp_rpl1_stride, pu_h, pu_w, step_x, step_y, weight);
            mc_c->gpm_weighted(dst.cr, dst.stride_c, (int16_t*) tmp_rpl0_cr, tmp_rpl0_stride,
                               (int16_t*) tmp_rpl1_cr, tmp_rpl1_stride, pu_h, pu_w, step_x, step_y, weight);
        } else {
            mc_c->rpr_sum(dst.cb, dst.stride_c, tmp_rpl0_cb, tmp_rpl0_stride, tmp_rpl1_cb, tmp_rpl1_stride,
                          pu_h, 0, 0, pu_w);
            mc_c->rpr_sum(dst.cr, dst.stride_c, tmp_rpl0_cr, tmp_rpl0_stride, tmp_rpl1_cr, tmp_rpl1_stride,
                          pu_h, 0, 0, pu_w);
        }
    } else {
        int wt1 = bcw_weights[mv0.bcw_idx_plus1-1];
        int wt0 = 8 - wt1;
        mc_c->rpr_w[log2_pb_w-2](dst.cb, dst.stride_c, tmp_rpl0_cb, tmp_rpl0_stride, tmp_rpl1_cb, tmp_rpl1_stride,
                                 pu_h, BCW_DENOM, wt0, wt1, 0, 0, pu_w);
        mc_c->rpr_w[log2_pb_w-2](dst.cr, dst.stride_c, tmp_rpl0_cr, tmp_rpl0_stride, tmp_rpl1_cr, tmp_rpl1_stride,
                                 pu_h, BCW_DENOM, wt0, wt1, 0, 0, pu_w);
    }

}

static void
rcn_mcp(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
        OVMV mv, uint8_t rpl_idx, uint8_t ref_idx)
{
    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx));

    if (no_rpr_l0) {
        mcp_l(ctudec, dst, x0, y0, log2_pu_w, log2_pu_h, mv, rpl_idx, ref_idx);
        mcp_c(ctudec, dst, x0, y0, log2_pu_w, log2_pu_h, mv, rpl_idx, ref_idx);
    } else {
        int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx][0];
        int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx][1];
        mcp_rpr(ctudec, dst, x0, y0, log2_pu_w, log2_pu_h, mv, rpl_idx, ref_idx, rpr_scale_h0, rpr_scale_v0);
    }

}

static void
rcn_mcp_b(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
          const OVPartInfo *const part_ctx,
          const OVMV mv0, const OVMV mv1,
          unsigned int x0, unsigned int y0,
          unsigned int log2_pb_w, unsigned int log2_pb_h,
          uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    uint8_t identical_motion = check_identical_motion(inter_ctx, inter_dir, mv0, mv1, ref_idx0, ref_idx1);

    if (inter_dir == 3 && !identical_motion) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        if (no_rpr_l0 && no_rpr_l1) {
            motion_compensation_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);
            rcn_motion_compensation_b_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];
            mc_rpr_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1,
                       rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, NULL);
            mc_rpr_b_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1,
                       rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, NULL);
        }

    } else if (inter_dir & 0x2 || identical_motion) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        if (no_rpr_l1) {
            mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);
            mcp_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);
        } else {
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];
            mcp_rpr(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1, rpr_scale_h1, rpr_scale_v1);
        }

    } else if (inter_dir & 0x1) {
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        if (no_rpr_l0) {
            mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);
            mcp_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            mcp_rpr(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0, rpr_scale_h0, rpr_scale_v0);
        }
    }
}

static void
rcn_mcp_b_l(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
            const OVPartInfo *const part_ctx,
            const OVMV mv0, const OVMV mv1,
            unsigned int x0, unsigned int y0,
            unsigned int log2_pb_w, unsigned int log2_pb_h,
            uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{

    uint8_t identical_motion = check_identical_motion(inter_ctx, inter_dir, mv0, mv1, ref_idx0, ref_idx1);

    if (inter_dir == 3 && !identical_motion) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        if (no_rpr_l0 && no_rpr_l1) {
            motion_compensation_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];
            mc_rpr_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1,
                            rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, NULL);
        }

    } else if (inter_dir & 0x2 || identical_motion) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        if (no_rpr_l1) {
            mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);
        } else {
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];
            mcp_rpr_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1, rpr_scale_h1, rpr_scale_v1);
        }

    } else if (inter_dir & 0x1) {
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        if (no_rpr_l0) {
            mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            mcp_rpr_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0, rpr_scale_h0, rpr_scale_v0);
        }
    }
}


static void
rcn_prof_mcp_b_l(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                 const OVPartInfo *const part_ctx,
                 const OVMV mv0, const OVMV mv1,
                 unsigned int x0, unsigned int y0,
                 unsigned int log2_pb_w, unsigned int log2_pb_h,
                 uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                 uint8_t prof_dir, const struct PROFInfo *const prof_info)
{
    if (inter_dir == 3) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        uint8_t apply_prof = (no_rpr_l0 && (prof_dir & 0x1)) || (no_rpr_l1 && (prof_dir & 0x2));
        if (no_rpr_l0 && no_rpr_l1) {
            prof_motion_compensation_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1,
                                         ref_idx0, ref_idx1, prof_dir, prof_info);
        } else if (apply_prof) {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];

            mc_rpr_prof_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1,
                            rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, prof_dir, prof_info);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];

            mc_rpr_b_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1,
                       rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, NULL);
        }

    } else if (inter_dir & 0x2) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        if (no_rpr_l1) {
            rcn_prof_mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1,
                           prof_info->dmv_scale_h_1,
                           prof_info->dmv_scale_v_1);
        } else {
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];

            mcp_rpr_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1,
                      rpr_scale_h1, rpr_scale_v1);
        }

    } else if (inter_dir & 0x1) {
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        if (no_rpr_l0) {
            rcn_prof_mcp_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0,
                           prof_info->dmv_scale_h_0,
                           prof_info->dmv_scale_v_0);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];

            mcp_rpr_l(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0,
                      rpr_scale_h0, rpr_scale_v0);
        }
    }
}

static void
rcn_mcp_b_c(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
            const OVPartInfo *const part_ctx,
            const OVMV mv0, const OVMV mv1,
            unsigned int x0, unsigned int y0,
            unsigned int log2_pb_w, unsigned int log2_pb_h,
            uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    uint8_t identical_motion = check_identical_motion(inter_ctx, inter_dir, mv0, mv1, ref_idx0, ref_idx1);
    if (inter_dir == 3 && !identical_motion) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));

        if (no_rpr_l0 && no_rpr_l1) {
            rcn_motion_compensation_b_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];

            mc_rpr_b_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, ref_idx0, ref_idx1,
                       rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, NULL);
        }

    } else if (inter_dir & 0x2 || identical_motion) {
        uint8_t no_rpr_l1 = !(inter_ctx->rpr_scale_msk1 & (1 << ref_idx1));
        if (no_rpr_l1) {
            mcp_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1);
        } else {
            int rpr_scale_h1 = inter_ctx->scale_fact_rpl1[ref_idx1][0];
            int rpr_scale_v1 = inter_ctx->scale_fact_rpl1[ref_idx1][1];
            mcp_rpr_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv1, 1, ref_idx1,
                      rpr_scale_h1, rpr_scale_v1);
        }

    } else if (inter_dir & 0x1) {
        uint8_t no_rpr_l0 = !(inter_ctx->rpr_scale_msk0 & (1 << ref_idx0));
        if (no_rpr_l0) {
            mcp_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0);
        } else {
            int rpr_scale_h0 = inter_ctx->scale_fact_rpl0[ref_idx0][0];
            int rpr_scale_v0 = inter_ctx->scale_fact_rpl0[ref_idx0][1];
            mcp_rpr_c(lc_ctx, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, 0, ref_idx0,
                      rpr_scale_h0, rpr_scale_v0);
        }
    }
}

static void
ciip_weighted_sum(OVCTUDec*const ctudec, struct OVBuffInfo* tmp_intra, struct OVBuffInfo* tmp_inter,
                  unsigned int x0, unsigned int y0,
                  unsigned int log2_pb_w, unsigned int log2_pb_h)
{
    struct CIIPFunctions *ciip = &ctudec->rcn_funcs.ciip;
    //Compute weight in function of neighboring coding modes
    int x_right  = x0 + (1 << log2_pb_w) - 1;
    int y_bottom = y0 + (1 << log2_pb_h) - 1;
    int mode_abv = ctudec->part_map.cu_mode_x[x_right  >> ctudec->part_ctx->log2_min_cb_s];
    int mode_lft = ctudec->part_map.cu_mode_y[y_bottom >> ctudec->part_ctx->log2_min_cb_s];
    uint8_t cu_intra_abv = mode_abv == OV_INTRA || mode_abv == OV_MIP;
    uint8_t cu_intra_lft = mode_lft == OV_INTRA || mode_lft == OV_MIP;
    int wt = 1 + (cu_intra_abv + cu_intra_lft);

    //Apply weighted sum to the final CIIP predicted block
    struct OVBuffInfo dst = ctudec->rcn_ctx.ctu_buff;
    dst.y  += x0 + y0 * dst.stride;
    tmp_intra->y  += x0 + y0 * tmp_intra->stride;
    tmp_inter->y  += x0 + y0 * tmp_inter->stride;
    ciip->weighted(dst.y, dst.stride, tmp_intra->y, tmp_inter->y, tmp_intra->stride, tmp_inter->stride,
                   1 << log2_pb_w, 1 << log2_pb_h, wt);

    dst.cb += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    tmp_intra->cb += (x0 >> 1) + (y0 >> 1) * tmp_intra->stride_c;
    tmp_inter->cb += (x0 >> 1) + (y0 >> 1) * tmp_inter->stride_c;

    dst.cr += (x0 >> 1) + (y0 >> 1) * dst.stride_c;
    tmp_intra->cr += (x0 >> 1) + (y0 >> 1) * tmp_intra->stride_c;
    tmp_inter->cr += (x0 >> 1) + (y0 >> 1) * tmp_inter->stride_c;

    struct MCFunctions *mc_c = &ctudec->rcn_funcs.mc_c;
    if (log2_pb_w <= 2) {
        mc_c->unidir[0][0](dst.cb, dst.stride_c, tmp_inter->cb, tmp_inter->stride_c, 1 << (log2_pb_h - 1), 0, 0, 1 << (log2_pb_w - 1));
        mc_c->unidir[0][0](dst.cr, dst.stride_c, tmp_inter->cr, tmp_inter->stride_c, 1 << (log2_pb_h - 1), 0, 0, 1 << (log2_pb_w - 1));
    } else {
        ciip->weighted(dst.cb, dst.stride_c, tmp_intra->cb, tmp_inter->cb, tmp_intra->stride_c, tmp_inter->stride_c,
                       1 << (log2_pb_w - 1), 1 << (log2_pb_h - 1), wt);
        ciip->weighted(dst.cr, dst.stride_c, tmp_intra->cr, tmp_inter->cr, tmp_intra->stride_c, tmp_inter->stride_c,
                       1 << (log2_pb_w - 1), 1 << (log2_pb_h - 1), wt);
    }
}

static void
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
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;

    tmp_inter.y  = (OVSample *)&rcn_ctx->data.tmp_inter_l[RCN_CTB_PADDING];
    tmp_inter.cb = (OVSample *)&rcn_ctx->data.tmp_inter_cb[RCN_CTB_PADDING];
    tmp_inter.cr = (OVSample *)&rcn_ctx->data.tmp_inter_cr[RCN_CTB_PADDING];
    tmp_inter.stride   = RCN_CTB_STRIDE;
    tmp_inter.stride_c = RCN_CTB_STRIDE;

    rcn_mcp_b(ctudec, tmp_inter, inter_ctx, part_ctx, mv0, mv1, x0, y0, log2_pb_w, log2_pb_h,
              inter_dir, ref_idx0, ref_idx1);

    //Intra Planar mode
    struct OVBuffInfo tmp_intra = ctudec->rcn_ctx.ctu_buff;
    ctudec->rcn_funcs.intra_pred(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0, y0, log2_pb_w, log2_pb_h, 0);
    ctudec->rcn_funcs.intra_pred_c(&ctudec->rcn_ctx, OVINTRA_PLANAR, x0 >> 1, y0 >> 1, log2_pb_w - 1, log2_pb_h - 1, 0);

    ciip_weighted_sum(ctudec, &tmp_intra, &tmp_inter, x0, y0, log2_pb_w, log2_pb_h);
}

static void
rcn_ciip(OVCTUDec *const ctudec,
         int x0, int y0, int log2_pb_w, int log2_pb_h,
         OVMV mv, uint8_t ref_idx)
{
    //Inter merge mode
    struct OVBuffInfo tmp_inter;
    struct OVRCNCtx    *const rcn_ctx   = &ctudec->rcn_ctx;
    OVSample * tmp_inter_l = (OVSample *)rcn_ctx->data.tmp_inter_l;
    OVSample *tmp_inter_cb = (OVSample *)rcn_ctx->data.tmp_inter_cb;
    OVSample *tmp_inter_cr = (OVSample *)rcn_ctx->data.tmp_inter_cr;

    tmp_inter.y  = &tmp_inter_l [RCN_CTB_PADDING];
    tmp_inter.cb = &tmp_inter_cb[RCN_CTB_PADDING];
    tmp_inter.cr = &tmp_inter_cr[RCN_CTB_PADDING];
    tmp_inter.stride   = RCN_CTB_STRIDE;
    tmp_inter.stride_c = RCN_CTB_STRIDE;

    rcn_mcp(ctudec, tmp_inter, x0, y0, log2_pb_w, log2_pb_h, mv, 0, ref_idx);

    //Intra Planar mode
    struct OVBuffInfo tmp_intra = ctudec->rcn_ctx.ctu_buff;
    ctudec->rcn_funcs.intra_pred(&ctudec->rcn_ctx, &tmp_intra, OVINTRA_PLANAR, x0, y0, log2_pb_w, log2_pb_h, 0);
    ctudec->rcn_funcs.intra_pred_c(&ctudec->rcn_ctx, OVINTRA_PLANAR, x0 >> 1, y0 >> 1, log2_pb_w - 1, log2_pb_h - 1, 0);

    ciip_weighted_sum(ctudec, &tmp_intra, &tmp_inter, x0, y0, log2_pb_w, log2_pb_h);

}

static const int8_t g_angle2mask[GEO_NUM_ANGLES] =
{
    0, -1,  1,  2,  3,  4, -1, -1,
    5, -1, -1,  4,  3,  2,  1, -1,

    0, -1,  1,  2,  3,  4, -1, -1,
    5, -1, -1,  4,  3,  2,  1, -1
};


static void
gpm_weights_and_steps(int split_dir, int log2_pb_w_l, int log2_pb_h_l, int* step_x, int* step_y,
                      int16_t** weight, int cr_scale)
{
static const int8_t g_angle2mirror[GEO_NUM_ANGLES] =
    {
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

    if (g_angle2mirror[angle] == 2) {
        int weight_offset = (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[split_dir][hIdx][wIdx][1])
            * GEO_WEIGHT_MASK_SIZE
            + g_weightOffset[split_dir][hIdx][wIdx][0];

        *step_y = -(int)((GEO_WEIGHT_MASK_SIZE << cr_scale) + (1 << log2_pb_w_l));
        *weight = &g_globalGeoWeights[g_angle2mask[angle]][weight_offset];
    } else if (g_angle2mirror[angle] == 1) {
        int weight_offset = g_weightOffset[split_dir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + (GEO_WEIGHT_MASK_SIZE - 1 - g_weightOffset[split_dir][hIdx][wIdx][0]);
        *step_x = -1u << cr_scale;
        *step_y = (GEO_WEIGHT_MASK_SIZE << cr_scale) + (1 << log2_pb_w_l);
        *weight = &g_globalGeoWeights[g_angle2mask[angle]][weight_offset];
    } else {
        int weight_offset = g_weightOffset[split_dir][hIdx][wIdx][1] * GEO_WEIGHT_MASK_SIZE + g_weightOffset[split_dir][hIdx][wIdx][0];
        *step_y = (GEO_WEIGHT_MASK_SIZE << cr_scale) - (1 << log2_pb_w_l);
        *weight = &g_globalGeoWeights[g_angle2mask[angle]][weight_offset];
    }
}


static void
rcn_gpm_b(OVCTUDec *const ctudec, struct VVCGPM* gpm_ctx,
          int x0, int y0, int log2_pb_w, int log2_pb_h)
{
    struct OVBuffInfo dst = ctudec->rcn_ctx.ctu_buff;
    const struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    OVMV mv0 = gpm_ctx->mv0;
    OVMV mv1 = gpm_ctx->mv1;

    const uint16_t* scale_fact_rpl ;
    scale_fact_rpl = gpm_ctx->inter_dir0 == 1 ? inter_ctx->scale_fact_rpl0[mv0.ref_idx]
                                              : inter_ctx->scale_fact_rpl1[mv0.ref_idx];
    int rpr_scale_h0 = scale_fact_rpl[0];
    int rpr_scale_v0 = scale_fact_rpl[1];

    scale_fact_rpl = gpm_ctx->inter_dir1 == 1 ? inter_ctx->scale_fact_rpl0[mv1.ref_idx]
                                              : inter_ctx->scale_fact_rpl1[mv1.ref_idx];
    int rpr_scale_h1 = scale_fact_rpl[0];
    int rpr_scale_v1 = scale_fact_rpl[1];

    /*The RPR functions apply the same bidir0 filters than non-RPR, therefore use only RPR function*/
    mc_rpr_b_l(ctudec, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, mv0.ref_idx, mv1.ref_idx,
               rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, gpm_ctx);
    mc_rpr_b_c(ctudec, dst, x0, y0, log2_pb_w, log2_pb_h, mv0, mv1, mv0.ref_idx, mv1.ref_idx,
               rpr_scale_h0, rpr_scale_v0, rpr_scale_h1, rpr_scale_v1, gpm_ctx);

}

void
BD_DECL(rcn_init_dmvr_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->dmvr.sad[0] = &rcn_dmvr_sad_8;
    rcn_funcs->dmvr.sad[1] = &rcn_dmvr_sad_16;

    rcn_funcs->dmvr.computeSB[0] = &dmvr_compute_sads_8;
    rcn_funcs->dmvr.computeSB[1] = &dmvr_compute_sads_16;
}

void
BD_DECL(rcn_init_inter_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->rcn_dmvr_mv_refine = &rcn_dmvr_mv_refine;
    rcn_funcs->rcn_bdof_mcp_l     = &rcn_bdof_mcp_l;
    rcn_funcs->rcn_mcp            = &rcn_mcp;
    rcn_funcs->rcn_mcp_b          = &rcn_mcp_b;
    rcn_funcs->rcn_mcp_b_l        = &rcn_mcp_b_l;
    rcn_funcs->rcn_prof_mcp_b_l   = &rcn_prof_mcp_b_l;
    rcn_funcs->rcn_mcp_b_c        = &rcn_mcp_b_c;
    rcn_funcs->rcn_ciip_b         = &rcn_ciip_b;
    rcn_funcs->rcn_ciip           = &rcn_ciip;
    rcn_funcs->rcn_gpm_b          = &rcn_gpm_b;
}

