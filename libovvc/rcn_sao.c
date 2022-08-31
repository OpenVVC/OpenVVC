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

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "ovutils.h"
#include "ovdefs.h"
#include "ovframe.h"
#include "ovdpb.h"
#include "slicedec.h"

#include "bitdepth.h"

static void
sao_band_filter(OVSample *_dst, OVSample *_src,
                ptrdiff_t stride_dst, ptrdiff_t stride_src,
                int width, int height,
                int8_t offset_val[],
                uint8_t band_pos)
{
    OVSample *dst = _dst;
    OVSample *src = _src;
    int offset_table[32] = { 0 };
    int k, y, x;
    int shift  = BITDEPTH - 5;

    for (k = 0; k < 4; k++) {
        offset_table[(k + band_pos) & 31] = offset_val[k];
    }

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(src[x] + offset_table[src[x] >> shift]);
        }
        dst += stride_dst;
        src += stride_src;
    }
}

#define CMP(a, b) ((a) > (b) ? 1 : ((a) == (b) ? 0 : -1))
static void
sao_edge_filter(OVSample *_dst, OVSample *_src,
                ptrdiff_t stride_dst, ptrdiff_t stride_src,
                int width, int height,
                int8_t offset_val[],
                uint8_t eo_dir)
{
    static const int8_t pos[4][2][2] = {
        { { -1,  0 }, {  1, 0 } }, // horizontal
        { {  0, -1 }, {  0, 1 } }, // vertical
        { { -1, -1 }, {  1, 1 } }, // 45 degree
        { {  1, -1 }, { -1, 1 } }, // 135 degree
    };

    const int16_t sao_offset_val[5] = {
        offset_val[0],
        offset_val[1],
        0,
        offset_val[2],
        offset_val[3]
    };

    OVSample *dst = _dst;
    OVSample *src = _src;

    int src_offset = 0;
    int dst_offset = 0;
    int x, y;

    int a_stride = pos[eo_dir][0][0] + pos[eo_dir][0][1] * stride_src;
    int b_stride = pos[eo_dir][1][0] + pos[eo_dir][1][1] * stride_src;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            int diff0 = CMP(src[x + src_offset], src[x + src_offset + a_stride]);
            int diff1 = CMP(src[x + src_offset], src[x + src_offset + b_stride]);
            int val   = 2 + diff0 + diff1;
            dst[x + dst_offset] = ov_bdclip(src[x + src_offset] + sao_offset_val[val]);
        }
        src_offset += stride_src;
        dst_offset += stride_dst;
    }
}

static void
sao_edge_filter2(OVSample *dst, OVSample *src_row, OVSample *src_col,
                ptrdiff_t stride_dst, ptrdiff_t stride_src,
                int width, int height,
                int8_t offset_val[],
                uint8_t eo_dir)
{
    const int16_t sao_offset_val[5] = {
        offset_val[0],
        offset_val[1],
        0,
        offset_val[2],
        offset_val[3]
    };

    int x, y;

    OVSample tmp_a[256];
    OVSample tmp_b[256];
    OVSample tmp_c[256];

    if (eo_dir == 0) {
        for (y = 0; y < height; y++) {
            memcpy(tmp_c    , dst      , sizeof(OVSample) * (width + 1));
            memcpy(tmp_a + 1, tmp_c    , sizeof(OVSample) * (width - 1));
            memcpy(tmp_b    , tmp_c + 1, sizeof(OVSample) * width);
            tmp_a[0] = src_col[0];
            for (x = 0; x < width; x++) {
                int diff0 = CMP(tmp_c[x], tmp_a[x]);
                int diff1 = CMP(tmp_c[x], tmp_b[x]);
                int val   = 2 + diff0 + diff1;
                dst[x] = ov_bdclip(tmp_c[x] + sao_offset_val[val]);
            }
            ++src_col;
            dst += stride_dst;
        }
    } else if (eo_dir == 1) {
        memcpy(tmp_a, src_row, sizeof(OVSample) * width);
        memcpy(tmp_c, dst    , sizeof(OVSample) * width);
        for (y = 0; y < height; y++) {
            memcpy(tmp_b, dst + stride_dst, sizeof(OVSample) * width);
            for (x = 0; x < width; x++) {
                int diff0 = CMP(tmp_c[x], tmp_a[x]);
                int diff1 = CMP(tmp_c[x], tmp_b[x]);
                int val   = 2 + diff0 + diff1;
                dst[x] = ov_bdclip(tmp_c[x] + sao_offset_val[val]);
            }
            memcpy(tmp_a, tmp_c, sizeof(OVSample) * width);
            memcpy(tmp_c, tmp_b, sizeof(OVSample) * width);
            dst += stride_dst;
        }
    } else if (eo_dir == 2) {
        OVSample tmp;
        memcpy(tmp_a + 1, src_row, sizeof(OVSample) * (width - 1));
        memcpy(tmp_c, dst        , sizeof(OVSample) * width);
        tmp_a[0] = src_col[-1];
        for (y = 0; y < height; y++) {
            memcpy(tmp_b, dst + stride_dst + 1, sizeof(OVSample) * width);
            tmp = dst[stride_dst];
            for (x = 0; x < width; x++) {
                int diff0 = CMP(tmp_c[x], tmp_a[x]);
                int diff1 = CMP(tmp_c[x], tmp_b[x]);
                int val   = 2 + diff0 + diff1;
                dst[x] = ov_bdclip(tmp_c[x] + sao_offset_val[val]);
            }
            memcpy(tmp_a + 1, tmp_c, sizeof(OVSample) * (width - 1));
            memcpy(tmp_c + 1, tmp_b, sizeof(OVSample) * (width - 1));
            tmp_a[0] = src_col[0];
            tmp_c[0] = tmp;
            ++src_col;
            dst += stride_dst;
        }
    } else {
        memcpy(tmp_a, src_row + 1, sizeof(OVSample) * width);
        memcpy(tmp_c, dst        , sizeof(OVSample) * (width + 1));
        for (y = 0; y < height; y++) {
            memcpy(tmp_b + 1, dst + stride_dst, sizeof(OVSample) * (width + 1));
            tmp_b[0] = src_col[1];
            for (x = 0; x < width; x++) {
                int diff0 = CMP(tmp_c[x], tmp_a[x]);
                int diff1 = CMP(tmp_c[x], tmp_b[x]);
                int val   = 2 + diff0 + diff1;
                dst[x] = ov_bdclip(tmp_c[x] + sao_offset_val[val]);
            }
            memcpy(tmp_a, tmp_c + 1, sizeof(OVSample) * width);
            memcpy(tmp_c, tmp_b + 1, sizeof(OVSample) * (width + 1));
            ++src_col;
            dst += stride_dst;
        }
    }
}

struct SAOBuff
{
    OVSample col[256];
    OVSample row[256];
};

static void
fill_sao_buff(struct SAOBuff *dst, const OVSample *src, int src_stride, int width, int height)
{
    const OVSample *src_row = src - src_stride;
    const OVSample *src_col = src - src_stride - 1;

    for (int y = 0; y < height + 2; y++) {
        dst->col[y] = src_col[0];
        src_col += src_stride;
    }

    memcpy(dst->row, src_row, sizeof(OVSample) * (width + 1));
}

static void
sao_call(const struct SAOFunctions *saofunc, struct SAOBuff *tmp, OVSample *dst, OVSample *src, int16_t dst_stride, int16_t src_stride,
         int width, int height, uint8_t is_border, int8_t *offsets, uint8_t mode_info, uint8_t type_idx)
{

    if (!offsets[0] && !offsets[1] && !offsets[2] && !offsets[3]) return;

    switch (type_idx) {
        case SAO_BAND:

            saofunc->band(dst, dst, dst_stride, dst_stride,
                          width, height, offsets, mode_info);

            break;

        case SAO_EDGE:
        {
            uint8_t eo_dir = mode_info;

            int x_start = (is_border & OV_BOUNDARY_LEFT_RECT)  && eo_dir != 1;
            int y_start = (is_border & OV_BOUNDARY_UPPER_RECT) && eo_dir != 0;

            int dst_offset = y_start * dst_stride + x_start;
            int src_offset = y_start * src_stride + x_start;

            /* Not vertical */
            width -= (is_border & OV_BOUNDARY_LEFT_RECT)  && eo_dir != 1;
            width -= (is_border & OV_BOUNDARY_RIGHT_RECT) && eo_dir != 1;

            /* Not horizontal */
            height -= (is_border & OV_BOUNDARY_UPPER_RECT)  && eo_dir != 0;
            height -= (is_border & OV_BOUNDARY_BOTTOM_RECT) && eo_dir != 0;

#if 0
            saofunc->edge[!(width % 8)](dst + dst_offset, src + src_offset, dst_stride,
                                        src_stride, width, height, offsets, eo_dir);
#else

            if (src_offset)
                fill_sao_buff(tmp, src + src_offset, src_stride, width, height);

            OVSample *src_col = &tmp->col[1];
            OVSample *src_row = &tmp->row[0];

            sao_edge_filter2(dst + dst_offset, src_row, src_col, dst_stride,
                             src_stride, width, height, offsets, eo_dir);
#endif

            break;
        }
    }
}

static void
rcn_sao_ctu(OVCTUDec *const ctudec, SAOParamsCtu *sao, int x_pic, int y_pic, int y_end_pic, int fb_offset, uint8_t is_border)
{

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    int ctb_w = (OVMIN(  (1 << log2_ctb_s), ctudec->pic_w - x_pic));
    int ctb_h = (OVMIN((y_end_pic - y_pic), ctudec->pic_h - y_pic));

    const OVFrame *frame = ctudec->rcn_ctx.frame_start;
    struct OVFilterBuffers* fb = &ctudec->rcn_ctx.filter_buffers;
    struct SAOBuff tmp;

    const struct SAOFunctions *sao_func = &ctudec->rcn_funcs.sao;

    int sao_l  = sao->sao_ctu_flag        & 0x3;
    int sao_cb = (sao->sao_ctu_flag >> 2) & 0x3;
    int sao_cr = sao_cb;

    if (sao_l) {
        ptrdiff_t dst_stride = frame->linesize[0] / sizeof(OVSample);
        ptrdiff_t src_stride = fb->filter_region_stride[0];

        int dst_offset = y_pic * dst_stride + x_pic;
        int src_offset = fb->filter_region_offset[0] + src_stride * fb_offset;

        int c_idx = 0;
        OVSample *dst = (OVSample *)frame->data[c_idx] + dst_offset;
        OVSample *src = (OVSample *)fb->filter_region[c_idx] + src_offset;

        uint8_t mode_info = sao->mode_info[c_idx];
        int8_t *offsets = sao->offset[c_idx];

        fill_sao_buff(&tmp, src, src_stride, ctb_w, ctb_h);

        sao_call(sao_func, &tmp, dst, src, dst_stride, src_stride, ctb_w, ctb_h, is_border,
                 offsets, mode_info, sao_l);

    }

    if (sao_cb || sao_cr) {
        ptrdiff_t dst_stride = frame->linesize[1] / sizeof(OVSample);
        ptrdiff_t src_stride = fb->filter_region_stride[1];

        int dst_offset = (y_pic >> 1) * dst_stride + (x_pic >> 1);;
        int src_offset = fb->filter_region_offset[1] + src_stride * (fb_offset >> 1);

        ctb_w >>= 1;
        ctb_h >>= 1;

        if (sao_cb) {
            int c_idx = 1;
            OVSample *dst = (OVSample *)frame->data[c_idx]       + dst_offset;
            OVSample *src = (OVSample *)fb->filter_region[c_idx] + src_offset;
            int8_t *offsets = sao->offset[c_idx];
            uint8_t mode_info = sao->mode_info[c_idx];

            fill_sao_buff(&tmp, src, src_stride, ctb_w, ctb_h);

            sao_call(sao_func, &tmp, dst, src, dst_stride, src_stride, ctb_w, ctb_h, is_border,
                     offsets, mode_info, sao_cb);
        }

        if (sao_cr) {
            int c_idx = 2;
            OVSample *dst = (OVSample *)frame->data[c_idx] + dst_offset;
            OVSample *src = (OVSample *)fb->filter_region[c_idx] + src_offset;
            int8_t *offsets = sao->offset[c_idx];
            uint8_t mode_info = sao->mode_info[c_idx];

            fill_sao_buff(&tmp, src, src_stride, ctb_w, ctb_h);

            sao_call(sao_func, &tmp, dst, src, dst_stride, src_stride, ctb_w, ctb_h, is_border,
                     offsets, mode_info, sao_cb);
        }
    }
}

static void
rcn_extend_filter_region2(struct OVRCNCtx *const rcn_ctx, int x_l,
                            int x_pic, int y_pic, uint8_t bnd_msk)
{
    struct OVFilterBuffers* fb = &rcn_ctx->filter_buffers;
    const OVFrame *f = rcn_ctx->frame_start;
    const int width_l  = (x_pic + fb->filter_region_w[0] > f->width) ? (f->width - x_pic)
                                                                     : fb->filter_region_w[0];
    const int height_l = (y_pic + fb->filter_region_h[0] > f->height) ? (f->height - y_pic)
                                                                      : fb->filter_region_h[0];

    for (int comp = 0; comp < 3; comp++) {
        const int width  = width_l >> (comp != 0);
        const int height = height_l >> (comp != 0);
        int stride_pic = f->linesize[comp] / sizeof(OVSample);
        const int pic_offset = (y_pic >> (comp != 0)) * stride_pic + (x_pic >> (comp != 0));

        const OVSample* src_pic = (OVSample*)f->data[comp] + pic_offset;

        int stride_filter = fb->filter_region_stride[comp];

        OVSample *const dst_0 = fb->filter_region[comp] + fb->filter_region_offset[comp];

        const OVSample *src_line = fb->saved_rows_sao[comp] + (x_l >> (comp != 0)) + 2 * fb->saved_rows_stride[comp];

        uint8_t not_bnd_rgt = !(bnd_msk & OV_BOUNDARY_RIGHT_RECT);
        uint8_t not_bnd_btm = !(bnd_msk & OV_BOUNDARY_BOTTOM_RECT);
        int cpy_s = sizeof(OVSample) * (width + not_bnd_rgt);

        if (!(bnd_msk & OV_BOUNDARY_LEFT_RECT)) {
            OVSample* filter_region = fb->filter_region[comp];
            int stride = fb->filter_region_stride[comp];

            const int width  = fb->filter_region_w[comp];
            const int height = height_l >> (comp != 0);

            OVSample *dst = dst_0 - stride - 1;
            for(int i = 0; i < height + 1; ++i) {
                dst[0] = dst[width];
                dst[1] = dst[width + 1];
                dst += stride;
            }
        }
        if (1) {
            OVSample *dst = dst_0;
            const OVSample *src = src_pic;

            memcpy(dst_0 - stride_filter, src_line, cpy_s);
            memcpy(dst_0                , src,      cpy_s);
            memcpy(dst_0 + (height - 1) * stride_filter, src + (height - 1)            * stride_pic, cpy_s);
            memcpy(dst_0 +  height      * stride_filter, src + (height - !not_bnd_btm) * stride_pic, cpy_s);

            for(int i = 0; i < height; ++i) {
                dst[width - 1] = src[width - 1];
                dst[width]     = src[width - !not_bnd_rgt];
                dst += stride_filter;
                src += stride_pic;
            }
        }

        if (bnd_msk & OV_BOUNDARY_LEFT_RECT) {
            OVSample *dst = dst_0 - 1;
            const OVSample *src = src_pic;
            for (int i = 0; i < height; ++i) {
                dst[0] = dst[1] = src[0];
                dst += stride_filter;
                src += stride_pic;
            }
        }
    }
}

static void
line_alloc(OVSample **dst, const OVFrame *f)
{
    dst[0] = ov_malloc(f->linesize[0] + 256);
    dst[1] = ov_malloc(f->linesize[1] + 256);
    dst[2] = ov_malloc(f->linesize[2] + 256);
}

static void
line_free(OVSample **dst)
{
    ov_free(dst[0]);
    ov_free(dst[1]);
    ov_free(dst[2]);
}

static void
line_to_saved(OVSample **lbck, struct OVFilterBuffers* fb)
{
    for (int comp = 0; comp < 3; comp++) {
        OVSample* saved_rows = fb->saved_rows_sao[comp];
        OVSample* src = lbck[comp];

        int stride_filter = fb->filter_region_stride[comp];

        const int width = fb->saved_rows_stride[comp];
        int stride_rows = fb->saved_rows_stride[comp];

        memcpy(&saved_rows[2 * stride_rows], src, width * sizeof(OVSample));
    }
}

static void
backup_line(OVSample **dst, const OVFrame *f, int32_t y)
{

    if (y + 5 < f->height) {
        uint8_t *src_y = f->data[0];
        ptrdiff_t offset_y = (y + 5) * f->linesize[0];

        memcpy(dst[0], src_y  + offset_y, f->linesize[0]);
    }

    if ((y >> 1) + 2 < (f->height >> 1)) {
        uint8_t *src_cb = f->data[1];
        uint8_t *src_cr = f->data[2];
        ptrdiff_t offset_c = ((y >> 1) + 2) * f->linesize[1];
        memcpy(dst[1], src_cb + offset_c, f->linesize[1]);
        memcpy(dst[2], src_cr + offset_c, f->linesize[2]);
    }
}

static void
backup_line2(OVSample **dst, const OVFrame *f, int32_t y)
{

    uint8_t *src_y = f->data[0];
    ptrdiff_t offset_y = y * f->linesize[0];

    memcpy(dst[0], src_y  + offset_y, f->linesize[0]);

    uint8_t *src_cb = f->data[1];
    uint8_t *src_cr = f->data[2];
    ptrdiff_t offset_c = (y >> 1) * f->linesize[1];
    memcpy(dst[1], src_cb + offset_c, f->linesize[1]);
    memcpy(dst[2], src_cr + offset_c, f->linesize[2]);
}

static void
backup_line3(OVSample **dst, const OVFrame *f, int32_t y)
{

    if (y > f->height) return;

    uint8_t *src_y = f->data[0];
    ptrdiff_t offset_y = (y - 1) * f->linesize[0];

    memcpy(dst[0], src_y  + offset_y, f->linesize[0]);

    uint8_t *src_cb = f->data[1];
    uint8_t *src_cr = f->data[2];
    ptrdiff_t offset_c = ((y >> 1) - 1) * f->linesize[1];
    memcpy(dst[1], src_cb + offset_c, f->linesize[1]);
    memcpy(dst[2], src_cr + offset_c, f->linesize[2]);
}

static void
rcn_sao_filter_line(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y)
{
    if (!ctudec->sao_info.sao_luma_flag && !ctudec->sao_info.sao_chroma_flag){
        return;
    }
    OVSample *lbck[3];
    OVSample *lbck2[3];

    line_alloc(lbck, ctudec->rcn_ctx.frame_start);
    line_alloc(lbck2, ctudec->rcn_ctx.frame_start);

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    int ctu_w  = 1 << log2_ctb_s;

    struct OVFilterBuffers* fb = &ctudec->rcn_ctx.filter_buffers;
    int margin = 2 * fb->margin;

    uint8_t border_init = -(ctb_y == einfo->nb_ctu_h - 1) & OV_BOUNDARY_BOTTOM_RECT;

    int y_pic = (ctb_y + einfo->ctb_y) << log2_ctb_s;

    int x_offset = einfo->ctb_x << log2_ctb_s;
    int x_pos = 0;

    int ctb_addr_rs = ctb_y * einfo->nb_ctu_w;

    if (!border_init)
        backup_line(lbck, ctudec->rcn_ctx.frame_start, y_pic + ctu_w);
    else if (!ctb_y) {
        backup_line2(lbck, ctudec->rcn_ctx.frame_start, y_pic);
        line_to_saved(lbck, &ctudec->rcn_ctx.filter_buffers);
    }

    backup_line3(lbck2, ctudec->rcn_ctx.frame_start, y_pic + ctu_w);

    for (int ctb_x = 0; ctb_x < einfo->nb_ctu_w; ++ctb_x) {
        uint8_t is_border = border_init;
        is_border |= -(ctb_x == 0)                   & OV_BOUNDARY_LEFT_RECT;
        is_border |= -(ctb_x == einfo->nb_ctu_w - 1) & OV_BOUNDARY_RIGHT_RECT;

        int x_pic = x_pos + x_offset;

        SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];

        rcn_extend_filter_region2(&ctudec->rcn_ctx, x_pos, x_pic,
                                  y_pic + margin, is_border);

        if (sao->sao_ctu_flag) {
            rcn_sao_ctu(ctudec, sao, x_pic, y_pic + margin, y_pic + ctu_w, 0, is_border);
        }

#if 0
        /* Apply partial SAO filter on next CTU line for ALF */
        if (!(is_border & OV_BOUNDARY_BOTTOM_RECT)) {
            sao = &ctudec->sao_info.sao_params[ctb_addr_rs + einfo->nb_ctu_w];
            if (sao->sao_ctu_flag) {
                rcn_sao_ctu(ctudec, sao, x_pic, y_pic + ctu_w, y_pic + ctu_w + margin,
                            ctu_w - margin, is_border);
            }
        }
#endif

        x_pos += ctu_w;
        ++ctb_addr_rs;
    }

    if (!(border_init & OV_BOUNDARY_BOTTOM_RECT)) {
        int x_pos = 0;

        ctb_y++;
        y_pic += ctu_w;

        line_to_saved(lbck2, &ctudec->rcn_ctx.filter_buffers);

        uint8_t border_init2 = -(ctb_y == einfo->nb_ctu_h - 1) & OV_BOUNDARY_BOTTOM_RECT;
        int ctb_addr_rs = ctb_y * einfo->nb_ctu_w;

        for (int ctb_x = 0; ctb_x < einfo->nb_ctu_w; ++ctb_x) {
            uint8_t is_border = border_init2;
            is_border |= -(ctb_x == 0)                   & OV_BOUNDARY_LEFT_RECT;
            is_border |= -(ctb_x == einfo->nb_ctu_w - 1) & OV_BOUNDARY_RIGHT_RECT;

            int x_pic = x_pos + x_offset;

            SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];

            rcn_extend_filter_region2(&ctudec->rcn_ctx, x_pos, x_pic,
                                      y_pic, is_border);

            if (sao->sao_ctu_flag) {
                is_border = 0;
            is_border |= -(ctb_x == 0)                   & OV_BOUNDARY_LEFT_RECT;
            is_border |= -(ctb_x == einfo->nb_ctu_w - 1) & OV_BOUNDARY_RIGHT_RECT;
                rcn_sao_ctu(ctudec, sao, x_pic, y_pic, y_pic + margin,
                            0, is_border);
            }

            x_pos += ctu_w;
            ++ctb_addr_rs;
        }

        line_to_saved(lbck, &ctudec->rcn_ctx.filter_buffers);
    }

    if (!border_init)
        line_to_saved(lbck, &ctudec->rcn_ctx.filter_buffers);

    line_free(lbck);
    line_free(lbck2);
}

static void
rcn_sao_first_pix_rows(OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y)
{
    if (!ctudec->sao_info.sao_luma_flag && !ctudec->sao_info.sao_chroma_flag){
        return;
    }

    const OVPartInfo *const pinfo = ctudec->part_ctx;
    uint8_t log2_ctb_s = pinfo->log2_ctu_s;
    int y_pic = (ctb_y + einfo->ctb_y) << log2_ctb_s;
    uint8_t border_init = -(ctb_y == einfo->nb_ctu_h - 1) & OV_BOUNDARY_BOTTOM_RECT;

    struct OVFilterBuffers* fb = &ctudec->rcn_ctx.filter_buffers;
    int margin = 2 * fb->margin;

    OVSample *lbck[3];

    line_alloc(lbck, ctudec->rcn_ctx.frame_start);

    if (!border_init)
        backup_line(lbck, ctudec->rcn_ctx.frame_start, y_pic);
    else {
        backup_line2(lbck, ctudec->rcn_ctx.frame_start, y_pic);
        line_to_saved(lbck, &ctudec->rcn_ctx.filter_buffers);
    }

    for (int ctb_x = 0; ctb_x < einfo->nb_ctu_w; ctb_x++) {
        int ctb_x_pic = ctb_x + einfo->ctb_x;
        int x_pos = ctb_x << log2_ctb_s;
        int x_pic = ctb_x_pic << log2_ctb_s;
        
        //left | right | up | down
        uint8_t is_border = 0; 
        is_border = (ctb_y == 0)                   ? is_border | OV_BOUNDARY_UPPER_RECT: is_border;
        is_border = (ctb_y == einfo->nb_ctu_h - 1) ? is_border | OV_BOUNDARY_BOTTOM_RECT: is_border;
        is_border = (ctb_x == 0)                   ? is_border | OV_BOUNDARY_LEFT_RECT: is_border;
        is_border = (ctb_x == einfo->nb_ctu_w - 1) ? is_border | OV_BOUNDARY_RIGHT_RECT: is_border;

        //Apply SAO of previous ctu line
        rcn_extend_filter_region2(&ctudec->rcn_ctx, x_pos, x_pic, y_pic, is_border);

        int fb_offset = 0;
        int ctb_addr_rs    = ctb_y * einfo->nb_ctu_w + ctb_x;
        SAOParamsCtu *sao  = &ctudec->sao_info.sao_params[ctb_addr_rs];
        if (sao->sao_ctu_flag) {
            rcn_sao_ctu(ctudec, sao, x_pos_pic, y_start_pic, y_start_pic + margin, fb_offset, is_border);
        }

    }

    if (!border_init)
        line_to_saved(lbck, &ctudec->rcn_ctx.filter_buffers);

    line_free(lbck);
}

void
BD_DECL(rcn_init_sao_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->sao.band= &sao_band_filter;
    rcn_funcs->sao.edge[0]= &sao_edge_filter;
    rcn_funcs->sao.edge[1]= &sao_edge_filter;
    rcn_funcs->sao.rcn_sao_filter_line    = &rcn_sao_filter_line;
    rcn_funcs->sao.rcn_sao_first_pix_rows = &rcn_sao_first_pix_rows;
}
