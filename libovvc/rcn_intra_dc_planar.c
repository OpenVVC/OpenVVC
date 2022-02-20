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

#include <stddef.h>
#include <stdint.h>

#include "ovutils.h"
#include "rcn_structures.h"

#include "bitdepth.h"

static const uint8_t vvc_pdpc_w[3][128] = {
        { 32, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 32, 16, 8, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 32, 32, 16, 16, 8, 8, 4, 4, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static void
intra_dc(const OVSample* const ref_abv, const OVSample* const ref_lft,
         OVSample* const dst, ptrdiff_t dst_stride,
         int log2_pb_w, int log2_pb_h)
{
    OVSample* _dst = dst;
    int dc_val = 0;
    int i, j;

    const int shift = OVMAX(log2_pb_w, log2_pb_h) + (log2_pb_w == log2_pb_h);
    const int offset = ((1 << shift) >> 1);
    const int pb_w = 1 << log2_pb_w;
    const int pb_h = 1 << log2_pb_h;

    if (log2_pb_w >= log2_pb_h) {
        for (i = 0; i < pb_w; i++) {
            dc_val += ref_abv[1 + i];
        }
    }
    if (log2_pb_w <= log2_pb_h) {
        for (i = 0; i < pb_h; i++) {
            dc_val += ref_lft[1 + i];
        }
    }

    dc_val = (dc_val + offset) >> shift;

    for (i = 0; i < pb_h; ++i) {
        for (j = 0; j < pb_w; j++) {
            _dst[j] = dc_val;
        }
        _dst += dst_stride;
    }
}

static void
intra_planar(const OVSample* const ref_abv,
             const OVSample* const ref_lft, OVSample* const dst,
             ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h)
{

    OVSample* _dst = dst;
    const uint32_t pb_w = 1 << log2_pb_w;
    const uint32_t pb_h = 1 << log2_pb_h;
    const uint32_t shift = 1 + log2_pb_w + log2_pb_h;
    const uint32_t offset = 1 << (log2_pb_w + log2_pb_h);
    int value;

    int top_row[128], bottom_row[128];
    int top_right = ref_abv[pb_w + 1];

    for (int j = 0; j < pb_w + 1; j++) {
        bottom_row[j] = ref_lft[pb_h + 1];
    }

    for (int k = 0; k < pb_w; k++) {
        value = ref_abv[k + 1];
        bottom_row[k] -= value; // bottom_left - val
        top_row[k] = value << log2_pb_h;
    }

    for (int y = 0; y < pb_h; y++) {
        int value = (int)ref_lft[y + 1]; // left_value y
        int hor_pred = value << log2_pb_w;
        int right_pred = top_right - value;

        for (int x = 0; x < pb_w; x++) {
            int vertPred;
            hor_pred += right_pred;
            top_row[x] += bottom_row[x];

            vertPred = top_row[x];

            _dst[x] = ((hor_pred << log2_pb_h) +
                       (vertPred << log2_pb_w) + offset) >> shift;
        }
        _dst += dst_stride;
    }
}

static void
intra_dc_pdpc(const OVSample* const ref_abv,
              const OVSample* const ref_lft, OVSample* const dst,
              ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h)
{
    OVSample* _dst = dst;
    int i;

    const int shift = OVMAX(log2_pb_w, log2_pb_h) + (log2_pb_w == log2_pb_h);
    const int offset = ((1 << shift) >> 1);

    const int pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;
    const uint8_t* pdpc_w = vvc_pdpc_w[pdpc_scale];
    const uint32_t pb_w = 1 << log2_pb_w;
    const uint32_t pb_h = 1 << log2_pb_h;
    uint32_t dc_val = 0;

    if (log2_pb_w >= log2_pb_h) {
        for (i = 0; i < pb_w; i++) {
            dc_val += ref_abv[1 + i];
        }
    }
    if (log2_pb_w <= log2_pb_h) {
        for (i = 0; i < pb_h; i++) {
            dc_val += ref_lft[1 + i];
        }
    }

    dc_val = (dc_val + offset) >> shift;

    for (int y = 0; y < pb_h; y++) {
        int x;
        int l_wgh = pdpc_w[y];
        const int16_t l_val = ref_lft[y + 1];
        for (x = 0; x < pb_w; x++) {
            const int16_t t_val = ref_abv[x + 1];
            int t_wgh = pdpc_w[x];
            int val = ((t_wgh * l_val) + (l_wgh * t_val) +
                       (64 - (t_wgh + l_wgh)) * dc_val + 32) >> 6;
            _dst[x] = ov_bdclip(val);
        }
        _dst += dst_stride;
    }
}

static void
intra_planar_pdpc(const OVSample* const ref_abv,
                  const OVSample* const ref_lft, OVSample* const dst,
                  ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h)
{

    OVSample* _dst = dst;
    const uint32_t pb_w = 1 << log2_pb_w;
    const uint32_t pb_h = 1 << log2_pb_h;
    const uint32_t w_scale = OVMAX(1, log2_pb_w);
    const uint32_t h_scale = OVMAX(1, log2_pb_h);
    const uint32_t s_shift = w_scale + h_scale + 1;
    const uint32_t offset = 1 << (w_scale + h_scale);

    const uint8_t pdpc_scale = (log2_pb_w + log2_pb_h - 2) >> 2;

    int32_t t_row[128], b_row[128];
    const uint8_t* pdpc_w = vvc_pdpc_w[pdpc_scale];
    const int16_t bl_val = ref_lft[pb_h + 1];
    const int16_t tr_val = ref_abv[pb_w + 1];
    int y;

    for (y = 0; y < pb_w; ++y) {
        int32_t t_val = ref_abv[y + 1];
        b_row[y] = bl_val - t_val;
        t_row[y] = t_val << h_scale;
    }

    for (y = 0; y < pb_h; ++y) {
        int x;
        int32_t l_val = ref_lft[y + 1];
        int y_wgh = pdpc_w[y];
        int32_t r_col_y = tr_val - l_val;
        int32_t l_col_y = l_val << w_scale;
        for (x = 0; x < pb_w; ++x) {
            int32_t val;
            int x_wgh = pdpc_w[x];
            int32_t t_val = ref_abv[x + 1];
            l_col_y += r_col_y;
            t_row[x] += b_row[x];
            val = ((l_col_y << h_scale) + (t_row[x] << w_scale) + offset) >> s_shift;
            val = ((x_wgh * l_val) + (y_wgh * t_val)
                   + (64 - (x_wgh + y_wgh)) * val + 32) >> 6;
            _dst[x] = ov_bdclip(val);
        }
        _dst += dst_stride;
    }

}

void
BD_DECL(rcn_init_dc_planar_functions)(struct RCNFunctions *const rcn_funcs)
{
    rcn_funcs->dc.func = &intra_dc;
    rcn_funcs->planar.func = &intra_planar;
    for (uint8_t i = 0; i < 5; i++)
    {
        for (uint8_t j = 0; j < 5; j++)
        {
            rcn_funcs->dc.pdpc[i][j] = &intra_dc_pdpc;
            rcn_funcs->planar.pdpc[i][j] = &intra_planar_pdpc;
        }
    }
    

}
