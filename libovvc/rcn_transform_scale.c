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

/* Operation performed on residual after residual coefficient have been
 * decoded and transform has been performed before adding them to
 * Prediction Block
 */

#include <stdlib.h>
#include "rcn_transform_scale.h"
#include "ovutils.h"
#include "ctudec.h"
#include "bitdepth.h"

static void
scale_add_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                   int log2_tb_w, int log2_tb_h,
                   int scale)
{
    int i, j;
    int32_t value;
    uint16_t sign;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value = _src[j];
            sign  = value & (1 << 15);
            value = (ov_bdclip(abs(value)) * scale + (1 << (11 - 1))) >> 11;
            value = (ov_clip(sign ? -value : value ,-(1 << 15),1 << 15));
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
scale_sub_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                   int log2_tb_w, int log2_tb_h,
                   int scale)
{
    int i, j;
    int32_t value;
    uint16_t sign;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value = -_src[j];
            sign  = value & (1 << 15);
            value = (ov_bdclip(abs(value)) * scale + (1 << (11 - 1))) >> 11;
            value = (ov_clip(sign ? -value : value ,-(1 << 15),1 << 15));
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
scale_add_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                        int log2_tb_w, int log2_tb_h,
                        int scale)
{
    int i, j;
    int32_t value;
    uint16_t sign;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value = _src[j] >> 1;
            sign  = value & (1 << 15);
            value = (ov_bdclip(abs(value)) * scale + (1 << (11 - 1))) >> 11;
            value = (ov_clip(sign ? -value : value ,-(1 << 15),1 << 15));
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
scale_sub_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                        int log2_tb_w, int log2_tb_h,
                        int scale)
{
    int i, j;
    int32_t value;
    uint16_t sign;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value = (-_src[j]) >> 1;
            sign  = value & (1 << 15);
            value = (abs(value) * scale + (1 << (11 - 1))) >> 11;
            value = (ov_clip(sign ? -value : value ,-(1 << 15),1 << 15));
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
add_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
             int log2_tb_w, int log2_tb_h,
             int scale)
{
    int i, j;
    int32_t value;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value   = _src[j];
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
sub_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
             int log2_tb_w, int log2_tb_h,
             int scale)
{
    int i, j;
    int32_t value;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value   = -_src[j];
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
add_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                  int log2_tb_w, int log2_tb_h,
                  int scale)
{
    int i, j;
    int32_t value;
    const int16_t *_src = src;
    uint16_t       *_dst = dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value   = _src[j] >> 1;
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}

static void
sub_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                  int log2_tb_w, int log2_tb_h,
                  int scale)
{
    int i, j;
    int32_t value;
    const int16_t *_src = src;
    int16_t       *_dst = (int16_t *)dst;
    const int tb_w = 1 << log2_tb_w;
    const int tb_h = 1 << log2_tb_h;
    for (i = 0; i < tb_h; ++i){
        for (j = 0; j < tb_w; ++j){
            value   = (-_src[j]) >> 1;
            _dst[j] = ov_bdclip((int32_t)_dst[j] + value);
        }
        _dst += dst_stride;
        _src += tb_w;
    }
}


BD_DECL(const struct ResidualScaleFunctions scale_add_resid) =
{
    .add_residual            = add_residual,
    .sub_residual            = sub_residual,
    .add_half_residual       = add_half_residual,
    .sub_half_residual       = sub_half_residual,
    .scale_add_residual      = scale_add_residual,
    .scale_sub_residual      = scale_sub_residual,
    .scale_add_half_residual = scale_add_half_residual,
    .scale_sub_half_residual = scale_sub_half_residual
};

