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
*   Ibrahim FARHAT
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

#include "rcn_neon.h"

#include <stdint.h>
#include <stddef.h>

#include "ctudec.h"
#include "rcn_transform.h"
#include "rcn.h"
#include "ovutils.h"

void ov_transform_add_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_add_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_add_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
// add half
void ov_transform_add_half_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_add_half_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_add_half_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
// sub assembly
void ov_transform_sub_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_sub_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_sub_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
// sub half
void ov_transform_sub_half_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_sub_half_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
void ov_transform_sub_half_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride);
// add scale
void ov_transform_scale_add_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride, int scale);
void ov_transform_scale_add_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
void ov_transform_scale_add_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
// add half scale
void ov_transform_scale_add_half_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride, int scale);
void ov_transform_scale_add_half_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
void ov_transform_scale_add_half_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
// sub scale
void ov_transform_scale_sub_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride, int scale);
void ov_transform_scale_sub_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
void ov_transform_scale_sub_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
// sub half scale
void ov_transform_scale_sub_half_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride, int scale);
void ov_transform_scale_sub_half_16_2_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);
void ov_transform_scale_sub_half_32_1_10_neon(uint16_t *dst, ptrdiff_t dst_stride,const int16_t *src, ptrdiff_t src_stride,int scale);

static void
ov_vvc_add_residual_8_4_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                int log2_tb_w, int log2_tb_h, int scale){
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_add_8_4_10_neon(_dst, dst_stride<<1,_src, tb_w<<1);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  }else {
    vvc_add_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }

}

static void
ov_vvc_add_residual_16_2_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                             int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_add_16_2_10_neon(_dst, dst_stride<<1,
                                     _src, tb_w<<1);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_add_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}
static void
ov_vvc_add_residual_32_1_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                             int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_add_32_1_10_neon(_dst, dst_stride<<1,
                                   _src, tb_w<<1);
    _dst += dst_stride;
    _src += tb_w;
  }
}

static void
ov_vvc_add_residual_64_1_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                             int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_add_32_1_10_neon(_dst, dst_stride,
                                   _src, tb_w);
    ov_transform_add_32_1_10_neon(_dst+32, dst_stride,
                                   _src+32, tb_w);
    _dst += dst_stride;
    _src += tb_w;
  }
}
// add half residu
static void
ov_vvc_add_half_residual_8_4_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                 int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_add_half_8_4_10_neon(_dst, dst_stride<<1,
                                         _src, tb_w<<1);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_add_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}
static void
ov_vvc_add_half_residual_16_2_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                  int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_add_half_16_2_10_neon(_dst, dst_stride<<1,
                                          _src, tb_w<<1);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_add_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}

static void
ov_vvc_add_half_residual_32_1_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                  int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_add_half_32_1_10_neon(_dst, dst_stride<<1,
                                        _src, tb_w<<1);
    _dst += dst_stride;
    _src += tb_w;
  }
}


// sub residu

static void
ov_vvc_sub_residual_8_4_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                            int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_sub_8_4_10_neon(_dst, dst_stride<<1,
                                    _src, tb_w<<1);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_sub_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}
static void
ov_vvc_sub_residual_16_2_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                             int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_sub_16_2_10_neon(_dst, dst_stride<<1,
                                     _src, tb_w<<1);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_sub_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}
static void
ov_vvc_sub_residual_32_1_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                             int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_sub_32_1_10_neon(_dst, dst_stride<<1,
                                   _src, tb_w<<1);
    _dst += dst_stride;
    _src += tb_w;
  }
}

// sub half
static void
ov_vvc_sub_half_residual_8_4_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                 int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_sub_half_8_4_10_neon(_dst, dst_stride<<1,
                                         _src, tb_w<<1);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_sub_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}

static void
ov_vvc_sub_half_residual_16_2_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                  int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_sub_half_16_2_10_neon(_dst, dst_stride<<1,
                                          _src, tb_w<<1);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_sub_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, 0);
  }
}

static void
ov_vvc_sub_half_residual_32_1_10_neon(const int16_t *const src, uint16_t *const dst, int16_t dst_stride,
                                  int log2_tb_w, int log2_tb_h, int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_sub_half_32_1_10_neon(_dst, dst_stride<<1,
                                        _src, tb_w<<1);
    _dst += dst_stride;
    _src += tb_w;
  }
}

// add scale
static void
ov_vvc_scale_add_residual_8_4_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                  int log2_tb_w, int log2_tb_h,
                                  int scale)
{

  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_scale_add_8_4_10_neon(_dst, dst_stride<<1,
                                          _src, tb_w<<1, scale);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_scale_add_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_add_residual_16_2_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                   int log2_tb_w, int log2_tb_h,
                                   int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_scale_add_16_2_10_neon(_dst, dst_stride<<1,
                                           _src, tb_w<<1, scale);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_scale_add_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_add_residual_32_1_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                   int log2_tb_w, int log2_tb_h,
                                   int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_scale_add_32_1_10_neon(_dst, dst_stride<<1,
                                         _src, tb_w<<1, scale);
    _dst += dst_stride;
    _src += tb_w;
  }
}

static void
ov_vvc_scale_add_half_residual_8_4_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                       int log2_tb_w, int log2_tb_h,
                                       int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_scale_add_half_8_4_10_neon(_dst, dst_stride<<1,
                                               _src, tb_w<<1, scale);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_scale_add_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}
static void
ov_vvc_scale_add_half_residual_16_2_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                        int log2_tb_w, int log2_tb_h,
                                        int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_scale_add_half_16_2_10_neon(_dst, dst_stride<<1,
                                                _src, tb_w<<1, scale);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_scale_add_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_add_half_residual_32_1_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                        int log2_tb_w, int log2_tb_h,
                                        int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_scale_add_half_32_1_10_neon(_dst, dst_stride<<1,
                                              _src, tb_w<<1, scale);
    _dst += dst_stride;
    _src += tb_w;
  }
}

// scale sub
static void
ov_vvc_scale_sub_residual_8_4_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                  int log2_tb_w, int log2_tb_h,
                                  int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_scale_sub_8_4_10_neon(_dst, dst_stride<<1,
                                          _src, tb_w<<1, scale);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_scale_sub_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_sub_residual_16_2_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                   int log2_tb_w, int log2_tb_h,
                                   int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_scale_sub_16_2_10_neon(_dst, dst_stride<<1,
                                           _src, tb_w<<1, scale);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_scale_sub_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_sub_residual_32_1_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                   int log2_tb_w, int log2_tb_h,
                                   int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_scale_sub_32_1_10_neon(_dst, dst_stride<<1,
                                         _src, tb_w<<1, scale);
    _dst += dst_stride;
    _src += tb_w;
  }
}


// scale sub half
static void
ov_vvc_scale_sub_half_residual_8_4_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                       int log2_tb_w, int log2_tb_h,
                                       int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 2; ++i){
      ov_transform_scale_sub_half_8_4_10_neon(_dst, dst_stride<<1,
                                               _src, tb_w<<1, scale);
      _dst += dst_stride << 2;
      _src += tb_w << 2;
    }
  } else {
    vvc_scale_sub_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_sub_half_residual_16_2_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                        int log2_tb_w, int log2_tb_h,
                                        int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  if (log2_tb_h > 1) {
    for (i = 0; i < tb_h >> 1; ++i){
      ov_transform_scale_sub_half_16_2_10_neon(_dst, dst_stride<<1,
                                                _src, tb_w<<1, scale);
      _dst += dst_stride << 1;
      _src += tb_w << 1;
    }
  } else {
    vvc_scale_sub_half_residual(src, dst, dst_stride, log2_tb_w, log2_tb_h, scale);
  }
}

static void
ov_vvc_scale_sub_half_residual_32_1_10_neon(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                                        int log2_tb_w, int log2_tb_h,
                                        int scale)
{
  int i;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  const int16_t *_src = (const int16_t *)src;
  uint16_t *_dst = dst;
  for (i = 0; i < tb_h; ++i){
    ov_transform_scale_sub_half_32_1_10_neon(_dst, dst_stride<<1,
                                              _src, tb_w<<1, scale);
    _dst += dst_stride;
    _src += tb_w;
  }
}
void
rcn_init_ict_functions_neon(struct RCNFunctions *rcn_func, uint8_t type)
{
 rcn_func->ict.add[3] = &ov_vvc_add_residual_8_4_10_neon;
 rcn_func->ict.add[4] = &ov_vvc_add_residual_16_2_10_neon;
 rcn_func->ict.add[5] = &ov_vvc_add_residual_32_1_10_neon;
 rcn_func->ict.add[6] = &ov_vvc_add_residual_64_1_10_neon;
 switch (type)
 {
   case 3:
     rcn_func->ict.ict[3][0] = &ov_vvc_scale_add_residual_8_4_10_neon;
     rcn_func->ict.ict[4][0] = &ov_vvc_scale_add_residual_16_2_10_neon;
     rcn_func->ict.ict[5][0] = &ov_vvc_scale_add_residual_32_1_10_neon;

     rcn_func->ict.ict[3][1] = &ov_vvc_scale_sub_residual_8_4_10_neon;
     rcn_func->ict.ict[4][1] = &ov_vvc_scale_sub_residual_16_2_10_neon;
     rcn_func->ict.ict[5][1] = &ov_vvc_scale_sub_residual_32_1_10_neon;

     rcn_func->ict.ict[3][2] = &ov_vvc_scale_sub_half_residual_8_4_10_neon;
     rcn_func->ict.ict[4][2] = &ov_vvc_scale_sub_half_residual_16_2_10_neon;
     rcn_func->ict.ict[5][2] = &ov_vvc_scale_sub_half_residual_32_1_10_neon;
     break;
   case 2:
     rcn_func->ict.ict[3][0] = &ov_vvc_add_residual_8_4_10_neon;
     rcn_func->ict.ict[4][0] = &ov_vvc_add_residual_16_2_10_neon;
     rcn_func->ict.ict[5][0] = &ov_vvc_add_residual_32_1_10_neon;

     rcn_func->ict.ict[3][1] = &ov_vvc_sub_residual_8_4_10_neon;
     rcn_func->ict.ict[4][1] = &ov_vvc_sub_residual_16_2_10_neon;
     rcn_func->ict.ict[5][1] = &ov_vvc_sub_residual_32_1_10_neon;

     rcn_func->ict.ict[3][2] = &ov_vvc_sub_half_residual_8_4_10_neon;
     rcn_func->ict.ict[4][2] = &ov_vvc_sub_half_residual_16_2_10_neon;
     rcn_func->ict.ict[5][2] = &ov_vvc_sub_half_residual_32_1_10_neon;
     break;
   case 1:
     rcn_func->ict.ict[3][0] = &ov_vvc_scale_add_residual_8_4_10_neon;
     rcn_func->ict.ict[4][0] = &ov_vvc_scale_add_residual_16_2_10_neon;
     rcn_func->ict.ict[5][0] = &ov_vvc_scale_add_residual_32_1_10_neon;

     rcn_func->ict.ict[3][1] = &ov_vvc_scale_add_residual_8_4_10_neon;
     rcn_func->ict.ict[4][1] = &ov_vvc_scale_add_residual_16_2_10_neon;
     rcn_func->ict.ict[5][1] = &ov_vvc_scale_add_residual_32_1_10_neon;

     rcn_func->ict.ict[3][2] = &ov_vvc_scale_add_half_residual_8_4_10_neon;
     rcn_func->ict.ict[4][2] = &ov_vvc_scale_add_half_residual_16_2_10_neon;
     rcn_func->ict.ict[5][2] = &ov_vvc_scale_add_half_residual_32_1_10_neon;
     break;
   default:

     rcn_func->ict.ict[3][0] = &ov_vvc_add_residual_8_4_10_neon;
     rcn_func->ict.ict[4][0] = &ov_vvc_add_residual_16_2_10_neon;
     rcn_func->ict.ict[5][0] = &ov_vvc_add_residual_32_1_10_neon;

     rcn_func->ict.ict[3][1] = &ov_vvc_add_residual_8_4_10_neon;
     rcn_func->ict.ict[4][1] = &ov_vvc_add_residual_16_2_10_neon;
     rcn_func->ict.ict[5][1] = &ov_vvc_add_residual_32_1_10_neon;

     rcn_func->ict.ict[3][2] = &ov_vvc_add_half_residual_8_4_10_neon;
     rcn_func->ict.ict[4][2] = &ov_vvc_add_half_residual_16_2_10_neon;
     rcn_func->ict.ict[5][2] = &ov_vvc_add_half_residual_32_1_10_neon;
     break;
 }
}