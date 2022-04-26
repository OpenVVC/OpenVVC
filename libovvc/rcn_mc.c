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
#include <string.h>

#include "ovutils.h"

#include "rcn_structures.h"

#include "bitdepth.h"

#define MAX_PB_SIZE 128

#define EPEL_EXTRA_BEFORE 1
#define EPEL_EXTRA_AFTER 2
#define EPEL_EXTRA EPEL_EXTRA_BEFORE + EPEL_EXTRA_AFTER

#define QPEL_EXTRA_BEFORE 3
#define QPEL_EXTRA_AFTER 4
#define QPEL_EXTRA QPEL_EXTRA_BEFORE + QPEL_EXTRA_AFTER

#define MCP_FILTER_L(src, stride, filter)                                      \
        (filter[0] * src[x - 3 * stride] + filter[1] * src[x - 2 * stride] +   \
         filter[2] * src[x - stride] + filter[3] * src[x] +                    \
         filter[4] * src[x + stride] + filter[5] * src[x + 2 * stride] +       \
         filter[6] * src[x + 3 * stride] + filter[7] * src[x + 4 * stride])

#define MCP_FILTER_C(src, stride, filter)                                      \
        (filter[0] * src[x - stride] + filter[1] * src[x] +                    \
         filter[2] * src[x + stride] + filter[3] * src[x + 2 * stride])

#define BLN_FILTER_L(src, stride, filter)                                      \
        (filter[0] * src[x - 1 * stride] + filter[1] * src[x - 0 * stride])

static const int8_t ov_mcp_filters_c[31][4] =
{
    //{  0, 64,  0,  0 },
    { -1, 63,  2,  0 },
    { -2, 62,  4,  0 },
    { -2, 60,  7, -1 },
    { -2, 58, 10, -2 },
    { -3, 57, 12, -2 },
    { -4, 56, 14, -2 },
    { -4, 55, 15, -2 },
    { -4, 54, 16, -2 },
    { -5, 53, 18, -2 },
    { -6, 52, 20, -2 },
    { -6, 49, 24, -3 },
    { -6, 46, 28, -4 },
    { -5, 44, 29, -4 },
    { -4, 42, 30, -4 },
    { -4, 39, 33, -4 },
    { -4, 36, 36, -4 },
    { -4, 33, 39, -4 },
    { -4, 30, 42, -4 },
    { -4, 29, 44, -5 },
    { -4, 28, 46, -6 },
    { -3, 24, 49, -6 },
    { -2, 20, 52, -6 },
    { -2, 18, 53, -5 },
    { -2, 16, 54, -4 },
    { -2, 15, 55, -4 },
    { -2, 14, 56, -4 },
    { -2, 12, 57, -3 },
    { -2, 10, 58, -2 },
    { -1,  7, 60, -2 },
    {  0,  4, 62, -2 },
    {  0,  2, 63, -1 },
};

static const int8_t ov_mc_filters[16][8] =
{
    {   0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  -1, 2,  -5, 62,  8,  -3,  1,  0 },
    {  -1, 3,  -8, 60, 13,  -4,  1,  0 },
    {  -1, 4, -10, 58, 17,  -5,  1,  0 },
    {  -1, 4, -11, 52, 26,  -8,  3, -1 },
    {  -1, 3,  -9, 47, 31, -10,  4, -1 },
    {  -1, 4, -11, 45, 34, -10,  4, -1 },
    {  -1, 4, -11, 40, 40, -11,  4, -1 },
    {  -1, 4, -10, 34, 45, -11,  4, -1 },
    {  -1, 4, -10, 31, 47,  -9,  3, -1 },
    {  -1, 3,  -8, 26, 52, -11,  4, -1 },
    {   0, 1,  -5, 17, 58, -10,  4, -1 },
    {   0, 1,  -4, 13, 60,  -8,  3, -1 },
    {   0, 1,  -3,  8, 62,  -5,  2, -1 },
    {   0, 1,  -2,  4, 63,  -3,  1,  0 },

    //Hpel for amvr
    {  0, 3, 9, 20, 20, 9, 3, 0 }
};

static const int8_t ov_mc_filters_4[15][8] =
{
    {  0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  0, 1,  -5, 62,  8,  -3,  1,  0 },
    {  0, 2,  -8, 60, 13,  -4,  1,  0 },
    {  0, 3, -10, 58, 17,  -5,  1,  0 }, //1/4
    {  0, 3, -11, 52, 26,  -8,  2,  0 },
    {  0, 2,  -9, 47, 31, -10,  3,  0 },
    {  0, 3, -11, 45, 34, -10,  3,  0 },
    {  0, 3, -11, 40, 40, -11,  3,  0 }, //1/2
    {  0, 3, -10, 34, 45, -11,  3,  0 },
    {  0, 3, -10, 31, 47,  -9,  2,  0 },
    {  0, 2,  -8, 26, 52, -11,  3,  0 },
    {  0, 1,  -5, 17, 58, -10,  3,  0 }, //3/4
    {  0, 1,  -4, 13, 60,  -8,  2,  0 },
    {  0, 1,  -3,  8, 62,  -5,  1,  0 },
    {  0, 1,  -2,  4, 63,  -3,  1,  0 }
};

static const int8_t ov_bilinear_filters_4[15][2] =
{
  /*{ 16,  0, },*/
  { 15,  1, },
  { 14,  2, },
  { 13, 3, },
  { 12, 4, },
  { 11, 5, },
  { 10, 6, },
  { 9, 7, },
  { 8, 8, },
  { 7, 9, },
  { 6, 10, },
  { 5, 11, },
  { 4, 12, },
  { 3, 13, },
  { 2, 14, },
  { 1, 15, }
};


static const int8_t ov_mc_filters_rpr[6][16][8] =
{
    {{  0, 0,   0, 64,  0,  0,  0,  0 },
    {   0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  -1, 2,  -5, 62,  8,  -3,  1,  0 },
    {  -1, 3,  -8, 60, 13,  -4,  1,  0 },
    {  -1, 4, -10, 58, 17,  -5,  1,  0 },
    {  -1, 4, -11, 52, 26,  -8,  3, -1 },
    {  -1, 3,  -9, 47, 31, -10,  4, -1 },
    {  -1, 4, -11, 45, 34, -10,  4, -1 },
    {  -1, 4, -11, 40, 40, -11,  4, -1 },
    {  -1, 4, -10, 34, 45, -11,  4, -1 },
    {  -1, 4, -10, 31, 47,  -9,  3, -1 },
    {  -1, 3,  -8, 26, 52, -11,  4, -1 },
    {   0, 1,  -5, 17, 58, -10,  4, -1 },
    {   0, 1,  -4, 13, 60,  -8,  3, -1 },
    {   0, 1,  -3,  8, 62,  -5,  2, -1 },
    {   0, 1,  -2,  4, 63,  -3,  1,  0 }},

    // 1.5x
    {{ -1, -5, 17, 42, 17, -5, -1,  0 },
    {  0, -5, 15, 41, 19, -5, -1,  0 },
    {  0, -5, 13, 40, 21, -4, -1,  0 },
    {  0, -5, 11, 39, 24, -4, -2,  1 },
    {  0, -5,  9, 38, 26, -3, -2,  1 },
    {  0, -5,  7, 38, 28, -2, -3,  1 },
    {  1, -5,  5, 36, 30, -1, -3,  1 },
    {  1, -4,  3, 35, 32,  0, -4,  1 },
    {  1, -4,  2, 33, 33,  2, -4,  1 },
    {  1, -4,  0, 32, 35,  3, -4,  1 },
    {  1, -3, -1, 30, 36,  5, -5,  1 },
    {  1, -3, -2, 28, 38,  7, -5,  0 },
    {  1, -2, -3, 26, 38,  9, -5,  0 },
    {  1, -2, -4, 24, 39, 11, -5,  0 },
    {  0, -1, -4, 21, 40, 13, -5,  0 },
    {  0, -1, -5, 19, 41, 15, -5,  0 }},

    // 2.0x
    {{ -4,  2, 20, 28, 20,  2, -4,  0 },
    { -4,  0, 19, 29, 21,  5, -4, -2 },
    { -4, -1, 18, 29, 22,  6, -4, -2 },
    { -4, -1, 16, 29, 23,  7, -4, -2 },
    { -4, -1, 16, 28, 24,  7, -4, -2 },
    { -4, -1, 14, 28, 25,  8, -4, -2 },
    { -3, -3, 14, 27, 26,  9, -3, -3 },
    { -3, -1, 12, 28, 25, 10, -4, -3 },
    { -3, -3, 11, 27, 27, 11, -3, -3 },
    { -3, -4, 10, 25, 28, 12, -1, -3 },
    { -3, -3,  9, 26, 27, 14, -3, -3 },
    { -2, -4,  8, 25, 28, 14, -1, -4 },
    { -2, -4,  7, 24, 28, 16, -1, -4 },
    { -2, -4,  7, 23, 29, 16, -1, -4 },
    { -2, -4,  6, 22, 29, 18, -1, -4 },
    { -2, -4,  5, 21, 29, 19,  0, -4 }},

    // Affine
    {{  0, 0,  0, 64,  0,   0,  0,  0 },
    {  0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  0, 1,  -5, 62,  8,  -3,  1,  0 },
    {  0, 2,  -8, 60, 13,  -4,  1,  0 },
    {  0, 3, -10, 58, 17,  -5,  1,  0 }, //1/4
    {  0, 3, -11, 52, 26,  -8,  2,  0 },
    {  0, 2,  -9, 47, 31, -10,  3,  0 },
    {  0, 3, -11, 45, 34, -10,  3,  0 },
    {  0, 3, -11, 40, 40, -11,  3,  0 }, //1/2
    {  0, 3, -10, 34, 45, -11,  3,  0 },
    {  0, 3, -10, 31, 47,  -9,  2,  0 },
    {  0, 2,  -8, 26, 52, -11,  3,  0 },
    {  0, 1,  -5, 17, 58, -10,  3,  0 }, //3/4
    {  0, 1,  -4, 13, 60,  -8,  2,  0 },
    {  0, 1,  -3,  8, 62,  -5,  1,  0 },
    {  0, 1,  -2,  4, 63,  -3,  1,  0 }},
    // 1.5x
    {{  0, -6, 17, 42, 17, -5, -1,  0 },
      {  0, -5, 15, 41, 19, -5, -1,  0 },
      {  0, -5, 13, 40, 21, -4, -1,  0 },
      {  0, -5, 11, 39, 24, -4, -1,  0 },
      {  0, -5,  9, 38, 26, -3, -1,  0 },
      {  0, -5,  7, 38, 28, -2, -2,  0 },
      {  0, -4,  5, 36, 30, -1, -2,  0 },
      {  0, -3,  3, 35, 32,  0, -3,  0 },
      {  0, -3,  2, 33, 33,  2, -3,  0 },
      {  0, -3,  0, 32, 35,  3, -3,  0 },
      {  0, -2, -1, 30, 36,  5, -4,  0 },
      {  0, -2, -2, 28, 38,  7, -5,  0 },
      {  0, -1, -3, 26, 38,  9, -5,  0 },
      {  0, -1, -4, 24, 39, 11, -5,  0 },
      {  0, -1, -4, 21, 40, 13, -5,  0 },
      {  0, -1, -5, 19, 41, 15, -5,  0 }},

    // 2x
    {{  0, -2, 20, 28, 20,  2, -4,  0 },
      {  0, -4, 19, 29, 21,  5, -6,  0 },
      {  0, -5, 18, 29, 22,  6, -6,  0 },
      {  0, -5, 16, 29, 23,  7, -6,  0 },
      {  0, -5, 16, 28, 24,  7, -6,  0 },
      {  0, -5, 14, 28, 25,  8, -6,  0 },
      {  0, -6, 14, 27, 26,  9, -6,  0 },
      {  0, -4, 12, 28, 25, 10, -7,  0 },
      {  0, -6, 11, 27, 27, 11, -6,  0 },
      {  0, -7, 10, 25, 28, 12, -4,  0 },
      {  0, -6,  9, 26, 27, 14, -6,  0 },
      {  0, -6,  8, 25, 28, 14, -5,  0 },
      {  0, -6,  7, 24, 28, 16, -5,  0 },
      {  0, -6,  7, 23, 29, 16, -5,  0 },
      {  0, -6,  6, 22, 29, 18, -5,  0 },
      {  0, -6,  5, 21, 29, 19, -4,  0 }}
};

static const int8_t ov_mc_filters_rpr_c[3][32][4] =
{
    {
        {  0, 64,  0,  0 },
        { -1, 63,  2,  0 },
        { -2, 62,  4,  0 },
        { -2, 60,  7, -1 },
        { -2, 58, 10, -2 },
        { -3, 57, 12, -2 },
        { -4, 56, 14, -2 },
        { -4, 55, 15, -2 },
        { -4, 54, 16, -2 },
        { -5, 53, 18, -2 },
        { -6, 52, 20, -2 },
        { -6, 49, 24, -3 },
        { -6, 46, 28, -4 },
        { -5, 44, 29, -4 },
        { -4, 42, 30, -4 },
        { -4, 39, 33, -4 },
        { -4, 36, 36, -4 },
        { -4, 33, 39, -4 },
        { -4, 30, 42, -4 },
        { -4, 29, 44, -5 },
        { -4, 28, 46, -6 },
        { -3, 24, 49, -6 },
        { -2, 20, 52, -6 },
        { -2, 18, 53, -5 },
        { -2, 16, 54, -4 },
        { -2, 15, 55, -4 },
        { -2, 14, 56, -4 },
        { -2, 12, 57, -3 },
        { -2, 10, 58, -2 },
        { -1,  7, 60, -2 },
        {  0,  4, 62, -2 },
        {  0,  2, 63, -1 }
    },
    {
        { 12, 40, 12,  0 },
        { 11, 40, 13,  0 },
        { 10, 40, 15, -1 },
        {  9, 40, 16, -1 },
        {  8, 40, 17, -1 },
        {  8, 39, 18, -1 },
        {  7, 39, 19, -1 },
        {  6, 38, 21, -1 },
        {  5, 38, 22, -1 },
        {  4, 38, 23, -1 },
        {  4, 37, 24, -1 },
        {  3, 36, 25,  0 },
        {  3, 35, 26,  0 },
        {  2, 34, 28,  0 },
        {  2, 33, 29,  0 },
        {  1, 33, 30,  0 },
        {  1, 31, 31,  1 },
        {  0, 30, 33,  1 },
        {  0, 29, 33,  2 },
        {  0, 28, 34,  2 },
        {  0, 26, 35,  3 },
        {  0, 25, 36,  3 },
        { -1, 24, 37,  4 },
        { -1, 23, 38,  4 },
        { -1, 22, 38,  5 },
        { -1, 21, 38,  6 },
        { -1, 19, 39,  7 },
        { -1, 18, 39,  8 },
        { -1, 17, 40,  8 },
        { -1, 16, 40,  9 },
        { -1, 15, 40, 10 },
        {  0, 13, 40, 11 }
    },

    {
        { 17, 30, 17,  0 },
        { 17, 30, 18, -1 },
        { 16, 30, 18,  0 },
        { 16, 30, 18,  0 },
        { 15, 30, 18,  1 },
        { 14, 30, 18,  2 },
        { 13, 29, 19,  3 },
        { 13, 29, 19,  3 },
        { 12, 29, 20,  3 },
        { 11, 28, 21,  4 },
        { 10, 28, 22,  4 },
        { 10, 27, 22,  5 },
        {  9, 27, 23,  5 },
        {  9, 26, 24,  5 },
        {  8, 26, 24,  6 },
        {  7, 26, 25,  6 },
        {  7, 25, 25,  7 },
        {  6, 25, 26,  7 },
        {  6, 24, 26,  8 },
        {  5, 24, 26,  9 },
        {  5, 23, 27,  9 },
        {  5, 22, 27, 10 },
        {  4, 22, 28, 10 },
        {  4, 21, 28, 11 },
        {  3, 20, 29, 12 },
        {  3, 19, 29, 13 },
        {  3, 19, 29, 13 },
        {  2, 18, 30, 14 },
        {  1, 18, 30, 15 },
        {  0, 18, 30, 16 },
        {  0, 18, 30, 16 },
        { -1, 18, 30, 17 }
    }
};

static void
put_vvc_pel_uni_pixels(OVSample* _dst, ptrdiff_t _dststride,
                       const OVSample* _src, ptrdiff_t _srcstride, int height,
                       intptr_t mx, intptr_t my, int width)
{
    int y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;

    for (y = 0; y < height; y++) {
        memcpy(dst, src, width * sizeof(OVSample));
        src += srcstride;
        dst += dststride;
    }
}

/*FIXME It might actually be better to use one function for bi pred instead of
 * 2*/
static void
put_vvc_pel_pixels(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = src[x] << (14 - BITDEPTH);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_pel_bi_pixels(OVSample* _dst, ptrdiff_t _dststride,
                      const OVSample* _src0, ptrdiff_t _srcstride,
                      const int16_t* _src1, int height, intptr_t mx,
                      intptr_t my, int width)
{
    int x, y;
    const OVSample* src0 = _src0;
    const int16_t* src1 = _src1;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 1;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip(((src0[x] << (14 - BITDEPTH)) + src1[x] + offset)
                               >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_uni_h(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_uni_v(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, srcstride, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_uni_hv(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width)
{
    int x, y;
    const int8_t* filter;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= QPEL_EXTRA_BEFORE * srcstride;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6)
                                + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_uni_pixels(OVSample* _dst, ptrdiff_t _dststride,
                        const OVSample* _src, ptrdiff_t _srcstride, int height,
                        int denom, int wx, int offset2,
                        intptr_t mx, intptr_t my, int width)
{
    int y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        int x;
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((src[x] << (14 - BITDEPTH)) * wx + offset) >> shift) + offset2);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_uni_h(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height,
                   int denom, int wx, int offset2,
                   intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8)) * wx + offset) >> shift) + offset2);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_uni_v(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height,
                   int denom, int wx, int offset2,
                   intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((MCP_FILTER_L(src, srcstride, filter) >> (BITDEPTH - 8)) * wx + offset) >> shift) + offset2);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_uni_hv(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                    ptrdiff_t _srcstride, int height,
                    int denom, int wx, int offset2,
                    intptr_t mx, intptr_t my,
                    int width)
{
    int x, y;
    const int8_t* filter;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= QPEL_EXTRA_BEFORE * srcstride;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6) * wx + offset) >> shift) + offset2);
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_uni_h_c(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height,
                   int denom, int wx, int offset2,
                   intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8)) * wx + offset) >> shift) + offset2);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_uni_v_c(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height,
                   int denom, int wx, int offset2,
                   intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[my - 1];
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8)) * wx + offset) >> shift) + offset2);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_uni_hv_c(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                    ptrdiff_t _srcstride, int height,
                    int denom, int wx, int offset2,
                    intptr_t mx, intptr_t my,
                    int width)
{
    int x, y;
    const int8_t* filter;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    int shift = denom + 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= EPEL_EXTRA_BEFORE * srcstride;
    filter =  ov_mcp_filters_c[mx - 1];
    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6) * wx + offset) >> shift) + offset2);
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_pel_rpr(uint16_t* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = src[x] << (14 - BITDEPTH);
        }
        src += srcstride;
        dst += dststride;
    }
}


static void
put_vvc_qpel_rpr_h(uint16_t* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const OVSample* src = (OVSample*)_src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mc_filters_rpr[filter_idx][mx];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_rpr_clip_v(OVSample* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const int16_t* src =  (int16_t*)_src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mc_filters_rpr[filter_idx][my];    
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, srcstride, filter) >> 6) +
                                offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}


static void
put_vvc_qpel_rpr_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const int16_t* src =  (int16_t*)_src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mc_filters_rpr[filter_idx][my];    

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
             dst[x] = MCP_FILTER_L(src, srcstride, filter) >> 6;
        }
        src += srcstride;
        dst += dststride;
    }
}


static void
put_vvc_pel_rpr_clip(OVSample* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip((src[x] + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}


static void
put_vvc_pel_rpr_bi(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const uint16_t* src = _src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = src[x] ;
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_rpr_bi_sum(OVSample* _dst, ptrdiff_t _dststride,
                      const uint16_t* _src0, ptrdiff_t _src0stride,
                      const uint16_t* _src1, ptrdiff_t _src1stride,
                      int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const int16_t* src0 = (int16_t*)_src0;
    const int16_t* src1 = (int16_t*)_src1;
    ptrdiff_t src0stride = _src0stride;
    ptrdiff_t src1stride = _src1stride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 1;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip((src0[x] + src1[x] + offset) >>
                                   shift);
        }
        src0 += src0stride;
        src1 += src1stride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_rpr_weighted(OVSample* _dst, ptrdiff_t _dststride,
                      const uint16_t* _src0, ptrdiff_t _src0stride,
                      const uint16_t* _src1, ptrdiff_t _src1stride,
                      int height, int denom, int wx0, int wx1,
                      intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const int16_t* src0 = (int16_t*)_src0;
    const int16_t* src1 = (int16_t*)_src1;
    ptrdiff_t src0stride = _src0stride;
    ptrdiff_t src1stride = _src1stride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + denom + 1;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip((src0[x]*wx0 + src1[x]*wx1 + offset) >>
                                   shift);
        }
        src0 += src0stride;
        src1 += src1stride;
        dst += dststride;
    }
}

static void
put_vvc_epel_rpr_h(uint16_t* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;
    ptrdiff_t srcstride = _srcstride;
    uint16_t* dst = (uint16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mc_filters_rpr_c[filter_idx][mx];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_epel_rpr_clip_v(OVSample* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const int16_t* src =  (int16_t*)_src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mc_filters_rpr_c[filter_idx][my];    
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, srcstride, filter) >> 6) +
                                offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_epel_rpr_bi_v(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width, uint8_t filter_idx)
{
    int x, y;
    const int16_t* src =  (int16_t*)_src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mc_filters_rpr_c[filter_idx][my];    

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(src, srcstride, filter) >> 6;
        }
        src += srcstride;
        dst += dststride;
    }
}


static void
put_vvc_pel_bilinear_pixels(uint16_t* _dst, ptrdiff_t _dststride,
                            const OVSample* _src, ptrdiff_t _srcstride, int height,
                            intptr_t mx, intptr_t my, int width)
{
    int y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;

    for (y = 0; y < height; y++) {
        int x;
        for (x = 0; x < width; ++x) {
            #if (BITDEPTH - 10) > 0
            dst [x] = (src[x] + (1 << (BITDEPTH - 10)) >> (BITDEPTH - 10);
            #else
            dst [x] = src[x] << (10 - BITDEPTH);
            #endif
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_h(uint16_t* _dst, ptrdiff_t _dststride, const OVSample* _src,
                        ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                        int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_bilinear_filters_4[mx - 1];
    int shift = 4 - (10 - BITDEPTH);
    int offset = 1 << (shift - 1);

    src += 1;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (BLN_FILTER_L(src, 1, filter) + offset) >> shift;
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_v(uint16_t* _dst, ptrdiff_t _dststride, const OVSample* _src,
                        ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                        int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_bilinear_filters_4[my - 1];
    int shift = 4 - (10 - BITDEPTH);
    int offset = 1 << (shift - 1);

    src += srcstride;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (BLN_FILTER_L(src, srcstride, filter) + offset) >> shift;
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bilinear_hv(uint16_t* _dst, ptrdiff_t _dststride, const OVSample* _src,
                         ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                         int width)
{
    int x, y;
    const int8_t* filter;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    int16_t* dst = (int16_t*)_dst;
    ptrdiff_t dststride = _dststride;
    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    int shift = 4 - (10 - BITDEPTH);
    int offset = 1 << (shift - 1);

    filter = ov_bilinear_filters_4[mx - 1];

    for (y = 0; y < height + 2; y++) {
        for (x = 0; x < width + 2; x++) {
            tmp[x] = (BLN_FILTER_L(src, 1, filter) + offset) >> shift;
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + MAX_PB_SIZE + 1;
    filter = ov_bilinear_filters_4[my - 1];

    shift = 4;
    offset = 1 << (shift - 1);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (BLN_FILTER_L(tmp, MAX_PB_SIZE, filter) + offset) >> shift;
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_h(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_qpel_v(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*)_src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(src, srcstride, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_qpel_hv(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const OVSample* src = (OVSample*)_src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    src -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6;
        }
        tmp += MAX_PB_SIZE;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_qpel_bi_h(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src0 = _src0;
    const int16_t* src1 = _src1;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src0, 1, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bi_v(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src0 = _src0;
    const int16_t* src1 = _src1;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src0, srcstride, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_qpel_bi_hv(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const OVSample* src0 = _src0;
    const int16_t* src1 = _src1;
    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mc_filters[mx - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src0 -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src0, 1, filter) >> (BITDEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mc_filters[my - 1];
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6)
                                + src1[x] + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_epel_uni_h(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_epel_uni_v(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                   ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                   int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[my - 1];
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8))
                                + offset) >> shift);
        }
        src += srcstride;
        dst += dststride;
    }
}

static void
put_vvc_epel_uni_hv(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src,
                    ptrdiff_t _srcstride, int height, intptr_t mx, intptr_t my,
                    int width)
{
    int x, y;
    const OVSample* src = _src;
    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;
    int shift = 14 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6)
                                + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_epel_h(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const OVSample* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_epel_v(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
               int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const OVSample* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_epel_hv(int16_t* _dst, const OVSample* _src, ptrdiff_t _srcstride,
                int height, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = _src;

    int16_t* dst = (int16_t*)_dst;

    ptrdiff_t srcstride = _srcstride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6;
        }
        tmp += MAX_PB_SIZE;
        dst += MAX_PB_SIZE;
    }
}

static void
put_vvc_epel_bi_h(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src0 = _src0;
    const int16_t* src1 = (int16_t*)_src1;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src0, 1, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_vvc_epel_bi_v(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src0,
                  ptrdiff_t _srcstride, const int16_t* _src1, int height,
                  intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src0 = _src0;
    const int16_t* src1 = _src1;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[my - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src0, srcstride, filter) >> (BITDEPTH - 8))
                                + src1[x] + offset) >> shift);
        }
        src0 += srcstride;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_vvc_epel_bi_hv(OVSample* _dst, ptrdiff_t _dststride, const OVSample* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src0 = _src0;
    const int16_t* src1 = _src1;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];

    int shift = 14 + 1 - BITDEPTH;
    int offset = 1 << (shift - 1);

    src0 -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src0, 1, filter) >> (BITDEPTH - 8);
        }
        src0 += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6)
                                + src1[x] + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src1 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_epel_bi_h(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;
    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;
    shift = log2Wd + 1;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_epel_bi_v(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;


    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter = ov_mcp_filters_c[my - 1];
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;
    shift = log2Wd + 1;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(src, srcstride, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_epel_bi_hv(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                   int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                   int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;


    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    int16_t tmp_array[(MAX_PB_SIZE + EPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter = ov_mcp_filters_c[mx - 1];
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;
    shift = log2Wd + 1;

    src -= EPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + EPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_C(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + EPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = ov_mcp_filters_c[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_C(tmp, MAX_PB_SIZE, filter) >> 6) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        tmp += MAX_PB_SIZE;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_pel_bi_pixels(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;

    ptrdiff_t srcstride = _srcstride;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;
    shift = log2Wd + 1;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip(((src2[x] * wx0 + (int32_t)((uint32_t)(src[x] * wx1) << (14 - BITDEPTH))
                                 + offset)  >> shift));
        }
        src2 += MAX_PB_SIZE;
        src += srcstride;
        dst += dststride;
    }
}

static void
put_weighted_qpel_bi_h(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;


    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;
    shift = log2Wd + 1;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_qpel_bi_v(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                  int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                  int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;
    const OVSample* src = (OVSample*) _src;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;
    shift = log2Wd + 1;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(src, srcstride, filter) >> (BITDEPTH - 8)) * wx1
                                + src2[x] * wx0 + offset) >> shift);
        }
        src += srcstride;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}

static void
put_weighted_qpel_bi_hv(uint8_t* _dst, ptrdiff_t _dststride, uint8_t* _src, ptrdiff_t _srcstride,
                   int16_t* src2, ptrdiff_t src2stride, int height, int denom,
                   int wx0, int wx1, int offset0, int offset1, intptr_t mx, intptr_t my, int width)
{
    int x, y;

    const OVSample* src = (OVSample*) _src;

    OVSample* dst = (OVSample*)_dst;

    ptrdiff_t srcstride = _srcstride;
    ptrdiff_t dststride = _dststride;

    int16_t tmp_array[(MAX_PB_SIZE + QPEL_EXTRA) * MAX_PB_SIZE];
    int16_t* tmp = tmp_array;

    const int8_t* filter;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[mx - 1] : ov_mc_filters[mx - 1];
    //int shift = 14 + denom -BITDEPTH;
    int shift = 14 + 1 - BITDEPTH;
    int log2Wd = denom + shift - 1;
    int offset = (offset0 + offset1 + 1) << log2Wd;

    src -= QPEL_EXTRA_BEFORE * srcstride;

    for (y = 0; y < height + QPEL_EXTRA; y++) {
        for (x = 0; x < width; x++) {
            tmp[x] = MCP_FILTER_L(src, 1, filter) >> (BITDEPTH - 8);
        }
        src += srcstride;
        tmp += MAX_PB_SIZE;
    }

    tmp = tmp_array + QPEL_EXTRA_BEFORE * MAX_PB_SIZE;
    filter = width == 4 && height == 4 ? ov_mc_filters_4[my - 1] : ov_mc_filters[my - 1];

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = ov_bdclip(((MCP_FILTER_L(tmp, MAX_PB_SIZE, filter) >> 6) * wx1
                                + src2[x] * wx0 + offset) >> (log2Wd + 1));
        }
        tmp += MAX_PB_SIZE;
        src2 += MAX_PB_SIZE;
        dst += dststride;
    }
}


static void
put_weighted_ciip_pixels(OVSample* dst, int dststride,
                      const OVSample* src_intra, const OVSample* src_inter, int intra_stride, int inter_stride,
                      int width, int height, int wt)
{
    int x, y;
    int shift  = 2;
    int offset = (1 << (shift - 1));
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            dst[x] = ov_bdclip(((src_intra[x] * wt + src_inter[x] * (4 - wt)) + offset) >> shift);
        }
        src_intra += intra_stride;
        src_inter += inter_stride;
        dst += dststride;
    }
}


static void
put_weighted_gpm_bi_pixels(OVSample* _dst, int _dststride, const int16_t* _src0,
                  int srcstride0, const int16_t* _src1, int srcstride1, int height,
                  int width, int step_x, int step_y, int16_t* weight)
{
    int x, y;
    const int16_t* src0 = _src0;
    const int16_t* src1 = _src1;
    OVSample* dst = (OVSample*)_dst;
    ptrdiff_t dststride = _dststride;
    int shift = 14 - BITDEPTH + 3;
    int offset = (1 << (shift - 1)) ;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; ++x) {
            int w0 = weight[0];
            // int w0 = 0;
            int w1 = 8 - w0;
            dst[x] = ov_bdclip(((src1[x] * w1 + src0[x] * w0 + offset) >> shift));
            weight += step_x;
        }
        src1 += srcstride1;
        src0 += srcstride0;
        dst += dststride;
        weight += step_y;
    }
}

void
BD_DECL(rcn_init_mc_functions)(struct RCNFunctions *const rcn_funcs)
{
    struct MCFunctions *const mc_l = &rcn_funcs->mc_l;
    struct MCFunctions *const mc_c = &rcn_funcs->mc_c;

    int i;
    mc_l->rpr_sum = &put_vvc_qpel_rpr_bi_sum;
    mc_c->rpr_sum = &put_vvc_qpel_rpr_bi_sum;

    for (i = 0; i < 8; ++i) {
        mc_l->rpr_w[i] = &put_vvc_qpel_rpr_weighted;
        mc_c->rpr_w[i] = &put_vvc_qpel_rpr_weighted;

        /* Luma functions */
        mc_l->unidir[0][i] = &put_vvc_pel_uni_pixels;
        mc_l->unidir[1][i] = &put_vvc_qpel_uni_h;
        mc_l->unidir[2][i] = &put_vvc_qpel_uni_v;
        mc_l->unidir[3][i] = &put_vvc_qpel_uni_hv;

        mc_l->bidir0[0][i] = &put_vvc_pel_pixels;
        mc_l->bidir0[1][i] = &put_vvc_qpel_h;
        mc_l->bidir0[2][i] = &put_vvc_qpel_v;
        mc_l->bidir0[3][i] = &put_vvc_qpel_hv;

        mc_l->bidir1[0][i] = &put_vvc_pel_bi_pixels;
        mc_l->bidir1[1][i] = &put_vvc_qpel_bi_h;
        mc_l->bidir1[2][i] = &put_vvc_qpel_bi_v;
        mc_l->bidir1[3][i] = &put_vvc_qpel_bi_hv;

        mc_l->unidir_w[0][i] = &put_weighted_uni_pixels;
        mc_l->unidir_w[1][i] = &put_weighted_uni_h;
        mc_l->unidir_w[2][i] = &put_weighted_uni_v;
        mc_l->unidir_w[3][i] = &put_weighted_uni_hv;

        mc_l->bidir_w[0][i] = &put_weighted_pel_bi_pixels;
        mc_l->bidir_w[1][i] = &put_weighted_qpel_bi_h;
        mc_l->bidir_w[2][i] = &put_weighted_qpel_bi_v;
        mc_l->bidir_w[3][i] = &put_weighted_qpel_bi_hv;

        mc_l->bilinear[0][i] = &put_vvc_pel_bilinear_pixels;
        mc_l->bilinear[1][i] = &put_vvc_qpel_bilinear_h;
        mc_l->bilinear[2][i] = &put_vvc_qpel_bilinear_v;
        mc_l->bilinear[3][i] = &put_vvc_qpel_bilinear_hv;

        mc_l->rpr_h[0][i] = &put_vvc_pel_rpr;
        mc_l->rpr_h[1][i] = &put_vvc_qpel_rpr_h;
        mc_l->rpr_v_uni[0][i] = &put_vvc_pel_rpr_clip;
        mc_l->rpr_v_uni[1][i] = &put_vvc_qpel_rpr_clip_v;
        mc_l->rpr_v_bi[0][i] = &put_vvc_pel_rpr_bi;
        mc_l->rpr_v_bi[1][i] = &put_vvc_qpel_rpr_bi_v;

        /* Chroma functions */
        mc_c->unidir[0][i] = &put_vvc_pel_uni_pixels;
        mc_c->unidir[1][i] = &put_vvc_epel_uni_h;
        mc_c->unidir[2][i] = &put_vvc_epel_uni_v;
        mc_c->unidir[3][i] = &put_vvc_epel_uni_hv;

        mc_c->bidir0[0][i] = &put_vvc_pel_pixels;
        mc_c->bidir0[1][i] = &put_vvc_epel_h;
        mc_c->bidir0[2][i] = &put_vvc_epel_v;
        mc_c->bidir0[3][i] = &put_vvc_epel_hv;

        mc_c->bidir1[0][i] = &put_vvc_pel_bi_pixels;
        mc_c->bidir1[1][i] = &put_vvc_epel_bi_h;
        mc_c->bidir1[2][i] = &put_vvc_epel_bi_v;
        mc_c->bidir1[3][i] = &put_vvc_epel_bi_hv;

        mc_c->unidir_w[0][i] = &put_weighted_uni_pixels;
        mc_c->unidir_w[1][i] = &put_weighted_uni_h_c;
        mc_c->unidir_w[2][i] = &put_weighted_uni_v_c;
        mc_c->unidir_w[3][i] = &put_weighted_uni_hv_c;

        mc_c->bidir_w[0][i] = &put_weighted_pel_bi_pixels;
        mc_c->bidir_w[1][i] = &put_weighted_epel_bi_h;
        mc_c->bidir_w[2][i] = &put_weighted_epel_bi_v;
        mc_c->bidir_w[3][i] = &put_weighted_epel_bi_hv;

        mc_c->rpr_h[0][i] = &put_vvc_pel_rpr;
        mc_c->rpr_h[1][i] = &put_vvc_epel_rpr_h;
        mc_c->rpr_v_uni[0][i] = &put_vvc_pel_rpr_clip;
        mc_c->rpr_v_uni[1][i] = &put_vvc_epel_rpr_clip_v;
        mc_c->rpr_v_bi[0][i] = &put_vvc_pel_rpr_bi;
        mc_c->rpr_v_bi[1][i] = &put_vvc_epel_rpr_bi_v;
    }
    mc_l->gpm_weighted = &put_weighted_gpm_bi_pixels;
    mc_c->gpm_weighted = &put_weighted_gpm_bi_pixels;
}

void
BD_DECL(rcn_init_ciip_functions)(struct RCNFunctions *const rcn_funcs)
{
    struct CIIPFunctions *const ciip = &rcn_funcs->ciip;
    ciip->weighted = &put_weighted_ciip_pixels;
}
