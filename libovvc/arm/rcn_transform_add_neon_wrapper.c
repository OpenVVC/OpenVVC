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

#include "simde/x86/sse4.2.h"
#include "simde/x86/mmx.h"
#include "simde/x86/sse.h"
#include "simde/x86/sse2.h"
#include "simde/x86/sse3.h"
#include "simde/x86/ssse3.h"
#include "simde/x86/sse4.2.h"
#include "simde/x86/avx.h"
#include "simde/x86/avx2.h"
#include "simde/x86/avx512.h"
#include <stdint.h>
#include <stddef.h>

#include "ctudec.h"
#include "rcn_transform.h"
#include "rcn.h"
#include "ovutils.h"

#define CLIP_10 ((1 << 10) - 1)
#define SIGN_16 (int16_t)(1 << 15)

#define bd_clip_10_4x128_epi16(a,b,c,d)\
  a = _mm_max_epi16(a, _mm_setzero_si128());\
  b = _mm_max_epi16(b, _mm_setzero_si128());\
  c = _mm_max_epi16(c, _mm_setzero_si128());\
  d = _mm_max_epi16(d, _mm_setzero_si128());\
\
  a = _mm_min_epi16(a, _mm_set1_epi16(CLIP_10));\
  b = _mm_min_epi16(b, _mm_set1_epi16(CLIP_10));\
  c = _mm_min_epi16(c, _mm_set1_epi16(CLIP_10));\
  d = _mm_min_epi16(d, _mm_set1_epi16(CLIP_10));


void
vvc_add_residual_8_4_10_neon(uint16_t *dst, ptrdiff_t dst_stride,
                             const int16_t *src, ptrdiff_t src_stride);



void
rcn_init_ict_functions_neon(struct RCNFunctions *rcn_func, uint8_t type)
{
 rcn_func->ict.add[3] = &vvc_add_residual_8_4_10_neon;
 /*rcn_func->ict.add[4] = &vvc_add_residual_16_2_10_sse;
 rcn_func->ict.add[5] = &vvc_add_residual_32_1_10_sse;
 rcn_func->ict.add[6] = &vvc_add_residual_64_1_10_sse;
 switch (type)
 {
   case 3:
     rcn_func->ict.ict[3][0] = &vvc_scale_add_residual_8_4_10_sse;
     rcn_func->ict.ict[4][0] = &vvc_scale_add_residual_16_2_10_sse;
     rcn_func->ict.ict[5][0] = &vvc_scale_add_residual_32_1_10_sse;

     rcn_func->ict.ict[3][1] = &vvc_scale_sub_residual_8_4_10_sse;
     rcn_func->ict.ict[4][1] = &vvc_scale_sub_residual_16_2_10_sse;
     rcn_func->ict.ict[5][1] = &vvc_scale_sub_residual_32_1_10_sse;

     rcn_func->ict.ict[3][2] = &vvc_scale_sub_half_residual_8_4_10_sse;
     rcn_func->ict.ict[4][2] = &vvc_scale_sub_half_residual_16_2_10_sse;
     rcn_func->ict.ict[5][2] = &vvc_scale_sub_half_residual_32_1_10_sse;
     break;
   case 2:
     rcn_func->ict.ict[3][0] = &vvc_add_residual_8_4_10_sse;
     rcn_func->ict.ict[4][0] = &vvc_add_residual_16_2_10_sse;
     rcn_func->ict.ict[5][0] = &vvc_add_residual_32_1_10_sse;

     rcn_func->ict.ict[3][1] = &vvc_sub_residual_8_4_10_sse;
     rcn_func->ict.ict[4][1] = &vvc_sub_residual_16_2_10_sse;
     rcn_func->ict.ict[5][1] = &vvc_sub_residual_32_1_10_sse;

     rcn_func->ict.ict[3][2] = &vvc_sub_half_residual_8_4_10_sse;
     rcn_func->ict.ict[4][2] = &vvc_sub_half_residual_16_2_10_sse;
     rcn_func->ict.ict[5][2] = &vvc_sub_half_residual_32_1_10_sse;
     break;
   case 1:
     rcn_func->ict.ict[3][0] = &vvc_scale_add_residual_8_4_10_sse;
     rcn_func->ict.ict[4][0] = &vvc_scale_add_residual_16_2_10_sse;
     rcn_func->ict.ict[5][0] = &vvc_scale_add_residual_32_1_10_sse;

     rcn_func->ict.ict[3][1] = &vvc_scale_add_residual_8_4_10_sse;
     rcn_func->ict.ict[4][1] = &vvc_scale_add_residual_16_2_10_sse;
     rcn_func->ict.ict[5][1] = &vvc_scale_add_residual_32_1_10_sse;

     rcn_func->ict.ict[3][2] = &vvc_scale_add_half_residual_8_4_10_sse;
     rcn_func->ict.ict[4][2] = &vvc_scale_add_half_residual_16_2_10_sse;
     rcn_func->ict.ict[5][2] = &vvc_scale_add_half_residual_32_1_10_sse;
     break;
   default:

     rcn_func->ict.ict[3][0] = &vvc_add_residual_8_4_10_sse;
     rcn_func->ict.ict[4][0] = &vvc_add_residual_16_2_10_sse;
     rcn_func->ict.ict[5][0] = &vvc_add_residual_32_1_10_sse;

     rcn_func->ict.ict[3][1] = &vvc_add_residual_8_4_10_sse;
     rcn_func->ict.ict[4][1] = &vvc_add_residual_16_2_10_sse;
     rcn_func->ict.ict[5][1] = &vvc_add_residual_32_1_10_sse;

     rcn_func->ict.ict[3][2] = &vvc_add_half_residual_8_4_10_sse;
     rcn_func->ict.ict[4][2] = &vvc_add_half_residual_16_2_10_sse;
     rcn_func->ict.ict[5][2] = &vvc_add_half_residual_32_1_10_sse;
     break;
 }*/
}