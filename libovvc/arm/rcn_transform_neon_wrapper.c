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

#include <stdio.h>
#include <stddef.h>

#include "ovutils.h"
#include "data_rcn_transform.h"
#include "rcn_transform.h"

//void ov_idct_x_4_1_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_4_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_4_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_4_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);


void vvc_inverse_dct_ii_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                         int num_lines, int line_brk, int shift)
{
  //ov_idct_ii_2_2_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_2);
  //ov_idct_ii_2_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_2);
  //ov_idct_ii_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_2);
}

void vvc_inverse_dct_ii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                         int num_lines, int line_brk, int shift)
{

  for (int j = 0; j < num_lines / 8; j++) {

      ov_idct_x_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4);
      src += 8;
      dst += 32;
  }

  if (!(num_lines & 0x7)) return;

  if (num_lines & 0x4){
    ov_idct_x_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4);
  }

  if (num_lines & 0x2){
    ov_idct_x_4_2_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4);
  }

  if (num_lines & 0x1){
    vvc_inverse_dct_ii_4(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
  }
}

void vvc_inverse_dst_vii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                         int num_lines, int line_brk, int shift)
{

  for (int j = 0; j < num_lines / 8; j++) {

      ov_idct_x_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_4);
      src += 8;
      dst += 32;
  }

  if (!(num_lines & 0x7)) return;

  if (num_lines & 0x4){
    ov_idct_x_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_4);
  }

  if (num_lines & 0x2){
    ov_idct_x_4_2_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_4);
  }

  if (num_lines & 0x1){
    //ov_idct_x_4_1_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_4);
  }
}

void vvc_inverse_dct_viii_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                         int num_lines, int line_brk, int shift)
{
  for (int j = 0; j < num_lines / 8; j++) {

      ov_idct_x_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_4);
      src += 8;
      dst += 32;
  }

  if (!(num_lines & 0x7)) return;

  if (num_lines & 0x4){
    ov_idct_x_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_4);
  }

  if (num_lines & 0x2){
    ov_idct_x_4_2_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_4);
  }

  if (num_lines & 0x1){
    //ov_idct_x_4_1_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_4);
  }
}


void rcn_init_tr_functions_neon(struct RCNFunctions *const rcn_funcs){
  rcn_funcs->tr.func[DST_VII][2] = &vvc_inverse_dst_vii_4_neon;
  /*rcn_funcs->tr.func[DST_VII][3] = &vvc_inverse_dst_vii_8_sse;
  rcn_funcs->tr.func[DST_VII][4] = &vvc_inverse_dst_vii_16_sse;
  rcn_funcs->tr.func[DST_VII][5] = &vvc_inverse_dst_vii_32_sse;
  */
  rcn_funcs->tr.func[DCT_VIII][2] = &vvc_inverse_dct_viii_4_neon;
  /*rcn_funcs->tr.func[DCT_VIII][3] = &vvc_inverse_dct_viii_8_sse;
  rcn_funcs->tr.func[DCT_VIII][4] = &vvc_inverse_dct_viii_16_sse;
  rcn_funcs->tr.func[DCT_VIII][5] = &vvc_inverse_dct_viii_32_sse;
   */
  //rcn_funcs->tr.func[DCT_II][1] = &vvc_inverse_dct_ii_2_neon;
  rcn_funcs->tr.func[DCT_II][2] = &vvc_inverse_dct_ii_4_neon;
  /*rcn_funcs->tr.func[DCT_II][3] = &vvc_inverse_dct_ii_8_sse;
  rcn_funcs->tr.func[DCT_II][4] = &vvc_inverse_dct_ii_16_sse;
  rcn_funcs->tr.func[DCT_II][5] = &vvc_inverse_dct_ii_32_sse;
  rcn_funcs->tr.func[DCT_II][6] = &vvc_inverse_dct_ii_64_sse;

  rcn_funcs->tr.dc = &vvc_inverse_dct_ii_dc_sse;*/
}
