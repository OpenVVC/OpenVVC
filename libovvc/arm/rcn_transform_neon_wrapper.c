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
// 4-pt DxT-II/VIII/VII
void ov_idct_x_4_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_4_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_4_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
// 8-pt DxT-II/VIII/VII
void ov_idct_x_8_2_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_x_8_2_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_x_8_4_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_x_8_4_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_x_8_8_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_x_8_8_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
// 8-pt DCT-II butterfly
void ov_idct_ii_8_2_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_ii_8_2_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_ii_8_4_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);
void ov_idct_ii_8_4_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatE, const int16_t *MatO);

void vvc_inverse_dct_ii_2_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{
  //ov_idct_x_4_2_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4);
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


void vvc_inverse_dct_ii_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                          int num_lines, int line_brk, int shift)
{

  const int16_t DCT_II_8_odd[4*4] = {
    89,  75,  50,  18,
    75, -18, -89, -50,
    50, -89,  18,  75,
    18, -50,  75, -89
  };
  // without butterfly
  const int16_t DCT_II_8_top[8*4] = {
    64,  64,   64,  64,  64,  64,  64,  64,
    89,  75,   50,  18, -18, -50, -75, -89,
    83,  36,  -36, -83, -83, -36,  36,  83,
    75, -18,  -89, -50,  50,  89,  18, -75 };
  const int16_t DCT_II_8_bot[8*4] = {
    64, -64,  -64,  64,  64, -64, -64,  64,
    50, -89,   18,  75, -75, -18,  89, -50,
    36, -83,   83, -36, -36,  83, -83,  36,
    18, -50,   75, -89,  89, -75,  50, -18 };

  for (int j = 0; j < num_lines / 8; j++) {
    if(line_brk <=4) {
      ov_idct_x_8_8_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_top, DCT_II_8_bot);
    }else{
      ov_idct_x_8_8_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_top, DCT_II_8_bot);
    }
    src += 8;
    dst += 64;
  }
  /* for (int j = 0; j < num_lines / 4; j++) {
     if(line_brk <=4) {
       ov_idct_ii_8_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4, DCT_II_8_odd);
     }else{
       ov_idct_ii_8_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4, DCT_II_8_odd);
     }
       src += 4;
       dst += 32;
   }*/
  if (num_lines & 0x4){
    if(line_brk <=4) {
      ov_idct_x_8_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_top, DCT_II_8_bot);
    }else{
      ov_idct_x_8_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_top, DCT_II_8_bot);
    }
  }
  if (num_lines & 0x2){
    if(line_brk <=4) {
      //ov_idct_ii_8_2_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4, DCT_II_8_odd);
      ov_idct_x_8_2_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_top, DCT_II_8_bot);
    }else{
      //ov_idct_ii_8_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_4, DCT_II_8_odd);
      ov_idct_x_8_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_top, DCT_II_8_bot);
    }
  }
  if (num_lines & 0x1){
    vvc_inverse_dct_ii_8(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
  }
}

void vvc_inverse_dst_vii_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                           int num_lines, int line_brk, int shift){
  const int16_t DCT_VII_8_top[8*4] = {
    17,  32,  46,  60,  71,  78,  85,  86,
    46,  78,  86,  71,  32, -17, -60, -85,
    71,  85,  32, -46, -86, -60,  17,  78,
    85,  46, -60, -78,  17,  86,  32, -71 };
  const int16_t DCT_VII_8_bot[8*4] = {
    86, -17, -85,  32,  78, -46, -71,  60,
    78, -71, -17,  85, -60, -32,  86, -46,
    60, -86,  71, -17, -46,  85, -78,  32,
    32, -60, 78,  -86, 85,  -71, 46,  -17};

  for (int j = 0; j < num_lines / 8; j++) {
    if(line_brk <=4) {
      ov_idct_x_8_8_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }else{
      ov_idct_x_8_8_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }
    src += 8;
    dst += 64;
  }
  if (num_lines & 0x4){
    if(line_brk <=4) {
      ov_idct_x_8_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }else{
      ov_idct_x_8_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }
  }
  if (num_lines & 0x2){
    if(line_brk <=4) {
      ov_idct_x_8_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }else{
      ov_idct_x_8_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }
  }
}

void vvc_inverse_dct_viii_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                            int num_lines, int line_brk, int shift){
  const int16_t DCT_VIII_8_top[8*4] = {
    86,  85,  78,  71,  60,  46,  32,  17,
    85,  60,  17, -32, -71, -86, -78, -46,
    78,  17, -60, -86, -46,  32,  85,  71,
    71, -32, -86, -17,  78,  60, -46, -85 };
  const int16_t DCT_VIII_8_bot[8*4] = {
    60, -71, -46,  78,  32, -85, -17,  86,
    46, -86,  32,  60, -85,  17,  71, -78,
    32, -78,  85, -46, -17,  71, -86,  60,
    17, -46,  71, -85,  86, -78,  60, -32 };

  for (int j = 0; j < num_lines / 8; j++) {
    if(line_brk <=4) {
      ov_idct_x_8_8_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }else{
      ov_idct_x_8_8_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }
    src += 8;
    dst += 64;
  }
  if (num_lines & 0x4){
    if(line_brk <=4) {
      ov_idct_x_8_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }else{
      ov_idct_x_8_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }
  }
  if (num_lines & 0x2){
    if(line_brk <=4) {
      ov_idct_x_8_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }else{
      ov_idct_x_8_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }
  }
}

void vvc_inverse_dct_ii_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                           int nb_lines, int line_brk, int shift)
{

  const int16_t DCT_II_8_eo[4*8] = {
    64, 64,  64,  64,
    83, 36,  -36, -83,
    64, -64, -64, 64,
    36, -83, 83,  -36,
    89,  75,  50,  18,
    75, -18, -89, -50,
    50, -89,  18,  75,
    18, -50,  75, -89
  };
  const int16_t DCT_II_16_odd1[8*4] = {
    90,  87,  80,  70,  57,  43,  25,   9,
    87,  57,   9, -43, -80, -90, -70, -25,
    80,   9, -70, -87, -25,  57,  90,  43,
    70, -43, -87,   9,  90,  25, -80, -57
  };
  const int16_t DCT_II_16_odd2[8*4] = {
    57, -80, -25,  90,  -9, -87,  43,  70,
    43, -90,  57,  25, -87,  70,   9, -80,
    25, -70,  90, -80,  43,   9, -57,  87,
    9 , -25,  43, -57,  70, -80,  87, -90
  };

  for (int j = 0; j < nb_lines / 2; j++) {
    if(line_brk <=4) {
      ov_idct_ii_16_2_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2);
    }else if(line_brk <= 8){
      ov_idct_ii_16_2_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2);
    }else {
      ov_idct_ii_16_2_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2);
    }

    src += 2;
    dst += 32;
  }

  if (nb_lines & 0x1) {
    vvc_inverse_dct_ii_16(src, dst, src_stride, nb_lines & 0x1, line_brk, shift);
  }
}


void rcn_init_tr_functions_neon(struct RCNFunctions *const rcn_funcs){
  rcn_funcs->tr.func[DST_VII][2] = &vvc_inverse_dst_vii_4_neon;
  //rcn_funcs->tr.func[DST_VII][3] = &vvc_inverse_dst_vii_8_neon;
  /*rcn_funcs->tr.func[DST_VII][4] = &vvc_inverse_dst_vii_16_sse;
  rcn_funcs->tr.func[DST_VII][5] = &vvc_inverse_dst_vii_32_sse;
  */
  rcn_funcs->tr.func[DCT_VIII][2] = &vvc_inverse_dct_viii_4_neon;
  //rcn_funcs->tr.func[DCT_VIII][3] = &vvc_inverse_dct_viii_8_neon;
  /*rcn_funcs->tr.func[DCT_VIII][4] = &vvc_inverse_dct_viii_16_sse;
  rcn_funcs->tr.func[DCT_VIII][5] = &vvc_inverse_dct_viii_32_sse;
   */
  //rcn_funcs->tr.func[DCT_II][1] = &vvc_inverse_dct_ii_2_neon;
  rcn_funcs->tr.func[DCT_II][2] = &vvc_inverse_dct_ii_4_neon;
  rcn_funcs->tr.func[DCT_II][3] = &vvc_inverse_dct_ii_8_neon;
  //rcn_funcs->tr.func[DCT_II][4] = &vvc_inverse_dct_ii_16_neon;
  /*rcn_funcs->tr.func[DCT_II][5] = &vvc_inverse_dct_ii_32_sse;
  rcn_funcs->tr.func[DCT_II][6] = &vvc_inverse_dct_ii_64_sse;

  rcn_funcs->tr.dc = &vvc_inverse_dct_ii_dc_sse;*/
}
