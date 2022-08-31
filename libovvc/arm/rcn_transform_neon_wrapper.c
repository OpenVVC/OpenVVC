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
// store functions for DC
void ov_store_8_neon(int16_t *dst, int value);
void ov_store_16_neon(int16_t *dst, int value);
void ov_store_32_neon(int16_t *dst, int value);
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
// 16-pt dct-8/dst-7
void ov_idct_x_16_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_16_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_16_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
// 16-pt DCT-II depth two butterfly
void ov_idct_ii_16_2_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2);
void ov_idct_ii_16_2_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2);
void ov_idct_ii_16_2_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2);
// 32-pt DCT-II depth two butterfly
void ov_idct_ii_32_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2, const int16_t *MatO32);
void ov_idct_ii_32_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2, const int16_t *MatO32);
void ov_idct_ii_32_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2, const int16_t *MatO32);
void ov_idct_ii_32_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16s1, const int16_t *MatO16s2, const int16_t *MatO32);
// 32-pt dct-8/dst-7
void ov_idct_x_32_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_32_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
void ov_idct_x_32_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, int shift_add, const int16_t *Mat);
// 64-pt dct-2 depth three butterfly
void ov_idct_ii_64_4_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16, const int16_t *MatO32, const int16_t *MatO64);
void ov_idct_ii_64_8_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16, const int16_t *MatO32, const int16_t *MatO64);
void ov_idct_ii_64_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16, const int16_t *MatO32, const int16_t *MatO64);
void ov_idct_ii_64_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride, int shift, const int16_t *MatEO, const int16_t *MatO16, const int16_t *MatO32, const int16_t *MatO64);


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
      ov_idct_x_8_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }else{
      ov_idct_x_8_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
    }
  }
  if (num_lines & 0x2){
    if(line_brk <=4) {
      ov_idct_x_8_2_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VII_8_top, DCT_VII_8_bot);
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
      ov_idct_x_8_4_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }else{
      ov_idct_x_8_4_8_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
    }
  }
  if (num_lines & 0x2){
    if(line_brk <=4) {
      ov_idct_x_8_2_4_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_8_top, DCT_VIII_8_bot);
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

void vvc_inverse_dst_vii_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                            int num_lines, int line_brk, int shift){

  for (int j = 0; j < num_lines; j++) {
    if(line_brk > 8 ) {
      ov_idct_x_16_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_16);
    }else if(line_brk > 4 ) {
      ov_idct_x_16_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_16);
    }else{
      ov_idct_x_16_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_16);
    }
    src += 1;
    dst += 16;
  }
}

void vvc_inverse_dct_viii_16_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                             int num_lines, int line_brk, int shift){

  for (int j = 0; j < num_lines; j++) {
    if(line_brk > 8 ) {
      ov_idct_x_16_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_16);
    }else if(line_brk > 4 ) {
      ov_idct_x_16_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_16);
    }else{
      ov_idct_x_16_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_16);
    }
    src += 1;
    dst += 16;
  }
}

void vvc_inverse_dct_ii_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                           int num_lines, int line_brk, int shift)
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
  const int16_t DCT_II_32_odd[16*16] = {
    90,  90,  88,  85,   82,  78,  73,  67,  61,  54,  46,  38,  31,  22,  13,   4,
    90,  82,  67,  46,   22,  -4, -31, -54, -73, -85, -90, -88, -78, -61, -38, -13,
    88,  67,  31,  -13, -54, -82, -90, -78, -46,  -4,  38,  73,  90,  85,  61,  22,
    85,  46, -13,  -67, -90, -73, -22,  38,  82,  88,  54,  -4, -61, -90, -78, -31,
    82,  22, -54,  -90, -61,  13,  78,  85,  31, -46, -90, -67,   4,  73,  88,  38,
    78,  -4, -82,  -73,  13,  85,  67, -22, -88, -61,  31,  90,  54, -38, -90, -46,
    73, -31, -90,  -22,  78,  67, -38, -90, -13,  82,  61, -46, -88,  -4,  85,  54,
    67, -54, -78,   38,  85, -22, -90,   4,  90,  13, -88, -31,  82,  46, -73, -61,
    61, -73, -46,   82,  31, -88, -13,  90,  -4, -90,  22,  85, -38, -78,  54,  67,
    54, -85,  -4,   88, -46, -61,  82,  13, -90,  38,  67, -78, -22,  90, -31, -73,
    46, -90,  38,   54, -90,  31,  61, -88,  22,  67, -85,  13,  73, -82,   4,  78,
    38, -88,  73,   -4, -67,  90, -46, -31,  85, -78,  13,  61, -90,  54,  22, -82,
    31, -78,  90,  -61,   4,  54, -88,  82, -38, -22,  73, -90,  67, -13, -46,  85,
    22, -61,  85,  -90,  73, -38,  -4,  46, -78,  90, -82,  54, -13, -31,  67, -88,
    13, -38,  61,  -78,  88, -90,  85, -73,  54, -31,   4,  22, -46,  67, -82,  90,
    4,  -13,  22,  -31,  38, -46,  54, -61,  67, -73,  78, -82,  85, -88,  90, -90
  };
  if(num_lines > 1){

    for (int j = 0; j < num_lines ; j++) {

      if (line_brk > 16)
        ov_idct_ii_32_32_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2, DCT_II_32_odd);
      else if (line_brk > 8)
        ov_idct_ii_32_16_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2, DCT_II_32_odd);
      else if (line_brk > 4)
        ov_idct_ii_32_8_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2, DCT_II_32_odd);
      else
        ov_idct_ii_32_4_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd1, DCT_II_16_odd2, DCT_II_32_odd);

      src += 1;
      dst += 32;
    }
  }

  if (num_lines & 0x1){
    vvc_inverse_dct_ii_32(src, dst, src_stride, num_lines & 0x1, line_brk, shift);
  }
}

void vvc_inverse_dst_vii_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                            int num_lines, int line_brk, int shift){

  for (int j = 0; j < num_lines; j++) {
    /*if(line_brk > 8 ) {
      ov_idct_x_32_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_32);
    }else if(line_brk > 4)  {
      ov_idct_x_32_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_32);
    }else{ */
    ov_idct_x_32_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DST_VII_32);

    //}
    src += 1;
    dst += 32;
  }



}

void vvc_inverse_dct_viii_32_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                             int num_lines, int line_brk, int shift){

  for (int j = 0; j < num_lines; j++) {
    /*if(line_brk > 8 ) {
      ov_idct_x_32_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_32);
    }else if(line_brk > 4)  {
      ov_idct_x_32_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_32);
    }else{ */
    ov_idct_x_32_16_neon(src, dst, src_stride<<1, -shift,(1<<shift-1), DCT_VIII_32);

    //}
    src += 1;
    dst += 32;
  }

}

void vvc_inverse_dct_ii_64_neon(const int16_t *src, int16_t *dst, ptrdiff_t src_stride,
                           int num_lines, int line_brk, int shift)
{
  const int16_t DCT_II_8_eo[4*8] = {
    64, 64,  64,  64,
    83, 36,  -36, -83,
    64, -64, -64, 64, // zeroed
    36, -83, 83,  -36, //
    89,  75,  50,  18,
    75, -18, -89, -50,
    50, -89,  18,  75, // zeroed
    18, -50,  75, -89 //
  };
  const int16_t DCT_II_16_odd[8*4] = {
    90,  87,  80,  70,  57,  43,  25,   9,
    87,  57,   9, -43, -80, -90, -70, -25,
    80,   9, -70, -87, -25,  57,  90,  43,
    70, -43, -87,   9,  90,  25, -80, -57
  };

  const int16_t DCT_II_32_odd[16*8] = {
    90,  90,  88,  85,   82,  78,  73,  67,  61,  54,  46,  38,  31,  22,  13,   4,
    90,  82,  67,  46,   22,  -4, -31, -54, -73, -85, -90, -88, -78, -61, -38, -13,
    88,  67,  31,  -13, -54, -82, -90, -78, -46,  -4,  38,  73,  90,  85,  61,  22,
    85,  46, -13,  -67, -90, -73, -22,  38,  82,  88,  54,  -4, -61, -90, -78, -31,
    82,  22, -54,  -90, -61,  13,  78,  85,  31, -46, -90, -67,   4,  73,  88,  38,
    78,  -4, -82,  -73,  13,  85,  67, -22, -88, -61,  31,  90,  54, -38, -90, -46,
    73, -31, -90,  -22,  78,  67, -38, -90, -13,  82,  61, -46, -88,  -4,  85,  54,
    67, -54, -78,   38,  85, -22, -90,   4,  90,  13, -88, -31,  82,  46, -73, -61
  };

  const int16_t DCT_II_64_OT_t[16 * 32] =  {
    91, 90, 90, 90, 88, 87, 86, 84, 83, 81, 79, 77, 73, 71, 69, 65, 62, 59, 56, 52, 48, 44, 41, 37, 33, 28, 24, 20, 15, 11, 7, 2,
    90, 88, 84, 79, 71, 62, 52, 41, 28, 15, 2, -11, -24, -37, -48, -59, -69, -77, -83, -87, -90, -91, -90, -86, -81, -73, -65, -56, -44, -33, -20, -7,
    90, 84, 73, 59, 41, 20, -2, -24, -44, -62, -77, -86, -90, -90, -83, -71, -56, -37, -15, 7, 28, 48, 65, 79, 87, 91, 88, 81, 69, 52, 33, 11,
    90, 79, 59, 33, 2, -28, -56, -77, -88, -90, -81, -62, -37, -7, 24, 52, 73, 87, 90, 83, 65, 41, 11, -20, -48, -71, -86, -91, -84, -69, -44, -15,
    88, 71, 41, 2, -37, -69, -87, -90, -73, -44, -7, 33, 65, 86, 90, 77, 48, 11, -28, -62, -84, -90, -79, -52, -15, 24, 59, 83, 91, 81, 56, 20,
    87, 62, 20, -28, -69, -90, -84, -56, -11, 37, 73, 90, 81, 48, 2, -44, -79, -91, -77, -41, 7, 52, 83, 90, 71, 33, -15, -59, -86, -88, -65, -24,
    86, 52, -2, -56, -87, -84, -48, 7, 59, 88, 83, 44, -11, -62, -90, -81, -41, 15, 65, 90, 79, 37, -20, -69, -90, -77, -33, 24, 71, 91, 73, 28,
    84, 41, -24, -77, -90, -56, 7, 65, 91, 69, 11, -52, -88, -79, -28, 37, 83, 86, 44, -20, -73, -90, -59, 2, 62, 90, 71, 15, -48, -87, -81, -33,
    83, 28, -44, -88, -73, -11, 59, 91, 62, -7, -71, -90, -48, 24, 81, 84, 33, -41, -87, -77, -15, 56, 90, 65, -2, -69, -90, -52, 20, 79, 86, 37,
    81, 15, -62, -90, -44, 37, 88, 69, -7, -77, -84, -24, 56, 91, 52, -28, -86, -73, -2, 71, 87, 33, -48, -90, -59, 20, 83, 79, 11, -65, -90, -41,
    79, 2, -77, -81, -7, 73, 83, 11, -71, -84, -15, 69, 86, 20, -65, -87, -24, 62, 88, 28, -59, -90, -33, 56, 90, 37, -52, -90, -41, 48, 91, 44,
    77, -11, -86, -62, 33, 90, 44, -52, -90, -24, 69, 83, 2, -81, -71, 20, 88, 56, -41, -91, -37, 59, 87, 15, -73, -79, 7, 84, 65, -28, -90, -48,
    73, -24, -90, -37, 65, 81, -11, -88, -48, 56, 86, 2, -84, -59, 44, 90, 15, -79, -69, 33, 91, 28, -71, -77, 20, 90, 41, -62, -83, 7, 87, 52,
    71, -37, -90, -7, 86, 48, -62, -79, 24, 91, 20, -81, -59, 52, 84, -11, -90, -33, 73, 69, -41, -88, -2, 87, 44, -65, -77, 28, 90, 15, -83, -56,
    69, -48, -83, 24, 90, 2, -90, -28, 81, 52, -65, -71, 44, 84, -20, -90, -7, 88, 33, -79, -56, 62, 73, -41, -86, 15, 91, 11, -87, -37, 77, 59,
    65, -59, -71, 52, 77, -44, -81, 37, 84, -28, -87, 20, 90, -11, -90, 2, 91, 7, -90, -15, 88, 24, -86, -33, 83, 41, -79, -48, 73, 56, -69, -62
  };

  if(num_lines > 3){

    for (int j = 0; j < num_lines ; j++) {

      if (line_brk > 16)
        ov_idct_ii_64_32_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd, DCT_II_32_odd, DCT_II_64_OT_t);
      else if (line_brk > 8)
        ov_idct_ii_64_16_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd, DCT_II_32_odd, DCT_II_64_OT_t);
      else if (line_brk > 4)
        ov_idct_ii_64_8_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd, DCT_II_32_odd, DCT_II_64_OT_t);
      else
        ov_idct_ii_64_4_neon(src, dst, src_stride<<1, shift, DCT_II_8_eo, DCT_II_16_odd, DCT_II_32_odd, DCT_II_64_OT_t);

      src += 1;
      dst += 64;
    }
  }

  if (num_lines & 0x3)
    vvc_inverse_dct_ii_64(src, dst, src_stride, num_lines & 0x3, line_brk, shift);

}

void vvc_inverse_dct_ii_dc_neon(int16_t *const dst, int log2_tb_w, int log2_tb_h,
                           int dc_val)
{

  int i, j;
  int tb_w = 1 << log2_tb_w;
  int tb_h = 1 << log2_tb_h;
  int clip_min = -(1 << 15);
  int clip_max = (1 << 15)-1;
  int16_t * _dst = (int16_t *)dst;
  int value = (((dc_val + 1) >> 1) + 8) >> 4;
  value = ov_clip(value, clip_min, clip_max);

  switch (log2_tb_w){
    case 3:{
      for (i = 0; i < tb_h; ++i){
        // assembly function to store 8 values
        ov_store_8_neon(_dst, value);
        _dst += tb_w;
      }
      break;
    }
#if 1
    case 4:{
      for (i = 0; i < tb_h; ++i){
        // assembly function to store data 16 values (two 128-bit reg)
        ov_store_16_neon(_dst, value);
        _dst += tb_w;
      }
      break;
    }
    case 5:{
      for (i = 0; i < tb_h; ++i){
        // assembly function to store data 32 values (four 128-bit reg)
        ov_store_32_neon(_dst, value);
        _dst += tb_w;
      }
      break;
    }
#endif
    default:
      vvc_inverse_dct_ii_dc(dst, log2_tb_w, log2_tb_h, dc_val);
  }
}

void rcn_init_tr_functions_neon(struct RCNFunctions *const rcn_funcs){
  rcn_funcs->tr.func[DST_VII][2] = &vvc_inverse_dst_vii_4_neon;
  rcn_funcs->tr.func[DST_VII][3] = &vvc_inverse_dst_vii_8_neon;
  rcn_funcs->tr.func[DST_VII][4] = &vvc_inverse_dst_vii_16_neon;
  rcn_funcs->tr.func[DST_VII][5] = &vvc_inverse_dst_vii_32_neon;

  rcn_funcs->tr.func[DCT_VIII][2] = &vvc_inverse_dct_viii_4_neon;
  rcn_funcs->tr.func[DCT_VIII][3] = &vvc_inverse_dct_viii_8_neon;
  rcn_funcs->tr.func[DCT_VIII][4] = &vvc_inverse_dct_viii_16_neon;
  rcn_funcs->tr.func[DCT_VIII][5] = &vvc_inverse_dct_viii_32_neon;

  rcn_funcs->tr.func[DCT_II][2] = &vvc_inverse_dct_ii_4_neon;
  rcn_funcs->tr.func[DCT_II][3] = &vvc_inverse_dct_ii_8_neon;
  rcn_funcs->tr.func[DCT_II][4] = &vvc_inverse_dct_ii_16_neon;
  rcn_funcs->tr.func[DCT_II][5] = &vvc_inverse_dct_ii_32_neon;
  rcn_funcs->tr.func[DCT_II][6] = &vvc_inverse_dct_ii_64_neon;

  rcn_funcs->tr.dc = &vvc_inverse_dct_ii_dc_neon;
}
