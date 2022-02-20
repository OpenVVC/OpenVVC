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

#ifndef RCN_TRANSFORM_H
#define RCN_TRANSFORM_H

#include <stddef.h>
#include <stdint.h>
#include "rcn_structures.h"


void rcn_init_tr_functions(struct RCNFunctions *const rcn_funcs);

void
vvc_inverse_dct_ii_2(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift);
void
vvc_inverse_dct_ii_4(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift);
void
vvc_inverse_dct_ii_8(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                     int num_lines, int num_columns, int shift);
void
vvc_inverse_dct_ii_16(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);
void
vvc_inverse_dct_ii_32(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);
void
vvc_inverse_dct_ii_64(const int16_t* src, int16_t* dst, ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);

void
vvc_inverse_dct_ii_dc(int16_t* const dst, int log2_tb_w, int log2_tb_h,
                      int dc_val);

#endif // RCN_TRANSFORM_H
