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

#ifndef RCN_INTRA_ANGULAR_H
#define RCN_INTRA_ANGULAR_H

#include <stddef.h>
#include <stdint.h>
#include "bitdepth.h"

struct IntraAngularFunctions {
    void (*pure)(const OVSample* ref_lft, OVSample* dst,
                 ptrdiff_t dst_stride, int8_t log2_pb_w,
                 int8_t log2_pb_h);

    void (*diagonal)(const OVSample* ref_lft, OVSample* dst,
                     ptrdiff_t dst_stride, int8_t log2_pb_w,
                     int8_t log2_pb_h);

    void (*angular)(const OVSample* ref_lft, OVSample* dst,
                    ptrdiff_t dst_stride, int8_t log2_pb_w,
                    int8_t log2_pb_h, int angle_val);

    void (*pure_pdpc)(const OVSample* ref_abv, const OVSample* ref_lft,
                      OVSample* const dst, ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h);

    void (*diagonal_pdpc)(const OVSample* ref_abv, const OVSample* ref_lft,
                         OVSample* const dst, ptrdiff_t dst_stride,
                         int8_t log2_pb_w, int8_t log2_pb_h);

    void (*angular_pdpc)(const OVSample* ref_abv, const OVSample* ref_lft,
                         OVSample* const dst, ptrdiff_t dst_stride,
                         int8_t log2_pb_w, int8_t log2_pb_h, int mode_idx);
};

struct IntraMRLFunctions
{
    void (*angular_h)(const OVSample* const ref_lft, OVSample* const dst,
                      ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h,
                      int angle_val, uint8_t multi_ref_idx);

    void (*angular_v)(const OVSample* const ref_lft, OVSample* const dst,
                      ptrdiff_t dst_stride,
                      int8_t log2_pb_w, int8_t log2_pb_h,
                      int angle_val, uint8_t multi_ref_idx);
};

#endif // RCN_INTRA_ANGULAR_H
