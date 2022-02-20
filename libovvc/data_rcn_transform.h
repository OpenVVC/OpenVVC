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

#ifndef RCN_TRANSFORM_DATA_H
#define RCN_TRANSFORM_DATA_H

#include <stdint.h>

extern const int16_t DCT_II_2[4];
extern const int16_t DCT_II_4[16];
extern const int16_t DCT_II_8[64];
extern const int16_t DCT_II_16[256];
extern const int16_t DCT_II_32[1024];
extern const int16_t DCT_II_64_OT[16*32];
extern const int16_t DCT_II_64_EOT[8*16];
extern const int16_t DCT_II_64_EEOT[4*8];
extern const int16_t DCT_II_64_EEEOT[2*4];
extern const int16_t DCT_II_64_EEEEO[2];
extern const int16_t DCT_II_64_EEEEE[2];

extern const int16_t DCT_VIII_4[16];
extern const int16_t DCT_VIII_8[64];
extern const int16_t DCT_VIII_16[256];
extern const int16_t DCT_VIII_32[1024];

extern const int16_t DST_VII_4[16];
extern const int16_t DST_VII_8[64];
extern const int16_t DST_VII_16[256];
extern const int16_t DST_VII_32[1024];

extern const int8_t* const lfnst[2][4][2];

#endif // RCN_TRANSFORM_DATA_H
