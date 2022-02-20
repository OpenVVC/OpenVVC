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

#ifndef DATA_SCAN_LUT_H
#define DATA_SCAN_LUT_H

#include <stdint.h>

extern const uint8_t ff_vvc_diag_scan_2x8_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_4x4_num_cg [16];
extern const uint8_t ff_vvc_diag_scan_8x2_num_cg [16];

extern const uint8_t ff_vvc_inv_diag_scan_2x8 [16];
extern const uint8_t ff_vvc_inv_diag_scan_4x4 [16];
extern const uint8_t ff_vvc_inv_diag_scan_8x2 [16];

extern const uint8_t *const ff_vvc_idx_2_num[6][6];

extern const uint8_t *const ff_vvc_scan_x[6][6];
extern const uint8_t *const ff_vvc_scan_y[6][6];

#endif
