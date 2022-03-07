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
#include "ovutils.h"
#include "rcn_structures.h"
#include "stdint.h"

#define DECALARE_DC_PDPC(w, h) void ov_dc_pdpc_w ## w ## _h ## h ##_neon(const uint16_t *const src_above,               \
                                                                         const uint16_t *const src_left,                \
                                                                         uint16_t *const dst, ptrdiff_t dst_stride,     \
                                                                         int log2_pb_w, int log2_pb_h);

DECALARE_DC_PDPC( 4,  4)
DECALARE_DC_PDPC( 4,  8)
DECALARE_DC_PDPC( 4, 16)
DECALARE_DC_PDPC( 4, 32)
DECALARE_DC_PDPC( 4, 64)

DECALARE_DC_PDPC( 8,  4)
DECALARE_DC_PDPC( 8,  8)
DECALARE_DC_PDPC( 8, 16)
DECALARE_DC_PDPC( 8, 32)
DECALARE_DC_PDPC( 8, 64)

DECALARE_DC_PDPC(16,  4)
DECALARE_DC_PDPC(16,  8)
DECALARE_DC_PDPC(16, 16)
DECALARE_DC_PDPC(16, 32)
DECALARE_DC_PDPC(16, 64)

DECALARE_DC_PDPC(32,  4)
DECALARE_DC_PDPC(32,  8)
DECALARE_DC_PDPC(32, 16)
DECALARE_DC_PDPC(32, 32)
DECALARE_DC_PDPC(32, 64)

DECALARE_DC_PDPC(64,  4)
DECALARE_DC_PDPC(64,  8)
DECALARE_DC_PDPC(64, 16)
DECALARE_DC_PDPC(64, 32)
DECALARE_DC_PDPC(64, 64)

void rcn_init_dc_planar_functions_neon(struct RCNFunctions *const rcn_funcs){
    //DC
    rcn_funcs->dc.pdpc[0][0] = &ov_dc_pdpc_w4_h4_neon;
    rcn_funcs->dc.pdpc[0][1] = &ov_dc_pdpc_w4_h8_neon;
    rcn_funcs->dc.pdpc[0][2] = &ov_dc_pdpc_w4_h16_neon;
    rcn_funcs->dc.pdpc[0][3] = &ov_dc_pdpc_w4_h32_neon;
    rcn_funcs->dc.pdpc[0][4] = &ov_dc_pdpc_w4_h64_neon;

    rcn_funcs->dc.pdpc[1][0] = &ov_dc_pdpc_w8_h4_neon;
    rcn_funcs->dc.pdpc[1][1] = &ov_dc_pdpc_w8_h8_neon;
    rcn_funcs->dc.pdpc[1][2] = &ov_dc_pdpc_w8_h16_neon;
    rcn_funcs->dc.pdpc[1][3] = &ov_dc_pdpc_w8_h32_neon;
    rcn_funcs->dc.pdpc[1][4] = &ov_dc_pdpc_w8_h64_neon;

    rcn_funcs->dc.pdpc[2][0] = &ov_dc_pdpc_w16_h4_neon;
    rcn_funcs->dc.pdpc[2][1] = &ov_dc_pdpc_w16_h8_neon;
    rcn_funcs->dc.pdpc[2][2] = &ov_dc_pdpc_w16_h16_neon;
    rcn_funcs->dc.pdpc[2][3] = &ov_dc_pdpc_w16_h32_neon;
    rcn_funcs->dc.pdpc[2][4] = &ov_dc_pdpc_w16_h64_neon;

    rcn_funcs->dc.pdpc[3][0] = &ov_dc_pdpc_w32_h4_neon;
    rcn_funcs->dc.pdpc[3][1] = &ov_dc_pdpc_w32_h8_neon;
    rcn_funcs->dc.pdpc[3][2] = &ov_dc_pdpc_w32_h16_neon;
    rcn_funcs->dc.pdpc[3][3] = &ov_dc_pdpc_w32_h32_neon;
    rcn_funcs->dc.pdpc[3][4] = &ov_dc_pdpc_w32_h64_neon;

    rcn_funcs->dc.pdpc[4][0] = &ov_dc_pdpc_w64_h4_neon;
    rcn_funcs->dc.pdpc[4][1] = &ov_dc_pdpc_w64_h8_neon;
    rcn_funcs->dc.pdpc[4][2] = &ov_dc_pdpc_w64_h16_neon;
    rcn_funcs->dc.pdpc[4][3] = &ov_dc_pdpc_w64_h32_neon;
    rcn_funcs->dc.pdpc[4][4] = &ov_dc_pdpc_w64_h64_neon;
}
