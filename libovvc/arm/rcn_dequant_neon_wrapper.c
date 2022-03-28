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

#include <stdint.h>
#include <string.h>
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

#define MAX_LOG2_TR_RANGE 15

void
ov_dequant_tb_4x4_neg_neon(int16_t *dst, const int16_t *src, int scale, int shift,
                   uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map);

void
ov_dequant_tb_4x4_neon(int16_t *dst, const int16_t *src, int scale, int shift,
                   uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map);

void
rcn_init_dequant_neon(struct RCNFunctions *rcn_funcs)
{
    rcn_funcs->tmp.dequant_tb_4x4 = &ov_dequant_tb_4x4_neon;
    rcn_funcs->tmp.dequant_tb_4x4_neg = &ov_dequant_tb_4x4_neg_neon;
}
