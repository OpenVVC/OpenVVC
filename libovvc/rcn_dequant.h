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

#ifndef RCN_DEQUANT_H
#define RCN_DEQUANT_H
#include <stdint.h>
#include "ovdefs.h"

struct IQScale
{
    int scale;
    int shift;
    void (*dequant_sb)(int16_t *const sb_coeffs, int scale, int shift);
};

struct ScalingMatrix;

struct VVCQPCTX;
void derive_dequant_ctx(OVCTUDec *const ctudec, const struct VVCQPCTX *const qp_ctx,
                        int cu_qp_delta);

//void dequant_sb_neg(int16_t *const sb_coeffs, int scale, int shift);

//void dequant_sb(int16_t *const sb_coeffs, int scale, int shift);

#endif
