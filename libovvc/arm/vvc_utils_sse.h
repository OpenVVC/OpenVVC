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

#ifndef VVC_UTILS_SSE_H
#define VVC_UTILS_SSE_H

#include <stdint.h>
#include <stddef.h>
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

typedef int16_t TCoeff;
typedef int16_t TMatrixCoeff;

void inverse_sse2_B4(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B8(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B16(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B32(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);

void afficherVecteur4SSE128(__m128i cible);
void afficherVecteur8SSE128(__m128i cible);
void afficherTableau2D(int v, int h, TCoeff * m);

#endif//VVC_UTILS_SSE_H
