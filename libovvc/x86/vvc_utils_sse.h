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
#include <emmintrin.h>

typedef int16_t TCoeff;
typedef int16_t TMatrixCoeff;

void inverse_sse2_B4(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B8(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B16(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);
void inverse_sse2_B32(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT);

void afficherVecteur4SSE128(__m128i cible);
void afficherVecteur8SSE128(__m128i cible);
void afficherTableau2D(int v, int h, TCoeff * m);

#define NORETURN_INSTRUCTION(instruction, arg1, ...) instruction(arg1, ## __VA_ARGS__);

#define RETURN_INSTRUCTION(dst, instruction, arg1, ...) dst = NORETURN_INSTRUCTION(instruction, arg1, ## __VA_ARGS__)

#define LOAD1(dst, instruction, src, jump) \
        RETURN_INSTRUCTION(dst ## 1, instruction, (__m128i *) (src))

#define LOAD2(dst, instruction, src, jump) \
        RETURN_INSTRUCTION(dst ## 1, instruction, (__m128i *) (src))\
        RETURN_INSTRUCTION(dst ## 2, instruction, (__m128i *) ((src) + 1 * jump))

#define LOAD4(dst, instruction, src, jump) \
        RETURN_INSTRUCTION(dst ## 1, instruction, (__m128i *) (src))\
        RETURN_INSTRUCTION(dst ## 2, instruction, (__m128i *) ((src) + 1 * jump))\
        RETURN_INSTRUCTION(dst ## 3, instruction, (__m128i *) ((src) + 2 * jump))\
        RETURN_INSTRUCTION(dst ## 4, instruction, (__m128i *) ((src) + 3 * jump))

#define LOAD8(dst, instruction, src, jump) \
        RETURN_INSTRUCTION(dst ## 1, instruction, (__m128i *) (src))\
        RETURN_INSTRUCTION(dst ## 2, instruction, (__m128i *) ((src) + 1 * jump))\
        RETURN_INSTRUCTION(dst ## 3, instruction, (__m128i *) ((src) + 2 * jump))\
        RETURN_INSTRUCTION(dst ## 4, instruction, (__m128i *) ((src) + 3 * jump))\
        RETURN_INSTRUCTION(dst ## 5, instruction, (__m128i *) ((src) + 4 * jump))\
        RETURN_INSTRUCTION(dst ## 6, instruction, (__m128i *) ((src) + 5 * jump))\
        RETURN_INSTRUCTION(dst ## 7, instruction, (__m128i *) ((src) + 6 * jump))\
        RETURN_INSTRUCTION(dst ## 8, instruction, (__m128i *) ((src) + 7 * jump))

#define STORE1(instruction, dst, src, jump) \
        NORETURN_INSTRUCTION(instruction, (__m128i *) (dst), src ## 1)

#define STORE2(instruction, dst, src, jump) \
        NORETURN_INSTRUCTION(instruction, (__m128i *) (dst), src ## 1)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 1 * jump), src ## 2)

#define STORE4(instruction, dst, src, jump) \
        NORETURN_INSTRUCTION(instruction, (__m128i *) (dst), src ## 1)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 1 * jump), src ## 2)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 2 * jump), src ## 3)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 3 * jump), src ## 4)

#define STORE8(instruction, dst, src, jump) \
        NORETURN_INSTRUCTION(instruction, (__m128i *) (dst), src ## 1)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 1 * jump), src ## 2)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 2 * jump), src ## 3)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 3 * jump), src ## 4)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 4 * jump), src ## 5)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 5 * jump), src ## 6)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 6 * jump), src ## 7)\
        NORETURN_INSTRUCTION(instruction, (__m128i *) ((dst) + 7 * jump), src ## 8)

#define U111(n, dst, instruction, arg1, arg2) RETURN_INSTRUCTION(dst ## n, instruction, arg1 ## n, arg2 ## n)
#define U101(n, dst, instruction, arg1, arg2) RETURN_INSTRUCTION(dst ## n, instruction, arg1, arg2 ## n)
#define U110(n, dst, instruction, arg1, arg2) RETURN_INSTRUCTION(dst ## n, instruction, arg1 ## n, arg2)

#define UNROLL1(type, dst, instruction, arg1, arg2)     type(1, dst, instruction, arg1, arg2)
#define UNROLL2(type, dst, instruction, arg1, arg2)     type(1, dst, instruction, arg1, arg2)\
                                                        type(2, dst, instruction, arg1, arg2)

#define UNROLL4(type, dst, instruction, arg1, arg2)     type(1, dst, instruction, arg1, arg2)\
                                                        type(2, dst, instruction, arg1, arg2)\
                                                        type(3, dst, instruction, arg1, arg2)\
                                                        type(4, dst, instruction, arg1, arg2)

#define UNROLL8(type, dst, instruction, arg1, arg2)     type(1, dst, instruction, arg1, arg2)\
                                                        type(2, dst, instruction, arg1, arg2)\
                                                        type(3, dst, instruction, arg1, arg2)\
                                                        type(4, dst, instruction, arg1, arg2)\
                                                        type(5, dst, instruction, arg1, arg2)\
                                                        type(6, dst, instruction, arg1, arg2)\
                                                        type(7, dst, instruction, arg1, arg2)\
                                                        type(8, dst, instruction, arg1, arg2)

#endif//VVC_UTILS_SSE_H
