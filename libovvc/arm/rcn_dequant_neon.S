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

#include "arm/asm.S"

.macro DERIVE_NB_COLS, out, in, tmp
    ORR \out, \in, 1
    CLZ \out, \out
    LSR \out, \out, #3
    MOV \tmp, #8
    SUB \out, \tmp, \out, LSL #2
.endm

.macro DERIVE_NB_ROWS, out, in, tmp
    MOV \out, #5
    CMP \in, \out
    CSEL \out, \in, \out, LT
    MOV \tmp, #1
    LSL \out, \tmp, \out
.endm

.macro MULLSHIFT, o1, o2, i1, i2, scale, shift, neg, s1=.4S, s2=.8H, s3=.4H
    SMULL2 \o1\s1, \i1\s2,\scale\s2
    SMULL \i1\s1, \i1\s3,\scale\s3
    SMULL2 \o2\s1, \i2\s2,\scale\s2
    SMULL \i2\s1, \i2\s3,\scale\s3
    .if \neg
        SSHL \o1\s1, \o1\s1, \shift\s1
        SSHL \i1\s1, \i1\s1, \shift\s1
        SSHL \o2\s1, \o2\s1, \shift\s1
        SSHL \i2\s1, \i2\s1, \shift\s1
    .else
        SRSHL \o1\s1, \o1\s1, \shift\s1
        SRSHL \i1\s1, \i1\s1, \shift\s1
        SRSHL \o2\s1, \o2\s1, \shift\s1
        SRSHL \i2\s1, \i2\s1, \shift\s1
    .endif
    XTN \o1\s3, \o1\s1
    XTN \i1\s3, \i1\s1
    XTN \o2\s3, \o2\s1
    XTN \i2\s3, \i2\s1
    UZP1 \o1\().2D, \i1\().2D, \o1\().2D
    UZP1 \o2\().2D, \i2\().2D, \o2\().2D
.endm

.macro DEQUANT_TB_4x4, neg
    .if \neg
        function dequant_tb_4x4_neg_neon, export=1
    .else
        function dequant_tb_4x4_neon, export=1
    .endif
    //DERIVE_NB_COLS  x9, x6, x7
    DERIVE_NB_ROWS  x9, x5, x7
    DERIVE_NB_ROWS x10, x4, x7
    // Force sig_sb_map to one in case of DC coefficient
    MOV x7, #1
    CMP x6, #0
    CSEL  x7, x7, xzr, EQ
    ORR x6, x6, x7
    // Compute stride
    LSL x4, x10, #1
    LSL x7, x4, #2
    DUP V2.8H, w2
    MOVI V31.8H, #0
    .if \neg == 0
        SUB w3, wzr, w3
    .endif
    DUP V3.4S, w3
    // LOOP1
    LSR x9, x9, #2
    MOV x5, #0
0:  CBZ x9, 5f
    SUB x9, x9, #1
    MOV x13, x0
    MOV x14, x1
    // uint8_t sig_sb_row = sig_sb_map >> (i << 3);
    LSL x11, x5, #3
    LSR x11, x6, x11
    //i++;
    ADD x5, x5, #1
    // LOOP2
    LSR x12, x10, #2
1:  CBZ x12, 4f
    SUB x12, x12, #1
    // int16_t *_dst = dst + (j << 2);
    MOV x15, x13
    ADD x13, x13, #8
    // if (sig_sb_row & 0x1)
    TST x11, #1
    B.NE 2f
    ADD x14, x14, #32
    ST1 {V31.D}[0], [x15], x4
    ST1 {V31.D}[1], [x15], x4
    ST1 {V31.D}[0], [x15], x4
    ST1 {V31.D}[1], [x15], x4
    B 3f
 2:
    LD1 {V12.8H}, [x14], #16
    LD1 {V13.8H}, [x14], #16
    MULLSHIFT V10, V11, V12, V13, V2, V3, \neg
    ST1 {V10.D}[0], [x15], x4
    ST1 {V10.D}[1], [x15], x4
    ST1 {V11.D}[0], [x15], x4
    ST1 {V11.D}[1], [x15], x4
3:
    LSR x11, x11, #1
    B 1b // LOOP2
4:    
    // dst += dst_stride << 2;
    ADD x0, x0, x7
    // src += src_stride << 2;
    ADD x1, x1, x7
    B 0b // LOOP1
5: 
    ret lr
endfunc
.endm

DEQUANT_TB_4x4 0
DEQUANT_TB_4x4 1