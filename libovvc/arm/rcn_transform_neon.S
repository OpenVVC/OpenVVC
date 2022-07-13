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
*   Ibrahim FARHAT
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


function idct_ii_4_neon, export=1


    mov w10, #0x03FF
    dup v31.8h, wzr
    dup v30.8h, w10
    dup v29.4s, w4
    dup v28.4s, w3

    ld1 {v1.8h}, [x0], x2
    ld1 {v2.8h}, [x0], x2
    ld1 {v3.8h}, [x0], x2
    ld1 {v4.8h}, [x0]

    ld1 {v5.8h-v6.8h}, [x5]


    //dct_4x4_neon v15.8h, v1, v2, v3, v4, v5, v6, v7, v8, #16
    //multiply and transpose
    SMULL  V10.4S, v5.4h, v1.H[0]
    SMLAL2 V10.4S, v5.8h, v2.H[0]
    SMLAL  V10.4S, v6.4h, V3.H[0]
    SMLAL2 V10.4S, v6.8h, v4.H[0]

    SMULL  V11.4S, v5.4h, v1.H[1]
    SMLAL2 V11.4S, v5.8h, v2.H[1]
    SMLAL  V11.4S, v6.4h, V3.H[1]
    SMLAL2 V11.4S, v6.8h, v4.H[1]

    SMULL  V12.4S, v5.4h, v1.H[2]
    SMLAL2 v12.4S, v5.8h, v2.H[2]
    SMLAL  V12.4S, v6.4h, V3.H[2]
    SMLAL2 V12.4S, v6.8h, v4.H[2]

    SMULL  V13.4S, v5.4h, v1.H[3]
    SMLAL2 v13.4S, v5.8h, v2.H[3]
    SMLAL  V13.4S, v6.4h, V3.H[3]
    SMLAL2 V13.4S, v6.8h, V4.H[3]

    ADD   V10.4S, V10.4S, V29.4S
    ADD   V11.4S, V11.4S, V29.4S
    ADD   V12.4S, V12.4S, V29.4S
    ADD   V13.4S, V13.4S, V29.4S

    SSHL  V10.4S, V10.4S, V28.4S // varialble shift right => negative shift left
    SSHL  V11.4S, V11.4S, V28.4S
    SSHL  V12.4S, V12.4S, V28.4S
    SSHL  v13.4S, V13.4S, V28.4S

    UZP1  v15.8h, V10.8H, V11.8H
    UZP1  v16.8h, V12.8H, V13.8H

    st1 {v15.8h-v16.8h}, [x1]

    ret lr
endfunc
