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

.macro inverse_neon, out, s1=.S

    mov v9.s[0], \out\s1[3]
    mov v9.s[1], \out\s1[2]
    mov v9.s[2], \out\s1[1]
    mov v9.s[3], \out\s1[0]

    mov \out\s1[0], v9.s[0]
    mov \out\s1[1], v9.s[1]
    mov \out\s1[2], v9.s[2]
    mov \out\s1[3], v9.s[3]

.endm

.macro idct_x_4, w

  function idct_x_4_\w\()_neon, export=1

    dup v29.4s, w4
    dup v28.4s, w3

    ld1 {v1.8h}, [x0], x2
    ld1 {v2.8h}, [x0], x2
    ld1 {v3.8h}, [x0], x2
    ld1 {v4.8h}, [x0]

    ld1 {v5.8h-v6.8h}, [x5]

    SMULL  V10.4S, v5.4h, v1.H[0]
    SMLAL2 V10.4S, v5.8h, v2.H[0]
    SMLAL  V10.4S, v6.4h, V3.H[0]
    SMLAL2 V10.4S, v6.8h, v4.H[0]
    .if \w == 2 || \w == 4 || \w == 8
        SMULL  V11.4S, v5.4h, v1.H[1]
        SMLAL2 V11.4S, v5.8h, v2.H[1]
        SMLAL  V11.4S, v6.4h, V3.H[1]
        SMLAL2 V11.4S, v6.8h, v4.H[1]
    .endif
    .if \w == 4 || \w == 8
        SMULL  V12.4S, v5.4h, v1.H[2]
        SMLAL2 v12.4S, v5.8h, v2.H[2]
        SMLAL  V12.4S, v6.4h, V3.H[2]
        SMLAL2 V12.4S, v6.8h, v4.H[2]

        SMULL  V13.4S, v5.4h, v1.H[3]
        SMLAL2 v13.4S, v5.8h, v2.H[3]
        SMLAL  V13.4S, v6.4h, V3.H[3]
        SMLAL2 V13.4S, v6.8h, V4.H[3]
    .endif
    .if \w == 8
        SMULL  V14.4S, v5.4h, v1.H[4]
        SMLAL2 V14.4S, v5.8h, v2.H[4]
        SMLAL  V14.4S, v6.4h, V3.H[4]
        SMLAL2 V14.4S, v6.8h, v4.H[4]

        SMULL  V15.4S, v5.4h, v1.H[5]
        SMLAL2 V15.4S, v5.8h, v2.H[5]
        SMLAL  V15.4S, v6.4h, V3.H[5]
        SMLAL2 V15.4S, v6.8h, v4.H[5]

        SMULL  V16.4S, v5.4h, v1.H[6]
        SMLAL2 v16.4S, v5.8h, v2.H[6]
        SMLAL  V16.4S, v6.4h, V3.H[6]
        SMLAL2 V16.4S, v6.8h, v4.H[6]

        SMULL  V17.4S, v5.4h, v1.H[7]
        SMLAL2 v17.4S, v5.8h, v2.H[7]
        SMLAL  V17.4S, v6.4h, V3.H[7]
        SMLAL2 V17.4S, v6.8h, V4.H[7]
    .endif

    ADD   V10.4S, V10.4S, V29.4S
    SSHL  V10.4S, V10.4S, V28.4S
    .if \w == 2 || \w == 4 || \w == 8
        ADD   V11.4S, V11.4S, V29.4S
        SSHL  V11.4S, V11.4S, V28.4S
    .endif
    .if \w == 4 || \w == 8
        ADD   V12.4S, V12.4S, V29.4S
        ADD   V13.4S, V13.4S, V29.4S
        SSHL  V12.4S, V12.4S, V28.4S
        SSHL  v13.4S, V13.4S, V28.4S
    .endif
    .if \w == 8
        ADD   V14.4S, V14.4S, V29.4S
        ADD   V15.4S, V15.4S, V29.4S
        ADD   V16.4S, V16.4S, V29.4S
        ADD   V17.4S, V17.4S, V29.4S
        SSHL  V14.4S, V14.4S, V28.4S
        SSHL  V15.4S, V15.4S, V28.4S
        SSHL  V16.4S, V16.4S, V28.4S
        SSHL  v17.4S, V17.4S, V28.4S
    .endif

    UZP1  v18.8h, V10.8H, V11.8H
    .if \w == 4 || \w == 8
        UZP1  v19.8h, V12.8H, V13.8H
    .endif
    .if \w == 8
        UZP1  v20.8h, V14.8H, V15.8H
        UZP1  v21.8h, V16.8H, V17.8H
    .endif

    .if \w == 1 || \w == 2
        st1 {v18.8h}, [x1]
    .endif
    .if \w == 4
        st1 {v18.8h-v19.8h}, [x1]
    .endif
    .if \w == 8
        st1 {v18.8h-v21.8h}, [x1]
    .endif
    ret lr
  endfunc
.endm

.macro idct_ii_8, w
    function idct_ii_8_\w\()_neon, export=1

        dup v29.4s, w4
        dup v28.4s, w3

        ld1 {v1.8h}, [x0], x2
        ld1 {v2.8h}, [x0], x2
        ld1 {v3.8h}, [x0], x2
        ld1 {v4.8h}, [x0], x2
        ld1 {v5.8h}, [x0], x2
        ld1 {v6.8h}, [x0], x2
        ld1 {v7.8h}, [x0], x2
        ld1 {v8.8h}, [x0]

        ld1 {v9.8h-v10.8h}, [x5]
        // even
        SMULL  V13.4S, v9.4h , v1.H[0]
        SMLAL2 V13.4S, v9.8h , v3.H[0]
        SMLAL  V13.4S, v10.4h, V5.H[0]
        SMLAL2 V13.4S, v10.8h, v7.H[0]

        SMULL  V14.4S, v9.4h , v1.H[1]
        SMLAL2 V14.4S, v9.8h , v3.H[1]
        SMLAL  V14.4S, v10.4h, V5.H[1]
        SMLAL2 V14.4S, v10.8h, v7.H[1]
        .if \w == 4
            SMULL  V15.4S, v9.4h , v1.H[2]
            SMLAL2 v15.4S, v9.8h , v3.H[2]
            SMLAL  V15.4S, v10.4h, V5.H[2]
            SMLAL2 V15.4S, v10.8h, v7.H[2]

            SMULL  V16.4S, v9.4h , v1.H[3]
            SMLAL2 v16.4S, v9.8h , v3.H[3]
            SMLAL  V16.4S, v10.4h, V5.H[3]
            SMLAL2 V16.4S, v10.8h, V7.H[3]
        .endif
        // odd
        ld1 {v9.8h-v10.8h}, [x6]
        SMULL  V17.4S, v9.4h , v2.H[0]
        SMLAL2 V17.4S, v9.8h , v4.H[0]
        SMLAL  V17.4S, v10.4h, V6.H[0]
        SMLAL2 V17.4S, v10.8h, v8.H[0]

        SMULL  V18.4S, v9.4h , v2.H[1]
        SMLAL2 V18.4S, v9.8h , v4.H[1]
        SMLAL  V18.4S, v10.4h, V6.H[1]
        SMLAL2 V18.4S, v10.8h, v8.H[1]
        .if \w == 4
            SMULL  V19.4S, v9.4h , v2.H[2]
            SMLAL2 v19.4S, v9.8h , v4.H[2]
            SMLAL  V19.4S, v10.4h, V6.H[2]
            SMLAL2 V19.4S, v10.8h, v8.H[2]

            SMULL  V20.4S, v9.4h , v2.H[3]
            SMLAL2 v20.4S, v9.8h , v4.H[3]
            SMLAL  V20.4S, v10.4h, V6.H[3]
            SMLAL2 V20.4S, v10.8h, V8.H[3]
        .endif
        // out butterfly (add sub)
        ADD v9.4s , v13.4s, v17.4s
        ADD v10.4s, v14.4s, v18.4s
        .if \w == 4
            ADD v11.4s, v15.4s, v19.4s
            ADD v12.4s, v16.4s, v20.4s
        .endif
        //
        sub v21.4s, v13.4s, v17.4s
        sub v22.4s, v14.4s, v18.4s
        .if \w == 4
            sub v23.4s, v15.4s, v19.4s
            sub v24.4s, v16.4s, v20.4s
        .endif

        // shift and round
        ADD   V13.4S, v9.4S , V29.4S
        ADD   V14.4S, V10.4S, V29.4S
        SSHL  V13.4S, V13.4S, V28.4S
        SSHL  V14.4S, V14.4S, V28.4S
        .if \w == 4
            ADD   V15.4S, V11.4S, V29.4S
            ADD   V16.4S, V12.4S, V29.4S
            SSHL  V15.4S, V15.4S, V28.4S
            SSHL  v16.4S, V16.4S, V28.4S
        .endif
        ADD   V17.4S, V21.4S, V29.4S
        ADD   V18.4S, V22.4S, V29.4S
        SSHL  V17.4S, V17.4S, V28.4S
        SSHL  V18.4S, V18.4S, V28.4S
        .if \w == 4
            ADD   V19.4S, V23.4S, V29.4S
            ADD   V20.4S, V24.4S, V29.4S
            SSHL  V19.4S, V19.4S, V28.4S
            SSHL  v20.4S, V20.4S, V28.4S
        .endif

        // inverse the order of the second vector before merge
        inverse_neon v17
        inverse_neon v18
        .if \w == 4
            inverse_neon v19
            inverse_neon v20
        .endif

        UZP1  v21.8h, V13.8H, V17.8H
        UZP1  v22.8h, V14.8H, V18.8H
        .if \w == 4
            UZP1  v23.8h, V15.8H, V19.8H
            UZP1  v24.8h, V16.8H, V20.8H
        .endif

        .if \w == 2
            st1 {v21.8h-v22.8h}, [x1]
        .endif
        .if \w == 4
            st1 {v21.8h-v24.8h}, [x1]
        .endif
        ret lr
    endfunc
.endm
// 4-pt dct-2/8 and dst-7
idct_x_4 2
idct_x_4 4
idct_x_4 8
// 8-pt dct-2
idct_ii_8 4
idct_ii_8 2

