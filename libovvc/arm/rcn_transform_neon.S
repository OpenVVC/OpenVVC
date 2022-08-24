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


.macro CLIP, val, min_val, max_val
       SMAX  \val, \val, \min_val
       SMIN  \val, \val, \max_val
.endm 

.macro inverse_neon, out, s1=.S

    mov v31.s[0], \out\s1[3]
    mov v31.s[1], \out\s1[2]
    mov v31.s[2], \out\s1[1]
    mov v31.s[3], \out\s1[0]

    mov \out\s1[0], v31.s[0]
    mov \out\s1[1], v31.s[1]
    mov \out\s1[2], v31.s[2]
    mov \out\s1[3], v31.s[3]
    
.endm

.macro store_fun, size
    function store_\size\()_neon, export=1
        dup v1.8h, w1
        dup v2.8h, w1
        dup v3.8h, w1
        dup v4.8h, w1

        .if \size == 8
            st1 {v1.8h}, [x0]
        .endif
        .if \size == 16
            st1 {v1.8h-v2.8h}, [x0]
        .endif 
        .if \size == 32
            st1 {v1.8h-v4.8h}, [x0]
        .endif
        ret lr
    endfunc
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

.macro idct_ii_8, w, red
    function idct_ii_8_\w\()_\red\()_neon, export=1
    
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
        .if \red == 8
            SMLAL  V13.4S, v10.4h, V5.H[0]
            SMLAL2 V13.4S, v10.8h, v7.H[0]
        .endif
        SMULL  V14.4S, v9.4h , v1.H[1]
        SMLAL2 V14.4S, v9.8h , v3.H[1]
        .if \red == 8
            SMLAL  V14.4S, v10.4h, V5.H[1]
            SMLAL2 V14.4S, v10.8h, v7.H[1]
        .endif
        .if \w == 4
            SMULL  V15.4S, v9.4h , v1.H[2]
            SMLAL2 v15.4S, v9.8h , v3.H[2]
            .if \red == 8
                SMLAL  V15.4S, v10.4h, V5.H[2]
                SMLAL2 V15.4S, v10.8h, v7.H[2]
            .endif    
            SMULL  V16.4S, v9.4h , v1.H[3]
            SMLAL2 v16.4S, v9.8h , v3.H[3]
            .if \red == 8
                SMLAL  V16.4S, v10.4h, V5.H[3]
                SMLAL2 V16.4S, v10.8h, V7.H[3]
            .endif
        .endif
        // odd 
        ld1 {v9.8h-v10.8h}, [x6]
        SMULL  V17.4S, v9.4h , v2.H[0]
        SMLAL2 V17.4S, v9.8h , v4.H[0]
        .if \red == 8
            SMLAL  V17.4S, v10.4h, V6.H[0]
            SMLAL2 V17.4S, v10.8h, v8.H[0]
        .endif
        SMULL  V18.4S, v9.4h , v2.H[1]
        SMLAL2 V18.4S, v9.8h , v4.H[1]
        .if \red == 8
            SMLAL  V18.4S, v10.4h, V6.H[1]
            SMLAL2 V18.4S, v10.8h, v8.H[1]
        .endif
        .if \w == 4 
            SMULL  V19.4S, v9.4h , v2.H[2]
            SMLAL2 v19.4S, v9.8h , v4.H[2]
            .if \red == 8
                SMLAL  V19.4S, v10.4h, V6.H[2]
                SMLAL2 V19.4S, v10.8h, v8.H[2]
            .endif
            SMULL  V20.4S, v9.4h , v2.H[3]
            SMLAL2 v20.4S, v9.8h , v4.H[3]
            .if \red ==8
                SMLAL  V20.4S, v10.4h, V6.H[3]
                SMLAL2 V20.4S, v10.8h, V8.H[3]
            .endif
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

.macro idct_x_8, w, red
    function idct_x_8_\w\()_\red\()_neon, export=1
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
        
        ld1 {v20.8h-v23.8h}, [x5]
        ld1 {v24.8h-v27.8h}, [x6]
        // without butterfly
        // top section of the first vector
        SMULL V10.4S, v20.4h, v1.H[0]
        SMLAL V10.4S, v21.4h, V2.H[0]
        SMLAL V10.4S, v22.4h, V3.H[0]
        SMLAL V10.4S, v23.4h, V4.H[0]
        .if \red == 8 
            SMLAL V10.4S, v24.4h, V5.H[0]
            SMLAL V10.4S, v25.4h, V6.H[0]
            SMLAL V10.4S, v26.4h, V7.H[0]
            SMLAL V10.4S, v27.4h, V8.H[0]
        .endif
        // bootom section of the first vector
        SMULL2 V11.4S, v20.8h, v1.H[0]
        SMLAL2 V11.4S, v21.8h, v2.H[0]
        SMLAL2 V11.4S, v22.8h, v3.H[0]
        SMLAL2 V11.4S, v23.8h, v4.H[0]
        .if \red == 8
            SMLAL2 V11.4S, v24.8h, v5.H[0]
            SMLAL2 V11.4S, v25.8h, v6.H[0]
            SMLAL2 V11.4S, v26.8h, v7.H[0]
            SMLAL2 V11.4S, v27.8h, v8.H[0]
        .endif
        // top section of the second vector
        SMULL V12.4S, v20.4h, v1.H[1]
        SMLAL V12.4S, v21.4h, V2.H[1]
        SMLAL V12.4S, v22.4h, V3.H[1]
        SMLAL V12.4S, v23.4h, V4.H[1]
        .if \red == 8
            SMLAL V12.4S, v24.4h, V5.H[1]
            SMLAL V12.4S, v25.4h, V6.H[1]
            SMLAL V12.4S, v26.4h, V7.H[1]
            SMLAL V12.4S, v27.4h, V8.H[1]
        .endif
        // bootom section of the second vector
        SMULL2 V13.4S, v20.8h, v1.H[1]
        SMLAL2 V13.4S, v21.8h, v2.H[1]
        SMLAL2 V13.4S, v22.8h, v3.H[1]
        SMLAL2 V13.4S, v23.8h, v4.H[1]
        .if \red == 8
            SMLAL2 V13.4S, v24.8h, v5.H[1]
            SMLAL2 V13.4S, v25.8h, v6.H[1]
            SMLAL2 V13.4S, v26.8h, v7.H[1]
            SMLAL2 V13.4S, v27.8h, v8.H[1]
        .endif
        // top section of the third ODD vector
        .if \w == 4 || \w == 8
            SMULL V14.4S, v20.4h, v1.H[2]
            SMLAL V14.4S, v21.4h, V2.H[2]
            SMLAL V14.4S, v22.4h, V3.H[2]
            SMLAL V14.4S, v23.4h, V4.H[2]
            .if \red == 8
                SMLAL V14.4S, v24.4h, V5.H[2]
                SMLAL V14.4S, v25.4h, V6.H[2]
                SMLAL V14.4S, v26.4h, V7.H[2]
                SMLAL V14.4S, v27.4h, V8.H[2]
            .endif
            // bootom section of the third vector
            SMULL2 V15.4S, v20.8h, v1.H[2]
            SMLAL2 V15.4S, v21.8h, v2.H[2]
            SMLAL2 V15.4S, v22.8h, v3.H[2]
            SMLAL2 V15.4S, v23.8h, v4.H[2]
            .if \red == 8
                SMLAL2 V15.4S, v24.8h, v5.H[2]
                SMLAL2 V15.4S, v25.8h, v6.H[2]
                SMLAL2 V15.4S, v26.8h, v7.H[2]
                SMLAL2 V15.4S, v27.8h, v8.H[2]
            .endif
            // top section of the fourth  vector
            SMULL V16.4S, v20.4h, v1.H[3]
            SMLAL V16.4S, v21.4h, V2.H[3]
            SMLAL V16.4S, v22.4h, V3.H[3]
            SMLAL V16.4S, v23.4h, V4.H[3]
            .if \red == 8
                SMLAL v16.4S, v24.4h, V5.H[3]
                SMLAL v16.4S, v25.4h, V6.H[3]
                SMLAL v16.4S, v26.4h, V7.H[3]
                SMLAL v16.4S, v27.4h, V8.H[3]
            .endif
            // bootom section of the fourth vector
            SMULL2 V17.4S, v20.8h, v1.H[3]
            SMLAL2 V17.4S, v21.8h, v2.H[3]
            SMLAL2 V17.4S, v22.8h, v3.H[3]
            SMLAL2 V17.4S, v23.8h, v4.H[3]
            .if \red == 8
                SMLAL2 V17.4S, v24.8h, v5.H[3]
                SMLAL2 V17.4S, v25.8h, v6.H[3]
                SMLAL2 V17.4S, v26.8h, v7.H[3]
                SMLAL2 V17.4S, v27.8h, v8.H[3]
            .endif
        .endif
        // shift and round 
        ADD   V10.4S, v10.4S, V29.4S
        ADD   V11.4S, V11.4S, V29.4S
        ADD   V12.4S, V12.4S, V29.4S
        ADD   V13.4S, V13.4S, V29.4S
        SSHL  V10.4S, V10.4S, V28.4S 
        SSHL  V11.4S, V11.4S, V28.4S
        SSHL  V12.4S, V12.4S, V28.4S 
        SSHL  V13.4S, V13.4S, V28.4S
        .if \w == 4 || \w == 8
            ADD   V14.4S, V14.4S, V29.4S
            ADD   V15.4S, v15.4S, V29.4S
            ADD   V16.4S, V16.4S, V29.4S
            ADD   V17.4S, V17.4S, V29.4S
            SSHL  V14.4S, V14.4S, V28.4S
            SSHL  V15.4S, V15.4S, V28.4S 
            SSHL  V16.4S, V16.4S, V28.4S
            SSHL  V17.4S, V17.4S, V28.4S 
        .endif

        UZP1  v18.8h , V10.8H, V11.8H
        UZP1  v19.8h, V12.8H, V13.8H
        .if \w == 4 || \w == 8
            UZP1  v30.8h, V14.8H, V15.8H
            UZP1  v31.8h, V16.8H, V17.8H
        .endif
        .if \w == 8
            // top section of the first vector
            SMULL V10.4S, v20.4h, v1.H[4]
            SMLAL V10.4S, v21.4h, V2.H[4]
            SMLAL V10.4S, v22.4h, V3.H[4]
            SMLAL V10.4S, v23.4h, V4.H[4]
            .if \red == 8
                SMLAL V10.4S, v24.4h, V5.H[4]
                SMLAL V10.4S, v25.4h, V6.H[4]
                SMLAL V10.4S, v26.4h, V7.H[4]
                SMLAL V10.4S, v27.4h, V8.H[4]
            .endif
            // bootom section of the first vector
            SMULL2 V11.4S, v20.8h, v1.H[4]
            SMLAL2 V11.4S, v21.8h, v2.H[4]
            SMLAL2 V11.4S, v22.8h, v3.H[4]
            SMLAL2 V11.4S, v23.8h, v4.H[4]
            .if \red == 8
                SMLAL2 V11.4S, v24.8h, v5.H[4]
                SMLAL2 V11.4S, v25.8h, v6.H[4]
                SMLAL2 V11.4S, v26.8h, v7.H[4]
                SMLAL2 V11.4S, v27.8h, v8.H[4]
            .endif
            // top section of the second vector
            SMULL V12.4S, v20.4h, v1.H[5]
            SMLAL V12.4S, v21.4h, V2.H[5]
            SMLAL V12.4S, v22.4h, V3.H[5]
            SMLAL V12.4S, v23.4h, V4.H[5]
            .if \red == 8
                SMLAL V12.4S, v24.4h, V5.H[5]
                SMLAL V12.4S, v25.4h, V6.H[5]
                SMLAL V12.4S, v26.4h, V7.H[5]
                SMLAL V12.4S, v27.4h, V8.H[5]
            .endif
            // bootom section of the second vector
            SMULL2 V13.4S, v20.8h, v1.H[5]
            SMLAL2 V13.4S, v21.8h, v2.H[5]
            SMLAL2 V13.4S, v22.8h, v3.H[5]
            SMLAL2 V13.4S, v23.8h, v4.H[5]
            .if \red == 8
                SMLAL2 V13.4S, v24.8h, v5.H[5]
                SMLAL2 V13.4S, v25.8h, v6.H[5]
                SMLAL2 V13.4S, v26.8h, v7.H[5]
                SMLAL2 V13.4S, v27.8h, v8.H[5]
            .endif
            // top section of the third ODD vector
            SMULL V14.4S, v20.4h, v1.H[6]
            SMLAL V14.4S, v21.4h, V2.H[6]
            SMLAL V14.4S, v22.4h, V3.H[6]
            SMLAL V14.4S, v23.4h, V4.H[6]
            .if \red == 8
                SMLAL V14.4S, v24.4h, V5.H[6]
                SMLAL V14.4S, v25.4h, V6.H[6]
                SMLAL V14.4S, v26.4h, V7.H[6]
                SMLAL V14.4S, v27.4h, V8.H[6]
            .endif
            // bootom section of the third vector
            SMULL2 V15.4S, v20.8h, v1.H[6]
            SMLAL2 V15.4S, v21.8h, v2.H[6]
            SMLAL2 V15.4S, v22.8h, v3.H[6]
            SMLAL2 V15.4S, v23.8h, v4.H[6]
            .if \red == 8
                SMLAL2 V15.4S, v24.8h, v5.H[6]
                SMLAL2 V15.4S, v25.8h, v6.H[6]
                SMLAL2 V15.4S, v26.8h, v7.H[6]
                SMLAL2 V15.4S, v27.8h, v8.H[6]
            .endif
            // top section of the fourth  vector
            SMULL V16.4S, v20.4h, v1.H[7]
            SMLAL V16.4S, v21.4h, V2.H[7]
            SMLAL V16.4S, v22.4h, V3.H[7]
            SMLAL V16.4S, v23.4h, V4.H[7]
            .if \red == 8
                SMLAL v16.4S, v24.4h, V5.H[7]
                SMLAL v16.4S, v25.4h, V6.H[7]
                SMLAL v16.4S, v26.4h, V7.H[7]
                SMLAL v16.4S, v27.4h, V8.H[7]
            .endif
            // bootom section of the fourth vector
            SMULL2 V17.4S, v20.8h, v1.H[7]
            SMLAL2 V17.4S, v21.8h, v2.H[7]
            SMLAL2 V17.4S, v22.8h, v3.H[7]
            SMLAL2 V17.4S, v23.8h, v4.H[7]
            .if \red == 8
                SMLAL2 V17.4S, v24.8h, v5.H[7]
                SMLAL2 V17.4S, v25.8h, v6.H[7]
                SMLAL2 V17.4S, v26.8h, v7.H[7]
                SMLAL2 V17.4S, v27.8h, v8.H[7]
            .endif

            ADD   V10.4S, v10.4S, V29.4S
            ADD   V11.4S, V11.4S, V29.4S
            ADD   V12.4S, V12.4S, V29.4S
            ADD   V13.4S, V13.4S, V29.4S
            ADD   V14.4S, V14.4S, V29.4S
            ADD   V15.4S, v15.4S, V29.4S
            ADD   V16.4S, V16.4S, V29.4S
            ADD   V17.4S, V17.4S, V29.4S

            SSHL  V10.4S, V10.4S, V28.4S 
            SSHL  V11.4S, V11.4S, V28.4S
            SSHL  V12.4S, V12.4S, V28.4S 
            SSHL  V13.4S, V13.4S, V28.4S
            SSHL  V14.4S, V14.4S, V28.4S
            SSHL  V15.4S, V15.4S, V28.4S 
            SSHL  V16.4S, V16.4S, V28.4S
            SSHL  V17.4S, V17.4S, V28.4S 

            UZP1  v9.8h , V10.8H, V11.8H
            UZP1  v10.8h, V12.8H, V13.8H
            UZP1  v11.8h, V14.8H, V15.8H
            UZP1  v12.8h, V16.8H, V17.8H
            
        .endif
        .if \w == 2 
            st1 {v18.8h-v19.8h}, [x1]
        .endif
        .if \w == 4
            st1 {v18.8h-v19.8h}, [x1], #32
            st1 {v30.8h-v31.8h}, [x1]
        .endif
        .if \w == 8
            st1 {v18.8h-v19.8h}, [x1], #32
            st1 {v30.8h-v31.8h}, [x1], #32
            st1 {v9.8h-v12.8h} , [x1]
        .endif
        ret lr
    endfunc
.endm

.macro idct_ii_16_2, w
    function idct_ii_16_2_\w\()_neon, export=1
    
        dup v29.4s, w4 
        dup v28.4s, w3

        ld1 {v0.8h}, [x0], x2
        ld1 {v1.8h}, [x0], x2
        ld1 {v2.8h}, [x0], x2
        ld1 {v3.8h}, [x0], x2
        ld1 {v4.8h}, [x0], x2
        ld1 {v5.8h}, [x0], x2
        ld1 {v6.8h}, [x0], x2
        ld1 {v7.8h}, [x0], x2
        ld1 {v8.8h}, [x0], x2
        ld1 {v9.8h}, [x0], x2
        ld1 {v10.8h}, [x0], x2
        ld1 {v11.8h}, [x0], x2
        ld1 {v12.8h}, [x0], x2
        ld1 {v13.8h}, [x0], x2
        ld1 {v14.8h}, [x0], x2
        ld1 {v15.8h}, [x0]

        ld1 {v20.8h-v23.8h}, [x5]
        // even 8x8 
        SMULL  V16.4S, v20.4h, v0.H[0]
        .if \w == 8 || \w == 16
            SMLAL2 V16.4S, v20.8h, v4.H[0]
        .endif
        .if \w == 16
            SMLAL  V16.4S, v21.4h, V8.H[0]
            SMLAL2 V16.4S, v21.8h, v12.H[0]
        .endif
        
        SMULL  V17.4S, v20.4h, v0.H[1]
        .if \w == 8 || \w == 16
            SMLAL2 V17.4S, v20.8h, v4.H[1]
        .endif
        .if \w == 16
            SMLAL  V17.4S, v21.4h, V8.H[1]
            SMLAL2 V17.4S, v21.8h, v12.H[1]
        .endif
        // odd 8x8
        SMULL  V18.4S, v22.4h, v2.H[0]
        .if \w == 8 || \w == 16
            SMLAL2 V18.4S, v22.8h, v6.H[0]
        .endif
        .if \w == 16
            SMLAL  V18.4S, v23.4h, V10.H[0]
            SMLAL2 V18.4S, v23.8h, v14.H[0]
        .endif
        SMULL  V19.4S, v22.4h , v2.H[1]
        .if \w == 8 || \w == 16
            SMLAL2 V19.4S, v22.8h , v6.H[1]
        .endif
        .if \w == 16
            SMLAL  V19.4S, v23.4h, V10.H[1]
            SMLAL2 V19.4S, v23.8h, v14.H[1]
        .endif

        // first out butterfly (add sub)
        ADD v2.4s, v16.4s, v18.4s // top section of the first even vector
        ADD v4.4s, v17.4s, v19.4s // top section of the second even vector 

        sub v6.4s, v16.4s, v18.4s // bottom section of the first even vector
        sub v8.4s, v17.4s, v19.4s // bottom section of the second even vector

        inverse_neon v6
        inverse_neon v8

        ld1 {v20.8h-v23.8h}, [x6]
        ld1 {v24.8h-v27.8h}, [x7]
        // top section of the first ODD vector
        SMULL V0.4S, v20.4h, v1.H[0]
        SMLAL V0.4S, v21.4h, V3.H[0]
        .if \w == 8 || \w == 16
            SMLAL V0.4S, v22.4h, V5.H[0]
            SMLAL V0.4S, v23.4h, V7.H[0]
        .endif
        .if \w == 16
            SMLAL V0.4S, v24.4h, V9.H[0]
            SMLAL V0.4S, v25.4h, V11.H[0]
            SMLAL V0.4S, v26.4h, V13.H[0]
            SMLAL V0.4S, v27.4h, V15.H[0]
        .endif
        // bootom section of the first ODD vector
        
        SMULL2 V10.4S, v20.8h, v1.H[0]
        SMLAL2 V10.4S, v21.8h, v3.H[0]
        .if \w == 8 || \w == 16
            SMLAL2 V10.4S, v22.8h, v5.H[0]
            SMLAL2 V10.4S, v23.8h, v7.H[0]
        .endif
        .if \w == 16
            SMLAL2 V10.4S, v24.8h, v9.H[0]
            SMLAL2 V10.4S, v25.8h, v11.H[0]
            SMLAL2 V10.4S, v26.8h, v13.H[0]
            SMLAL2 V10.4S, v27.8h, v15.H[0]
        .endif
        // top section of the second ODD vector
        SMULL V12.4S, v20.4h, v1.H[1]
        SMLAL V12.4S, v21.4h, V3.H[1]
        .if \w == 8 || \w == 16
            SMLAL V12.4S, v22.4h, V5.H[1]
            SMLAL V12.4S, v23.4h, V7.H[1]
        .endif
        .if \w == 16
            SMLAL V12.4S, v24.4h, V9.H[1]
            SMLAL V12.4S, v25.4h, V11.H[1]
            SMLAL V12.4S, v26.4h, V13.H[1]
            SMLAL V12.4S, v27.4h, V15.H[1]
        .endif
        // bootom section of the second ODD vector
        SMULL2 V14.4S, v20.8h, v1.H[1]
        SMLAL2 V14.4S, v21.8h, v3.H[1]
        .if \w == 8 || \w == 16
            SMLAL2 V14.4S, v22.8h, v5.H[1]
            SMLAL2 V14.4S, v23.8h, v7.H[1]
        .endif
        .if \w == 16
            SMLAL2 V14.4S, v24.8h, v9.H[1]
            SMLAL2 V14.4S, v25.8h, v11.H[1]
            SMLAL2 V14.4S, v26.8h, v13.H[1]
            SMLAL2 V14.4S, v27.8h, v15.H[1]
        .endif

        // second out butterfly (add sub): first vector
        ADD v1.4s, v2.4s, v0.4s
        ADD v3.4s, v6.4s, v10.4s
        SUB v5.4s, v2.4s, v0.4s
        SUB v7.4s, v6.4s, v10.4s
        // second out butterfly (add sub): second vector
        ADD v9.4s , v4.4s, v12.4s
        ADD v11.4s, v8.4s, v14.4s
        SUB v13.4s, v4.4s, v12.4s
        SUB v15.4s, v8.4s, v14.4s

        // shift and round 
        ADD   V1.4S , v1.4S , V29.4S
        ADD   V3.4S , V3.4S , V29.4S
        ADD   V5.4S , V5.4S , V29.4S
        ADD   V7.4S , V7.4S , V29.4S

        ADD   V9.4S , V9.4S , V29.4S
        ADD   V11.4S, V11.4S, V29.4S
        ADD   V13.4S, V13.4S , V29.4S
        ADD   V15.4S, V15.4S, V29.4S

        SSHL  V1.4S , V1.4S , V28.4S 
        SSHL  V3.4S , V3.4S , V28.4S
        SSHL  V5.4S , V5.4S , V28.4S 
        SSHL  V7.4S , V7.4S , V28.4S

        SSHL  V9.4S , V9.4S , V28.4S 
        SSHL  V11.4S, V11.4S, V28.4S
        SSHL  V13.4S, V13.4S, V28.4S 
        SSHL  V15.4S, V15.4S, V28.4S
        
        // inverse the order of the second vector before merge : complete
        inverse_neon v5
        inverse_neon v7 
        inverse_neon v13
        inverse_neon v15 
        //
        UZP1  v21.8h, V1.8H , V3.8H
        UZP1  v22.8h, V7.8H , V5.8H
        UZP1  v23.8h, V9.8H , V11.8H
        UZP1  v24.8h, V15.8H, V13.8H

        st1 {v21.8h-v24.8h}, [x1]
        ret lr
    endfunc
.endm

.macro idct_ii_32, w
    function idct_ii_32_\w\()_neon, export=1

        ld1 {v0.8h}, [x0], x2
        ld1 {v1.8h}, [x0], x2
        ld1 {v2.8h}, [x0], x2
        ld1 {v3.8h}, [x0], x2
        .if \w == 8 || \w == 16 || \w == 32
            ld1 {v4.8h}, [x0], x2
            ld1 {v5.8h}, [x0], x2
            ld1 {v6.8h}, [x0], x2
            ld1 {v7.8h}, [x0], x2
        .endif
        .if \w == 16 || \w == 32
            ld1 {v8.8h}, [x0], x2
            ld1 {v9.8h}, [x0], x2
            ld1 {v10.8h}, [x0], x2
            ld1 {v11.8h}, [x0], x2
            ld1 {v12.8h}, [x0], x2
            ld1 {v13.8h}, [x0], x2
            ld1 {v14.8h}, [x0], x2
            ld1 {v15.8h}, [x0], x2
        .endif
        .if \w == 32
            ld1 {v16.8h}, [x0], x2
            ld1 {v17.8h}, [x0], x2
            ld1 {v18.8h}, [x0], x2
            ld1 {v19.8h}, [x0], x2
            ld1 {v20.8h}, [x0], x2
            ld1 {v21.8h}, [x0], x2
            ld1 {v22.8h}, [x0], x2
            ld1 {v23.8h}, [x0], x2
            ld1 {v24.8h}, [x0], x2
            ld1 {v25.8h}, [x0], x2
            ld1 {v26.8h}, [x0], x2
            ld1 {v27.8h}, [x0], x2
            ld1 {v28.8h}, [x0], x2
            ld1 {v29.8h}, [x0], x2
            ld1 {v30.8h}, [x0], x2
            ld1 {v31.8h}, [x0]
        .endif

        mov v0.H[1], v1.H[0]
        mov v0.H[2], v2.H[0]
        mov v0.H[3], v3.H[0]
        .if \w == 8 || \w == 16 || \w == 32
            mov v0.H[4], v4.H[0]
            mov v0.H[5], v5.H[0]
            mov v0.H[6], v6.H[0]
            mov v0.H[7], v7.H[0]
        .endif
        .if \w == 16 || \w == 32
            mov v1.H[0], v8.H[0]
            mov v1.H[1], v9.H[0]
            mov v1.H[2], v10.H[0]
            mov v1.H[3], v11.H[0]
            mov v1.H[4], v12.H[0]
            mov v1.H[5], v13.H[0]
            mov v1.H[6], v14.H[0]
            mov v1.H[7], v15.H[0]
        .endif 
        .if \w == 32
            mov v2.H[0], v16.H[0]
            mov v2.H[1], v17.H[0]
            mov v2.H[2], v18.H[0]
            mov v2.H[3], v19.H[0]
            mov v2.H[4], v20.H[0]
            mov v2.H[5], v21.H[0]
            mov v2.H[6], v22.H[0]
            mov v2.H[7], v23.H[0]
            //
            mov v3.H[0], v24.H[0]
            mov v3.H[1], v25.H[0]
            mov v3.H[2], v26.H[0]
            mov v3.H[3], v27.H[0]
            mov v3.H[4], v28.H[0]
            mov v3.H[5], v29.H[0]
            mov v3.H[6], v30.H[0]
            mov v3.H[7], v31.H[0]
        .endif

        sub  w10, wzr , w3 // -shift
        sub  w11, w3  , #1
        mov  w9, #1
        lsl  w9, w9 , w11   // 1<< (shift-1)
        dup v28.4s, w10
        dup v29.4s, w9
        
        ld1 {v20.8h-v23.8h}, [x4]
        // even 8x8 
        SMULL  V16.4S, v20.4h, v0.H[0]
        .if \w == 16 || \w == 32
            SMLAL2 V16.4S, v20.8h, v1.H[0]
        .endif
        .if \w == 32
            SMLAL  V16.4S, v21.4h, V2.H[0]
            SMLAL2 V16.4S, v21.8h, v3.H[0]
        .endif
        // odd 8x8
        .if \w == 4
            dup v17.4S, wzr
        .endif
        .if \w == 8 || \w == 16 || \w == 32
            SMULL  V17.4S, v22.4h, v0.H[4]
        .endif
        .if \w == 16 || \w == 32
            SMLAL2 V17.4S, v22.8h, v1.H[4]
        .endif
        .if \w == 32
            SMLAL  V17.4S, v23.4h, V2.H[4]
            SMLAL2 V17.4S, v23.8h, v3.H[4]
        .endif
        // first out butterfly (add sub)
        ADD v5.4s, v16.4s, v17.4s // top section of the first even vector
        sub v6.4s, v16.4s, v17.4s // bottom section of the first even vector

        inverse_neon v6

        ld1 {v20.8h-v23.8h}, [x5]
        ld1 {v24.8h-v27.8h}, [x6]
        // top section of the first ODD vector
        SMULL V16.4S, v20.4h, v0.H[2]
        .if \w == 8 || \w == 16 || \w == 32
            SMLAL V16.4S, v21.4h, V0.H[6]
        .endif
        .if \w == 16 || \w == 32
            SMLAL V16.4S, v22.4h, V1.H[2]
            SMLAL V16.4S, v23.4h, V1.H[6]
        .endif
        .if \w == 32
            SMLAL V16.4S, v24.4h, V2.H[2]
            SMLAL V16.4S, v25.4h, V2.H[6]
            SMLAL V16.4S, v26.4h, V3.H[2]
            SMLAL V16.4S, v27.4h, V3.H[6]
        .endif
        // bootom section of the first ODD vector
        SMULL2 V17.4S, v20.8h, v0.H[2]
        .if \w == 8 || \w == 16 || \w == 32
            SMLAL2 V17.4S, v21.8h, v0.H[6]
        .endif
        .if \w == 16 || \w == 32
            SMLAL2 V17.4S, v22.8h, v1.H[2]
            SMLAL2 V17.4S, v23.8h, v1.H[6]
        .endif
        .if \w == 32
            SMLAL2 V17.4S, v24.8h, v2.H[2]
            SMLAL2 V17.4S, v25.8h, v2.H[6]
            SMLAL2 V17.4S, v26.8h, v3.H[2]
            SMLAL2 V17.4S, v27.8h, v3.H[6]
        .endif
        // second out butterfly (add sub): first vector
        ADD v7.4s , v5.4s, v16.4s
        ADD v8.4s , v6.4s, v17.4s
        SUB v9.4s , v5.4s, v16.4s
        SUB v10.4s, v6.4s, v17.4s
        
        inverse_neon v9 
        inverse_neon v10 

        ld1 {v20.8h-v23.8h}, [x7], #64

        SMULL  V11.4S, v20.4h, v0.H[1]
        SMULL2 V12.4S, v20.8h, v0.H[1]
        SMULL  V13.4S, v21.4h, v0.H[1]
        SMULL2 V14.4S, v21.8h, v0.H[1]
        
        SMLAL  V11.4S, v22.4h, V0.H[3]
        SMLAL2 V12.4S, v22.8h, V0.H[3]
        SMLAL  V13.4S, v23.4h, V0.H[3]
        SMLAL2 V14.4S, v23.8h, V0.H[3]
        .if \w == 8 || \w == 16 || \w == 32
            ld1 {v20.8h-v23.8h}, [x7], #64
            
            SMLAL  V11.4S, v20.4h, V0.H[5]
            SMLAL2 V12.4S, v20.8h, V0.H[5]
            SMLAL  V13.4S, v21.4h, V0.H[5]
            SMLAL2 V14.4S, v21.8h, V0.H[5]

            SMLAL  V11.4S, v22.4h, V0.H[7]
            SMLAL2 V12.4S, v22.8h, V0.H[7]
            SMLAL  V13.4S, v23.4h, V0.H[7]
            SMLAL2 V14.4S, v23.8h, V0.H[7]
        .endif
        .if \w == 16 || \w == 32
            ld1 {v20.8h-v23.8h}, [x7], #64 

            SMLAL  V11.4S, v20.4h, V1.H[1]
            SMLAL2 V12.4S, v20.8h, V1.H[1]
            SMLAL  V13.4S, v21.4h, V1.H[1]
            SMLAL2 V14.4S, v21.8h, V1.H[1]

            SMLAL  V11.4S, v22.4h, V1.H[3]
            SMLAL2 V12.4S, v22.8h, V1.H[3]
            SMLAL  V13.4S, v23.4h, V1.H[3]
            SMLAL2 V14.4S, v23.8h, V1.H[3]

            ld1 {v20.8h-v23.8h}, [x7], #64
            
            SMLAL  V11.4S, v20.4h, V1.H[5]
            SMLAL2 V12.4S, v20.8h, V1.H[5]
            SMLAL  V13.4S, v21.4h, V1.H[5]
            SMLAL2 V14.4S, v21.8h, V1.H[5]

            SMLAL  V11.4S, v22.4h, V1.H[7]
            SMLAL2 V12.4S, v22.8h, V1.H[7]
            SMLAL  V13.4S, v23.4h, V1.H[7]
            SMLAL2 V14.4S, v23.8h, V1.H[7]
        .endif
        .if \w == 32
            ld1 {v20.8h-v23.8h}, [x7], #64
            
            SMLAL  V11.4S, v20.4h, V2.H[1]
            SMLAL2 V12.4S, v20.8h, V2.H[1]
            SMLAL  V13.4S, v21.4h, V2.H[1]
            SMLAL2 V14.4S, v21.8h, V2.H[1]

            SMLAL  V11.4S, v22.4h, V2.H[3]
            SMLAL2 V12.4S, v22.8h, V2.H[3]
            SMLAL  V13.4S, v23.4h, V2.H[3]
            SMLAL2 V14.4S, v23.8h, V2.H[3]

            ld1 {v20.8h-v23.8h}, [x7], #64
            
            SMLAL  V11.4S, v20.4h, V2.H[5]
            SMLAL2 V12.4S, v20.8h, V2.H[5]
            SMLAL  V13.4S, v21.4h, V2.H[5]
            SMLAL2 V14.4S, v21.8h, V2.H[5]

            SMLAL  V11.4S, v22.4h, V2.H[7]
            SMLAL2 V12.4S, v22.8h, V2.H[7]
            SMLAL  V13.4S, v23.4h, V2.H[7]
            SMLAL2 V14.4S, v23.8h, V2.H[7]

            ld1 {v20.8h-v23.8h}, [x7], #64
            
            SMLAL  V11.4S, v20.4h, V3.H[1]
            SMLAL2 V12.4S, v20.8h, V3.H[1]
            SMLAL  V13.4S, v21.4h, V3.H[1]
            SMLAL2 V14.4S, v21.8h, V3.H[1]

            SMLAL  V11.4S, v22.4h, V3.H[3]
            SMLAL2 V12.4S, v22.8h, V3.H[3]
            SMLAL  V13.4S, v23.4h, V3.H[3]
            SMLAL2 V14.4S, v23.8h, V3.H[3]

            ld1 {v20.8h-v23.8h}, [x7]
            
            SMLAL  V11.4S, v20.4h, V3.H[5]
            SMLAL2 V12.4S, v20.8h, V3.H[5]
            SMLAL  V13.4S, v21.4h, V3.H[5]
            SMLAL2 V14.4S, v21.8h, V3.H[5]

            SMLAL  V11.4S, v22.4h, V3.H[7]
            SMLAL2 V12.4S, v22.8h, V3.H[7]
            SMLAL  V13.4S, v23.4h, V3.H[7]
            SMLAL2 V14.4S, v23.8h, V3.H[7]
        .endif    
        // thid out butterfly (add sub)
        ADD v20.4s , v7.4s , v11.4s
        ADD v21.4s , v8.4s , v12.4s
        ADD v22.4s , v10.4s, v13.4s
        ADD v23.4s , v9.4s , v14.4s

        SUB v27.4s , v7.4s , v11.4s
        SUB v26.4s , v8.4s , v12.4s
        SUB v25.4s , v10.4s, v13.4s
        SUB v24.4s , v9.4s , v14.4s

        // shift and round 
        ADD   V20.4S , v20.4S , V29.4S
        ADD   V21.4S , V21.4S , V29.4S
        ADD   V22.4S , V22.4S , V29.4S
        ADD   V23.4S , V23.4S , V29.4S

        ADD   V24.4S, V24.4S , V29.4S
        ADD   V25.4S, V25.4S , V29.4S
        ADD   V26.4S, V26.4S , V29.4S
        ADD   V27.4S, V27.4S , V29.4S

        SSHL  V20.4S , V20.4S , V28.4S 
        SSHL  V21.4S , V21.4S , V28.4S
        SSHL  V22.4S , V22.4S , V28.4S 
        SSHL  V23.4S , V23.4S , V28.4S

        SSHL  V24.4S , V24.4S , V28.4S 
        SSHL  V25.4S , V25.4S , V28.4S
        SSHL  V26.4S , V26.4S , V28.4S 
        SSHL  V27.4S , V27.4S , V28.4S
        
        inverse_neon v24
        inverse_neon V25
        inverse_neon V26
        inverse_neon V27

        UZP1  v0.8h, V20.8H , V21.8H
        UZP1  v1.8h, V22.8H , V23.8H
        UZP1  v2.8h, V24.8H , V25.8H
        UZP1  v3.8h, V26.8H , V27.8H

        st1 {v0.8h-v3.8h}, [x1]
        ret lr
    endfunc
.endm

.macro idct_x_16, red
    function idct_x_16_\red\()_neon, export=1

        dup v28.4s, w3
        dup v29.4s, w4 


        ld1 {v0.8h}, [x0], x2
        ld1 {v1.8h}, [x0], x2
        ld1 {v2.8h}, [x0], x2
        ld1 {v3.8h}, [x0], x2
        .if \red == 8 || \red == 16
            ld1 {v4.8h}, [x0], x2
            ld1 {v5.8h}, [x0], x2
            ld1 {v6.8h}, [x0], x2
            ld1 {v7.8h}, [x0], x2
        .endif
        .if \red == 16
            ld1 {v8.8h}, [x0], x2
            ld1 {v9.8h}, [x0], x2
            ld1 {v10.8h}, [x0], x2
            ld1 {v11.8h}, [x0], x2
            ld1 {v12.8h}, [x0], x2
            ld1 {v13.8h}, [x0], x2
            ld1 {v14.8h}, [x0], x2
            ld1 {v15.8h}, [x0]
        .endif

        ld1 {v21.8h-v24.8h}, [x5], #64

        SMULL  V17.4S, v21.4h, v0.H[0]
        SMULL2 V18.4S, v21.8h, v0.H[0]
        SMULL  V19.4S, v22.4h, v0.H[0]
        SMULL2 V20.4S, v22.8h, v0.H[0]

        SMLAL  V17.4S, v23.4h, V1.H[0]
        SMLAL2 V18.4S, v23.8h, V1.H[0]
        SMLAL  V19.4S, v24.4h, V1.H[0]
        SMLAL2 V20.4S, v24.8h, V1.H[0]

        ld1 {v21.8h-v24.8h}, [x5], #64

        SMLAL  V17.4S, v21.4h, V2.H[0]
        SMLAL2 V18.4S, v21.8h, V2.H[0]
        SMLAL  V19.4S, v22.4h, V2.H[0]
        SMLAL2 V20.4S, v22.8h, V2.H[0]

        SMLAL  V17.4S, v23.4h, V3.H[0]
        SMLAL2 V18.4S, v23.8h, V3.H[0]
        SMLAL  V19.4S, v24.4h, V3.H[0]
        SMLAL2 V20.4S, v24.8h, V3.H[0]

        .if \red == 8 || \red == 16
            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V4.H[0]
            SMLAL2 V18.4S, v21.8h, V4.H[0]
            SMLAL  V19.4S, v22.4h, V4.H[0]
            SMLAL2 V20.4S, v22.8h, V4.H[0]

            SMLAL  V17.4S, v23.4h, V5.H[0]
            SMLAL2 V18.4S, v23.8h, V5.H[0]
            SMLAL  V19.4S, v24.4h, V5.H[0]
            SMLAL2 V20.4S, v24.8h, V5.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V6.H[0]
            SMLAL2 V18.4S, v21.8h, V6.H[0]
            SMLAL  V19.4S, v22.4h, V6.H[0]
            SMLAL2 V20.4S, v22.8h, V6.H[0]

            SMLAL  V17.4S, v23.4h, V7.H[0]
            SMLAL2 V18.4S, v23.8h, V7.H[0]
            SMLAL  V19.4S, v24.4h, V7.H[0]
            SMLAL2 V20.4S, v24.8h, V7.H[0]
        .endif
        .if \red == 16
            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V8.H[0]
            SMLAL2 V18.4S, v21.8h, V8.H[0]
            SMLAL  V19.4S, v22.4h, V8.H[0]
            SMLAL2 V20.4S, v22.8h, V8.H[0]

            SMLAL  V17.4S, v23.4h, V9.H[0]
            SMLAL2 V18.4S, v23.8h, V9.H[0]
            SMLAL  V19.4S, v24.4h, V9.H[0]
            SMLAL2 V20.4S, v24.8h, V9.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V10.H[0]
            SMLAL2 V18.4S, v21.8h, V10.H[0]
            SMLAL  V19.4S, v22.4h, V10.H[0]
            SMLAL2 V20.4S, v22.8h, V10.H[0]

            SMLAL  V17.4S, v23.4h, V11.H[0]
            SMLAL2 V18.4S, v23.8h, V11.H[0]
            SMLAL  V19.4S, v24.4h, V11.H[0]
            SMLAL2 V20.4S, v24.8h, V11.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V12.H[0]
            SMLAL2 V18.4S, v21.8h, V12.H[0]
            SMLAL  V19.4S, v22.4h, V12.H[0]
            SMLAL2 V20.4S, v22.8h, V12.H[0]

            SMLAL  V17.4S, v23.4h, V13.H[0]
            SMLAL2 V18.4S, v23.8h, V13.H[0]
            SMLAL  V19.4S, v24.4h, V13.H[0]
            SMLAL2 V20.4S, v24.8h, V13.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V14.H[0]
            SMLAL2 V18.4S, v21.8h, V14.H[0]
            SMLAL  V19.4S, v22.4h, V14.H[0]
            SMLAL2 V20.4S, v22.8h, V14.H[0]

            SMLAL  V17.4S, v23.4h, V15.H[0]
            SMLAL2 V18.4S, v23.8h, V15.H[0]
            SMLAL  V19.4S, v24.4h, V15.H[0]
            SMLAL2 V20.4S, v24.8h, V15.H[0]
        .endif
        // shift and round 

        ADD   V17.4S, v17.4S, V29.4S
        ADD   V18.4S, V18.4S, V29.4S
        ADD   V19.4S, V19.4S, V29.4S
        ADD   V20.4S, V20.4S, V29.4S


        SSHL  V17.4S, V17.4S, V28.4S 
        SSHL  V18.4S, V18.4S, V28.4S
        SSHL  V19.4S, V19.4S, V28.4S 
        SSHL  V20.4S, V20.4S, V28.4S

        UZP1  v9.8h , V17.8H, V18.8H
        UZP1  v10.8h, V19.8H, V20.8H

        st1 {v9.8h-v10.8h}, [x1]

        ret lr
    endfunc
.endm

.macro idct_x_32, red
    function idct_x_32_\red\()_neon, export=1

        dup v30.4s, w3
        dup v29.4s, w4 


        ld1 {v0.8h}, [x0], x2
        ld1 {v1.8h}, [x0], x2
        ld1 {v2.8h}, [x0], x2
        ld1 {v3.8h}, [x0], x2
        .if \red == 8 || \red == 16
            ld1 {v4.8h}, [x0], x2
            ld1 {v5.8h}, [x0], x2
            ld1 {v6.8h}, [x0], x2
            ld1 {v7.8h}, [x0], x2
        .endif
        .if \red == 16
            ld1 {v8.8h}, [x0], x2
            ld1 {v9.8h}, [x0], x2
            ld1 {v10.8h}, [x0], x2
            ld1 {v11.8h}, [x0], x2
            ld1 {v12.8h}, [x0], x2
            ld1 {v13.8h}, [x0], x2
            ld1 {v14.8h}, [x0], x2
            ld1 {v15.8h}, [x0]
        .endif

        ld1 {v21.8h-v24.8h}, [x5], #64

        SMULL  V17.4S, v21.4h, v0.H[0]
        SMULL2 V18.4S, v21.8h, v0.H[0]
        SMULL  V19.4S, v22.4h, v0.H[0]
        SMULL2 V20.4S, v22.8h, v0.H[0]
        SMLAL  V25.4S, v23.4h, V0.H[0]
        SMLAL2 V26.4S, v23.8h, V0.H[0]
        SMLAL  V27.4S, v24.4h, V0.H[0]
        SMLAL2 V28.4S, v24.8h, V0.H[0]

        ld1 {v21.8h-v24.8h}, [x5], #64

        SMLAL  V17.4S, v21.4h, V1.H[0]
        SMLAL2 V18.4S, v21.8h, V1.H[0]
        SMLAL  V19.4S, v22.4h, V1.H[0]
        SMLAL2 V20.4S, v22.8h, V1.H[0]
        SMLAL  V25.4S, v23.4h, V1.H[0]
        SMLAL2 V26.4S, v23.8h, V1.H[0]
        SMLAL  V27.4S, v24.4h, V1.H[0]
        SMLAL2 V28.4S, v24.8h, V1.H[0]

        ld1 {v21.8h-v24.8h}, [x5], #64

        SMLAL  V17.4S, v21.4h, V2.H[0]
        SMLAL2 V18.4S, v21.8h, V2.H[0]
        SMLAL  V19.4S, v22.4h, V2.H[0]
        SMLAL2 V20.4S, v22.8h, V2.H[0]
        SMLAL  V25.4S, v23.4h, V2.H[0]
        SMLAL2 V26.4S, v23.8h, V2.H[0]
        SMLAL  V27.4S, v24.4h, V2.H[0]
        SMLAL2 V28.4S, v24.8h, V2.H[0]

        ld1 {v21.8h-v24.8h}, [x5], #64

        SMLAL  V17.4S, v21.4h, V3.H[0]
        SMLAL2 V18.4S, v21.8h, V3.H[0]
        SMLAL  V19.4S, v22.4h, V3.H[0]
        SMLAL2 V20.4S, v22.8h, V3.H[0]
        SMLAL  V25.4S, v23.4h, V3.H[0]
        SMLAL2 V26.4S, v23.8h, V3.H[0]
        SMLAL  V27.4S, v24.4h, V3.H[0]
        SMLAL2 V28.4S, v24.8h, V3.H[0]

        .if \red == 8 || \red == 16
            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V4.H[0]
            SMLAL2 V18.4S, v21.8h, V4.H[0]
            SMLAL  V19.4S, v22.4h, V4.H[0]
            SMLAL2 V20.4S, v22.8h, V4.H[0]
            SMLAL  V25.4S, v23.4h, V4.H[0]
            SMLAL2 V26.4S, v23.8h, V4.H[0]
            SMLAL  V27.4S, v24.4h, V4.H[0]
            SMLAL2 V28.4S, v24.8h, V4.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V5.H[0]
            SMLAL2 V18.4S, v21.8h, V5.H[0]
            SMLAL  V19.4S, v22.4h, V5.H[0]
            SMLAL2 V20.4S, v22.8h, V5.H[0]
            SMLAL  V25.4S, v23.4h, V5.H[0]
            SMLAL2 V26.4S, v23.8h, V5.H[0]
            SMLAL  V27.4S, v24.4h, V5.H[0]
            SMLAL2 V28.4S, v24.8h, V5.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V6.H[0]
            SMLAL2 V18.4S, v21.8h, V6.H[0]
            SMLAL  V19.4S, v22.4h, V6.H[0]
            SMLAL2 V20.4S, v22.8h, V6.H[0]
            SMLAL  V25.4S, v23.4h, V6.H[0]
            SMLAL2 V26.4S, v23.8h, V6.H[0]
            SMLAL  V27.4S, v24.4h, V6.H[0]
            SMLAL2 V28.4S, v24.8h, V6.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V7.H[0]
            SMLAL2 V18.4S, v21.8h, V7.H[0]
            SMLAL  V19.4S, v22.4h, V7.H[0]
            SMLAL2 V20.4S, v22.8h, V7.H[0]
            SMLAL  V25.4S, v23.4h, V7.H[0]
            SMLAL2 V26.4S, v23.8h, V7.H[0]
            SMLAL  V27.4S, v24.4h, V7.H[0]
            SMLAL2 V28.4S, v24.8h, V7.H[0]
        .endif
        .if \red == 16
            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V8.H[0]
            SMLAL2 V18.4S, v21.8h, V8.H[0]
            SMLAL  V19.4S, v22.4h, V8.H[0]
            SMLAL2 V20.4S, v22.8h, V8.H[0]
            SMLAL  V25.4S, v23.4h, V8.H[0]
            SMLAL2 V26.4S, v23.8h, V8.H[0]
            SMLAL  V27.4S, v24.4h, V8.H[0]
            SMLAL2 V28.4S, v24.8h, V8.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V9.H[0]
            SMLAL2 V18.4S, v21.8h, V9.H[0]
            SMLAL  V19.4S, v22.4h, V9.H[0]
            SMLAL2 V20.4S, v22.8h, V9.H[0]
            SMLAL  V25.4S, v23.4h, V9.H[0]
            SMLAL2 V26.4S, v23.8h, V9.H[0]
            SMLAL  V27.4S, v24.4h, V9.H[0]
            SMLAL2 V28.4S, v24.8h, V9.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V10.H[0]
            SMLAL2 V18.4S, v21.8h, V10.H[0]
            SMLAL  V19.4S, v22.4h, V10.H[0]
            SMLAL2 V20.4S, v22.8h, V10.H[0]
            SMLAL  V25.4S, v23.4h, V10.H[0]
            SMLAL2 V26.4S, v23.8h, V10.H[0]
            SMLAL  V27.4S, v24.4h, V10.H[0]
            SMLAL2 V28.4S, v24.8h, V10.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V11.H[0]
            SMLAL2 V18.4S, v21.8h, V11.H[0]
            SMLAL  V19.4S, v22.4h, V11.H[0]
            SMLAL2 V20.4S, v22.8h, V11.H[0]
            SMLAL  V25.4S, v23.4h, V11.H[0]
            SMLAL2 V26.4S, v23.8h, V11.H[0]
            SMLAL  V27.4S, v24.4h, V11.H[0]
            SMLAL2 V28.4S, v24.8h, V11.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V12.H[0]
            SMLAL2 V18.4S, v21.8h, V12.H[0]
            SMLAL  V19.4S, v22.4h, V12.H[0]
            SMLAL2 V20.4S, v22.8h, V12.H[0]
            SMLAL  V25.4S, v23.4h, V12.H[0]
            SMLAL2 V26.4S, v23.8h, V12.H[0]
            SMLAL  V27.4S, v24.4h, V12.H[0]
            SMLAL2 V28.4S, v24.8h, V12.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V13.H[0]
            SMLAL2 V18.4S, v21.8h, V13.H[0]
            SMLAL  V19.4S, v22.4h, V13.H[0]
            SMLAL2 V20.4S, v22.8h, V13.H[0]
            SMLAL  V25.4S, v23.4h, V13.H[0]
            SMLAL2 V26.4S, v23.8h, V13.H[0]
            SMLAL  V27.4S, v24.4h, V13.H[0]
            SMLAL2 V28.4S, v24.8h, V13.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V14.H[0]
            SMLAL2 V18.4S, v21.8h, V14.H[0]
            SMLAL  V19.4S, v22.4h, V14.H[0]
            SMLAL2 V20.4S, v22.8h, V14.H[0]
            SMLAL  V25.4S, v23.4h, V14.H[0]
            SMLAL2 V26.4S, v23.8h, V14.H[0]
            SMLAL  V27.4S, v24.4h, V14.H[0]
            SMLAL2 V28.4S, v24.8h, V14.H[0]

            ld1 {v21.8h-v24.8h}, [x5], #64

            SMLAL  V17.4S, v21.4h, V15.H[0]
            SMLAL2 V18.4S, v21.8h, V15.H[0]
            SMLAL  V19.4S, v22.4h, V15.H[0]
            SMLAL2 V20.4S, v22.8h, V15.H[0]
            SMLAL  V25.4S, v23.4h, V15.H[0]
            SMLAL2 V26.4S, v23.8h, V15.H[0]
            SMLAL  V27.4S, v24.4h, V15.H[0]
            SMLAL2 V28.4S, v24.8h, V15.H[0]
        .endif
        // shift and round 

        ADD   V17.4S, v17.4S, V29.4S
        ADD   V18.4S, V18.4S, V29.4S
        ADD   V19.4S, V19.4S, V29.4S
        ADD   V20.4S, V20.4S, V29.4S
        ADD   V25.4S, v25.4S, V29.4S
        ADD   V26.4S, V26.4S, V29.4S
        ADD   V27.4S, V27.4S, V29.4S
        ADD   V28.4S, V28.4S, V29.4S

        SSHL  V17.4S, V17.4S, v30.4S
        SSHL  V18.4S, V18.4S, v30.4S
        SSHL  V19.4S, V19.4S, v30.4S
        SSHL  V20.4S, V20.4S, v30.4S
        SSHL  V25.4S, V25.4S, v30.4S
        SSHL  V26.4S, V26.4S, v30.4S
        SSHL  V27.4S, V27.4S, v30.4S
        SSHL  V28.4S, V28.4S, v30.4S


        UZP1  v9.8h , V17.8H, V18.8H
        UZP1  v10.8h, V19.8H, V20.8H
        UZP1  v11.8h, V25.8H, V26.8H
        UZP1  v12.8h, V27.8H, V28.8H

        st1 {v9.8h-v10.8h}, [x1]

        ret lr
    endfunc
.endm

// store data for DC
store_fun 8
store_fun 16 
store_fun 32

// 4-pt dct-2/8 and dst-7
idct_x_4 2
idct_x_4 4
idct_x_4 8

// idct-ii butterfly 
idct_ii_8 2 8
idct_ii_8 2 4
idct_ii_8 4 8
idct_ii_8 4 4

// 8-pt dxt-2
idct_x_8 2 4
idct_x_8 2 8
idct_x_8 4 4
idct_x_8 4 8
idct_x_8 8 4
idct_x_8 8 8

// 16-pt dct-2
idct_ii_16_2 4
idct_ii_16_2 8
idct_ii_16_2 16

// 16-pt dxt 7/8 
idct_x_16 4 
idct_x_16 8 
idct_x_16 16

// 32-pt DCT2
idct_ii_32 4
idct_ii_32 8
idct_ii_32 16
idct_ii_32 32

// 32-pt dxt 7/8 
idct_x_32 4 
idct_x_32 8 
idct_x_32 16

