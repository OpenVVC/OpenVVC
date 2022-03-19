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

.macro  movrel rd, val, offset=0
#if defined(__APPLE__)
  .if \offset < 0
        adrp            \rd, \val@PAGE
        add             \rd, \rd, \val@PAGEOFF
        sub             \rd, \rd, -(\offset)
  .else
        adrp            \rd, \val+(\offset)@PAGE
        add             \rd, \rd, \val+(\offset)@PAGEOFF
  .endif
#elif defined(PIC) && defined(_WIN32)
  .if \offset < 0
        adrp            \rd, \val
        add             \rd, \rd, :lo12:\val
        sub             \rd, \rd, -(\offset)
  .else
        adrp            \rd, \val+(\offset)
        add             \rd, \rd, :lo12:\val+(\offset)
  .endif
#elif defined(PIC)
        adrp            \rd, \val+(\offset)
        add             \rd, \rd, :lo12:\val+(\offset)
#else
        ldr             \rd, =\val+\offset
#endif
.endm


.macro UCLIP, val, min_val, max_val
    UMAX  \val, \val, \min_val
    UMIN  \val, \val, \max_val
.endm

.macro SCLIP, val, min_val, max_val
    SMAX  \val, \val, \min_val
    SMIN  \val, \val, \max_val
.endm

.macro INIT_CLIP_VECTORS, suffix
    MOV   w9, #1023
    DUP   V30\suffix, w9
    DUP   V31\suffix, wzr
.endm

.macro INPUT_SCALING, w
    // shift stride
    LSL x3, x3, #1
    .if \w == 64
        SUB x3, x3, #64
    .endif
    // src_above + 1
    ADD x0, x0, #2
    // src_left + 1
    ADD x1, x1, #2
.endm

.macro LOAD_REF_4,  ref, v, s1=.4H
    LD1   {\v\()0\s1}, [\ref]
.endm

.macro LOAD_REF_8,  ref, v, s1=.8H
    LD1   {\v\()0\s1}, [\ref]
.endm

.macro LOAD_REF_16, ref, v, s1=.8H
    LD1   {\v\()0\s1-\v\()1\s1}, [\ref]
.endm

.macro LOAD_REF_32, ref, v, s1=.8H
    LD1   {\v\()0\s1-\v\()3\s1}, [\ref]
.endm

.macro LOAD_REF_64, ref, v, s1=.8H, stride=#64
    LD1   {\v\()0\s1-\v\()3\s1}, [\ref], \stride
    LD1   {\v\()4\s1-\v\()7\s1}, [\ref]
    SUB \ref, \ref, \stride
.endm

.macro STORE_4,  v, ref, stride, s1=.4H, pdpc=1, scale=0, enable_stride=1
    .if \enable_stride
        ST1   {\v\()0\s1}, [\ref], \stride
    .else
        ST1   {\v\()0\s1}, [\ref]
    .endif
.endm

.macro STORE_8,  v, ref, stride, s1=.8H, pdpc=1, scale=0, enable_stride=1
    .if \enable_stride
        ST1   {\v\()0\s1}, [\ref], \stride
    .else
        ST1   {\v\()0\s1}, [\ref]
    .endif
.endm

.macro STORE_16, v, ref, stride, s1=.8H, pdpc=1, scale=0, internstride=#16, enable_stride=1
    .if \pdpc == 0 && \scale < 2
        ST1   {\v\()0\s1}, [\ref], \internstride
            .if \enable_stride
                ST1   {\v\()2\s1}, [\ref], \stride
            .else
                ST1   {\v\()2\s1}, [\ref]
            .endif
    .else
        .if \enable_stride
            ST1   {\v\()0\s1-\v\()1\s1}, [\ref], \stride
        .else
            ST1   {\v\()0\s1-\v\()1\s1}, [\ref]
        .endif
    .endif
.endm

.macro STORE_32, v, ref, stride, s1=.8H, pdpc=1, scale=0, internstride=#16, enable_stride=1
    .if \pdpc == 0
        .if \scale == 2
            ST1   {\v\()0\s1}, [\ref], \internstride
            ST1   {\v\()1\s1}, [\ref], \internstride
            ST1   {\v\()2\s1}, [\ref], \internstride
            .if \enable_stride
                ST1   {\v\()2\s1}, [\ref], \stride
            .else
                ST1   {\v\()2\s1}, [\ref]
            .endif
        .else
            ST1   {\v\()0\s1}, [\ref], \internstride
            ST1   {\v\()2\s1}, [\ref], \internstride
            ST1   {\v\()2\s1}, [\ref], \internstride
            .if \enable_stride
                ST1   {\v\()2\s1}, [\ref], \stride
            .else
                ST1   {\v\()2\s1}, [\ref]
            .endif
        .endif
    .else
        .if \enable_stride
            ST1   {\v\()0\s1-\v\()3\s1}, [\ref], \stride
        .else
            ST1   {\v\()0\s1-\v\()3\s1}, [\ref]
        .endif
    .endif
.endm

.macro STORE_64, v, ref, stride, s1=.8H, pdpc=1, scale=0, internstride=#64, internstride2=#16, enable_stride=1
    .if \pdpc == 0
        .if \scale == 2
            ST1   {\v\()0\s1}, [\ref], \internstride2
            ST1   {\v\()1\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            .if \enable_stride
                ST1   {\v\()2\s1}, [\ref], \stride
            .else
                ST1   {\v\()2\s1}, [\ref]
            .endif
        .else
            ST1   {\v\()0\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            ST1   {\v\()2\s1}, [\ref], \internstride2
            .if \enable_stride
                ST1   {\v\()2\s1}, [\ref], \stride
            .else
                ST1   {\v\()2\s1}, [\ref]
            .endif
        .endif
    .else
        ST1   {\v\()0\s1-\v\()3\s1}, [\ref], \internstride
        .if \enable_stride
            ST1   {\v\()4\s1-\v\()7\s1}, [\ref], \stride
        .else
            ST1   {\v\()4\s1-\v\()7\s1}, [\ref]
        .endif
    .endif
.endm

.macro ADD_REF_16, out, in, s1=.8H
    ADD   \out\s1, \in\()0\s1, \in\()1\s1
.endm

.macro ADD_REF_32, out, in, s1=.8H
    ADD   V0\s1, \in\()0\s1, \in\()2\s1
    ADD   V1\s1, \in\()1\s1, \in\()3\s1
    ADD   \out\s1, V0\s1, V1\s1
.endm

.macro ADD_REF_64, out, in, s1=.8H
    ADD   V0\s1, \in\()0\s1, \in\()4\s1
    ADD   V1\s1, \in\()1\s1, \in\()5\s1
    ADD   V2\s1, \in\()2\s1, \in\()6\s1
    ADD   V3\s1, \in\()3\s1, \in\()7\s1
    ADD   V0\s1, V0\s1, V2\s1
    ADD   V1\s1, V1\s1, V3\s1
    ADD   \out\s1, V0\s1, V1\s1
.endm

.macro COMPUTE_DC, w, h, out, offset, shift, s1=.8H
    LOAD_REF_\w x0, V2
    .if \w >= \h
        .if \w > 8
            ADD_REF_\w V12, V2
        .else
            MOV V12\s1, V20\s1
        .endif
    .endif
    .if \w <= \h
            LOAD_REF_\h x1, V
        .if \h > 8
            ADD_REF_\h V0, V
        .endif
    .endif
    .if \w == \h 
        ADD   V0\s1, V0\s1, V12\s1
    .endif
    .if \w > \h
        UADDLV  S0, V12\s1
    .else
        UADDLV  S0, V0\s1
    .endif
    UMOV  w\out, V0.S[0]
    ADD   w\out, w\out, \offset
    LSR   w\out, w\out, \shift
    DUP   V\out\s1, w\out
.endm
    

.macro INIT_SCALE, scale, suffix
    //Load PDPC
    movrel x9, PDPC_SCALE\scale
    .if \scale < 2
        LD1   {V28\suffix}, [x9]
    .else
        LD1   {V28\suffix-V29\suffix}, [x9]
    .endif
    USHL  V14\suffix, V12\suffix, V28\suffix
    .if \scale == 2
        USHL  V15\suffix, V12\suffix, V29\suffix
    .endif
    //add_v = 32
    //    __m128i tst = _mm_slli_epi16(dc_v, 6);               
    //    tst = _mm_adds_epu16(tst, add_v);                    
    LSL w10, w12, #6
    ADD w10, w10, #32
    DUP V17\suffix, w10
    //    __m128i tst0 = _mm_subs_epu16(tst, w_y);
    UQSUB V18\suffix, V17\suffix, V14\suffix
    .if \scale == 2
        UQSUB  V19\suffix, V17\suffix, V15\suffix
    .endif
.endm


.macro LOAD_PDPC, s1
    // int y_wgh = pdpc_sh[y];
    LDRSH  w14, [x9], #2
    DUP   V11\s1, w14
.endm

.macro LOAD_REF_LEFT, s1
    // const int16_t l_val = src_left[y + 1];              
    // l_v  = _mm_set1_epi16(l_val);                       
    LD1R {V16\s1}, [x1], #2
.endm

.macro COMPUTE_WX, s1
        // w_x = _mm_slli_epi16(dc_v, y_wgh);                  
        LSL  x13, x12, x14
        DUP  V13\s1, w13
.endm

.macro dc_pdpc, w, h, scale, offset, shift, loop1, loop2=0, s1=.8H, s2=.8H
    function dc_pdpc_w\w\()_h\h\()_neon, export=1
        INPUT_SCALING \w
        // pb_h compute
        INIT_CLIP_VECTORS \s1
        //compute DC_VAL
        COMPUTE_DC \w, \h, 12, \offset, \shift, \s2
        INIT_SCALE \scale, \s1
    .rept \loop1
        LOAD_PDPC \s1
        LOAD_REF_LEFT \s1
        COMPUTE_WX \s1
        // tst = _mm_subs_epu16(tst0, w_x);                    
        UQSUB V8\s1, V18\s1, V13\s1
        .if \scale == 2
            UQSUB V9\s1, V19\s1, V13\s1
        .endif
        .if (\scale == 2 && \w > 16) || (\scale < 2 && \w > 8)
            UQSUB V10\s1, V17\s1, V13\s1
        .endif
        // xl = _mm_mullo_epi16(l_v, x_v);                     
        USHL  V14\s1, V16\s1, V28\s1
        .if \scale == 2
            USHL  V15\s1, V16\s1, V29\s1
        .endif
        // yt = _mm_slli_epi16(a1, y_wgh);
        USHL V0\s1, V20\s1, V11\s1
        .if \w >= 16
            USHL V1\s1, V21\s1, V11\s1
        .endif
        .if \w >= 32
            USHL V2\s1, V22\s1, V11\s1
            USHL V3\s1, V23\s1, V11\s1
        .endif
        .if \w >= 64
            USHL V4\s1, V24\s1, V11\s1
            USHL V5\s1, V25\s1, V11\s1
            USHL V6\s1, V26\s1, V11\s1
            USHL V7\s1, V27\s1, V11\s1
        .endif
        // pdpc_rnd = _mm_adds_epu16(xl, yt);                  
        UQADD V0\s1, V0\s1,V14\s1
        .if \scale == 2
            UQADD V1\s1, V1\s1,V15\s1
        .endif
        // tst = _mm_adds_epu16(tst, pdpc_rnd);                
        UQADD V0\s1, V0\s1,V8\s1
        .if \w >= 16
            .if \scale == 2
                UQADD V1\s1, V1\s1, V9\s1
            .else
                UQADD V1\s1, V1\s1, V10\s1
            .endif
        .endif
        .if \w >= 32
            UQADD V2\s1, V2\s1, V10\s1
            UQADD V3\s1, V3\s1, V10\s1
        .endif
        .if \w >= 64
            UQADD V4\s1, V4\s1, V10\s1
            UQADD V5\s1, V5\s1, V10\s1
            UQADD V6\s1, V6\s1, V10\s1
            UQADD V7\s1, V7\s1, V10\s1
        .endif

        USHR V0\s1, V0\s1, #6
        .if \w >= 16
            USHR V1\s1, V1\s1, #6
        .endif
        .if \w >= 32
            USHR V2\s1, V2\s1, #6
            USHR V3\s1, V3\s1, #6
        .endif
        .if \w >= 64
            USHR V4\s1, V4\s1, #6
            USHR V5\s1, V5\s1, #6
            USHR V6\s1, V6\s1, #6
            USHR V7\s1, V7\s1, #6
        .endif

        UCLIP  V0\s1, V31\s1, V30\s1
        .if \w >= 16
            UCLIP V1\s1, V31\s1, V30\s1
        .endif
        .if \w >= 32
            UCLIP V2\s1, V31\s1, V30\s1
            UCLIP V3\s1, V31\s1, V30\s1
        .endif
        .if \w >= 64
            UCLIP V4\s1, V31\s1, V30\s1
            UCLIP V5\s1, V31\s1, V30\s1
            UCLIP V6\s1, V31\s1, V30\s1
            UCLIP V7\s1, V31\s1, V30\s1
        .endif

        STORE_\w  V, x2, x3, \s1
    .endr

    .if \loop2 > 0
        .if \w == 16 && \scale < 2
            SUB x3, x3, #16
        .endif
        .if \w >= 32
            SUB x3, x3, #48
        .endif
        .rept \loop2
            LOAD_REF_LEFT \s1
            USHL  V14\s1, V16\s1, V28\s1
            .if \scale == 2
                USHL  V15\s1, V16\s1, V29\s1
            .endif               
            UQADD V0\s1, V14\s1, V18\s1
            .if \scale == 2
                UQADD V1\s1, V15\s1,V19\s1
            .endif

            USHR V0\s1, V0\s1, #6
            .if \w >= 16 && \scale == 2
                USHR V1\s1, V1\s1, #6
            .endif
            .if (\scale == 2 && \w > 16) || (\scale < 2 && \w > 8)
                USHR V2\s1, V17\s1, #6
            .endif

            UCLIP  V0\s1, V31\s1, V30\s1
            .if \w >= 16 && \scale == 2
                UCLIP V1\s1, V31\s1, V30\s1
            .endif
            .if (\scale == 2 && \w > 16) || (\scale < 2 && \w > 8)
                UCLIP V2\s1, V31\s1, V30\s1
            .endif

            STORE_\w  V, x2, x3, \s1, 0, \scale
        .endr
    .endif
        RET   lr
    endfunc
.endm

.macro COMPUTE_PLANAR_VECTOR w, max_scale rc_scale, tr_scale, s1, s2, s3, s4=.4H
    // Load above and tr_v
    LD1 {V8\s1}, [x6], jump
    LD1 {V23\s1}, [x7]
    .if ! first_it
        //rcol_mul_v = _mm_add_epi16(rcol_mul_v, _mm_slli_epi16(rcol_v, 3));
        ADD V9\s1, V9\s1, V10\s1
    .endif
    .set first_it, 0
    //__m128i rcol_x1 = _mm_add_epi16(lcol_v, rcol_mul_v);
    ADD V11\s1, V22\s1, V9\s1
    // tr_v1 = _mm_sub_epi16(tr_v1, a1);
    SUB V23\s1, V23\s1, V8\s1 
    // tr_v1 = _mm_add_epi16(tr_v1, bl_val_v);
    ADD V23\s1, V23\s1, V27\s1
    ST1 {V23\s1}, [x7], jump
    .if is64
        // __m128i tr_v_hi1 = _mm_unpackhi_epi16(tr_v1, _mm_setzero_si128());
        .if \w > 4
            UXTL2  V13\s3, V23\s1
        .endif
        // __m128i tr_v_lo1 = _mm_unpacklo_epi16(tr_v1, _mm_setzero_si128());
        UXTL  V23\s3, V23\s4
        // __m128i rc_v_hi1 = _mm_unpackhi_epi16(rcol_x1, _mm_setzero_si128());
        .if \w > 4
            UXTL2  V12\s3, V11\s1
        .endif
        // __m128i rc_v_lo1 = _mm_unpacklo_epi16(rcol_x1, _mm_setzero_si128());
        UXTL  V11\s3, V11\s4
    .endif
    .if \rc_scale
        // __m128i src_v_lo1 = _mm_slli_epi32(rc_v_lo1, rc_scale);
        USHL V11\s3, V11\s3, V16\s3
        .if is64
            // __m128i src_v_hi1 = _mm_slli_epi32(rc_v_hi1, rc_scale);
            USHL V12\s3, V12\s3, V16\s3
        .endif
    .endif
    .if \tr_scale
        // __m128i str_v_lo1 = _mm_slli_epi32(tr_v_lo1, tr_scale);
        USHL V23\s3, V23\s3, V17\s3
        .if is64
            // __m128i str_v_hi1 = _mm_slli_epi32(tr_v_hi1, tr_scale);
            USHL V13\s3, V13\s3, V17\s3
        .endif
    .endif
    // __m128i out_lo1 = _mm_add_epi32(src_v_lo1, str_v_lo1);
    ADD V11\s3, V11\s3, V23\s3
    // out_lo1 = _mm_add_epi32(out_lo1, rnd_v);
    // out_lo1 = _mm_srli_epi32(out_lo1, max_scale + 1);
    URSHR V11\s3, V11\s3, #\max_scale+1
    .if is64
        // __m128i out_hi1 = _mm_add_epi32(src_v_hi1, str_v_hi1);
        ADD V12\s3, V12\s3, V13\s3
        // out_hi1 = _mm_add_epi32(out_hi1, rnd_v);
        // out_hi1 = _mm_srli_epi32(out_hi1, max_scale + 1);
        URSHR V12\s3, V12\s3, #\max_scale+1
        // out_lo1 = _mm_packs_epi32(out_lo1, out_hi1);
        XTN V11\s4, V11\s3
        .if \w > 4
            XTN V12\s4, V12\s3
            MOV V11.D[1], V12.D[0]
        .endif
    .endif
.endm

.macro PLANAR_LOOP_PDPC out, w_pdpc=1, h_pdpc=1, s1, scale_vect=V28
    // __m128i tst1 = _mm_slli_epi16(out_lo1, 6);
    MOVI V7\s1, #6
    USHL \out\s1, V11\s1, V7\s1
    // __m128i w_x1 = _mm_mullo_epi16(out_lo1, y_v);
    .if \h_pdpc
        USHL V14\s1, V11\s1, V19\s1
    .endif
    // __m128i w_y1 = _mm_mullo_epi16(out_lo1, x_v1);
    .if \w_pdpc
        USHL V15\s1, V11\s1, \scale_vect\s1
    .endif
    // tst1 = _mm_subs_epu16(tst1, w_x1);
    .if \h_pdpc
        UQSUB \out\s1, \out\s1, V14\s1
    .endif
    // tst1 = _mm_subs_epu16(tst1, w_y1);
    .if \w_pdpc
        UQSUB \out\s1, \out\s1, V15\s1
    .endif 
    // tst1 = _mm_adds_epu16(tst1, add_v);
    UQADD \out\s1, \out\s1, V26\s1
    // __m128i xl1 = _mm_mullo_epi16(x_v1, l_v);
    .if \w_pdpc
        USHL V14\s1, V18\s1, \scale_vect\s1
    .endif
    // __m128i yt1 = _mm_mullo_epi16(y_v, a1);
    .if \h_pdpc
        USHL V15\s1, V8\s1, V19\s1
    .endif
    // tst1 = _mm_adds_epu16(tst1, xl1);
    .if \w_pdpc
        UQADD \out\s1, \out\s1, V14\s1
    .endif
    // tst1 = _mm_adds_epu16(tst1, yt1);
    .if \h_pdpc
        UQADD \out\s1, \out\s1, V15\s1
    .endif

    // __m128i out_v1 = _mm_srli_epi16(tst1, 6);
    USHR \out\s1, \out\s1, #6 
    // out_v1 = _mm_min_epi16(out_v1, _mm_set1_epi16(1023));
    // out_v1 = _mm_max_epi16(out_v1, _mm_setzero_si128());
    UCLIP \out\s1, V31\s1, V30\s1
.endm

.macro planar_pdpc, w, h, scale, max_scale, rc_scale, tr_scale, loop1, loop2=0, s1=.8H, s2=.8H, s3=.8H
    function planar_pdpc_w\w\()_h\h\()_neon, export=1
        .set is64,\w == 64 || \h == 64
        .if \w > 4
            .set xloop, \w/8
            .set jump, 16
        .else
            .set xloop, 1
            .set jump, 8
        .endif
        INPUT_SCALING \w
        // pb_h compute
        INIT_CLIP_VECTORS \s1
        //Load PDPC
        movrel x9, PDPC_SCALE\scale
        .if \scale < 2
            LD1   {V28\s1}, [x9]
        .else
            LD1   {V28\s1-V29\s1}, [x9]
        .endif
        movrel x10, MUL_COL_VECT
        LD1   {V24\s1}, [x10]
        // __m128i tr_val_v = _mm_set1_epi16(src_above[width + 1]);
        ADD x11, x0, \w*2
        LD1R {V23\s2}, [x11]
        // __m128i bl_val_v = _mm_set1_epi16(src_left[height + 1]);
        ADD x12, x1, \h*2
        LD1R {V27\s1}, [x12]
        //__m128i add_v = _mm_set1_epi16(32);
        MOVI V26\s1, #32
        // __m128i a1 = _mm_loadu_si128((__m128i *) (src_above + 1));
        LOAD_REF_\w x0, V
        DUP V8\s1, w5
        // __m128i tr_v1 = _mm_slli_epi16(a1, log2_pb_h);
        USHL V0\s1, V0\s1, V8\s1
        .if \w >= 16
            USHL V1\s1, V1\s1, V8\s1
        .endif
        .if \w >= 32
            USHL V2\s1, V2\s1, V8\s1
            USHL V3\s1, V3\s1, V8\s1
        .endif
        .if \w >= 64
            USHL V4\s1, V4\s1, V8\s1
            USHL V5\s1, V5\s1, V8\s1
            USHL V6\s1, V6\s1, V8\s1
            USHL V7\s1, V7\s1, V8\s1
        .endif
        //save SP
        MOV x14, SP
        //reserve stack space
        SUB SP, x14, 3*128
        SUB x10, x14, 1*128
        SUB x11, x14, 2*128
        MOV x12, SP
        //store in stack
        MOV x13, x10
        STORE_\w  V, x13, , \s1, enable_stride=0
        // __m128i l1 = _mm_loadu_si128((__m128i *) (src_left + 1));
        LOAD_REF_\h x1, V
        // __m128i rc_v1 = _mm_sub_epi16(tr_val_v, l1);
        SUB V10\s2, V23\s2, V0\s2
        .if \h >= 16
            SUB V11\s2, V23\s2, V1\s2
        .endif
        .if \h >= 32
            SUB V12\s2, V23\s2, V2\s2
            SUB V13\s2, V23\s2, V3\s2
        .endif
        .if \h >= 64
            SUB V14\s2, V23\s2, V4\s2
            SUB V15\s2, V23\s2, V5\s2
            SUB V16\s2, V23\s2, V6\s2
            SUB V17\s2, V23\s2, V7\s2
        .endif
        // __m128i lc_v1 = _mm_slli_epi16(l1, log2_pb_w);
        DUP V8\s2, w4
        USHL V0\s2, V0\s2, V8\s2
        .if \h >= 16
            USHL V1\s2, V1\s2, V8\s2
        .endif
        .if \h >= 32
            USHL V2\s2, V2\s2, V8\s2
            USHL V3\s2, V3\s2, V8\s2
        .endif
        .if \h >= 64
            USHL V4\s2, V4\s2, V8\s2
            USHL V5\s2, V5\s2, V8\s2
            USHL V6\s2, V6\s2, V8\s2
            USHL V7\s2, V7\s2, V8\s2
        .endif
        MOV x13, x11
        STORE_\h  V1, x13, , \s2, enable_stride=0
        MOV x13, x12
        STORE_\h  V, x13, , \s2, enable_stride=0
        .if \rc_scale
            MOVI V16\s3, #\rc_scale
        .endif
        .if \tr_scale
            MOVI V17\s3, #\tr_scale
        .endif
        .rept \loop1
        // __m128i rcol_v = _mm_set1_epi16(r_col[y]);
            LD1R {V21\s1}, [x11], #2
        // __m128i lcol_v = _mm_set1_epi16(l_col[y]);
            LD1R {V22\s1}, [x12], #2
        // int32_t l_val = (int32_t)src_left[y + 1];
        // __m128i l_v = _mm_set1_epi16(l_val);
            LD1R {V18\s1}, [x1], #2
        // int y_wgh = pdpc_w[y];
        // __m128i y_v = _mm_set1_epi16(y_wgh);
            LDRSH  w6, [x9], #2
            DUP   V19\s1, w6
        // __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);
            MUL V9\s1, V24\s1, V21\s1
        // _mm_slli_epi16(rcol_v, 3)
            MOVI V0\s1, #3
            USHL V10\s1, V21\s1, V0\s1
            .set first_it,1
            MOV x6, x0
            MOV x7, x10
            COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
            PLANAR_LOOP_PDPC V0, 1, 1, \s1, scale_vect=V28
            .if \w >= 16
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                .if \scale == 2
                    PLANAR_LOOP_PDPC V1, 1, 1, \s1, scale_vect=V29
                .else
                    PLANAR_LOOP_PDPC V1, 0, 1, \s1
                .endif
            .endif
            .if \w >= 32
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V2, 0, 1, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V3, 0, 1, \s1
            .endif
            .if \w >= 64
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V4, 0, 1, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V5, 0, 1, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V6, 0, 1, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V7, 0, 1, \s1
            .endif
            STORE_\w  V, x2, x3, \s1, 1, \scale
        .endr

    .if \loop2 > 0
        .rept \loop2
            // __m128i rcol_v = _mm_set1_epi16(r_col[y]);
            LD1R {V21\s1}, [x11], #2
            // __m128i lcol_v = _mm_set1_epi16(l_col[y]);
            LD1R {V22\s1}, [x12], #2
            // int32_t l_val = (int32_t)src_left[y + 1];
            // __m128i l_v = _mm_set1_epi16(l_val);
            LD1R {V18\s1}, [x1], #2
            // __m128i rcol_mul_v = _mm_mullo_epi16(rcol_v, mul_col);
            MUL V9\s1, V24\s1, V21\s1
            // _mm_slli_epi16(rcol_v, 3)
            MOVI V0\s1, #3
            USHL V10\s1, V21\s1, V0\s1
            .set first_it,1
            MOV x6, x0
            MOV x7, x10
            COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
            PLANAR_LOOP_PDPC V0, 1, 0, \s1, scale_vect=V28
            .if \w >= 16
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                .if \scale == 2
                    PLANAR_LOOP_PDPC V1, 1, 0, \s1, scale_vect=V29
                .else
                    PLANAR_LOOP_PDPC V1, 0, 0, \s1
                .endif
            .endif
            .if \w >= 32
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V2, 0, 0, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V3, 0, 0, \s1
            .endif
            .if \w >= 64
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V4, 0, 0, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V5, 0, 0, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V6, 0, 0, \s1
                COMPUTE_PLANAR_VECTOR \w, \max_scale, \rc_scale, \tr_scale, \s1, \s2, \s3
                PLANAR_LOOP_PDPC V7, 0, 0, \s1
            .endif
            STORE_\w  V, x2, x3, \s1, 1, \scale
        .endr
    .endif
        MOV SP, x14
        RET   lr
    endfunc
.endm


// static void
// vvc_intra_dc_pdpc_wX_hX_neon(const uint16_t *const src_above,
//                       const uint16_t *const src_left,
//                       uint16_t *const dst, ptrdiff_t dst_stride,
//                       int log2_pb_w, int log2_pb_h)
// dc_pdpc, w, h, scale, offset, shift, loop1, loop2=0, s1=.8H, s2=.8H
dc_pdpc  4,  4, 0,  #4, #3,  3,  1, .4H, .4h
dc_pdpc  4,  8, 0,  #4, #3,  3,  5, .4H
dc_pdpc  4, 16, 1,  #8, #4,  6, 10, .4H
dc_pdpc  4, 32, 1, #16, #5,  6, 26, .4H
dc_pdpc  4, 64, 1, #32, #6,  6, 58, .4H

dc_pdpc  8,  4, 0,  #4, #3,  3,  1
dc_pdpc  8,  8, 1,  #8, #4,  6,  2
dc_pdpc  8, 16, 1,  #8, #4,  6, 10
dc_pdpc  8, 32, 1, #16, #5,  6, 26
dc_pdpc  8, 64, 1, #32, #6,  6, 58

dc_pdpc 16,  4, 1,  #8, #4,  4,  0
dc_pdpc 16,  8, 1,  #8, #4,  6,  2
dc_pdpc 16, 16, 1, #16, #5,  6, 10
dc_pdpc 16, 32, 1, #16, #5,  6, 26
dc_pdpc 16, 64, 2, #32, #6, 12, 52

dc_pdpc 32,  4, 1, #16, #5,  4,  0
dc_pdpc 32,  8, 1, #16, #5,  6,  2
dc_pdpc 32, 16, 1, #16, #5,  6, 10
dc_pdpc 32, 32, 2, #32, #6, 12, 20
dc_pdpc 32, 64, 2, #32, #6, 12, 52

dc_pdpc 64,  4, 1, #32, #6,  4,  0
dc_pdpc 64,  8, 1, #32, #6,  6,  2
dc_pdpc 64, 16, 2, #32, #6, 12,  4
dc_pdpc 64, 32, 2, #32, #6, 12, 20
dc_pdpc 64, 64, 2, #64, #7, 12, 52

planar_pdpc  4,  4, 0, 2, 0, 0, 3,  1, .4H, .4H, .4H
planar_pdpc  4,  8, 0, 3, 1, 0, 3,  5, .4H, .8H, .4H
planar_pdpc  4, 16, 1, 4, 2, 0, 6, 10, .4H, .8H, .4H
planar_pdpc  4, 32, 1, 5, 3, 0, 6, 26, .4H, .8H, .4H
planar_pdpc  4, 64, 1, 6, 4, 0, 6, 58, .4H, .8H, .4S
// 
planar_pdpc  8,  4, 0, 3, 0, 1,  3,  1, s2=.4H
planar_pdpc  8,  8, 1, 3, 0, 0,  6,  2
planar_pdpc  8, 16, 1, 4, 1, 0,  6, 10
planar_pdpc  8, 32, 1, 5, 2, 0,  6, 26
planar_pdpc  8, 64, 1, 6, 3, 0,  6, 58, s3=.4S
// 
planar_pdpc 16,  4, 1, 4, 0, 2,  4,  0, s2=.4H
planar_pdpc 16,  8, 1, 4, 0, 1,  6,  2
planar_pdpc 16, 16, 1, 4, 0, 0,  6, 10
planar_pdpc 16, 32, 1, 5, 1, 0,  6, 26
planar_pdpc 16, 64, 2, 6, 2, 0, 12, 52, s3=.4S
// 
planar_pdpc 32,  4, 1, 5, 0, 3,  4,  0, s2=.4H
planar_pdpc 32,  8, 1, 5, 0, 2,  6,  2
planar_pdpc 32, 16, 1, 5, 0, 1,  6, 10
planar_pdpc 32, 32, 2, 5, 0, 0, 12, 20
planar_pdpc 32, 64, 2, 6, 1, 0, 12, 52, s3=.4S
// 
planar_pdpc 64,  4, 1, 6, 0, 4,  4,  0, s2=.4H, s3=.4S
planar_pdpc 64,  8, 1, 6, 0, 3,  6,  2, s3=.4S
planar_pdpc 64, 16, 2, 6, 0, 2, 12,  4, s3=.4S
planar_pdpc 64, 32, 2, 6, 0, 1, 12, 20, s3=.4S
planar_pdpc 64, 64, 2, 6, 0, 0, 12, 52, s3=.4S

const PDPC_SCALE0
    .2byte  5,  3,  1,  0x2F, 0x2F, 0x2F, 0x2F, 0x2F
endconst

const PDPC_SCALE1
    .2byte  5,  4,  3,     2,    1,    0, 0x2F, 0x2F
endconst

const PDPC_SCALE2
    .2byte  5,  5,  4,     4,    3,    3,    2,    2,    1,    1,    0,    0, 0x2F, 0x2F, 0x2F, 0x2F
endconst

const MUL_COL_VECT
    .2byte  1, 2, 3, 4, 5, 6, 7, 8
endconst

