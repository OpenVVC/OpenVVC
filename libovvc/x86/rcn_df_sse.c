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
#include <tmmintrin.h>
#include <smmintrin.h>

#include "ovutils.h"
#include "rcn_structures.h"

static void
filter_luma_weak_h(OVSample* src, const int stride, const int tc, const uint8_t extend_p, const uint8_t extend_q)
{
    const int th_cut = tc * 10;
    const int tc2_p = -extend_p & (tc >> 1);
    const int tc2_q = -extend_q & (tc >> 1);

    for (int i = 0; i < 4; i++) {
        const int16_t p2  = src[-3];
        const int16_t p1  = src[-2];
        const int16_t p0  = src[-1];
        const int16_t q0  = src[0];
        const int16_t q1  = src[1];
        const int16_t q2  = src[2];

        /* Weak filter */
        int delta = (((q0 - p0) << 3) + (q0 - p0) - (((q1 - p1) << 1) + (q1 - p1)) + 8) >> 4;

        if (abs(delta) < th_cut) {
            delta = ov_clip(delta, -tc, tc);
            const int delta1 = ov_clip(((((p2 + p0 + 1) >> 1) - p1 + delta) >> 1), -tc2_p, tc2_p);
            const int delta2 = ov_clip(((((q2 + q0 + 1) >> 1) - q1 - delta) >> 1), -tc2_q, tc2_q);
            src[-2] = ov_bdclip(p1 + delta1);
            src[-1] = ov_bdclip(p0 + delta);
            src[0]  = ov_bdclip(q0 - delta);
            src[1]  = ov_bdclip(q1 + delta2);
        }
        src += stride;
    }
}

static void
filter_luma_strong_small_h(OVSample* src, const int stride, const int tc)
{
    static const int16_t tc_c[8] = { 0, 1, 2 ,3, 3, 2, 1, 0 };
    static const int16_t msk[8*8] = {
        0, 2, 0, 0, 0, 0, 0, 0,
        0, 3, 2, 1, 0, 0, 0, 0,
        0, 1, 2, 2, 1, 0, 0, 0,
        0, 1, 2, 2, 2, 2, 1, 0,
        0, 1, 2, 2, 2, 2, 1, 0,
        0, 0, 0, 1, 2, 2, 1, 0,
        0, 0, 0, 0, 1, 2, 3, 0,
        0, 0, 0, 0, 0, 0, 2, 0
    };

    const __m128i msk0  = _mm_loadu_si128((__m128i *)msk + 0);
    const __m128i msk1  = _mm_loadu_si128((__m128i *)msk + 1);
    const __m128i msk2  = _mm_loadu_si128((__m128i *)msk + 2);
    const __m128i msk3  = _mm_loadu_si128((__m128i *)msk + 3);
    const __m128i msk4  = _mm_loadu_si128((__m128i *)msk + 4);
    const __m128i msk5  = _mm_loadu_si128((__m128i *)msk + 5);
    const __m128i msk6  = _mm_loadu_si128((__m128i *)msk + 6);
    const __m128i msk7  = _mm_loadu_si128((__m128i *)msk + 7);

    const __m128i clp_val = _mm_mullo_epi16(_mm_set1_epi16(tc), _mm_loadu_si128((__m128i *)tc_c));
    const __m128i add_4   = _mm_set1_epi16(4);

    for (int i = 0; i < 4; i++) {
        __m128i dval;
        __m128i line = _mm_loadu_si128((__m128i *) &src[-4]);
        __m128i clp_min = _mm_sub_epi16(line, clp_val);
        __m128i clp_max = _mm_add_epi16(line, clp_val);

        __m128i p3  = _mm_set1_epi16(src[-4]);
        __m128i p2  = _mm_set1_epi16(src[-3]);
        __m128i p1  = _mm_set1_epi16(src[-2]);
        __m128i p0  = _mm_set1_epi16(src[-1]);
        __m128i q0  = _mm_set1_epi16(src[ 0]);
        __m128i q1  = _mm_set1_epi16(src[ 1]);
        __m128i q2  = _mm_set1_epi16(src[ 2]);
        __m128i q3  = _mm_set1_epi16(src[ 3]);

        p3  = _mm_mullo_epi16(p3, msk0);
        p2  = _mm_mullo_epi16(p2, msk1);
        p1  = _mm_mullo_epi16(p1, msk2);
        p0  = _mm_mullo_epi16(p0, msk3);
        q0  = _mm_mullo_epi16(q0, msk4);
        q1  = _mm_mullo_epi16(q1, msk5);
        q2  = _mm_mullo_epi16(q2, msk6);
        q3  = _mm_mullo_epi16(q3, msk7);

        dval = _mm_add_epi16(p2, p3);
        dval = _mm_add_epi16(dval, p1);
        dval = _mm_add_epi16(dval, p0);
        dval = _mm_add_epi16(dval, q0);
        dval = _mm_add_epi16(dval, q1);
        dval = _mm_add_epi16(dval, q2);
        dval = _mm_add_epi16(dval, q3);

        dval = _mm_add_epi16(dval, add_4);

        dval = _mm_srli_epi16(dval, 3);

        dval = _mm_max_epi16(dval, clp_min);
        dval = _mm_min_epi16(dval, clp_max);

        _mm_storeu_si128((__m128i *) &src[-4], dval);
        src += stride;
    }
}

static void
filter_h_7_7(OVSample *src, const int stride, const int tc)
{
    static const int16_t db7_p[8] = { 0, 5, 14, 23, 32, 41, 50, 59 };
    static const int16_t tc7_p[8] = { 0, 1, 1, 2, 3, 4, 5, 6, };
    static const int16_t db7_q[8] = { 59, 50, 41, 32, 23, 14, 5, 0 };
    static const int16_t tc7_q[8] = { 6, 5, 4, 3, 2, 1, 1, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc7_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc7_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db7_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db7_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcP[0] + srcQ[0])
                     + srcP[-1] + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1] + srcP[-6 * 1]
                     + srcQ[ 1] + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + srcQ[ 6 * 1] + 8) >> 4;
        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_h_5_7(OVSample *src, const int stride, const int tc)
{
    static const int16_t db5_p[8] = { 0, 0, 0, 6, 19, 32, 45, 58 };
    static const int16_t tc5_p[8] = { 0, 0, 0, 2, 3, 4, 5, 6 };
    static const int16_t db7_q[8] = { 59, 50, 41, 32, 23, 14, 5, 0 };
    static const int16_t tc7_q[8] = { 6, 5, 4, 3, 2, 1, 1, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc5_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc7_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db5_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db7_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-4 * 1] + srcP[-5 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcP[0] + srcP[-1] + srcQ[0] + srcQ[ 1])
                         + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1]
                         + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + 8) >> 4;
        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_h_7_5(OVSample *src, const int stride, const int tc)
{

    static const int16_t db7_p[8] = { 0, 5, 14, 23, 32, 41, 50, 59 };
    static const int16_t tc7_p[8] = { 0, 1, 1, 2, 3, 4, 5, 6 };
    static const int16_t db5_q[8] = { 58, 45, 32, 19, 6 , 0, 0, 0 };
    static const int16_t tc5_q[8] = { 6, 5, 4, 3, 2 , 0, 0, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc7_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc5_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db7_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db5_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 4 * 1] + srcQ[ 5 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcP[0] + srcP[-1] + srcQ[0] + srcQ[ 1])
                + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1]
                + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + 8) >> 4;
        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_h_3_5(OVSample *src, const int stride, const int tc)
{

    static const int16_t db3_p[8] = { 0, 0, 0, 0, 0, 11, 32, 53 };
    static const int16_t tc3_p[8] = { 0, 0, 0, 0, 0, 2, 4, 6 };
    static const int16_t db5_q[8] = { 58, 45, 32, 19, 6 , 0, 0, 0 };
    static const int16_t tc5_q[8] = { 6, 5, 4, 3, 2 , 0, 0, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc3_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc5_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db3_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db5_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-2 * 1] + srcP[-3 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 4 * 1] + srcQ[ 5 * 1] + 1) >> 1;

        db_ref[i] = (srcP[0] + srcP[-1] + srcP[-2] + srcP[-3]
                    + srcQ[0] + srcQ[1] + srcQ[2] + srcQ[3] + 4) >> 3;


        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_h_5_5(OVSample *src, const int stride, const int tc)
{

    static const int16_t db5_p[8] = { 0, 0, 0, 6, 19, 32, 45, 58 };
    static const int16_t tc5_p[8] = { 0, 0, 0, 2, 3, 4, 5, 6 };
    static const int16_t db5_q[8] = { 58, 45, 32, 19, 6 , 0, 0, 0 };
    static const int16_t tc5_q[8] = { 6, 5, 4, 3, 2 , 0, 0, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc5_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc5_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db5_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db5_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-4 * 1] + srcP[-5 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 4 * 1] + srcQ[ 5 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcP[0] + srcP[-1] + srcP[-2 * 1]
                        + srcQ[0] + srcQ[ 1] + srcQ[ 2 * 1])
                     + srcP[-3 * 1] + srcP[-4 * 1]
                     + srcQ[ 3 * 1] + srcQ[ 4 * 1] + 8) >> 4;

        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_h_7_3(OVSample *src, const int stride, const int tc)
{

    static const int16_t db7_p[8] = { 0, 5, 14, 23, 32, 41, 50, 59 };
    static const int16_t tc7_p[8] = { 0, 1, 1, 2, 3, 4, 5, 6 };
    static const int16_t db3_q[8] = { 53, 32, 11, 0, 0, 0, 0, 0 };
    static const int16_t tc3_q[8] = { 6, 4, 2, 0, 0, 0, 0, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc7_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc3_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db7_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db3_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 2 * 1] + srcQ[ 3 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcP[0] + srcQ[0])
            + srcP[     -1] + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1] + srcP[-6 * 1]
            + srcQ[      0] + srcQ[     1] + srcQ[     1] + srcQ[2 *  1] + srcQ[ 2 * 1] + srcQ[     1] + 8) >> 4;

        _src += stride;
    }

    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_h_3_7(OVSample *src, const int stride, const int tc)
{

    static const int16_t db3_p[8] = { 0, 0, 0, 0, 0, 11, 32, 53 };
    static const int16_t tc3_p[8] = { 0, 0, 0, 0, 0, 2, 4, 6 };
    static const int16_t db7_q[8] = { 59, 50, 41, 32, 23, 14, 5, 0 };
    static const int16_t tc7_q[8] = { 6, 5, 4, 3, 2, 1, 1, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc3_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc7_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db3_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db7_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-2 * 1] + srcP[-3 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

        db_ref[i] = (2 * (srcQ[0] + srcP[0])
            + srcP[      0] + srcP[    -1] + srcP[    -1] + srcP[-2 * 1] + srcP[-2 * 1] + srcP[    -1]
            + srcQ[      1] + srcQ[2 *  1] + srcQ[3 *  1] + srcQ[4 *  1] + srcQ[5 *  1] + srcQ[6 *  1] + 8) >> 4;

        _src += stride;
    }

    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_luma_weak_v(OVSample* src, const int stride, const int tc, const uint8_t extend_p, const uint8_t extend_q)
{
    const int th_cut = tc * 10;
    const int tc2_p = -extend_p & (tc >> 1);
    const int tc2_q = -extend_q & (tc >> 1);

    for (int i = 0; i < 4; i++) {
        const int16_t p2  = src[-3*stride];
        const int16_t p1  = src[-2*stride];
        const int16_t p0  = src[-1*stride];
        const int16_t q0  = src[0*stride];
        const int16_t q1  = src[1*stride];
        const int16_t q2  = src[2*stride];

        /* Weak filter */
        int delta = (9 * (q0 - p0) - 3 * (q1 - p1) + 8) >> 4;

        if (abs(delta) < th_cut) {
            delta = ov_clip(delta, -tc, tc);
            const int delta1 = ov_clip(((((p2 + p0 + 1) >> 1) - p1 + delta) >> 1), -tc2_p, tc2_p);
            const int delta2 = ov_clip(((((q2 + q0 + 1) >> 1) - q1 - delta) >> 1), -tc2_q, tc2_q);
            src[-2*stride] = ov_bdclip(p1 + delta1);
            src[-1*stride] = ov_bdclip(p0 + delta);
            src[0*stride]  = ov_bdclip(q0 - delta);
            src[1*stride]  = ov_bdclip(q1 + delta2);
        }
        ++src;
    }
}

static void
filter_v_7_7(OVSample *src, const int stride, const int tc)
{
    static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

    int16_t* srcP = src - stride;
    int16_t* srcQ = src;

    __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
    __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
    __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
    __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
    __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
    __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );
    __m128i srcP6  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -6] );
    __m128i srcP7  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -7] );

    __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
    __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
    __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
    __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
    __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
    __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );
    __m128i srcQ6  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 6] );
    __m128i srcQ7  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 7] );

    __m128i refPx4 = _mm_add_epi16(srcP6, srcP7);
    __m128i refQx4 = _mm_add_epi16(srcQ6, srcQ7);
    refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
    refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
    refPx4 = _mm_srai_epi16(refPx4, 1);
    refQx4 = _mm_srai_epi16(refQx4, 1);

    srcP0  = _mm_cvtepi16_epi32(srcP0);
    srcP1  = _mm_cvtepi16_epi32(srcP1);
    srcP2  = _mm_cvtepi16_epi32(srcP2);
    srcP3  = _mm_cvtepi16_epi32(srcP3);
    srcP4  = _mm_cvtepi16_epi32(srcP4);
    srcP5  = _mm_cvtepi16_epi32(srcP5);
    srcP6  = _mm_cvtepi16_epi32(srcP6);

    srcQ0  = _mm_cvtepi16_epi32(srcQ0);
    srcQ1  = _mm_cvtepi16_epi32(srcQ1);
    srcQ2  = _mm_cvtepi16_epi32(srcQ2);
    srcQ3  = _mm_cvtepi16_epi32(srcQ3);
    srcQ4  = _mm_cvtepi16_epi32(srcQ4);
    srcQ5  = _mm_cvtepi16_epi32(srcQ5);
    srcQ6  = _mm_cvtepi16_epi32(srcQ6);

    srcP0 = _mm_add_epi32(srcP0, srcQ0);
    srcP1 = _mm_add_epi32(srcP1, srcQ1);
    srcP2 = _mm_add_epi32(srcP2, srcQ2);
    srcP3 = _mm_add_epi32(srcP3, srcQ3);
    srcP4 = _mm_add_epi32(srcP4, srcQ4);
    srcP5 = _mm_add_epi32(srcP5, srcQ5);
    srcP6 = _mm_add_epi32(srcP6, srcQ6);

    srcP0 = _mm_slli_epi32(srcP0, 1);
    srcP1 = _mm_add_epi32(srcP1, srcP2);
    srcP3 = _mm_add_epi32(srcP3, srcP4);
    srcP5 = _mm_add_epi32(srcP5, srcP6);

    srcP0 = _mm_add_epi32(srcP0, srcP1);
    srcP3 = _mm_add_epi32(srcP3, srcP5);

    srcP0 = _mm_add_epi32(srcP0, srcP3);
    __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
    refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
    refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

   for( int pos = 0; pos < 7; pos++ )
   {
     __m128i vref1 = refPx4;
     __m128i vref0 = refMiddlex4;
     __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
     __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
     __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
     vmax          = _mm_add_epi16( vsrc, vmax );
     vref0         = _mm_unpacklo_epi16( vref0, vref1 );
     __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
     vtmp          = _mm_madd_epi16( vref0, vtmp );
     vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
     vtmp          = _mm_srli_epi32( vtmp, 6 );
     vtmp          = _mm_packs_epi32( vtmp, vtmp );
     vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
     _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
   }

   for( int pos = 0; pos < 7; pos++ )
   {
     __m128i vref1 = refQx4;
     __m128i vref0 = refMiddlex4;
     __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
     __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
     __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
     vmax          = _mm_add_epi16( vsrc, vmax );
     vref0         = _mm_unpacklo_epi16( vref0, vref1 );
     __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
     vtmp          = _mm_madd_epi16( vref0, vtmp );
     vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
     vtmp          = _mm_srli_epi32( vtmp, 6 );
     vtmp          = _mm_packs_epi32( vtmp, vtmp );
     vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
     _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
   }
}

static void
filter_v_7_5(OVSample *src, const int stride, const int tc)
{
      static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
      static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
      static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

      int16_t* srcP = src - stride;
      int16_t* srcQ = src;

      __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
      __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
      __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
      __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
      __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
      __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );
      __m128i srcP6  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -6] );
      __m128i srcP7  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -7] );

      __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
      __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
      __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
      __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
      __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
      __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );

      __m128i refPx4 = _mm_add_epi16(srcP6, srcP7);
      __m128i refQx4 = _mm_add_epi16(srcQ4, srcQ5);
      refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
      refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
      refPx4 = _mm_srai_epi16(refPx4, 1);
      refQx4 = _mm_srai_epi16(refQx4, 1);

      srcP0  = _mm_cvtepi16_epi32(srcP0);
      srcP1  = _mm_cvtepi16_epi32(srcP1);
      srcP2  = _mm_cvtepi16_epi32(srcP2);
      srcP3  = _mm_cvtepi16_epi32(srcP3);
      srcP4  = _mm_cvtepi16_epi32(srcP4);
      srcP5  = _mm_cvtepi16_epi32(srcP5);

      srcQ0  = _mm_cvtepi16_epi32(srcQ0);
      srcQ1  = _mm_cvtepi16_epi32(srcQ1);
      srcQ2  = _mm_cvtepi16_epi32(srcQ2);
      srcQ3  = _mm_cvtepi16_epi32(srcQ3);
      srcQ4  = _mm_cvtepi16_epi32(srcQ4);
      srcQ5  = _mm_cvtepi16_epi32(srcQ5);

      srcP0 = _mm_add_epi32(srcP0, srcQ0);
      srcP1 = _mm_add_epi32(srcP1, srcQ1);
      srcP2 = _mm_add_epi32(srcP2, srcQ2);
      srcP3 = _mm_add_epi32(srcP3, srcQ3);
      srcP4 = _mm_add_epi32(srcP4, srcQ4);
      srcP5 = _mm_add_epi32(srcP5, srcQ5);

      srcP0 = _mm_add_epi32(srcP0, srcP1);
      srcP2 = _mm_add_epi32(srcP2, srcP3);
      srcP4 = _mm_add_epi32(srcP4, srcP5);
      srcP0 = _mm_slli_epi32(srcP0, 1);

      srcP0 = _mm_add_epi32(srcP0, srcP2);

      srcP0 = _mm_add_epi32(srcP0, srcP4);
      __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
      refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
      refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

     for( int pos = 0; pos < 7; pos++ )
     {
       __m128i vref1 = refPx4;
       __m128i vref0 = refMiddlex4;
       __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
       __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
       __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
       vmax          = _mm_add_epi16( vsrc, vmax );
       vref0         = _mm_unpacklo_epi16( vref0, vref1 );
       __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
       vtmp          = _mm_madd_epi16( vref0, vtmp );
       vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
       vtmp          = _mm_srli_epi32( vtmp, 6 );
       vtmp          = _mm_packs_epi32( vtmp, vtmp );
       vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
       _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
     }

     for( int pos = 0; pos < 5; pos++ )
     {
       __m128i vref1 = refQx4;
       __m128i vref0 = refMiddlex4;
       __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
       __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
       __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
       vmax          = _mm_add_epi16( vsrc, vmax );
       vref0         = _mm_unpacklo_epi16( vref0, vref1 );
       __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
       vtmp          = _mm_madd_epi16( vref0, vtmp );
       vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
       vtmp          = _mm_srli_epi32( vtmp, 6 );
       vtmp          = _mm_packs_epi32( vtmp, vtmp );
       vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
       _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
     }
}

static void
filter_v_5_7(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );
  __m128i srcQ6  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 6] );
  __m128i srcQ7  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 7] );

  __m128i refPx4 = _mm_add_epi16(srcP4, srcP5);
  __m128i refQx4 = _mm_add_epi16(srcQ6, srcQ7);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);
  srcP4  = _mm_cvtepi16_epi32(srcP4);
  srcP5  = _mm_cvtepi16_epi32(srcP5);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);
  srcQ4  = _mm_cvtepi16_epi32(srcQ4);
  srcQ5  = _mm_cvtepi16_epi32(srcQ5);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);
  srcP4 = _mm_add_epi32(srcP4, srcQ4);
  srcP5 = _mm_add_epi32(srcP5, srcQ5);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP2 = _mm_add_epi32(srcP2, srcP3);
  srcP4 = _mm_add_epi32(srcP4, srcP5);
  srcP0 = _mm_slli_epi32(srcP0, 1);

  srcP0 = _mm_add_epi32(srcP0, srcP2);

  srcP0 = _mm_add_epi32(srcP0, srcP4);
  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 7; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }

}

static void
filter_v_5_5(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );

  __m128i refPx4 = _mm_add_epi16(srcP4, srcP5);
  __m128i refQx4 = _mm_add_epi16(srcQ4, srcQ5);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);
  srcP4  = _mm_cvtepi16_epi32(srcP4);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);
  srcQ4  = _mm_cvtepi16_epi32(srcQ4);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);
  srcP4 = _mm_add_epi32(srcP4, srcQ4);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP0 = _mm_add_epi32(srcP0, srcP2);
  srcP3 = _mm_add_epi32(srcP3, srcP4);
  srcP0 = _mm_slli_epi32(srcP0, 1);

  srcP0 = _mm_add_epi32(srcP0, srcP3);

  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_7_3(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };
  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );
  __m128i srcP6  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -6] );
  __m128i srcP7  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -7] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );

  __m128i refPx4 = _mm_add_epi16(srcP6, srcP7);
  __m128i refQx4 = _mm_add_epi16(srcQ2, srcQ3);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);
  srcP4  = _mm_cvtepi16_epi32(srcP4);
  srcP5  = _mm_cvtepi16_epi32(srcP5);
  srcP6  = _mm_cvtepi16_epi32(srcP6);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);


  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ0);
  srcP2 = _mm_add_epi32(srcP2, srcQ1);
  srcP3 = _mm_add_epi32(srcP3, srcQ1);
  srcP4 = _mm_add_epi32(srcP4, srcQ1);
  srcP5 = _mm_add_epi32(srcP5, srcQ2);
  srcP6 = _mm_add_epi32(srcP6, srcQ2);

  srcP0 = _mm_slli_epi32(srcP0, 1);
  srcP1 = _mm_add_epi32(srcP1, srcP2);
  srcP3 = _mm_add_epi32(srcP3, srcP4);
  srcP5 = _mm_add_epi32(srcP5, srcP6);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP3 = _mm_add_epi32(srcP3, srcP5);

  srcP0 = _mm_add_epi32(srcP0, srcP3);
  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 7; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_3_7(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs7[7] = { 59, 50, 41, 32, 23, 14, 5 };
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );
  __m128i srcQ6  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 6] );
  __m128i srcQ7  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 7] );

  __m128i refPx4 = _mm_add_epi16(srcP2, srcP3);
  __m128i refQx4 = _mm_add_epi16(srcQ6, srcQ7);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);
  srcQ4  = _mm_cvtepi16_epi32(srcQ4);
  srcQ5  = _mm_cvtepi16_epi32(srcQ5);
  srcQ6  = _mm_cvtepi16_epi32(srcQ6);

  srcQ0 = _mm_add_epi32(srcP0, srcQ0);
  srcQ1 = _mm_add_epi32(srcP0, srcQ1);
  srcQ2 = _mm_add_epi32(srcP1, srcQ2);
  srcQ3 = _mm_add_epi32(srcP1, srcQ3);
  srcQ4 = _mm_add_epi32(srcP1, srcQ4);
  srcQ5 = _mm_add_epi32(srcP2, srcQ5);
  srcQ6 = _mm_add_epi32(srcP2, srcQ6);

  srcQ0 = _mm_slli_epi32(srcQ0, 1);
  srcQ1 = _mm_add_epi32(srcQ1, srcQ2);
  srcQ3 = _mm_add_epi32(srcQ3, srcQ4);
  srcQ5 = _mm_add_epi32(srcQ5, srcQ6);

  srcQ0 = _mm_add_epi32(srcQ0, srcQ1);
  srcQ3 = _mm_add_epi32(srcQ3, srcQ5);

  srcQ0 = _mm_add_epi32(srcQ0, srcQ3);
  __m128i refMiddlex4 = _mm_add_epi32(srcQ0, _mm_set1_epi32(8));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 4);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 7; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs7[pos] | ( ( 64 - dbCoeffs7[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_v_5_3(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );
  __m128i srcP4  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -4] );
  __m128i srcP5  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -5] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );


  __m128i refPx4 = _mm_add_epi16(srcP4, srcP5);
  __m128i refQx4 = _mm_add_epi16(srcQ2, srcQ3);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP2 = _mm_add_epi32(srcP2, srcP3);

  srcP0 = _mm_add_epi32(srcP0, srcP2);

  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(4));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 3);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_h_5_3(OVSample *src, const int stride, const int tc)
{
    static const int16_t db5_p[8] = { 0, 0, 0, 6, 19, 32, 45, 58 };
    static const int16_t tc5_p[8] = { 0, 0, 0, 2, 3, 4, 5, 6 };
    static const int16_t db3_q[8] = { 53, 32, 11, 0, 0, 0, 0, 0 };
    static const int16_t tc3_q[8] = { 6, 4, 2, 0, 0, 0, 0, 0 };

    __m128i tmp = _mm_set1_epi16(tc);
    __m128i clp_p = _mm_loadu_si128((__m128i *)tc5_p);
    clp_p = _mm_mullo_epi16(clp_p, tmp);
    clp_p = _mm_srli_epi16(clp_p, 1);

    __m128i clp_q = _mm_loadu_si128((__m128i *)tc3_q);
    clp_q = _mm_mullo_epi16(clp_q, tmp);
    clp_q = _mm_srli_epi16(clp_q, 1);

    __m128i db7p = _mm_loadu_si128((__m128i *)db5_p);
    __m128i db7q = _mm_loadu_si128((__m128i *)db3_q);


    int16_t ref_p[4];
    int16_t ref_q[4];
    int16_t db_ref[4];

    OVSample* _src = src;
    for (int i = 0; i < 4; ++i) {
        OVSample* srcP = _src - 1;
        OVSample* srcQ = _src;
        ref_p[i] = (srcP[-4 * 1] + srcP[-5 * 1] + 1) >> 1;
        ref_q[i] = (srcQ[ 2 * 1] + srcQ[ 3 * 1] + 1) >> 1;

        db_ref[i] = (srcP[0] + srcP[-1] + srcP[-2 * 1] + srcP[-3 * 1]
                   + srcQ[0] + srcQ[ 1] + srcQ[ 2 * 1] + srcQ[ 3 * 1] + 4) >> 3;

        _src += stride;
    }
    __m128i add_32 = _mm_set1_epi16(32);

    for (int i = 0; i < 4; i++) {
        OVSample* _src = src - 8;
        __m128i p = _mm_load_si128((__m128i *) _src);
        __m128i q = _mm_load_si128((__m128i *)  src);
        __m128i clp_min_p = _mm_sub_epi16(p, clp_p);
        __m128i clp_max_p = _mm_add_epi16(p, clp_p);
        __m128i clp_min_q = _mm_sub_epi16(q, clp_q);
        __m128i clp_max_q = _mm_add_epi16(q, clp_q);
        __m128i dref = _mm_set1_epi16(db_ref[i]);
        __m128i rf_p = _mm_set1_epi16(ref_p[i]);
        __m128i rf_q = _mm_set1_epi16(ref_q[i]);

        __m128i dref_p = _mm_mullo_epi16(dref, db7p);
        __m128i dref_q = _mm_mullo_epi16(dref, db7q);
        __m128i x1 = _mm_mullo_epi16(rf_p, db7p);
        __m128i x2 = _mm_mullo_epi16(rf_q, db7q);

        __m128i x0 = _mm_subs_epu16(_mm_slli_epi16(rf_p, 6), x1);
        __m128i x3 = _mm_subs_epu16(_mm_slli_epi16(rf_q, 6), x2);

        x0 = _mm_adds_epu16(x0, dref_p);
        x0 = _mm_adds_epu16(x0, add_32);
        x3 = _mm_adds_epu16(x3, dref_q);
        x3 = _mm_adds_epu16(x3, add_32);

        x3 = _mm_srli_epi16(x3, 6);
        x0 = _mm_srli_epi16(x0, 6);

        x0 = _mm_max_epi16(x0, clp_min_p);
        x3 = _mm_max_epi16(x3, clp_min_q);
        x0 = _mm_min_epi16(x0, clp_max_p);
        x3 = _mm_min_epi16(x3, clp_max_q);

        _mm_store_si128((__m128i *) _src, x0);
        _mm_store_si128((__m128i *) src, x3);

        src += stride;
    }
}

static void
filter_luma_strong_small_v(OVSample* src, const int stride, const int tc)
{
    int16_t tc_c[6] = { 1 * tc, 2 * tc, 3 * tc , 3 * tc, 2 * tc, 1 * tc};

    for (int i = 0; i < 4; i++) {
        const int16_t p3  = src[-stride * 4];
        const int16_t p2  = src[-stride * 3];
        const int16_t p1  = src[-stride * 2];
        const int16_t p0  = src[-stride    ];
        const int16_t q0  = src[ 0         ];
        const int16_t q1  = src[ stride    ];
        const int16_t q2  = src[ stride * 2];
        const int16_t q3  = src[ stride * 3];

        src[-stride * 3] = ov_clip((((p2 + p2) + (p3 + p3) + (p2 + p1) + (p0 + q0) + 4) >> 3), p2 - tc_c[0], p2 + tc_c[0]);
        src[-stride * 2] = ov_clip((((p1 + p2) + (p2 + p1) + (p0 + p0) + (q0 + q0) + 4) >> 3), p1 - tc_c[1], p1 + tc_c[1]);
        src[-stride]     = ov_clip((((p0 + p1) + (p2 + p1) + (p0 + q0) + (q0 + q1) + 4) >> 3), p0 - tc_c[2], p0 + tc_c[2]);
        src[0]           = ov_clip((((q0 + q1) + (p1 + p0) + (p0 + q0) + (q1 + q2) + 4) >> 3), q0 - tc_c[3], q0 + tc_c[3]);
        src[stride]      = ov_clip((((q1 + q2) + (p0 + q0) + (p0 + q0) + (q1 + q2) + 4) >> 3), q1 - tc_c[4], q1 + tc_c[4]);
        src[stride * 2]  = ov_clip((((q2 + q2) + (p0 + q0) + (q1 + q2) + (q3 + q3) + 4) >> 3), q2 - tc_c[5], q2 + tc_c[5]);
        src++;
    }
}

static void
filter_v_3_5(OVSample *src, const int stride, const int tc)
{
  static const int dbCoeffs3[3] = { 53, 32, 11 };
  static const int dbCoeffs5[5] = { 58, 45, 32, 19, 6};
  static const int8_t tc7[7] = { 6, 5, 4, 3, 2, 1, 1};
  static const int8_t tc3[3] = { 6, 4, 2 };

  int16_t* srcP = src - stride;
  int16_t* srcQ = src;

  __m128i srcP0  = _mm_loadl_epi64( ( const __m128i* ) &srcP[0] );
  __m128i srcP1  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride] );
  __m128i srcP2  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -2] );
  __m128i srcP3  = _mm_loadl_epi64( ( const __m128i* ) &srcP[stride * -3] );

  __m128i srcQ0  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[0] );
  __m128i srcQ1  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride] );
  __m128i srcQ2  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 2] );
  __m128i srcQ3  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 3] );
  __m128i srcQ4  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 4] );
  __m128i srcQ5  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * 5] );


  __m128i refPx4 = _mm_add_epi16(srcP2, srcP3);
  __m128i refQx4 = _mm_add_epi16(srcQ4, srcQ5);
  refPx4 = _mm_add_epi16(refPx4, _mm_set1_epi16(1));
  refQx4 = _mm_add_epi16(refQx4, _mm_set1_epi16(1));
  refPx4 = _mm_srai_epi16(refPx4, 1);
  refQx4 = _mm_srai_epi16(refQx4, 1);

  srcP0  = _mm_cvtepi16_epi32(srcP0);
  srcP1  = _mm_cvtepi16_epi32(srcP1);
  srcP2  = _mm_cvtepi16_epi32(srcP2);
  srcP3  = _mm_cvtepi16_epi32(srcP3);

  srcQ0  = _mm_cvtepi16_epi32(srcQ0);
  srcQ1  = _mm_cvtepi16_epi32(srcQ1);
  srcQ2  = _mm_cvtepi16_epi32(srcQ2);
  srcQ3  = _mm_cvtepi16_epi32(srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcQ0);
  srcP1 = _mm_add_epi32(srcP1, srcQ1);
  srcP2 = _mm_add_epi32(srcP2, srcQ2);
  srcP3 = _mm_add_epi32(srcP3, srcQ3);

  srcP0 = _mm_add_epi32(srcP0, srcP1);
  srcP2 = _mm_add_epi32(srcP2, srcP3);

  srcP0 = _mm_add_epi32(srcP0, srcP2);

  __m128i refMiddlex4 = _mm_add_epi32(srcP0, _mm_set1_epi32(4));
  refMiddlex4 = _mm_srai_epi32(refMiddlex4, 3);
  refMiddlex4 = _mm_packs_epi32(refMiddlex4, _mm_setzero_si128());

 for( int pos = 0; pos < 3; pos++ )
 {
   __m128i vref1 = refPx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcP[-stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc3[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs3[pos] | ( ( 64 - dbCoeffs3[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcP[-stride * pos], vtmp );
 }

 for( int pos = 0; pos < 5; pos++ )
 {
   __m128i vref1 = refQx4;
   __m128i vref0 = refMiddlex4;
   __m128i vsrc  = _mm_loadl_epi64( ( const __m128i* ) &srcQ[stride * pos] );
   __m128i vmax  = _mm_set1_epi16( ( tc * tc7[pos] ) >> 1 );
   __m128i vmin  = _mm_sub_epi16( vsrc, vmax );
   vmax          = _mm_add_epi16( vsrc, vmax );
   vref0         = _mm_unpacklo_epi16( vref0, vref1 );
   __m128i vtmp  = _mm_set1_epi32( dbCoeffs5[pos] | ( ( 64 - dbCoeffs5[pos] ) << 16 ) );
   vtmp          = _mm_madd_epi16( vref0, vtmp );
   vtmp          = _mm_add_epi32( vtmp, _mm_set1_epi32( 32 ) );
   vtmp          = _mm_srli_epi32( vtmp, 6 );
   vtmp          = _mm_packs_epi32( vtmp, vtmp );
   vtmp          = _mm_min_epi16( _mm_max_epi16( vtmp, vmin ), vmax );
   _mm_storel_epi64( ( __m128i* ) &srcQ[stride * pos], vtmp );
 }
}

static void
filter_chroma_strong_c_v(OVSample* src, const int stride, const int tc, uint8_t is_ctb_b)
{
    int j;
    if (is_ctb_b) {
        for (j = 0; j < 2; ++j) {
            const int16_t p3 = src[-stride * 4];
            const int16_t p2 = src[-stride * 3];
            const int16_t p1 = src[-stride * 2];
            const int16_t p0 = src[-stride    ];
            const int16_t q0 = src[0          ];
            const int16_t q1 = src[ stride    ];
            const int16_t q2 = src[ stride * 2];
            const int16_t q3 = src[ stride * 3];

            src[-1 * stride] = ov_clip((((p0 + p1) + (p1 + p1) + (p0 + q0) + (q1 + q2) + 4) >> 3), p0 - tc, p0 + tc);
            src[ 0 * stride] = ov_clip((((q0 + q3) + (p1 + p1) + (p0 + q0) + (q1 + q2) + 4) >> 3), q0 - tc, q0 + tc);
            src[ 1 * stride] = ov_clip((((q1 + p1) + (p0 + q0) + (q1 + q2) + (q3 + q3) + 4) >> 3), q1 - tc, q1 + tc);
            src[ 2 * stride] = ov_clip((((q2 + q3) + (p0 + q0) + (q1 + q2) + (q3 + q3) + 4) >> 3), q2 - tc, q2 + tc);
            src++;
        }
    } else {
        for (j = 0; j < 2; ++j) {
            const int16_t p3 = src[-stride * 4];
            const int16_t p2 = src[-stride * 3];
            const int16_t p1 = src[-stride * 2];
            const int16_t p0 = src[-stride    ];
            const int16_t q0 = src[0          ];
            const int16_t q1 = src[ stride    ];
            const int16_t q2 = src[ stride * 2];
            const int16_t q3 = src[ stride * 3];

            src[-3 * stride] = ov_clip((((p2 + p3) + (p3 + p3) + (p2 + p1) + (p0 + q0) + 4) >> 3), p2 - tc, p2 + tc);
            src[-2 * stride] = ov_clip((((p1 + q1) + (p3 + p3) + (p2 + p1) + (p0 + q0) + 4) >> 3), p1 - tc, p1 + tc);
            src[-1 * stride] = ov_clip((((p0 + p3) + (p2 + p1) + (p0 + q0) + (q1 + q2) + 4) >> 3), p0 - tc, p0 + tc);
            src[ 0 * stride] = ov_clip((((q0 + q3) + (p2 + p1) + (p0 + q0) + (q1 + q2) + 4) >> 3), q0 - tc, q0 + tc);
            src[ 1 * stride] = ov_clip((((q1 + p1) + (p0 + q0) + (q1 + q2) + (q3 + q3) + 4) >> 3), q1 - tc, q1 + tc);
            src[ 2 * stride] = ov_clip((((q2 + q3) + (p0 + q0) + (q1 + q2) + (q3 + q3) + 4) >> 3), q2 - tc, q2 + tc);
            src++;
        }
    }
}

static void
filter_chroma_strong_c_h(OVSample* src, const int stride, const int tc)
{

    const __m128i add_4 = _mm_set1_epi16(4);
    const __m128i tc_c = _mm_set1_epi16(tc);

    __m128i line0 = _mm_loadu_si128((__m128i*)(src - 4 +      0));
    __m128i line1 = _mm_loadu_si128((__m128i*)(src - 4 + stride));

    __m128i p0 = _mm_bslli_si128(line0, 6);
    __m128i q0 = _mm_bsrli_si128(line0, 6);
    __m128i p1 = _mm_bslli_si128(line1, 6);
    __m128i q1 = _mm_bsrli_si128(line1, 6);

    p0 = _mm_shufflelo_epi16(p0, 0xFF);
    q0 = _mm_shufflehi_epi16(q0, 0x00);
    p1 = _mm_shufflelo_epi16(p1, 0xFF);
    q1 = _mm_shufflehi_epi16(q1, 0x00);

    __m128i y0 = _mm_blend_epi16(p0, q0, 0x55);
    __m128i y1 = _mm_blend_epi16(p1, q1, 0x55);

    q0 = _mm_bsrli_si128(q0, 4);
    q1 = _mm_bsrli_si128(q1, 4);

    __m128i h0 = _mm_hadd_epi16(p0, q0);
    __m128i h1 = _mm_hadd_epi16(p1, q1);

    __m128i v0 = _mm_add_epi16(h0, _mm_bsrli_si128(h0, 2));
    __m128i v1 = _mm_add_epi16(h1, _mm_bsrli_si128(h1, 2));

    __m128i w0 = _mm_add_epi16(v0, _mm_bsrli_si128(h0, 4));
    __m128i w1 = _mm_add_epi16(v1, _mm_bsrli_si128(h1, 4));

    __m128i x0 = _mm_bsrli_si128(_mm_unpacklo_epi16(w0, w0), 2);
    __m128i x1 = _mm_bsrli_si128(_mm_unpacklo_epi16(w1, w1), 2);

    y0 = _mm_add_epi16(line0, y0);
    y1 = _mm_add_epi16(line1, y1);
    x0 = _mm_add_epi16(x0, y0);
    x1 = _mm_add_epi16(x1, y1);
    x0 = _mm_add_epi16(x0, add_4);
    x1 = _mm_add_epi16(x1, add_4);
    x0 = _mm_srli_epi16(x0, 3);
    x1 = _mm_srli_epi16(x1, 3);

    __m128i clp_min0 = _mm_subs_epu16(line0, tc_c);
    __m128i clp_max0 = _mm_adds_epu16(line0, tc_c);
    __m128i clp_min1 = _mm_subs_epu16(line1, tc_c);
    __m128i clp_max1 = _mm_adds_epu16(line1, tc_c);

    x0 = _mm_max_epi16(x0, clp_min0);
    x0 = _mm_min_epi16(x0, clp_max0);
    x1 = _mm_max_epi16(x1, clp_min1);
    x1 = _mm_min_epi16(x1, clp_max1);

    x0 = _mm_blend_epi16(x0, line0, 0x81);
    x1 = _mm_blend_epi16(x1, line1, 0x81);

    _mm_storeu_si128((__m128i *) &src[-4 +      0], x0);
    _mm_storeu_si128((__m128i *) &src[-4 + stride], x1);
}

static void
filter_chroma_weak_h(OVSample* src, const int stride, const int tc)
{
    int delta;
    int j;

    for (j = 0; j < 2; ++j) {
        const int16_t p1 = src[-2];
        const int16_t p0 = src[-1];
        const int16_t q0 = src[ 0];
        const int16_t q1 = src[ 1];

        delta = ov_clip((((q0 << 2) - (p0 << 2) + p1 - q1 + 4) >> 3), -tc, tc);
        src[-1] = ov_bdclip(p0 + delta);
        src[0]       = ov_bdclip(q0 - delta);
        src += stride;
    }
}

static void
filter_chroma_weak_v(OVSample* src, const int stride, const int tc)
{
    int delta;
    int j;

    for (j = 0; j < 2; ++j) {
        const int16_t p1 = src[-stride * 2];
        const int16_t p0 = src[-stride    ];
        const int16_t q0 = src[0          ];
        const int16_t q1 = src[ stride    ];

        delta = ov_clip((((q0 << 2) - (p0 << 2) + p1 - q1 + 4) >> 3), -tc, tc);
        src[-stride] = ov_bdclip(p0 + delta);
        src[0]       = ov_bdclip(q0 - delta);
        src++;
    }
}

void
rcn_init_df_functions_sse(struct RCNFunctions *const rcn_funcs)
{
  rcn_funcs->df.filter_v[0] = &filter_luma_strong_small_v;
  rcn_funcs->df.filter_v[1] = &filter_v_3_5;
  rcn_funcs->df.filter_v[2] = &filter_v_3_7;
  rcn_funcs->df.filter_v[4] = &filter_v_5_3;
  rcn_funcs->df.filter_v[5] = &filter_v_5_5;
  rcn_funcs->df.filter_v[6] = &filter_v_5_7;
  rcn_funcs->df.filter_v[8] = &filter_v_7_3;
  rcn_funcs->df.filter_v[9] = &filter_v_7_5;
  rcn_funcs->df.filter_v[10]= &filter_v_7_7;

  rcn_funcs->df.filter_h[0] = &filter_luma_strong_small_h;
  rcn_funcs->df.filter_h[1] = &filter_h_3_5;
  rcn_funcs->df.filter_h[2] = &filter_h_3_7;

  rcn_funcs->df.filter_h[4] = &filter_h_5_3;
  rcn_funcs->df.filter_h[5] = &filter_h_5_5;
  rcn_funcs->df.filter_h[6] = &filter_h_5_7;

  rcn_funcs->df.filter_h[8] = &filter_h_7_3;
  rcn_funcs->df.filter_h[9] = &filter_h_7_5;
  rcn_funcs->df.filter_h[10]= &filter_h_7_7;
  rcn_funcs->df.filter_weak_h = filter_luma_weak_h;
  rcn_funcs->df.filter_weak_v = filter_luma_weak_v;
  rcn_funcs->df.filter_weak_h_c = filter_chroma_weak_h;
  rcn_funcs->df.filter_weak_v_c = filter_chroma_weak_v;
  rcn_funcs->df.filter_strong_h_c = filter_chroma_strong_c_h;
  rcn_funcs->df.filter_strong_v_c = filter_chroma_strong_c_v;
}
