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
#include <stdlib.h>
#include "ovutils.h"
#include "dec_structures.h"
#include "ctudec.h"
#include "rcn_structures.h"

#define DEFAULT_INTRA_TC_OFFSET 2 ///< Default intra TC offset
#define MAX_QP 64

#include "bitdepth.h"

static inline uint16_t ov_abs(uint16_t x)
{
     uint16_t msk = -((int16_t)x < 0);
     return (x + msk) ^ msk;
}

static const uint16_t tc_lut[MAX_QP + DEFAULT_INTRA_TC_OFFSET] =
{
      0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
      3,   4,   4,   4,   4,   5,   5,   5,
      5,   7,   7,   8,   9,  10,  10,  11,
     13,  14,  15,  17,  19,  21,  24,  25,
     29,  33,  36,  41,  45,  51,  57,  64,
     71,  80,  89, 100, 112, 125, 141, 157,
    177, 198, 222, 250, 280, 314, 352, 395
};

static const uint8_t beta_lut[MAX_QP] =
{
     0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,
     6,  7,  8,  9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 20, 22, 24,
    26, 28, 30, 32, 34, 36, 38, 40,
    42, 44, 46, 48, 50, 52, 54, 56,
    58, 60, 62, 64, 66, 68, 70, 72,
    74, 76, 78, 80, 82, 84, 86, 88
};

static inline uint8_t
use_strong_filter_l0(const OVSample* src, const int stride, const int beta, const int tc, int max_l_p, int max_l_q)
{
    const int16_t p3 = src[-stride * 4];
    const int16_t p0 = src[-stride    ];
    const int16_t q0 = src[ 0         ];
    const int16_t q3 = src[ stride * 3];

    int sp3 = ov_abs(p3 - p0);
    int sq3 = ov_abs(q3 - q0);

    /*FIXME we could use branchless to derive eiher mxmy coordinate */
    if (max_l_p == 7) {
        const int16_t p7 = src[-stride * 8];
        const int16_t p6 = src[-stride * 7];
        const int16_t p5 = src[-stride * 6];
        const int16_t p4 = src[-stride * 5];
        sp3 += ov_abs((p4 - p5) - p6 + p7);
        sp3 += ov_abs(p3 - p7) + 1;
        sp3 >>= 1;
    } else if (max_l_p == 5) {
        /* max_l_p == 5 */
        const int16_t p5 = src[-stride * 6];
        sp3 += ov_abs(p3 - p5) + 1;
        sp3 >>= 1;
    }

    if (max_l_q == 7) {
        const int16_t q4  = src[stride * 4];
        const int16_t q5  = src[stride * 5];
        const int16_t q6 = src[stride * 6];
        const int16_t q7 = src[stride * 7];
        sq3 += ov_abs((q4 - q5) - q6 + q7);
        sq3 += ov_abs(q7 - q3) + 1;
        sq3 >>= 1;
    } else if (max_l_q == 5) {
        const int16_t q5  = src[stride * 5];
        sq3 += ov_abs(q5 - q3) + 1;
        sq3 >>= 1;
    }

    return ((sp3 + sq3) < (beta * 3 >> 5)) && (ov_abs(p0 - q0) < ((tc * 5 + 1) >> 1));
}

static inline uint8_t
use_strong_filter_l1(const OVSample* src, const int stride, const int beta, const int tc)
{
    const int16_t p3 = src[-stride * 4];
    const int16_t p0 = src[-stride    ];
    const int16_t q0 = src[ 0         ];
    const int16_t q3 = src[ stride * 3];

    int sp3 = ov_abs(p3 - p0);
    int sq3 = ov_abs(q3 - q0);

    const int d_strong = sp3 + sq3;

    return ((d_strong < (beta >> 3)) && (ov_abs(p0 - q0) < ((tc * 5 + 1) >> 1)));
}

/* FIXME Macros ? */
static inline uint16_t
compute_dp_c(const OVSample* src, const int stride , const uint8_t is_ctb_b)
{
    return ov_abs((int16_t)src[-stride * (3 - is_ctb_b)] - 2 * (int16_t)src[-stride * 2] + (int16_t)src[-stride]);
}

static inline uint16_t
compute_dp(const OVSample* src, const int stride)
{
    const int16_t p2 = src[-stride * 3];
    const int16_t p1 = src[-stride * 2];
    const int16_t p0 = src[-stride    ];
    return ov_abs((p2 - p1) + (p0 -p1));
}

static inline uint16_t
compute_dq(const OVSample* src, const int stride)
{
    const int16_t q0 = src[0         ];
    const int16_t q1 = src[stride * 1];
    const int16_t q2 = src[stride * 2];
    return ov_abs((q0 - q1) + (q2 - q1));
}

struct DBFParams{
    int beta;
    int tc;
};

/* FIXME This could be done once at init if not if delta QP is disabled */
/* Note QP is average QP between the two blocks instead of storing qp
 * we might be able to store beta and tc indexes directly
 */
static struct DBFParams
compute_dbf_limits(const struct DBFInfo *const dbf_info, int qp, int bs, int8_t ladf_qp_offfset, uint8_t comp_idx)
{
    int beta_offset  = dbf_info->beta_offset;
    int tc_offset    = dbf_info->tc_offset;

    if (comp_idx == 1) {
        beta_offset  = dbf_info->beta_offset_cb;
        tc_offset    = dbf_info->tc_offset_cb;
    } else if (comp_idx == 2) {
        beta_offset  = dbf_info->beta_offset_cr;
        tc_offset    = dbf_info->tc_offset_cr;
    }
    /*FIXME add LADF handling */
    qp += ladf_qp_offfset;

    const int tc_idx   = ov_clip((qp + DEFAULT_INTRA_TC_OFFSET * (bs - 1) + tc_offset), 0, 63 + DEFAULT_INTRA_TC_OFFSET);
    const int beta_idx = ov_clip(qp + beta_offset, 0, 63);

    #if (BITDEPTH - 10) >= 0
    const int tc = tc_lut[tc_idx] << (BITDEPTH - 10);
    #else
    const int tc = (tc_lut[tc_idx] + (1 << (9 - BITDEPTH))) >> (10 - BITDEPTH);
    #endif

    const int beta = beta_lut[beta_idx] << (BITDEPTH - 8);
    struct DBFParams dbf_params;
    dbf_params.beta = beta;
    dbf_params.tc = tc;

    return dbf_params;
}

static inline uint64_t
derive_size_3_map(const uint64_t *edge_map)
{
    uint64_t is_large_map = edge_map[0];
    is_large_map |= edge_map[1];
    is_large_map |= edge_map[2];
    is_large_map |= edge_map[3];
    is_large_map |= edge_map[4];
    is_large_map |= edge_map[5];
    is_large_map |= edge_map[6];
    return ~is_large_map;
}

static inline uint8_t
derive_filter_idx(int filter_l_p, int filter_l_q)
{
    uint8_t p_idx = (filter_l_p >> 1) - 1;
    uint8_t q_idx = (filter_l_q >> 1) - 1;
    return (((p_idx & 0x3) << 2) | (q_idx & 0x3));
}

static void
filter_h_7_7(OVSample *src, const int stride, const int tc)
{
    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db7_p[7] = { 5, 14, 23, 32, 41, 50, 59 };
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1 };
    static const int8_t tc7_p[7] = { 1, 1, 2, 3, 4, 5, 6 };

    int16_t clp_q[7];
    int16_t clp_p[7];

    clp_q[0] = (tc * tc7_q[0]) >> 1;
    clp_q[1] = (tc * tc7_q[1]) >> 1;
    clp_q[2] = (tc * tc7_q[2]) >> 1;
    clp_q[3] = (tc * tc7_q[3]) >> 1;
    clp_q[4] = (tc * tc7_q[4]) >> 1;
    clp_q[5] = (tc * tc7_q[5]) >> 1;
    clp_q[6] = (tc * tc7_q[6]) >> 1;

    clp_p[0] = clp_q[6];
    clp_p[1] = clp_q[5];
    clp_p[2] = clp_q[4];
    clp_p[3] = clp_q[3];
    clp_p[4] = clp_q[2];
    clp_p[5] = clp_q[1];
    clp_p[6] = clp_q[0];

    for (int i = 0; i < 4; i++) {
        OVSample* srcP = src - 1;
        OVSample* srcQ = src;
        OVSample* _src = src - 7;

        int ref_p = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
        int ref_q = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

        int db_ref = (2 * (srcP[0] + srcQ[0])
                      + srcP[-1] + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1] + srcP[-6 * 1]
                      + srcQ[ 1] + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + srcQ[ 6 * 1] + 8) >> 4;

        _src[0] = ov_clip(((db_ref * db7_p[0] + ref_p * (64 - db7_p[0]) + 32) >> 6), _src[0] - clp_p[0], _src[0] + clp_p[0]);
        _src[1] = ov_clip(((db_ref * db7_p[1] + ref_p * (64 - db7_p[1]) + 32) >> 6), _src[1] - clp_p[1], _src[1] + clp_p[1]);
        _src[2] = ov_clip(((db_ref * db7_p[2] + ref_p * (64 - db7_p[2]) + 32) >> 6), _src[2] - clp_p[2], _src[2] + clp_p[2]);
        _src[3] = ov_clip(((db_ref * db7_p[3] + ref_p * (64 - db7_p[3]) + 32) >> 6), _src[3] - clp_p[3], _src[3] + clp_p[3]);
        _src[4] = ov_clip(((db_ref * db7_p[4] + ref_p * (64 - db7_p[4]) + 32) >> 6), _src[4] - clp_p[4], _src[4] + clp_p[4]);
        _src[5] = ov_clip(((db_ref * db7_p[5] + ref_p * (64 - db7_p[5]) + 32) >> 6), _src[5] - clp_p[5], _src[5] + clp_p[5]);
        _src[6] = ov_clip(((db_ref * db7_p[6] + ref_p * (64 - db7_p[6]) + 32) >> 6), _src[6] - clp_p[6], _src[6] + clp_p[6]);

        src[0] = ov_clip(((db_ref * db7_q[0] + ref_q * (64 - db7_q[0]) + 32) >> 6), src[0] - clp_q[0], src[0] + clp_q[0]);
        src[1] = ov_clip(((db_ref * db7_q[1] + ref_q * (64 - db7_q[1]) + 32) >> 6), src[1] - clp_q[1], src[1] + clp_q[1]);
        src[2] = ov_clip(((db_ref * db7_q[2] + ref_q * (64 - db7_q[2]) + 32) >> 6), src[2] - clp_q[2], src[2] + clp_q[2]);
        src[3] = ov_clip(((db_ref * db7_q[3] + ref_q * (64 - db7_q[3]) + 32) >> 6), src[3] - clp_q[3], src[3] + clp_q[3]);
        src[4] = ov_clip(((db_ref * db7_q[4] + ref_q * (64 - db7_q[4]) + 32) >> 6), src[4] - clp_q[4], src[4] + clp_q[4]);
        src[5] = ov_clip(((db_ref * db7_q[5] + ref_q * (64 - db7_q[5]) + 32) >> 6), src[5] - clp_q[5], src[5] + clp_q[5]);
        src[6] = ov_clip(((db_ref * db7_q[6] + ref_q * (64 - db7_q[6]) + 32) >> 6), src[6] - clp_q[6], src[6] + clp_q[6]);
        src += stride;
    }
}

static void
filter_h_7_5(OVSample *src, const int stride, const int tc)
{
    static const int8_t db7_p[7] = { 5, 14, 23, 32, 41, 50, 59 };
    static const int8_t tc7_p[7] = { 1, 1, 2, 3, 4, 5, 6 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6 };
    static const int8_t tc5_q[7] = { 6, 5, 4, 3, 2 };

    int16_t clp_p[7];
    int16_t clp_q[5];

    clp_p[0] = (tc * tc7_p[0]) >> 1;
    clp_p[1] = (tc * tc7_p[1]) >> 1;
    clp_p[2] = (tc * tc7_p[2]) >> 1;
    clp_p[3] = (tc * tc7_p[3]) >> 1;
    clp_p[4] = (tc * tc7_p[4]) >> 1;
    clp_p[5] = (tc * tc7_p[5]) >> 1;
    clp_p[6] = (tc * tc7_p[6]) >> 1;

    clp_q[0] = clp_p[6];
    clp_q[1] = clp_p[5];
    clp_q[2] = clp_p[4];
    clp_q[3] = clp_p[3];
    clp_q[4] = clp_p[2];

  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - 1;
    OVSample* srcQ = src;
    OVSample* _src = src - 7;

    int ref_p = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
    int ref_q = (srcQ[ 4 * 1] + srcQ[ 5 * 1] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcP[-1] + srcQ[0] + srcQ[ 1])
                + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1]
                + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + 8) >> 4;

        _src[0] = ov_clip(((db_ref * db7_p[0] + ref_p * (64 - db7_p[0]) + 32) >> 6), _src[0] - clp_p[0], _src[0] + clp_p[0]);
        _src[1] = ov_clip(((db_ref * db7_p[1] + ref_p * (64 - db7_p[1]) + 32) >> 6), _src[1] - clp_p[1], _src[1] + clp_p[1]);
        _src[2] = ov_clip(((db_ref * db7_p[2] + ref_p * (64 - db7_p[2]) + 32) >> 6), _src[2] - clp_p[2], _src[2] + clp_p[2]);
        _src[3] = ov_clip(((db_ref * db7_p[3] + ref_p * (64 - db7_p[3]) + 32) >> 6), _src[3] - clp_p[3], _src[3] + clp_p[3]);
        _src[4] = ov_clip(((db_ref * db7_p[4] + ref_p * (64 - db7_p[4]) + 32) >> 6), _src[4] - clp_p[4], _src[4] + clp_p[4]);
        _src[5] = ov_clip(((db_ref * db7_p[5] + ref_p * (64 - db7_p[5]) + 32) >> 6), _src[5] - clp_p[5], _src[5] + clp_p[5]);
        _src[6] = ov_clip(((db_ref * db7_p[6] + ref_p * (64 - db7_p[6]) + 32) >> 6), _src[6] - clp_p[6], _src[6] + clp_p[6]);

        src[0] = ov_clip(((db_ref * db5_q[0] + ref_q * (64 - db5_q[0]) + 32) >> 6), src[0] - clp_q[0], src[0] + clp_q[0]);
        src[1] = ov_clip(((db_ref * db5_q[1] + ref_q * (64 - db5_q[1]) + 32) >> 6), src[1] - clp_q[1], src[1] + clp_q[1]);
        src[2] = ov_clip(((db_ref * db5_q[2] + ref_q * (64 - db5_q[2]) + 32) >> 6), src[2] - clp_q[2], src[2] + clp_q[2]);
        src[3] = ov_clip(((db_ref * db5_q[3] + ref_q * (64 - db5_q[3]) + 32) >> 6), src[3] - clp_q[3], src[3] + clp_q[3]);
        src[4] = ov_clip(((db_ref * db5_q[4] + ref_q * (64 - db5_q[4]) + 32) >> 6), src[4] - clp_q[4], src[4] + clp_q[4]);
    src += stride;
  }
}

static void
filter_h_5_7(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - 1;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};

    int ref_p = (srcP[-4 * 1] + srcP[-5 * 1] + 1) >> 1;
    int ref_q = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcP[-1] + srcQ[0] + srcQ[ 1])
            + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1]
            + srcQ[ 2 * 1] + srcQ[ 3 * 1] + srcQ[ 4 * 1] + srcQ[ 5 * 1] + 8) >> 4;

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-1 * pos];
        srcP[-1 * pos] = ov_clip(((db_ref * db5_q[pos] + ref_p * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ 1 * pos];
        srcQ[ 1 * pos] = ov_clip(((db_ref * db7_q[pos] + ref_q * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src += stride;
  }
}

static void
filter_h_5_5(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - 1;
    OVSample* srcQ = src;

    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t db_c5_p[5] = { 6, 19, 32, 45, 58 };
    static const int8_t tc5[5] = { 6, 5, 4, 3, 2 };
    static const int8_t tc5_p[5] = { 2, 3, 4, 5, 6 };

    int ref_p = (srcP[-4 * 1] + srcP[-5 * 1] + 1) >> 1;
    int ref_q = (srcQ[ 4 * 1] + srcQ[ 5 * 1] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcP[-1] + srcP[-2 * 1]
                        + srcQ[0] + srcQ[ 1] + srcQ[ 2 * 1])
              + srcP[-3 * 1] + srcP[-4 * 1]
              + srcQ[ 3 * 1] + srcQ[ 4 * 1] + 8) >> 4;


    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc5[pos]) >> 1;
        int32_t val = srcP[-1 * pos];
        srcP[-1 * pos] = ov_clip(((db_ref * db5_q[pos] + ref_p * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc5[pos]) >> 1;
        int32_t val = srcQ[ 1 * pos];
        srcQ[ 1 * pos] = ov_clip(((db_ref * db5_q[pos] + ref_q * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src += stride;
  }
}

static void
filter_h_7_3(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - 1;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };

    int ref_p = (srcP[-6 * 1] + srcP[-7 * 1] + 1) >> 1;
    int ref_q = (srcQ[ 2 * 1] + srcQ[ 3 * 1] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcQ[0])
            + srcP[     -1] + srcP[-2 * 1] + srcP[-3 * 1] + srcP[-4 * 1] + srcP[-5 * 1] + srcP[-6 * 1]
            + srcQ[      0] + srcQ[     1] + srcQ[     1] + srcQ[2 *  1] + srcQ[ 2 * 1] + srcQ[     1] + 8) >> 4;

    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-1 * pos];
        srcP[-1 * pos] = ov_clip(((db_ref * db7_q[pos] + ref_p * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcQ[ 1 * pos];
        srcQ[ 1 * pos] = ov_clip(((db_ref * db3_q[pos] + ref_q * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src += stride;
  }
}

static void
filter_h_3_7(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - 1;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };


    int ref_p = (srcP[-2 * 1] + srcP[-3 * 1] + 1) >> 1;
    int ref_q = (srcQ[ 6 * 1] + srcQ[ 7 * 1] + 1) >> 1;

    int db_ref = (2 * (srcQ[0] + srcP[0])
            + srcP[      0] + srcP[    -1] + srcP[    -1] + srcP[-2 * 1] + srcP[-2 * 1] + srcP[    -1]
            + srcQ[      1] + srcQ[2 *  1] + srcQ[3 *  1] + srcQ[4 *  1] + srcQ[5 *  1] + srcQ[6 *  1] + 8) >> 4;

    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcP[-1 * pos];
        srcP[-1 * pos] = ov_clip(((db_ref * db3_q[pos] + ref_p * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ 1 * pos];
        srcQ[ 1 * pos] = ov_clip(((db_ref * db7_q[pos] + ref_q * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src += stride;
  }
}

static void
filter_h_5_3(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - 1;
    OVSample* srcQ = src;

    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };

    int ref_p = (srcP[-4 * 1] + srcP[-5 * 1] + 1) >> 1;
    int ref_q = (srcQ[ 2 * 1] + srcQ[ 3 * 1] + 1) >> 1;

    int db_ref = (srcP[0] + srcP[-1] + srcP[-2 * 1] + srcP[-3 * 1]
                   + srcQ[0] + srcQ[ 1] + srcQ[ 2 * 1] + srcQ[ 3 * 1] + 4) >> 3;

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-1 * pos];
        srcP[-1 * pos] = ov_clip(((db_ref * db5_q[pos] + ref_p * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcQ[ 1 * pos];
        srcQ[ 1 * pos] = ov_clip(((db_ref * db3_q[pos] + ref_q * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src += stride;
  }
}

static void
filter_h_3_5(OVSample *src, const int stride, const int tc)
{
    static const int8_t db_c[8] = { 11, 32, 53, 58, 45, 32, 19, 6};
    static const int8_t tc_c[8] = {  2,  4,  6,  6,  5,  4,  3, 2};

    int16_t tc_val[8];

    tc_val[0] = (tc * tc_c[0]) >> 1;
    tc_val[1] = (tc * tc_c[1]) >> 1;
    tc_val[2] = (tc * tc_c[2]) >> 1;
    tc_val[3] = (tc * tc_c[3]) >> 1;
    tc_val[4] = (tc * tc_c[4]) >> 1;
    tc_val[5] = (tc * tc_c[5]) >> 1;
    tc_val[6] = (tc * tc_c[6]) >> 1;
    tc_val[7] = (tc * tc_c[7]) >> 1;

    src -= 3;

    for (int i = 0; i < 4; i++) {
        int ref_p = (src[0] + src[-1] + 1) >> 1;
        int ref_q = (src[7] + src[ 8] + 1) >> 1;

        int db_ref = (src[2] + src[1] + src[0] + src[-1]
                    + src[3] + src[4] + src[5] + src[ 6] + 4) >> 3;

        int16_t clp_min[8];
        int16_t clp_max[8];
        int16_t db_val[8];

        clp_min[0] = (int)src[0] - tc_val[0];
        clp_min[1] = (int)src[1] - tc_val[1];
        clp_min[2] = (int)src[2] - tc_val[2];
        clp_min[3] = (int)src[3] - tc_val[3];
        clp_min[4] = (int)src[4] - tc_val[4];
        clp_min[5] = (int)src[5] - tc_val[5];
        clp_min[6] = (int)src[6] - tc_val[6];
        clp_min[7] = (int)src[7] - tc_val[7];

        clp_max[0] = (int)src[0] + tc_val[0];
        clp_max[1] = (int)src[1] + tc_val[1];
        clp_max[2] = (int)src[2] + tc_val[2];
        clp_max[3] = (int)src[3] + tc_val[3];
        clp_max[4] = (int)src[4] + tc_val[4];
        clp_max[5] = (int)src[5] + tc_val[5];
        clp_max[6] = (int)src[6] + tc_val[6];
        clp_max[7] = (int)src[7] + tc_val[7];

        db_val[0] = (db_ref * db_c[0] + ref_p * (64 - db_c[0]) + 32) >> 6;
        db_val[1] = (db_ref * db_c[1] + ref_p * (64 - db_c[1]) + 32) >> 6;
        db_val[2] = (db_ref * db_c[2] + ref_p * (64 - db_c[2]) + 32) >> 6;
        db_val[3] = (db_ref * db_c[3] + ref_q * (64 - db_c[3]) + 32) >> 6;
        db_val[4] = (db_ref * db_c[4] + ref_q * (64 - db_c[4]) + 32) >> 6;
        db_val[5] = (db_ref * db_c[5] + ref_q * (64 - db_c[5]) + 32) >> 6;
        db_val[6] = (db_ref * db_c[6] + ref_q * (64 - db_c[6]) + 32) >> 6;
        db_val[7] = (db_ref * db_c[7] + ref_q * (64 - db_c[7]) + 32) >> 6;

        src[0] = ov_clip(db_val[0], clp_min[0], clp_max[0]);
        src[1] = ov_clip(db_val[1], clp_min[1], clp_max[1]);
        src[2] = ov_clip(db_val[2], clp_min[2], clp_max[2]);
        src[3] = ov_clip(db_val[3], clp_min[3], clp_max[3]);
        src[4] = ov_clip(db_val[4], clp_min[4], clp_max[4]);
        src[5] = ov_clip(db_val[5], clp_min[5], clp_max[5]);
        src[6] = ov_clip(db_val[6], clp_min[6], clp_max[6]);
        src[7] = ov_clip(db_val[7], clp_min[7], clp_max[7]);

        src += stride;
    }
}

static void
filter_v_7_7(OVSample *src, const int stride, const int tc)
{
    for (int i = 0; i < 4; i++) {
      OVSample* srcP = src - stride;
      OVSample* srcQ = src;

      static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
      static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};

      int ref_p = (srcP[-6 * stride] + srcP[-7 * stride] + 1) >> 1;
      int ref_q = (srcQ[ 6 * stride] + srcQ[ 7 * stride] + 1) >> 1;

      int db_ref = (2 * (srcP[0] + srcQ[0])
                      + srcP[-stride] + srcP[-2 * stride] + srcP[-3 * stride] + srcP[-4 * stride] + srcP[-5 * stride] + srcP[-6 * stride]
                      + srcQ[ stride] + srcQ[ 2 * stride] + srcQ[ 3 * stride] + srcQ[ 4 * stride] + srcQ[ 5 * stride] + srcQ[ 6 * stride] + 8) >> 4;

      for (int pos = 0; pos < 7; pos++) {
          int cvalue = (tc * tc7_q[pos]) >> 1;
          int32_t val = srcP[-stride * pos];
          srcP[-stride * pos] = ov_clip(((db_ref * db7_q[pos] + ref_p * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
      }

      for (int pos = 0; pos < 7; pos++) {
          int cvalue = (tc * tc7_q[pos]) >> 1;
          int32_t val = srcQ[ stride * pos];
          srcQ[ stride * pos] = ov_clip(((db_ref * db7_q[pos] + ref_q * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
      }
      src+=1;
    }
}

static void
filter_v_7_5(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {
    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};

    int ref_p = (srcP[-6 * stride] + srcP[-7 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 4 * stride] + srcQ[ 5 * stride] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcP[-stride] + srcQ[0] + srcQ[ stride])
                + srcP[-2 * stride] + srcP[-3 * stride] + srcP[-4 * stride] + srcP[-5 * stride]
                + srcQ[ 2 * stride] + srcQ[ 3 * stride] + srcQ[ 4 * stride] + srcQ[ 5 * stride] + 8) >> 4;


    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db7_q[pos] + ref_p * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db5_q[pos] + ref_q * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
  }
}

static void
filter_v_5_7(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {

    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};

    int ref_p = (srcP[-4 * stride] + srcP[-5 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 6 * stride] + srcQ[ 7 * stride] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcP[-stride] + srcQ[0] + srcQ[ stride])
            + srcP[-2 * stride] + srcP[-3 * stride] + srcP[-4 * stride] + srcP[-5 * stride]
            + srcQ[ 2 * stride] + srcQ[ 3 * stride] + srcQ[ 4 * stride] + srcQ[ 5 * stride] + 8) >> 4;

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db5_q[pos] + ref_p * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db7_q[pos] + ref_q * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
  }
}

static void
filter_v_5_5(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {

    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};

    int ref_p = (srcP[-4 * stride] + srcP[-5 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 4 * stride] + srcQ[ 5 * stride] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcP[-stride] + srcP[-2 * stride]
                        + srcQ[0] + srcQ[ stride] + srcQ[ 2 * stride])
                        + srcP[-3 * stride] + srcP[-4 * stride]
                        + srcQ[ 3 * stride] + srcQ[ 4 * stride] + 8) >> 4;


    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db5_q[pos] + ref_p * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db5_q[pos] + ref_q * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
  }
}

static void
filter_v_7_3(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {

    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };


    int ref_p = (srcP[-6 * stride] + srcP[-7 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 2 * stride] + srcQ[ 3 * stride] + 1) >> 1;

    int db_ref = (2 * (srcP[0] + srcQ[0])
            + srcP[-stride] + srcP[-2 * stride] + srcP[-3 * stride] + srcP[-4 * stride] + srcP[-5 * stride] + srcP[-6 * stride]
            + srcQ[      0] + srcQ[     stride] + srcQ[     stride] + srcQ[2 *  stride] + srcQ[ 2 * stride] + srcQ[     stride] + 8) >> 4;

    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db7_q[pos] + ref_p * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db3_q[pos] + ref_q * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
  }
}

static void
filter_v_3_7(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {

    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db7_q[7] = { 59, 50, 41, 32, 23, 14, 5 };
    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };


    int ref_p = (srcP[-2 * stride] + srcP[-3 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 6 * stride] + srcQ[ 7 * stride] + 1) >> 1;

    int db_ref = (2 * (srcQ[0] + srcP[0])
            + srcP[      0] + srcP[    -stride] + srcP[    -stride] + srcP[-2 * stride] + srcP[-2 * stride] + srcP[    -stride]
            + srcQ[ stride] + srcQ[2 *  stride] + srcQ[3 *  stride] + srcQ[4 *  stride] + srcQ[5 *  stride] + srcQ[6 *  stride] + 8) >> 4;

    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db3_q[pos] + ref_p * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 7; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db7_q[pos] + ref_q * (64 - db7_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
  }
}

static void
filter_v_5_3(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {

    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };


    int ref_p = (srcP[-4 * stride] + srcP[-5 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 2 * stride] + srcQ[ 3 * stride] + 1) >> 1;

    int db_ref = (srcP[0] + srcP[-stride] + srcP[-2 * stride] + srcP[-3 * stride]
                   + srcQ[0] + srcQ[ stride] + srcQ[ 2 * stride] + srcQ[ 3 * stride] + 4) >> 3;

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db5_q[pos] + ref_p * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db3_q[pos] + ref_q * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
  }
}

static void
filter_v_3_5(OVSample *src, const int stride, const int tc)
{
  for (int i = 0; i < 4; i++) {

    OVSample* srcP = src - stride;
    OVSample* srcQ = src;

    static const int8_t db3_q[3] = { 53, 32, 11 };
    static const int8_t db3_p[3] = { 11, 32, 53 };
    static const int8_t db5_q[5] = { 58, 45, 32, 19, 6};
    static const int8_t tc7_q[7] = { 6, 5, 4, 3, 2, 1, 1};
    static const int8_t tc3_q[3] = { 6, 4, 2 };
    static const int8_t tc3_p[3] = { 2, 4, 6 };

    int ref_p = (srcP[-2 * stride] + srcP[-3 * stride] + 1) >> 1;
    int ref_q = (srcQ[ 4 * stride] + srcQ[ 5 * stride] + 1) >> 1;

    int db_ref = (srcP[0] + srcP[-stride] + srcP[-2 * stride] + srcP[-3 * stride]
                   + srcQ[0] + srcQ[ stride] + srcQ[ 2 * stride] + srcQ[ 3 * stride] + 4) >> 3;


    for (int pos = 0; pos < 3; pos++) {
        int cvalue = (tc * tc3_q[pos]) >> 1;
        int32_t val = srcP[-stride * pos];
        srcP[-stride * pos] = ov_clip(((db_ref * db3_q[pos] + ref_p * (64 - db3_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }

    for (int pos = 0; pos < 5; pos++) {
        int cvalue = (tc * tc7_q[pos]) >> 1;
        int32_t val = srcQ[ stride * pos];
        srcQ[ stride * pos] = ov_clip(((db_ref * db5_q[pos] + ref_q * (64 - db5_q[pos]) + 32) >> 6),val - cvalue, val + cvalue);
    }
    src+=1;
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

        src[-stride * 3] = ov_clip(((2 * p3 + 3 * p2 + p1 + p0 + q0 + 4) >> 3),     p2 - tc_c[0], p2 + tc_c[0]);
        src[-stride * 2] = ov_clip(((p2 + p1 + p0 + q0 + 2) >> 2),                  p1 - tc_c[1], p1 + tc_c[1]);
        src[-stride]     = ov_clip(((p2 + 2 * p1 + 2 * p0 + 2 * q0 + q1 + 4) >> 3), p0 - tc_c[2], p0 + tc_c[2]);
        src[0]           = ov_clip(((p1 + 2 * p0 + 2 * q0 + 2 * q1 + q2 + 4) >> 3), q0 - tc_c[3], q0 + tc_c[3]);
        src[stride]      = ov_clip(((p0 + q0 + q1 + q2 + 2) >> 2),                  q1 - tc_c[4], q1 + tc_c[4]);
        src[stride * 2]  = ov_clip(((p0 + q0 + q1 + 3 * q2 + 2 * q3 + 4) >> 3),     q2 - tc_c[5], q2 + tc_c[5]);
        src++;
    }
}

static void
filter_luma_strong_small_h(OVSample* src, const int stride, const int tc)
{

    int16_t tc_c[6] = { 1 * tc, 2 * tc, 3 * tc , 3 * tc, 2 * tc, 1 * tc};

    for (int i = 0; i < 4; i++) {
        const int16_t p3  = src[-4];
        const int16_t p2  = src[-3];
        const int16_t p1  = src[-2];
        const int16_t p0  = src[-1];
        const int16_t q0  = src[ 0];
        const int16_t q1  = src[ 1];
        const int16_t q2  = src[ 2];
        const int16_t q3  = src[ 3];

        src[-3] = ov_clip(((2 * p3 + 3 * p2 + p1 + p0 + q0 + 4) >> 3),     p2 - tc_c[0], p2 + tc_c[0]);
        src[-2] = ov_clip(((p2 + p1 + p0 + q0 + 2) >> 2),                  p1 - tc_c[1], p1 + tc_c[1]);
        src[-1] = ov_clip(((p2 + 2 * p1 + 2 * p0 + 2 * q0 + q1 + 4) >> 3), p0 - tc_c[2], p0 + tc_c[2]);
        src[0]  = ov_clip(((p1 + 2 * p0 + 2 * q0 + 2 * q1 + q2 + 4) >> 3), q0 - tc_c[3], q0 + tc_c[3]);
        src[1]  = ov_clip(((p0 + q0 + q1 + q2 + 2) >> 2),                  q1 - tc_c[4], q1 + tc_c[4]);
        src[2]  = ov_clip(((p0 + q0 + q1 + 3 * q2 + 2 * q3 + 4) >> 3),     q2 - tc_c[5], q2 + tc_c[5]);
        src += stride;
    }
}

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
        int delta = (9 * (q0 - p0) - 3 * (q1 - p1) + 8) >> 4;

        if (ov_abs(delta) < th_cut) {
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
filter_luma_weak_v(OVSample* src, const int stride, const int tc, const uint8_t extend_p, const uint8_t extend_q)
{
    const int th_cut = tc * 10;
    const int tc2_p = -extend_p & (tc >> 1);
    const int tc2_q = -extend_q & (tc >> 1);
    for (int i = 0; i < 4; i++) {
        const int16_t p2  = src[-stride * 3];
        const int16_t p1  = src[-stride * 2];
        const int16_t p0  = src[-stride    ];
        const int16_t q0  = src[ 0         ];
        const int16_t q1  = src[ stride    ];
        const int16_t q2  = src[ stride * 2];

        /* Weak filter */
        int delta = (9 * (q0 - p0) - 3 * (q1 - p1) + 8) >> 4;

        if (ov_abs(delta) < th_cut) {
            delta = ov_clip(delta, -tc, tc);
            const int delta1 = ov_clip(((((p2 + p0 + 1) >> 1) - p1 + delta) >> 1), -tc2_p, tc2_p);
            const int delta2 = ov_clip(((((q2 + q0 + 1) >> 1) - q1 - delta) >> 1), -tc2_q, tc2_q);
            src[stride * -2] = ov_bdclip(p1 + delta1);
            src[stride * -1] = ov_bdclip(p0 + delta);
            src[stride * 0]  = ov_bdclip(q0 - delta);
            src[stride * 1]  = ov_bdclip(q1 + delta2);
        }
        ++src;
    }
}


static inline uint8_t
use_strong_filter_c2(const OVSample* src, const int stride, const int beta, const int tc, uint8_t is_ctb_b)
{
    const int16_t p3 = src[(-stride * 4) >> is_ctb_b];
    const int16_t p0 = src[-stride    ];
    const int16_t q0 = src[ 0         ];
    const int16_t q3 = src[ stride * 3];

    int sp3 = ov_abs(p3 - p0);
    int sq3 = ov_abs(q3 - q0);

    const int d_strong = sp3 + sq3;

    return ((d_strong < (beta >> 3)) && (ov_abs(p0 - q0) < ((tc * 5 + 1) >> 1)));
}

static inline uint8_t
use_strong_filter_c(const OVSample* src, const int stride, const int beta, const int tc)
{
    const int16_t p3 = src[-stride * 4];
    const int16_t p0 = src[-stride    ];
    const int16_t q0 = src[ 0         ];
    const int16_t q3 = src[ stride * 3];

    int sp3 = ov_abs(p3 - p0);
    int sq3 = ov_abs(q3 - q0);

    const int d_strong = sp3 + sq3;

    return ((d_strong < (beta >> 3)) && (ov_abs(p0 - q0) < ((tc * 5 + 1) >> 1)));
}

static void
filter_chroma_strong_c_v(OVSample* src, const int stride, const int tc, uint8_t is_ctb_b)
{
    int j;
    for (j = 0; j < 2; ++j) {
        const int16_t p3 = src[-stride * 4];
        const int16_t p2 = src[-stride * 3];
        const int16_t p1 = src[-stride * 2];
        const int16_t p0 = src[-stride    ];
        const int16_t q0 = src[0          ];
        const int16_t q1 = src[ stride    ];
        const int16_t q2 = src[ stride * 2];
        const int16_t q3 = src[ stride * 3];

        if (is_ctb_b) {
            src[-stride * 1] = ov_clip(((3 * p1 + 2 * p0 + q0 + q1 + q2 + 4) >> 3),p0 - tc, p0 + tc); // p)0
            src[0          ] = ov_clip(((2 * p1 + p0 + 2 * q0 + q1 + q2 + q3 + 4) >> 3),q0 - tc, q0 + tc); // q)0
            src[ stride * 1] = ov_clip(((p1 + p0 + q0 + 2 * q1 + q2 + 2 * q3 + 4) >> 3),q1 - tc, q1 + tc); // q)1
            src[ stride * 2] = ov_clip(((p0 + q0 + q1 + 2 * q2 + 3 * q3 + 4) >> 3),q2 - tc, q2 + tc);      // q)2
        } else {
            src[-stride * 3] = ov_clip(((3 * p3 + 2 * p2 + p1 + p0 + q0 + 4) >> 3),p2 - tc, p2 + tc);      // p)2
            src[-stride * 2] = ov_clip(((2 * p3 + p2 + 2 * p1 + p0 + q0 + q1 + 4) >> 3),p1 - tc, p1 + tc); // p)1
            src[-stride * 1] = ov_clip(((p3 + p2 + p1 + 2 * p0 + q0 + q1 + q2 + 4) >> 3),p0 - tc, p0 + tc); // p)0
            src[ 0         ] = ov_clip(((p2 + p1 + p0 + 2 * q0 + q1 + q2 + q3 + 4) >> 3),q0 - tc, q0 + tc); // q)0
            src[ stride * 1] = ov_clip(((p1 + p0 + q0 + 2 * q1 + q2 + 2 * q3 + 4) >> 3),q1 - tc, q1 + tc);  // q)1
            src[ stride * 2] = ov_clip(((p0 + q0 + q1 + 2 * q2 + 3 * q3 + 4) >> 3),q2 - tc, q2 + tc);       // q)2
        }
        src++;
    }
}

static void
filter_chroma_strong_c_h(OVSample* src, const int stride, const int tc)
{
    int j;
    for (j = 0; j < 2; ++j) {
        const int16_t p3 = src[-1 * 4];
        const int16_t p2 = src[-1 * 3];
        const int16_t p1 = src[-1 * 2];
        const int16_t p0 = src[-1    ];
        const int16_t q0 = src[0          ];
        const int16_t q1 = src[ 1    ];
        const int16_t q2 = src[ 1 * 2];
        const int16_t q3 = src[ 1 * 3];

        src[-1 * 3] = ov_clip(((3 * p3 + 2 * p2 + p1 + p0 + q0 + 4) >> 3),p2 - tc, p2 + tc);      // p)2
        src[-1 * 2] = ov_clip(((2 * p3 + p2 + 2 * p1 + p0 + q0 + q1 + 4) >> 3),p1 - tc, p1 + tc); // p)1
        src[-1 * 1] = ov_clip(((p3 + p2 + p1 + 2 * p0 + q0 + q1 + q2 + 4) >> 3),p0 - tc, p0 + tc); // p)0
        src[ 0    ] = ov_clip(((p2 + p1 + p0 + 2 * q0 + q1 + q2 + q3 + 4) >> 3),q0 - tc, q0 + tc); // q)0
        src[ 1 * 1] = ov_clip(((p1 + p0 + q0 + 2 * q1 + q2 + 2 * q3 + 4) >> 3),q1 - tc, q1 + tc);  // q)1
        src[ 1 * 2] = ov_clip(((p0 + q0 + q1 + 2 * q2 + 3 * q3 + 4) >> 3),q2 - tc, q2 + tc);       // q)2
        src += stride;
    }
}

static void
filter_chroma_weak_h(OVSample* src, const int stride, const int tc)
{
    int delta;
    int j;

    for (j = 0; j < 2; ++j) {
        const int16_t p1 = src[-1 * 2];
        const int16_t p0 = src[-1    ];
        const int16_t q0 = src[0          ];
        const int16_t q1 = src[ 1    ];

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

/* Check if filter is 3 or 1 sample large based on other left edges */
static uint64_t
derive_large_map_from_ngh(const uint64_t *src_map)
{
    const uint64_t *fwd = src_map + 1;
    const uint64_t *bwd = src_map - 1;
    uint64_t dst_map = 0;

    /* If an edge is detected on either of those maps
     * the corresponding bit in dst_map is set to one
     */
    dst_map |= bwd[ 0] | fwd[0];
    dst_map |= bwd[-1] | fwd[1];
    dst_map |= bwd[-2] | fwd[2];

    /* Revert the bits on the map so the that output
     * bits set to one correspond to large blocks edges
     */
    return ~dst_map;
}

static void
filter_veritcal_edge_c(const struct DFFunctions *df, const struct DBFInfo *const dbf_info, OVSample *src, ptrdiff_t stride,
                       int8_t qp, uint64_t bs2_map, uint64_t large_map_q, uint8_t comp_idx)
{
    const uint8_t is_large = large_map_q & 0x1;
    const uint8_t is_bs2   = bs2_map     & 0x1;

    /* Note there is no need to check anything here since condition
     * is already checked in edge_map
     */
    OVSample *src0 = src;
    OVSample *src1 = src + stride;

    const struct DBFParams dbf_params = compute_dbf_limits(dbf_info, qp, 1 + is_bs2, 0, comp_idx);
    if (dbf_params.tc == 0 || dbf_params.beta == 0) return;
    uint8_t is_strong = 0;

    if (is_large) {
        const int dp0 = compute_dp(src0, 1);
        const int dq0 = compute_dq(src0, 1);
        const int dp3 = compute_dp(src1, 1);
        const int dq3 = compute_dq(src1, 1);

        const int d0 = dp0 + dq0;
        const int d3 = dp3 + dq3;

        const int d = d0 + d3;

        is_strong = (d < dbf_params.beta) &&
            (2 * d0 < (dbf_params.beta >> 2)) &&
            (2 * d3 < (dbf_params.beta >> 2)) &&
            use_strong_filter_c(src0, 1, dbf_params.beta, dbf_params.tc) &&
            use_strong_filter_c(src1, 1, dbf_params.beta, dbf_params.tc);

    }

    if (is_strong) {
        df->filter_strong_h_c(src, stride, dbf_params.tc);
    } else {
        df->filter_weak_h_c(src, stride, dbf_params.tc);
    }
}

/* Filter vertical edges */
static void
vvc_dbf_chroma_hor(const struct DFFunctions *df, OVSample *src_cb, OVSample *src_cr, int stride,
                   const struct DBFInfo *const dbf_info,
                   uint8_t nb_unit_h, int is_last_h, uint8_t nb_unit_w,
                   uint8_t ctu_lft)
{

    const int blk_stride = stride << 1;
    /* Mask applied to edge_mask based on CTU height */
    const uint64_t vedge_mask = ((uint64_t)1 << nb_unit_h) - 1;
    const uint64_t *const edg_map_tab = &dbf_info->ctb_bound_ver_c[8];
    int i;

    /* The number of edges to process (should be ((1 << log2_ctu_s) >> 4) since we start
     * finish from CTU left border to the last edge before the right CTU border
     * (Note we know the CTU right border is an implicit edge so we can set it to 0xFF)
     */
    const uint8_t nb_vedge = ((nb_unit_w + 3) >> 2);
    const uint8_t skip_first = !ctu_lft;

    src_cb += skip_first << 3;
    src_cr += skip_first << 3;

    for (i = skip_first; i < nb_vedge; i++) {
        /* FIXME chroma_edges could be stored on a smaller grid */
        uint8_t edge_idx = i << 2;

        uint64_t bs2_map = dbf_info->bs2_map_c.ver  [edge_idx];
        uint64_t bs1_map = dbf_info->bs1_map_cb.ver [edge_idx];

        uint64_t edge_map = edg_map_tab[edge_idx];

        /* FIXME Use directly CTU start and modify maps storage and rotation
         */
        /* Discard first edge corresponding to upper CTU and mask
         * out of picture edge
         */
        edge_map &= vedge_mask;

        /* Discard non filtered edges from edge_map */
        edge_map &= bs2_map | bs1_map;

        if (edge_map) {
            uint64_t large_map_q = derive_large_map_from_ngh(&edg_map_tab[edge_idx]);

            /* FIXME use absolute QP maps */
            const int8_t *qp_col = &dbf_info->qp_map_cb.hor[36 + edge_idx];
            OVSample *src = src_cb;

            /* Discard non filtered edges from edge_map */
            edge_map &= bs2_map | (bs1_map & large_map_q);

            /* Note while is expected instead of do since edge_map can be zero
             * if edge_map consists only in bs1 with non large block sizes
             */
            while (edge_map) {
                uint8_t nb_skipped_blk = ov_ctz64(edge_map);
                int8_t qp;

                /* Skip non filtered edges */
                large_map_q >>= nb_skipped_blk;
                bs2_map     >>= nb_skipped_blk;
                qp_col       += nb_skipped_blk * 34;
                src          += nb_skipped_blk * blk_stride;

                qp = (qp_col[-1] + qp_col[0] + 1) >> 1;

                filter_veritcal_edge_c(df, dbf_info, src, stride, qp, bs2_map, large_map_q, 1);

                edge_map    >>= nb_skipped_blk + 1;
                large_map_q >>= 1;
                bs2_map     >>= 1;

                src += blk_stride;
                qp_col += 34;
            }
        }
        src_cb += 1 << 3;
    }

    for (i = skip_first; i < nb_vedge; i++) {
        uint8_t edge_idx = i << 2;

        uint64_t bs2_map = dbf_info->bs2_map_c.ver  [edge_idx];
        uint64_t bs1_map = dbf_info->bs1_map_cr.ver [edge_idx];

        uint64_t edge_map = edg_map_tab[edge_idx];

        edge_map &= vedge_mask;
        edge_map &= bs2_map | bs1_map;

        if (edge_map) {
            uint64_t large_map_q = derive_large_map_from_ngh(&edg_map_tab[edge_idx]);

            const int8_t *qp_col = &dbf_info->qp_map_cr.hor[36 + edge_idx];

            OVSample *src = src_cr;

            edge_map &= bs2_map | (bs1_map & large_map_q);

             while (edge_map) {
                uint8_t nb_skipped_blk = ov_ctz64(edge_map);
                int8_t qp;

                /* Skip non filtered edges */
                large_map_q >>= nb_skipped_blk;
                bs2_map     >>= nb_skipped_blk;
                qp_col       += nb_skipped_blk * 34;
                src          += nb_skipped_blk * blk_stride;

                qp = (qp_col[-1] + qp_col[0] + 1) >> 1;

                filter_veritcal_edge_c(df, dbf_info, src, stride, qp, bs2_map, large_map_q, 2);

                edge_map    >>= nb_skipped_blk + 1;
                large_map_q >>= 1;
                bs2_map     >>= 1;

                src += blk_stride;
                qp_col += 34;
            }
        }
        src_cr += 1 << 3;
    }
}

static void
filter_horizontal_edge_c(const struct DFFunctions *df, const struct DBFInfo *const dbf_info, OVSample *src, ptrdiff_t stride,
                         int8_t qp, uint64_t bs2_map, uint64_t large_map_q, uint8_t is_ctb_b, uint8_t comp_idx)
{
    const uint8_t is_large = large_map_q & 0x1;
    const uint8_t is_bs2   = bs2_map     & 0x1;

    /* Note there is no need to check anything here since condition
     * is already checked in edge_map
     */
    OVSample *src0 = src;
    OVSample *src1 = src + 1;

    const struct DBFParams dbf_params = compute_dbf_limits(dbf_info, qp, 1 + is_bs2, 0, comp_idx);
    if (dbf_params.tc == 0 || dbf_params.beta == 0) return;
    uint8_t is_strong = 0;

    if (is_large) {
        const int dp0 = compute_dp_c(src0, stride, is_ctb_b);
        const int dq0 = compute_dq(src0, stride);
        const int dp3 = compute_dp_c(src1, stride, is_ctb_b);
        const int dq3 = compute_dq(src1, stride);

        const int d0 = dp0 + dq0;
        const int d3 = dp3 + dq3;

        const int d = d0 + d3;

        is_strong = (d < dbf_params.beta) &&
            (2 * d0 < (dbf_params.beta >> 2)) &&
            (2 * d3 < (dbf_params.beta >> 2)) &&
            use_strong_filter_c2(src0, stride, dbf_params.beta, dbf_params.tc, is_ctb_b) &&
            use_strong_filter_c2(src1, stride, dbf_params.beta, dbf_params.tc, is_ctb_b);

        if (is_strong) {
            df->filter_strong_v_c(src, stride, dbf_params.tc, is_ctb_b);
        }
    }

    if (!is_strong) {
        df->filter_weak_v_c(src, stride, dbf_params.tc);
    }
}

static void
vvc_dbf_chroma_ver(const struct DFFunctions *df, OVSample *src_cb, OVSample *src_cr, int stride,
                   const struct DBFInfo *const dbf_info,
                   uint8_t nb_unit_w, int is_last_w, uint8_t nb_unit_h, uint8_t is_last_h,
                   uint8_t ctu_abv)
{
    const int blk_stride = 1 << 1;
    const uint64_t hedge_mask = (((uint64_t)1 << (nb_unit_w + (!!is_last_w << 1))) - 1);
    const uint8_t nb_hedge = ((nb_unit_h + 3) >> 2);
    const uint8_t skip_first = !ctu_abv;
    const uint64_t *const edg_map_tab = &dbf_info->ctb_bound_hor_c[8];
    int i;

    src_cb -= blk_stride << 1;
    src_cr -= blk_stride << 1;

    src_cb += (skip_first * stride) << 3;
    src_cr += (skip_first * stride) << 3;

    for (i = skip_first; i < nb_hedge; i++) {
        uint8_t edge_idx = i << 2;

        uint64_t edge_map = edg_map_tab[edge_idx];

        uint64_t bs2_map  = dbf_info->bs2_map_c.hor[edge_idx];
        uint64_t bs1_map  = dbf_info->bs1_map_cb.hor[edge_idx];

        edge_map &= hedge_mask;
        edge_map &= bs2_map | bs1_map;

        if (edge_map) {
            const int8_t *qp_row = &dbf_info->qp_map_cb.hor[edge_idx * 34];

            uint64_t large_map_q = derive_large_map_from_ngh(&edg_map_tab[edge_idx]);

            uint8_t is_ctb_b = i == 0;
            OVSample *src = src_cb;

            edge_map &= bs2_map | (bs1_map & large_map_q);

            while(edge_map) {
                uint8_t nb_skipped_blk = ov_ctz64(edge_map);
                int8_t qp;

                /* Skip non filtered edges */
                large_map_q >>= nb_skipped_blk;
                bs2_map     >>= nb_skipped_blk;
                qp_row       += nb_skipped_blk;
                src          += nb_skipped_blk * blk_stride;

                qp = (qp_row[0] + qp_row[34] + 1) >> 1;

                filter_horizontal_edge_c(df, dbf_info, src, stride, qp, bs2_map,
                                         large_map_q, is_ctb_b, 1);

                edge_map    >>= nb_skipped_blk + 1;
                large_map_q >>= 1;
                bs2_map     >>= 1;

                src += blk_stride;
                qp_row++;
            }
        }
        src_cb += stride << 3;
    }

    for (i = skip_first; i < nb_hedge; i++) {
        uint8_t edge_idx = i << 2;

        uint64_t edge_map = edg_map_tab[edge_idx];

        uint64_t bs2_map  = dbf_info->bs2_map_c.hor[edge_idx];
        uint64_t bs1_map  = dbf_info->bs1_map_cr.hor[edge_idx];

        edge_map &= hedge_mask;
        edge_map &= bs2_map | bs1_map;

        if (edge_map) {
            uint64_t large_map_q = derive_large_map_from_ngh(&edg_map_tab[edge_idx]);
            const int8_t *qp_row = &dbf_info->qp_map_cr.hor[edge_idx * 34];
            OVSample *src = src_cr;
            uint8_t is_ctb_b = i == 0;

            edge_map &= bs2_map | (bs1_map & large_map_q);

            while(edge_map) {
                uint8_t nb_skipped_blk = ov_ctz64(edge_map);
                int8_t qp;

                /* Skip non filtered edges */
                large_map_q >>= nb_skipped_blk;
                bs2_map     >>= nb_skipped_blk;
                qp_row       += nb_skipped_blk;
                src          += nb_skipped_blk * blk_stride;

                qp = (qp_row[0] + qp_row[34] + 1) >> 1;

                filter_horizontal_edge_c(df, dbf_info, src, stride, qp, bs2_map,
                                         large_map_q, is_ctb_b, 2);

                edge_map    >>= nb_skipped_blk + 1;
                large_map_q >>= 1;
                bs2_map     >>= 1;

                src += blk_stride;
                qp_row++;
            }
        }
        src_cr += stride << 3;
    }
}

static void
filter_vertical_edge(const struct DFFunctions *df, const struct DBFParams *const dbf_params,
                     OVSample *src, ptrdiff_t stride,
                     uint8_t max_l_p, uint8_t max_l_q)
{


    OVSample* src0 = src;
    OVSample* src3 = src + stride * 3;

    const int dp0 = compute_dp(src0, 1);
    const int dq0 = compute_dq(src0, 1);
    const int dp3 = compute_dp(src3, 1);
    const int dq3 = compute_dq(src3, 1);

    const int d0 = dp0 + dq0;
    const int d3 = dp3 + dq3;
    const int d  = d0  + d3;

    if (d < dbf_params->beta) {
        uint8_t use_strong_large = 0;

        if (max_l_p > 3 || max_l_q > 3) {
            int dp0L = dp0;
            int dq0L = dq0;
            int dp3L = dp3;
            int dq3L = dq3;

            if (max_l_p > 3) {
                dp0L += compute_dp(src0 - 3, 1) + 1;
                dp3L += compute_dp(src3 - 3, 1) + 1;
                dp0L >>= 1;
                dp3L >>= 1;
            }

            if (max_l_q > 3) {
                dq0L += compute_dq(src0 + 3, 1) + 1;
                dq3L += compute_dq(src3 + 3, 1) + 1;
                dq0L >>= 1;
                dq3L >>= 1;
            }

            int d0L = dp0L + dq0L;
            int d3L = dp3L + dq3L;

            int dL = d0L + d3L;

            use_strong_large = (dL < dbf_params->beta) &&
                (d0L < (dbf_params->beta + 0x10 >> 5)) &&
                (d3L < (dbf_params->beta + 0x10 >> 5)) &&
                use_strong_filter_l0(src0, 1, dbf_params->beta, dbf_params->tc, max_l_p, max_l_q) &&
                use_strong_filter_l0(src3, 1, dbf_params->beta, dbf_params->tc, max_l_p, max_l_q);
        }

        if (use_strong_large) {
            const int filter_idx = derive_filter_idx(max_l_p, max_l_q);
            df->filter_h[filter_idx](src0, stride, dbf_params->tc);
        } else {
            uint8_t sw = max_l_p > 2;

            sw = sw && (d0 < (dbf_params->beta + 0x4 >> 3))
                && (d3 < (dbf_params->beta + 0x4 >> 3))
                && use_strong_filter_l1(src0, 1, dbf_params->beta, dbf_params->tc)
                && use_strong_filter_l1(src3, 1, dbf_params->beta, dbf_params->tc);

            if (sw){
                df->filter_h[0](src0, stride, dbf_params->tc);
            } else {
                const int dp = dp0 + dp3;
                const int dq = dq0 + dq3;
                const int side_thd = (dbf_params->beta + (dbf_params->beta >> 1)) >> 3;
                uint8_t extend_p = dp < side_thd && max_l_p > 1;
                uint8_t extend_q = dq < side_thd && max_l_p > 1;
                df->filter_weak_h(src0, stride, dbf_params->tc, extend_p, extend_q);
            }
        }
    }
}

#define LF_MV_THRESHOLD 8
static inline uint8_t
mv_threshold_check(OVMV a, OVMV b)
{
    uint32_t abs_delta_x = OVABS(a.x - b.x);
    uint32_t abs_delta_y = OVABS(a.y - b.y);

    uint8_t chk = (abs_delta_x >= LF_MV_THRESHOLD) || (abs_delta_y >= LF_MV_THRESHOLD);

    return chk;
}

#define PB_POS_IN_BUF(x,y) (35 + (x) + ((y) * 34))

static uint64_t
check_dbf_enabled_p(const int16_t *dist_ref_p, const int16_t *dist_ref_q, OVMV mv_p0, OVMV mv_q0)
{
    int16_t ref0_p = dist_ref_p[mv_p0.ref_idx];

    int16_t ref0_q = dist_ref_q[mv_q0.ref_idx];
    uint8_t bs = 1;

    if (ref0_p == ref0_q) {
        bs  = mv_threshold_check(mv_q0, mv_p0);
    }

    return (uint64_t)bs;
}

static uint64_t
check_dbf_enabled(const struct InterDRVCtx *const inter_ctx,
                  OVMV mv_p0, OVMV mv_p1, OVMV mv_q0, OVMV mv_q1)
{
    const int16_t *dist_0 = inter_ctx->dist_ref_0;
    const int16_t *dist_1 = inter_ctx->dist_ref_1;

    int16_t ref0_p = dist_0[mv_p0.ref_idx];
    int16_t ref1_p = dist_1[mv_p1.ref_idx];

    int16_t ref0_q = dist_0[mv_q0.ref_idx];
    int16_t ref1_q = dist_1[mv_q1.ref_idx];

    uint8_t paired_ref_pq  = (ref0_p == ref0_q) && (ref1_p == ref1_q);
    uint8_t swapped_ref_pq = (ref0_p == ref1_q) && (ref1_p == ref0_q);

    /* FIXME ref check can be done on q */
    uint8_t coupled_l0_l1 = ref0_p == ref1_p; // Same L0 & L1
    uint8_t bs = 1;

    /* No need to check for both paired and swapped since coupled L0 L1 implies
     * paired_ref_pq == swapped_ref_pq
     */
    if ((coupled_l0_l1) && (paired_ref_pq)) {
        bs  = mv_threshold_check(mv_q0, mv_p0) || mv_threshold_check(mv_q1, mv_p1);
        bs &= mv_threshold_check(mv_q1, mv_p0) || mv_threshold_check(mv_q0, mv_p1);
    } else if (paired_ref_pq){
        bs  = mv_threshold_check(mv_q0, mv_p0);
        bs |= mv_threshold_check(mv_q1, mv_p1);
    } else if (swapped_ref_pq) {
        bs  = mv_threshold_check(mv_q1, mv_p0);
        bs |= mv_threshold_check(mv_q0, mv_p1);
    }

    return (uint64_t)bs;
}

static uint64_t
dbf_mv_set_hedges(const struct InterDRVCtx *const inter_ctx,
                  const struct IBCMVCtx *const ibc_ctx,
                  int x0_unit, int y0_unit,
                  int nb_unit_w, int nb_unit_h, uint64_t msk)
{
    const struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    const struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t unit_msk_h = (uint64_t)((uint64_t)1 << nb_unit_w) - 1llu;

    uint64_t abv0_p_msk = (mv_ctx0->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_h;
    uint64_t abv1_p_msk = (mv_ctx1->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_h;
    uint64_t ibc_p_msk = (ibc_ctx->ctu_map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_h;

    uint64_t abv0_q_msk = (mv_ctx0->map.hfield[y0_unit + 1] >> (x0_unit + 1)) & unit_msk_h;
    uint64_t abv1_q_msk = (mv_ctx1->map.hfield[y0_unit + 1] >> (x0_unit + 1)) & unit_msk_h;
    uint64_t ibc_q_msk = (ibc_ctx->ctu_map.hfield[y0_unit + 1] >> (x0_unit + 1)) & unit_msk_h;

    uint64_t bs1_map_h = ((~msk) >> x0_unit + 2) & unit_msk_h;
    bs1_map_h |= ibc_p_msk & (abv1_q_msk | abv0_q_msk);
    bs1_map_h |= ibc_q_msk & (abv1_p_msk | abv0_p_msk);

    uint64_t mv_q_b  = (abv0_q_msk &  abv1_q_msk);
    uint64_t mv_q_p0 = (abv0_q_msk & ~abv1_q_msk);
    uint64_t mv_q_p1 = (abv1_q_msk & ~abv0_q_msk);

    uint64_t mv_p_b  = (abv0_p_msk &  abv1_p_msk);
    uint64_t mv_p_p0 = (abv0_p_msk & ~abv1_p_msk);
    uint64_t mv_p_p1 = (abv1_p_msk & ~abv0_p_msk);

    uint64_t chk_b  = mv_q_b & mv_p_b;

    uint64_t chk_p0 = mv_q_p0 & (mv_p_p0 | mv_p_p1);
    uint64_t chk_p1 = mv_q_p1 & (mv_p_p0 | mv_p_p1);

    const int16_t *dist_ref0 = inter_ctx->dist_ref_0;
    const int16_t *dist_ref1 = inter_ctx->dist_ref_1;

    chk_b  &= (~(bs1_map_h | (ibc_p_msk & ibc_q_msk))) & unit_msk_h;
    chk_p0 &= (~(bs1_map_h | (ibc_p_msk & ibc_q_msk))) & unit_msk_h;
    chk_p1 &= (~(bs1_map_h | (ibc_p_msk & ibc_q_msk))) & unit_msk_h;

    uint64_t dst_map_h = (ibc_p_msk & ibc_q_msk) ^ (~(chk_p0 | chk_p1 | chk_b)) & unit_msk_h;

    if (chk_b) {
        const OVMV *mv0_p = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv1_p = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv0_q = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        const OVMV *mv1_q = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];

        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_b);
            mv0_p += nb_skipped_blk;
            mv1_p += nb_skipped_blk;
            mv0_q += nb_skipped_blk;
            mv1_q += nb_skipped_blk;
            pos_shift += nb_skipped_blk;

            uint64_t abv_th = check_dbf_enabled(inter_ctx, *mv0_p, *mv1_p, *mv0_q, *mv1_q);

            dst_map_h |= abv_th << pos_shift;

            mv0_p++;
            mv1_p++;
            mv0_q++;
            mv1_q++;
            pos_shift++;

            chk_b >>= nb_skipped_blk + 1;

        } while (chk_b);
    }

    if (chk_p0 | chk_p1) {
        const OVMV *mv0_p = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv1_p = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv0_q = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        const OVMV *mv1_q = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];

        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_p0 | chk_p1);
            mv0_p += nb_skipped_blk;
            mv1_p += nb_skipped_blk;
            mv0_q += nb_skipped_blk;
            mv1_q += nb_skipped_blk;

            chk_p0 >>= nb_skipped_blk;
            chk_p1 >>= nb_skipped_blk;
            mv_p_p0 >>= nb_skipped_blk;

            pos_shift += nb_skipped_blk;
            uint8_t is_l0_p = mv_p_p0 & 0x1;
            uint8_t is_l0_q =  chk_p0 & 0x1;
            const int16_t *dist_p = is_l0_p ? dist_ref0 : dist_ref1;
            const int16_t *dist_q = is_l0_q ? dist_ref0 : dist_ref1;
            OVMV mv_p = is_l0_p ? *mv0_p : *mv1_p;
            OVMV mv_q = is_l0_q ? *mv0_q : *mv1_q;

            uint64_t abv_th = check_dbf_enabled_p(dist_p, dist_q, mv_p, mv_q);

            dst_map_h |= abv_th << pos_shift;

            mv0_p++;
            mv1_p++;
            mv0_q++;
            mv1_q++;
            pos_shift++;

            chk_p0 >>= 1;
            chk_p1 >>= 1;
            mv_p_p0 >>= 1;

        } while (chk_p0 | chk_p1);
    }

    return ((dst_map_h | bs1_map_h) << (x0_unit + 2) & msk);
}

static uint64_t
dbf_mv_set_vedges(const struct InterDRVCtx *const inter_ctx,
                  const struct IBCMVCtx *const ibc_ctx,
                  int x0_unit, int y0_unit,
                  int nb_unit_w, int nb_unit_h, uint64_t msk)
{
    const struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    const struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t unit_msk_v = (uint64_t)((uint64_t)1 << nb_unit_h) - 1llu;

    uint64_t lft0_p_msk = (mv_ctx0->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_v;
    uint64_t lft1_p_msk = (mv_ctx1->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_v;
    uint64_t ibc_p_msk = (ibc_ctx->ctu_map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_v;

    uint64_t lft0_q_msk = (mv_ctx0->map.vfield[x0_unit + 1] >> (y0_unit + 1)) & unit_msk_v;
    uint64_t lft1_q_msk = (mv_ctx1->map.vfield[x0_unit + 1] >> (y0_unit + 1)) & unit_msk_v;
    uint64_t ibc_q_msk = (ibc_ctx->ctu_map.vfield[x0_unit + 1] >> (y0_unit + 1)) & unit_msk_v;

    uint64_t bs1_map_v = ((~msk) >> y0_unit )& unit_msk_v;
    bs1_map_v |= ibc_p_msk & (lft1_q_msk | lft0_q_msk);
    bs1_map_v |= ibc_q_msk & (lft1_p_msk | lft0_p_msk);

    uint64_t mv_q_b  = (lft0_q_msk &  lft1_q_msk);
    uint64_t mv_q_p0 = (lft0_q_msk & ~lft1_q_msk);
    uint64_t mv_q_p1 = (lft1_q_msk & ~lft0_q_msk);

    uint64_t mv_p_b  = (lft0_p_msk &  lft1_p_msk);
    uint64_t mv_p_p0 = (lft0_p_msk & ~lft1_p_msk);
    uint64_t mv_p_p1 = (lft1_p_msk & ~lft0_p_msk);

    uint64_t chk_b  = mv_q_b & mv_p_b;

    uint64_t chk_p0 = mv_q_p0 & (mv_p_p0 | mv_p_p1);
    uint64_t chk_p1 = mv_q_p1 & (mv_p_p0 | mv_p_p1);

    const int16_t *dist_ref0 = inter_ctx->dist_ref_0;
    const int16_t *dist_ref1 = inter_ctx->dist_ref_1;

    chk_b  &= (~(bs1_map_v | (ibc_p_msk & ibc_q_msk))) & unit_msk_v;

    chk_p0 &= (~(bs1_map_v | (ibc_p_msk & ibc_q_msk))) & unit_msk_v;
    chk_p1 &= (~(bs1_map_v | (ibc_p_msk & ibc_q_msk))) & unit_msk_v;

    uint64_t dst_map_v = (ibc_p_msk & ibc_q_msk) ^ (~(chk_p0 | chk_p1 | chk_b)) & unit_msk_v;

    if (chk_b) {
        const OVMV *mv0_p = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv1_p = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv0_q = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        const OVMV *mv1_q = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];

        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_b);
            mv0_p += 34 * nb_skipped_blk;
            mv1_p += 34 * nb_skipped_blk;
            mv0_q += 34 * nb_skipped_blk;
            mv1_q += 34 * nb_skipped_blk;
            pos_shift += nb_skipped_blk;

            uint64_t lft_th = check_dbf_enabled(inter_ctx, *mv0_p, *mv1_p, *mv0_q, *mv1_q);

            dst_map_v |= lft_th << pos_shift;

            mv0_p += 34;
            mv1_p += 34;
            mv0_q += 34;
            mv1_q += 34;
            pos_shift++;

            chk_b >>= nb_skipped_blk + 1;

        } while (chk_b);
    }

    if (chk_p0 | chk_p1) {
        const OVMV *mv0_p = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv1_p = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv0_q = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        const OVMV *mv1_q = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];

        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_p0 | chk_p1);
            mv0_p += 34 * nb_skipped_blk;
            mv1_p += 34 * nb_skipped_blk;
            mv0_q += 34 * nb_skipped_blk;
            mv1_q += 34 * nb_skipped_blk;

            chk_p0 >>= nb_skipped_blk;
            chk_p1 >>= nb_skipped_blk;
            mv_p_p0 >>= nb_skipped_blk;

            pos_shift += nb_skipped_blk;
            uint8_t is_l0_p = mv_p_p0 & 0x1;
            uint8_t is_l0_q =  chk_p0 & 0x1;
            const int16_t *dist_p = is_l0_p ? dist_ref0 : dist_ref1;
            const int16_t *dist_q = is_l0_q ? dist_ref0 : dist_ref1;
            OVMV mv_p = is_l0_p ? *mv0_p : *mv1_p;
            OVMV mv_q = is_l0_q ? *mv0_q : *mv1_q;

            uint64_t lft_th = check_dbf_enabled_p(dist_p, dist_q, mv_p, mv_q);

            dst_map_v |= lft_th << pos_shift;

            mv0_p += 34;
            mv1_p += 34;
            mv0_q += 34;
            mv1_q += 34;
            pos_shift++;

            chk_p0 >>= 1;
            chk_p1 >>= 1;
            mv_p_p0 >>= 1;

        } while (chk_p0 | chk_p1);
    }

    return ((dst_map_v  | bs1_map_v) << y0_unit) & msk;
}

static void
dbf_ctu_preproc_h(const struct InterDRVCtx *const inter_ctx, struct DBFInfo *const dbf_info,
                  uint8_t nb_unit_h, uint8_t nb_unit_w)
{
    int i;
    const uint64_t *cu_edg_map = &dbf_info->cu_edge.hor[0];
    const uint64_t *sb_edg_map = &dbf_info->aff_edg_hor[8];

    for (i = 0; i < nb_unit_h; ++i) {
        uint64_t edg_msk = cu_edg_map[i] | sb_edg_map[i];
        uint64_t bs_map  = dbf_info->bs2_map.hor[i];
                 bs_map |= dbf_info->bs1_map.hor[i];
                 bs_map &= edg_msk;
        if (edg_msk ^ bs_map) {
            /* If some edges are not set yet check for MV condition */
            uint64_t no_filter = edg_msk ^ bs_map;

            if (no_filter) {
                uint64_t tmp = bs_map;
                tmp = dbf_mv_set_hedges(inter_ctx, dbf_info->ibc_ctx, 0, i, nb_unit_w, nb_unit_h, no_filter);
                dbf_info->bs1_map.hor[i] |= tmp;
            }
        }
    }
}

static void
dbf_ctu_preproc_v(const struct InterDRVCtx *const inter_ctx, struct DBFInfo *const dbf_info,
                  uint8_t nb_unit_h, uint8_t nb_unit_w)
{
    const uint64_t vedge_mask = ((uint64_t)1 << nb_unit_h) - (uint64_t)1;

    const uint64_t *cu_edg_map = &dbf_info->cu_edge.ver[0];
    const uint64_t *sb_edg_map = &dbf_info->aff_edg_ver[8];

    int i;

    for (i = 0; i < nb_unit_w; ++i) {
        uint64_t edg_msk = cu_edg_map[i] | sb_edg_map[i];
        uint64_t bs_map  = dbf_info->bs2_map.ver[i];
                 bs_map |= dbf_info->bs1_map.ver[i];
                 bs_map &= edg_msk;
        if (edg_msk ^ bs_map) {
            /* If some edges are not set yet check for MV condition */
            uint64_t no_filter = edg_msk ^ bs_map;

            if (no_filter) {
                uint64_t tmp = bs_map;
                tmp = dbf_mv_set_vedges(inter_ctx, dbf_info->ibc_ctx, i, 0, nb_unit_w, nb_unit_h, no_filter);
                dbf_info->bs1_map.ver[i] |= tmp;
            }
        }
    }
}

struct EdgeCtx
{
    uint64_t large_p_map;
    uint64_t large_q_map;
    uint64_t small_map;
    uint64_t aff_edg_1;
};

struct DBFLength
{
    uint8_t lgth_p;
    uint8_t lgth_q;
};

static struct DBFLength
derive_filter_length(const struct EdgeCtx *edg_ctx, uint64_t affine_p, uint64_t affine_q, uint64_t pos_msk)
{
    struct DBFLength lgth_info;
    uint8_t max_l_p, max_l_q;

    if (edg_ctx->small_map & pos_msk) {
        max_l_p = 1;
        max_l_q = 1;
    } else if (edg_ctx->aff_edg_1 & pos_msk) {
        max_l_p = 2;
        max_l_q = 2;
    } else {
        max_l_p = 3;
        max_l_q = 3;

        if (edg_ctx->large_p_map & pos_msk) {
            max_l_p = 7;
            if (affine_p & pos_msk) {
                max_l_p = 5;
            }
        }

        if (edg_ctx->large_q_map & pos_msk) {
            max_l_q = 7;
            if (affine_q & pos_msk) {
                max_l_q = 5;
            }
        }
    }

    lgth_info.lgth_p = max_l_p;
    lgth_info.lgth_q = max_l_q;

    return lgth_info;
}

static void
set_edge_context(struct EdgeCtx *edg_ctx, const uint64_t *edg_map, const uint64_t *sb_edg_map, uint8_t i)
{
    edg_ctx->large_q_map = i % 4 ? 0 : derive_size_3_map(&edg_map[i + 1]);

    edg_ctx->small_map = edg_map[i - 1] | edg_map[i + 1] | sb_edg_map[i - 1] | sb_edg_map[i + 1];
    edg_ctx->aff_edg_1 = sb_edg_map[i] & (edg_map[i - 2] | edg_map[i + 2])
                                       & ~edg_map[i];

    edg_ctx->large_p_map  &= ~(sb_edg_map[i] & (~edg_map[i]));
    edg_ctx->large_q_map  &= ~(sb_edg_map[i] & (~edg_map[i]));
}

static int8_t
derive_ladf_offset(const struct LADFParams *const ladf, int16_t luma_lvl)
{
     int8_t ladf_qp_offset = ladf->qp_offset[0];
     for (int i = 1; i < ladf->nb_intervals; ++i) {
         int16_t ladf_threshold  = ladf->threshold[i];
         if (luma_lvl > ladf_threshold) {
             ladf_qp_offset = ladf->qp_offset[i];
         } else {
           break;
         }
     }

     return ladf_qp_offset;
  
}

static void
vvc_dbf_ctu_hor(const struct DFFunctions * df, OVSample *src, int stride, const struct DBFInfo *const dbf_info,
                uint8_t nb_unit_h, int is_last_h, uint8_t nb_unit_w, uint8_t ctu_lft)
{
    const int blk_stride = stride << 2;
    const uint64_t vedge_mask = ((uint64_t)1 << nb_unit_h) - 1;

    const uint64_t *edg_map = &dbf_info->ctb_bound_ver[8];
    const uint64_t *sb_edg_map = &dbf_info->aff_edg_ver[8];
    const uint8_t skip_first = !ctu_lft;

    int i;

    src += skip_first << 2;

    for (i = skip_first; i < nb_unit_w; ++i) {
        OVSample* src_tmp = src;

        uint64_t edg_msk = edg_map[i] | sb_edg_map[i];
        uint64_t bs1_map  = dbf_info->bs1_map.ver[i];
        uint64_t bs2_map  = dbf_info->bs2_map.ver[i];

        edg_msk &= vedge_mask;
        edg_msk &= bs2_map | bs1_map;

        if (edg_msk) {
            uint64_t pos_msk = 1;
            uint8_t once = 1;
            struct EdgeCtx edg_ctx={{0}};
            const int8_t *qp_col = &dbf_info->qp_map_y.hor[36 + i];

            do {
                int8_t qp, bs;
                /* Skip non filtered edges */
                uint8_t nb_skipped_blk = ov_ctz64(edg_msk);
                int8_t ladf_qp_offset = 0;

                pos_msk <<= nb_skipped_blk;

                qp_col       += nb_skipped_blk * 34;
                src_tmp      += nb_skipped_blk * blk_stride;

                if (dbf_info->ladf_prms.nb_intervals) {
                    uint16_t luma_lvl = (src_tmp[0] + src_tmp[3*stride] + src_tmp[-1] + src_tmp[3*stride - 1]) >> 2;
                    ladf_qp_offset = derive_ladf_offset(&dbf_info->ladf_prms, luma_lvl);
                }

                bs = 1 + !!(bs2_map & pos_msk);

                qp = (qp_col[-1] + qp_col[0] + 1) >> 1;

                const struct DBFParams dbf_params = compute_dbf_limits(dbf_info, qp, bs, ladf_qp_offset, 0);
                if (dbf_params.tc || dbf_params.beta) {
                    const uint64_t affine_p = dbf_info->affine_map.ver[i    ];
                    const uint64_t affine_q = dbf_info->affine_map.ver[i + 1];

                    if (once) {
                        edg_ctx.large_p_map = i % 4 ? 0 : derive_size_3_map(&edg_map[i - 7]);
                        set_edge_context(&edg_ctx, edg_map, sb_edg_map, i);
                        once = 0;
                    }

                    struct DBFLength max_lgth_info = derive_filter_length(&edg_ctx, affine_p, affine_q,
                                                                          pos_msk);

                    filter_vertical_edge(df, &dbf_params, src_tmp, stride,
                                         max_lgth_info.lgth_p, max_lgth_info.lgth_q);
                }

                edg_msk >>= nb_skipped_blk + 1;
                pos_msk <<= 1;

                src_tmp += blk_stride;
                qp_col  += 34;

            } while (edg_msk);
        }
        src += 1 << 2;
    }
}

static void
filter_horizontal_edge(const struct DFFunctions *df, const struct DBFParams *const dbf_params,
                       OVSample *src, ptrdiff_t stride,
                       uint8_t max_l_p, uint8_t max_l_q)
{
    OVSample *src0 = (OVSample *)src;
    OVSample *src3 = (OVSample *)src + 3;

    const int dp0 = compute_dp(src0, stride);
    const int dq0 = compute_dq(src0, stride);
    const int dp3 = compute_dp(src3, stride);
    const int dq3 = compute_dq(src3, stride);

    const int d0 = dp0 + dq0;
    const int d3 = dp3 + dq3;
    const int d  = d0 + d3;

    if (d < dbf_params->beta) {
        uint8_t use_strong_large = 0;

        if (max_l_p > 3 || max_l_q > 3) {
            int dp0L = dp0;
            int dq0L = dq0;
            int dp3L = dp3;
            int dq3L = dq3;

            if (max_l_p > 3) {
                dp0L = (dp0L + compute_dp(src0 - 3 * stride, stride) + 1) >> 1;
                dp3L = (dp3L + compute_dp(src3 - 3 * stride, stride) + 1) >> 1;
            }

            if (max_l_q > 3) {
                dq0L = (dq0L + compute_dq(src0 + 3 * stride, stride) + 1) >> 1;
                dq3L = (dq3L + compute_dq(src3 + 3 * stride, stride) + 1) >> 1;
            }

            int d0L = dp0L + dq0L;
            int d3L = dp3L + dq3L;

            int dL = d0L + d3L;
            use_strong_large = (dL < dbf_params->beta) &&
                (d0L < (dbf_params->beta + 0x10 >> 5)) &&
                (d3L < (dbf_params->beta + 0x10 >> 5)) &&
                use_strong_filter_l0(src0, stride, dbf_params->beta, dbf_params->tc, max_l_p, max_l_q) &&
                use_strong_filter_l0(src3, stride, dbf_params->beta, dbf_params->tc, max_l_p, max_l_q);
        }

        if (use_strong_large) {
            const int filter_idx = derive_filter_idx(max_l_p, max_l_q);
            df->filter_v[filter_idx](src0, stride, dbf_params->tc);
        } else {
            uint8_t sw = max_l_p > 2;

            sw = sw && (d0 < (dbf_params->beta + 0x4 >> 3))
                && (d3 < (dbf_params->beta + 0x4 >> 3))
                && use_strong_filter_l1(src0, stride, dbf_params->beta, dbf_params->tc)
                && use_strong_filter_l1(src3, stride, dbf_params->beta, dbf_params->tc);

            if (sw){
                df->filter_v[0](src0, stride, dbf_params->tc);
            } else {
                const int dp = dp0 + dp3;
                const int dq = dq0 + dq3;
                const int side_thd = (dbf_params->beta + (dbf_params->beta >> 1)) >> 3;
                uint8_t extend_p = dp < side_thd && max_l_p > 1;
                uint8_t extend_q = dq < side_thd && max_l_p > 1;
                df->filter_weak_v(src0, stride, dbf_params->tc, extend_p, extend_q);
            }
        }
    }
}

static void
vvc_dbf_ctu_ver(const struct DFFunctions *df, OVSample *src, int stride, const struct DBFInfo *const dbf_info,
                uint8_t nb_unit_w, int is_last_w, uint8_t nb_unit_h, uint8_t ctu_abv)
{
    const int blk_stride = 1 << 2;
    const uint64_t hedge_mask = ((uint64_t)1 << (nb_unit_w + (!!is_last_w << 1))) - 1llu;
    int i;

    const uint64_t *edg_map = &dbf_info->ctb_bound_hor[8];
    const uint64_t *sb_edg_map = &dbf_info->aff_edg_hor[8];
    uint8_t skip_first = !ctu_abv;

    /* Filtering vertical edges on the whole would overlap with next CTU first
     * vertical edge.
     * Since max filter length is 7 we stop the horizontal process 2 units before
     * current CTU end and process the two previous units in order to finish filtering
     * the horizontal edges of previous CTU.
     */
    src -= blk_stride << 1;
    src += (skip_first * stride) << 2;

    for (i = skip_first; i < nb_unit_h; ++i) {
        OVSample *src_tmp = src;

        uint64_t edg_msk = edg_map[i] | sb_edg_map[i];
        uint64_t bs2_map = dbf_info->bs2_map.hor[i];
        uint64_t bs1_map = dbf_info->bs1_map.hor[i];

        edg_msk &= hedge_mask;
        edg_msk &= bs2_map | bs1_map;

        if (edg_msk) {
            uint64_t pos_msk = 1;
            uint8_t once = 1;
            struct EdgeCtx edg_ctx ={{0}};
            const int8_t *qp_row = &dbf_info->qp_map_y.hor[34 * i];

            do {
                int8_t qp, bs;
                uint8_t nb_skipped_blk = ov_ctz64(edg_msk);
                int8_t ladf_qp_offset = 0;

                /* Skip non filtered edges */
                pos_msk <<= nb_skipped_blk;

                qp_row       += nb_skipped_blk;
                src_tmp      += nb_skipped_blk * blk_stride;

                if (dbf_info->ladf_prms.nb_intervals) {
                    uint16_t luma_lvl = (src_tmp[0] + src_tmp[3] + src_tmp[-stride] + src_tmp[-stride + 3]) >> 2;
                    ladf_qp_offset = derive_ladf_offset(&dbf_info->ladf_prms, luma_lvl);
                }

                bs = 1 + !!(bs2_map & pos_msk);

                qp = (qp_row[0] + qp_row[34] + 1) >> 1;

                const struct DBFParams dbf_params = compute_dbf_limits(dbf_info, qp, bs, ladf_qp_offset, 0);
                if (dbf_params.tc || dbf_params.beta) {
                    const uint64_t affine_p = dbf_info->affine_map.hor[i    ];
                    const uint64_t affine_q = dbf_info->affine_map.hor[i + 1];

                    if (once) {
                        edg_ctx.large_p_map = i % 4 || i < 7 ? 0 : derive_size_3_map(&edg_map[i - 7]);
                        set_edge_context(&edg_ctx, edg_map, sb_edg_map, i);
                        once = 0;
                    }

                    struct DBFLength max_lgth_info = derive_filter_length(&edg_ctx, affine_p, affine_q,
                                                                          pos_msk);

                    filter_horizontal_edge(df, &dbf_params, src_tmp, stride,
                                           max_lgth_info.lgth_p, max_lgth_info.lgth_q);

                }

                edg_msk  >>= nb_skipped_blk + 1;
                pos_msk <<= 1;

                src_tmp += blk_stride;
                qp_row++;

            } while (edg_msk);
        }
        src += stride << 2;
    }
}

static void
rcn_dbf_ctu(const struct OVRCNCtx  *const rcn_ctx, struct DBFInfo *const dbf_info,
            uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y)
{
    const struct OVBuffInfo *const fbuff = &rcn_ctx->frame_buff;
    const struct DFFunctions *df = &rcn_ctx->ctudec->rcn_funcs.df;
    uint8_t nb_unit = (1 << log2_ctu_s) >> 2;
    /* FIXME give as argument */
    uint8_t ctu_lft = rcn_ctx->ctudec->ctu_ngh_flags & CTU_LFT_FLG;
    uint8_t ctu_abv = rcn_ctx->ctudec->ctu_ngh_flags & CTU_UP_FLG;

    if (rcn_ctx->ctudec->tmp_slice_type !=2){
        dbf_ctu_preproc_v(&rcn_ctx->ctudec->drv_ctx.inter_ctx, dbf_info, nb_unit, nb_unit);
        dbf_ctu_preproc_h(&rcn_ctx->ctudec->drv_ctx.inter_ctx, dbf_info, nb_unit, nb_unit);
    }

    if (!dbf_info->disable_h)
    vvc_dbf_ctu_hor(df, fbuff->y, fbuff->stride, dbf_info, nb_unit, !!last_y, nb_unit, ctu_lft);
    if (!dbf_info->disable_v)
    vvc_dbf_ctu_ver(df, fbuff->y, fbuff->stride, dbf_info, nb_unit, !!last_x, nb_unit, ctu_abv);

    if (!dbf_info->disable_h)
    vvc_dbf_chroma_hor(df, fbuff->cb, fbuff->cr, fbuff->stride_c, dbf_info,
                       nb_unit, !!last_y, nb_unit, ctu_lft);

    if (!dbf_info->disable_v)
    vvc_dbf_chroma_ver(df, fbuff->cb, fbuff->cr, fbuff->stride_c, dbf_info,
                       nb_unit, !!last_x, nb_unit, !!last_y, ctu_abv);

}

static void
rcn_dbf_truncated_ctu(const struct OVRCNCtx  *const rcn_ctx, struct DBFInfo *const dbf_info,
                      uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y, uint8_t ctu_w, uint8_t ctu_h)
{
    const struct OVBuffInfo *const fbuff = &rcn_ctx->frame_buff;
    const struct DFFunctions *df = &rcn_ctx->ctudec->rcn_funcs.df;

    uint8_t nb_unit_w = (ctu_w) >> 2;
    uint8_t nb_unit_h = (ctu_h) >> 2;
    /* FIXME give as argument */
    uint8_t ctu_lft = rcn_ctx->ctudec->ctu_ngh_flags & CTU_LFT_FLG;
    uint8_t ctu_abv = rcn_ctx->ctudec->ctu_ngh_flags & CTU_UP_FLG;

    if (rcn_ctx->ctudec->tmp_slice_type !=2){
        dbf_ctu_preproc_v(&rcn_ctx->ctudec->drv_ctx.inter_ctx, dbf_info, nb_unit_h, nb_unit_w);
        dbf_ctu_preproc_h(&rcn_ctx->ctudec->drv_ctx.inter_ctx, dbf_info, nb_unit_h, nb_unit_w);
    }

    if (!dbf_info->disable_h)
    vvc_dbf_ctu_hor(df, fbuff->y, fbuff->stride, dbf_info, nb_unit_h, !!last_y, nb_unit_w, ctu_lft);
    if (!dbf_info->disable_v)
    vvc_dbf_ctu_ver(df, fbuff->y, fbuff->stride, dbf_info, nb_unit_w, !!last_x, nb_unit_h, ctu_abv);

    if (!dbf_info->disable_h)
    vvc_dbf_chroma_hor(df, fbuff->cb, fbuff->cr, fbuff->stride_c, dbf_info,
                       nb_unit_h, !!last_y, nb_unit_w, ctu_lft);

    if (!dbf_info->disable_v)
    vvc_dbf_chroma_ver(df, fbuff->cb, fbuff->cr, fbuff->stride_c, dbf_info,
                       nb_unit_w, !!last_x, nb_unit_h, !!last_y, ctu_abv);

}

void
BD_DECL(rcn_init_df_functions)(struct RCNFunctions *const rcn_funcs)
{
  rcn_funcs->df.filter_h[0] = &filter_luma_strong_small_h;
  rcn_funcs->df.filter_h[1] = &filter_h_3_5;
  rcn_funcs->df.filter_h[2] = &filter_h_3_7;
  rcn_funcs->df.filter_h[3] = NULL;
  rcn_funcs->df.filter_h[4] = &filter_h_5_3;
  rcn_funcs->df.filter_h[5] = &filter_h_5_5;
  rcn_funcs->df.filter_h[6] = &filter_h_5_7;
  rcn_funcs->df.filter_h[7] = NULL;
  rcn_funcs->df.filter_h[8] = &filter_h_7_3;
  rcn_funcs->df.filter_h[9] = &filter_h_7_5;
  rcn_funcs->df.filter_h[10]= &filter_h_7_7;

  rcn_funcs->df.filter_v[0] = &filter_luma_strong_small_v;
  rcn_funcs->df.filter_v[1] = &filter_v_3_5;
  rcn_funcs->df.filter_v[2] = &filter_v_3_7;
  rcn_funcs->df.filter_v[3] = NULL;
  rcn_funcs->df.filter_v[4] = &filter_v_5_3;
  rcn_funcs->df.filter_v[5] = &filter_v_5_5;
  rcn_funcs->df.filter_v[6] = &filter_v_5_7;
  rcn_funcs->df.filter_v[7] = NULL;
  rcn_funcs->df.filter_v[8] = &filter_v_7_3;
  rcn_funcs->df.filter_v[9] = &filter_v_7_5;
  rcn_funcs->df.filter_v[10]= &filter_v_7_7;
  rcn_funcs->df.filter_weak_h = filter_luma_weak_h;
  rcn_funcs->df.filter_weak_v = filter_luma_weak_v;
  rcn_funcs->df.filter_weak_h_c = filter_chroma_weak_h;
  rcn_funcs->df.filter_weak_v_c = filter_chroma_weak_v;
  rcn_funcs->df.filter_strong_h_c = filter_chroma_strong_c_h;
  rcn_funcs->df.filter_strong_v_c = filter_chroma_strong_c_v;

  rcn_funcs->df.rcn_dbf_ctu = &rcn_dbf_ctu;
  rcn_funcs->df.rcn_dbf_truncated_ctu = &rcn_dbf_truncated_ctu;
}
