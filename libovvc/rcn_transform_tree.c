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
#include <string.h>

#include "ovutils.h"
#include "ovmem.h"
#include "ctudec.h"
#include "dbf_utils.h"
#include "drv_utils.h"
#include "drv.h"
#include "vcl.h"
#include "rcn_dequant.h"
#include "bitdepth.h"

#define TR_SHIFT_V (6 + 1)
#define TR_SHIFT_H ((6 + 15 - 1) - BITDEPTH)
#define LOG2_MIN_CU_S 2

struct TBInfo {
   uint16_t last_pos;
   uint64_t sig_sb_map;
};

struct TUInfo {
   uint8_t is_sbt;
   uint8_t cbf_mask;
   uint16_t pos_offset;
   uint8_t tr_skip_mask;
   uint8_t cu_mts_flag;
   uint8_t cu_mts_idx;
   uint8_t lfnst_flag;
   uint8_t lfnst_idx;
   struct TBInfo tb_info[3];
};

struct ISPTUInfo {
   uint8_t cbf_mask;
   uint8_t tr_skip_mask;
   uint8_t cu_mts_flag;
   uint8_t cu_mts_idx;
   uint8_t lfnst_flag;
   uint8_t lfnst_idx;
   struct TBInfo tb_info[4];
};

static int
derive_nb_rows(uint64_t sig_sb_map)
{
    uint8_t col = sig_sb_map | 1;

    col |= sig_sb_map >> 8;
    col |= sig_sb_map >> 16;
    col |= sig_sb_map >> 24;
    col |= sig_sb_map >> 32;
    col |= sig_sb_map >> 40;
    col |= sig_sb_map >> 48;
    col |= sig_sb_map >> 56;

    return (32 - ov_clz(col)) << 2;
}

static int
derive_nb_cols(uint64_t sig_sb_map)
{
    uint8_t num_z = ov_clz64(sig_sb_map | 1);
    uint8_t z_div8 = (num_z >> 3);

    return (8 - z_div8) << 2;
}

#define MAX_LOG2_TR_RANGE 15
static const uint16_t tb_lut_offset[64] =
{      0,    0,    2,    6,   14,   30,   62,   94,
      94,   96,  100,  108,  124,  156,  220,  284,
     284,  288,  296,  312,  344,  408,  536,  664,
     664,  672,  688,  720,  784,  912, 1168, 1424,
    1424, 1440, 1472, 1536, 1664, 1920, 2432, 2944,
    2944, 2976, 3040, 3168, 3424, 3936, 4960, 5984,
    5984, 6016, 6080, 6208, 6464, 6976, 8000, 9024,
    9024, 9024, 9024, 9024, 9024, 9024, 9024, 9024
};

static void
scale_dequant_tb(const int16_t *const luts, int16_t *dst, const int16_t *src, int scale, int shift, uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t nb_rows, uint8_t nb_cols)
{
    uint8_t lut_id = (log2_tb_w << 3) | log2_tb_h;
    const int16_t *sl = &luts[tb_lut_offset[lut_id]];

    int add = (1 << shift) >> 1;
    uint8_t stride = 1 << (OVMIN(5, log2_tb_w));
    scale /= 16;
    for (int i = 0; i < nb_rows ; i++) {
        for (int j = 0; j < nb_cols ; j++) {
            dst[j] = ov_clip_intp2((int32_t)(src[j] * sl[j + (i << log2_tb_w)] * scale + add) >> shift,
                                   MAX_LOG2_TR_RANGE + 1);
        }
        src += stride;
        dst += stride;
    }
}

static void
dequant_tb(int16_t *dst, const int16_t *src, int scale, int shift, uint8_t log2_tb_w, uint8_t nb_rows, uint8_t nb_cols)
{
    int add = (1 << shift) >> 1;
    uint8_t stride = 1 << (OVMIN(5, log2_tb_w));
    for (int i = 0; i < nb_rows ; i++) {
        for (int j = 0; j < nb_cols ; j++) {
            dst[j] = ov_clip_intp2((int32_t)(src[j] * scale + add) >> shift,
                                   MAX_LOG2_TR_RANGE + 1);
        }
        src += stride;
        dst += stride;
    }
}

static void
reorder_tb_4x4(int16_t *dst, const int16_t *src, int scale, int shift,
               uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map)
{
    int nb_rows = 1 << log2_tb_h; //derive_nb_cols(sig_sb_map);
    int nb_cols = 1 << log2_tb_w; //derive_nb_rows(sig_sb_map);
    uint8_t src_stride = 1 << (OVMIN(5, log2_tb_w));
    uint8_t dst_stride = 1 << (OVMIN(5, log2_tb_w));

    for (int i = 0; i < nb_rows/4 ; i++) {
        uint8_t sig_sb_row = sig_sb_map >> (i << 3);
        for (int j = 0; j < nb_cols/4 ; j++) {
            int16_t *_dst = dst + (j << 2);
            const int16_t *_src = src + (j << 4);
            if (sig_sb_row & 0x1) {
                _dst[0] = _src[ 0];
                _dst[1] = _src[ 1];
                _dst[2] = _src[ 2];
                _dst[3] = _src[ 3];

                _dst += dst_stride;

                _dst[0] = _src[ 4];
                _dst[1] = _src[ 5];
                _dst[2] = _src[ 6];
                _dst[3] = _src[ 7];

                _dst += dst_stride;

                _dst[0] = _src[ 8];
                _dst[1] = _src[ 9];
                _dst[2] = _src[10];
                _dst[3] = _src[11];

                _dst += dst_stride;

                _dst[0] = _src[12];
                _dst[1] = _src[13];
                _dst[2] = _src[14];
                _dst[3] = _src[15];
            } else {
                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;

                _dst += dst_stride;

                _dst[0] = 0;
                _dst[1] = 0;
                _dst[2] = 0;
                _dst[3] = 0;
            }
            sig_sb_row >>= 1;
        }
        src += src_stride << 2;
        dst += dst_stride << 2;
    }
}

static void
dequant_4x4_ts(OVCTUDec *const ctudec, int16_t *dst, const int16_t *src, uint64_t sig_sb_map, uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t qp)
{
    struct IQScale deq_prms;
    int nb_cols = 1 << log2_tb_w; //derive_nb_cols(sig_sb_map);
    int nb_rows = 1 << log2_tb_h; //derive_nb_rows(sig_sb_map);

    deq_prms = ctudec->rcn_funcs.tmp.derive_dequant_ts(qp, log2_tb_w, log2_tb_h);

    ctudec->rcn_funcs.tmp.dequant_tb_4x4(dst, src, deq_prms.scale, deq_prms.shift, log2_tb_w, log2_tb_h, sig_sb_map);
}

/*TODO at decode init */
static inline struct IQScale
derive_dequant(const OVCTUDec *const ctudec, uint8_t qp, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    return ctudec->rcn_funcs.tmp.derive_dequant(qp, log2_tb_w, log2_tb_h);
}

static void
dequant_4x4_sb(OVCTUDec *const ctudec, int16_t *dst, const int16_t *src, uint64_t sig_sb_map, uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t qp, uint8_t is_intra, uint8_t comp_idx, uint8_t is_lfnst)
{
    int nb_cols = derive_nb_cols(sig_sb_map);
    int nb_rows = derive_nb_rows(sig_sb_map);

    struct IQScale deq_prms = derive_dequant(ctudec, qp, log2_tb_w, log2_tb_h);

    if (!ctudec->scaling_list_enabled || is_lfnst) {
        ctudec->rcn_funcs.tmp.dequant_tb_4x4(dst, src, deq_prms.scale, deq_prms.shift, log2_tb_w, log2_tb_h, sig_sb_map);
    } else {
        const uint16_t *lut;
        if (is_intra) {
          lut = ctudec->tb_scaling_luts.intra_luts;
          if (comp_idx == 1) {
              lut = ctudec->tb_scaling_luts.intra_luts_cb;
          }
          if (comp_idx == 2) {
              lut = ctudec->tb_scaling_luts.intra_luts_cr;
          }
        } else {
          lut = ctudec->tb_scaling_luts.inter_luts;
          if (comp_idx == 1) {
              lut = ctudec->tb_scaling_luts.inter_luts_cb;
          }
          if (comp_idx == 2) {
              lut = ctudec->tb_scaling_luts.inter_luts_cr;
          }
        }
        ctudec->rcn_funcs.tmp.scale_dequant_tb_4x4(lut, dst, src, deq_prms.scale, deq_prms.shift, log2_tb_w, log2_tb_h, sig_sb_map);
    }
}

static void
rcn_residual(OVCTUDec *const ctudec,
             int16_t *const dst, int16_t *src,
             uint8_t x0, uint8_t y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t cu_mts_flag, uint8_t cu_mts_idx,
             uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx, uint64_t sig_sb_map)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_funcs.tr;
    fill_bs_map(&ctudec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);

    DECLARE_ALIGNED(32, int16_t, tmp)[64*64];
    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[32*32];

    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    int qp = ctudec->dequant_luma.qp;
    dequant_4x4_sb(ctudec, dequant_coeffs, src, sig_sb_map, log2_tb_w, log2_tb_h, qp, 1, 0, lfnst_flag);


    if (!is_mip && !cu_mts_flag && ctudec->mts_implicit && (log2_tb_w <= 4 || log2_tb_h <= 4) && !lfnst_flag) {

        enum DCTType tr_h_idx = log2_tb_w <= 4 ? DST_VII : DCT_II;
        enum DCTType tr_v_idx = log2_tb_h <= 4 ? DST_VII : DCT_II;

        int nb_col = derive_nb_cols(sig_sb_map);
        int nb_row = derive_nb_rows(sig_sb_map);

        int cb_w = OVMIN(32, tb_w);

        /* FIXME use coefficient zeroing in MTS */
        memset(&tmp[nb_row << log2_tb_h], 0, sizeof(int16_t) * ((1 << (log2_tb_w + log2_tb_h)) - (nb_row << log2_tb_h)));

        TRFunc->func[tr_v_idx][log2_tb_h](dequant_coeffs, tmp, cb_w, nb_row, nb_col, TR_SHIFT_V);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);

    } else if (!cu_mts_flag) {

        int nb_row = derive_nb_rows(sig_sb_map);
        int nb_col = derive_nb_cols(sig_sb_map);

        if (lfnst_flag) {
            int16_t lfnst_sb[16];

            int tmp_shift = OVMIN(5,log2_tb_w);

            int8_t intra_mode = is_mip ? OVINTRA_PLANAR : ctudec->intra_mode;

            int8_t lfnst_intra_mode = drv_lfnst_mode_l(log2_tb_w, log2_tb_h, intra_mode);
            uint8_t is_8x8 = log2_tb_w >= 3 && log2_tb_h >= 3;

            memcpy(lfnst_sb     , &dequant_coeffs[0], sizeof(int16_t) * 4);
            memcpy(lfnst_sb +  4, &dequant_coeffs[1 << tmp_shift], sizeof(int16_t) * 4);
            memcpy(lfnst_sb +  8, &dequant_coeffs[2 << tmp_shift], sizeof(int16_t) * 4);
            memcpy(lfnst_sb + 12, &dequant_coeffs[3 << tmp_shift], sizeof(int16_t) * 4);

            process_lfnst_luma(ctudec, dequant_coeffs, lfnst_sb, OVMIN(5, log2_tb_w), OVMIN(log2_tb_h, 5),
                               lfnst_idx, lfnst_intra_mode);

            nb_row = 4 << is_8x8;
            nb_col = 4 << is_8x8;

            is_dc = 0;
        }

        if (is_dc) {
            TRFunc->dc(dst, log2_tb_w, log2_tb_h, dequant_coeffs[0]);
        } else {
            int cb_w = OVMIN(32, tb_w);

            /* FIXME use coefficient zeroing in MTS */
            memset(&tmp[nb_row << log2_tb_h], 0, sizeof(int16_t) * ((1 << (log2_tb_w + log2_tb_h)) - (nb_row << log2_tb_h)));

            TRFunc->func[DCT_II][log2_tb_h](dequant_coeffs, tmp, cb_w, nb_row, nb_col, TR_SHIFT_V);
            TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);
        }
    } else {
        enum DCTType tr_h_idx = cu_mts_idx  & 1;
        enum DCTType tr_v_idx = cu_mts_idx >> 1;

        int nb_row = derive_nb_rows(sig_sb_map);
        int nb_col = derive_nb_cols(sig_sb_map);

        int cb_w = OVMIN(32, tb_w);

        /* FIXME use coefficient zeroing in MTS */
        memset(&tmp[nb_row << log2_tb_h], 0, sizeof(int16_t) * ((1 << (log2_tb_w + log2_tb_h)) - (nb_row << log2_tb_h)));

        TRFunc->func[tr_v_idx][log2_tb_h](dequant_coeffs, tmp, cb_w, nb_row, nb_col, TR_SHIFT_V);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);
    }
}

static void
dequant_non_4x4_sb(const OVCTUDec *const ctudec, int16_t *dst, int16_t *src, uint64_t sig_sb_map,
                   uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t qp)
{
    /* Note in case of small blocks TBs are supposed to be small so processing the
     * whole block should not be an issue regarding performance
     */
    uint8_t tb_w = 1 << log2_tb_w;
    uint8_t tb_h = 1 << log2_tb_h;

    struct IQScale deq_prms = derive_dequant(ctudec, qp, log2_tb_w, log2_tb_h);

    if (!ctudec->scaling_list_enabled) {
        dequant_tb(dst, src, deq_prms.scale, deq_prms.shift, log2_tb_w, tb_h, tb_w);
    } else {
        const uint16_t *lut;
        lut = ctudec->tb_scaling_luts.intra_luts;
        scale_dequant_tb(lut, dst, src, deq_prms.scale, deq_prms.shift, log2_tb_w, log2_tb_h, tb_h, tb_w);
    }
}

static void
rcn_residual_c(OVCTUDec *const ctudec,
               int16_t *const dst, int16_t *src,
               uint8_t x0, uint8_t y0,
               uint8_t log2_tb_w, uint8_t log2_tb_h,
               uint16_t last_pos,
               uint8_t lfnst_flag, uint8_t lfnst_idx, uint64_t sig_sb_map, uint8_t qp)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_funcs.tr;

    DECLARE_ALIGNED(32, int16_t, tmp)[32*32];
    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[32*32];

    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    if (log2_tb_h > 1 && log2_tb_w > 1) {
        dequant_4x4_sb(ctudec, dequant_coeffs, src, sig_sb_map, log2_tb_w, log2_tb_h, qp, 1, 1, lfnst_flag);
    } else {
        int log2_sb_w = 2;
        int log2_sb_h = 2;

        if (log2_tb_w < 2 || log2_tb_h < 2) {
            log2_sb_w = 1;
            log2_sb_h = 1;
            if (log2_tb_h > 2) log2_sb_h = 3;
            if (log2_tb_w > 2) log2_sb_w = 3;
        }

        int nb_col = (derive_nb_cols(sig_sb_map) >> 2) << log2_sb_h;
        int nb_row = (derive_nb_rows(sig_sb_map) >> 2) << log2_sb_w;

        if (!ctudec->scaling_list_enabled) {
            dequant_non_4x4_sb(ctudec, dequant_coeffs, src, sig_sb_map, log2_tb_w, log2_tb_h, qp);
        } else {
            dequant_non_4x4_sb(ctudec, dequant_coeffs, src, sig_sb_map, log2_tb_w, log2_tb_h, qp);
        }
    }

    if (!last_pos && !lfnst_flag) {

        TRFunc->dc(dst, log2_tb_w, log2_tb_h, dequant_coeffs[0]);

    } else {
        int log2_sb_w = 2;
        int log2_sb_h = 2;

        if (log2_tb_w < 2 || log2_tb_h < 2) {
            log2_sb_w = 1;
            log2_sb_h = 1;
            if (log2_tb_h > 2) log2_sb_h = 3;
            if (log2_tb_w > 2) log2_sb_w = 3;
        }

        int nb_col = (derive_nb_cols(sig_sb_map) >> 2) << log2_sb_h;
        int nb_row = (derive_nb_rows(sig_sb_map) >> 2) << log2_sb_w;

        if (lfnst_flag) {
            /* FIXME separate lfnst mode derivation from lfnst reconstruction */
            int16_t lfnst_sb[16];
            uint8_t is_8x8 = log2_tb_w >= 3 && log2_tb_h >= 3;

            memcpy(lfnst_sb     , &dequant_coeffs[0], sizeof(int16_t) * 4);
            memcpy(lfnst_sb +  4, &dequant_coeffs[1 << log2_tb_w], sizeof(int16_t) * 4);
            memcpy(lfnst_sb +  8, &dequant_coeffs[2 << log2_tb_w], sizeof(int16_t) * 4);
            memcpy(lfnst_sb + 12, &dequant_coeffs[3 << log2_tb_w], sizeof(int16_t) * 4);

            process_lfnst(ctudec, dequant_coeffs, lfnst_sb, log2_tb_w, log2_tb_h,
                          x0, y0, lfnst_idx);


            nb_row = 4 << is_8x8;
            nb_col = 4 << is_8x8;
        }

        memset(&tmp[nb_row << log2_tb_h], 0, sizeof(int16_t) * ((1 << (log2_tb_w + log2_tb_h)) - (nb_row << log2_tb_h)));

        TRFunc->func[DCT_II][log2_tb_h](dequant_coeffs, tmp, tb_w, nb_row, nb_col, TR_SHIFT_V);
        TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);
    }
}

static void
apply_bdpcm_1(int16_t *dst, const int16_t *src, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint8_t tb_w = 1 << log2_tb_w;
    uint8_t tb_h = 1 << log2_tb_h;

    for(int y = 0; y < tb_h; y++) {
      dst[0] = src[0];

      for (int x = 1; x < tb_w; x++) {
        dst[x] = ov_clip((int32_t)dst[x - 1] + (int32_t)src[x], -((1<<15)), (1 << 15) - 1);
      }
      src += tb_w;
      dst += tb_w;
    }
}

static void
apply_bdpcm_2(int16_t *dst, const int16_t *src, uint8_t log2_tb_w, uint8_t log2_tb_h)
{
    uint8_t tb_w = 1 << log2_tb_w;
    uint8_t tb_h = 1 << log2_tb_h;
    for (int x = 0; x < tb_w; x++) {
      dst[x] = src[x];
    }

    for(int y = 0; y < tb_h - 1; y++) {
      for (int x = 0; x < tb_w; x++) {
        dst[x + tb_w] = ov_clip((int32_t)dst[x] + (int32_t)src[x + tb_w], -((1<<15)), (1 << 15) - 1);
      }
      src += tb_w;
      dst += tb_w;
    }
}

static void
rcn_bdpcm_tb(OVCTUDec *const ctu_dec, int16_t *dst, int16_t *src, const struct TBInfo *tb_info,
             uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t bdpcm_dir, uint8_t qp)
{
    DECLARE_ALIGNED(32, int16_t, tmp)[32*32];
    const struct IQScale deq_prms = ctu_dec->rcn_funcs.tmp.derive_dequant_ts(qp, log2_tb_w, log2_tb_h);
    int16_t *bdpcm_src = src;

    if (ctu_dec->sh_ts_disabled && log2_tb_w > 1 && log2_tb_h > 1) {

        reorder_tb_4x4(tmp, src, deq_prms.scale, deq_prms.shift, log2_tb_w, log2_tb_h, tb_info->sig_sb_map);
        bdpcm_src = tmp;
    }

    if (bdpcm_dir) {
        apply_bdpcm_2(dst, bdpcm_src, log2_tb_w, log2_tb_h);
    } else {
        apply_bdpcm_1(dst, bdpcm_src, log2_tb_w, log2_tb_h);
    }

    for (int i = 0; i < ((1 << (log2_tb_w + log2_tb_h)) >> 4); i++) {
        deq_prms.dequant_sb(&dst[16*i], deq_prms.scale, deq_prms.shift);
    }
}

static void
rcn_transform_skip_tb_c(OVCTUDec *const ctu_dec, int16_t *dst, int16_t *src,
                        const struct TBInfo *const tb_info,
                        uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t qp, CUFlags cu_flags)
{
    if (cu_flags & flg_intra_bdpcm_chroma_flag) {
        uint8_t bdpcm_dir = !!(cu_flags & flg_intra_bdpcm_chroma_dir);

        rcn_bdpcm_tb(ctu_dec, dst, src, tb_info, log2_tb_w, log2_tb_h,
                     bdpcm_dir, qp);

    } else {
        if (ctu_dec->sh_ts_disabled) {
            if (log2_tb_w > 1 && log2_tb_h > 1) {
                dequant_4x4_ts(ctu_dec, dst, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp);
            } else {
                memcpy(dst, src, sizeof(int16_t) << (log2_tb_w + log2_tb_h));
                const struct IQScale deq_prms = ctu_dec->rcn_funcs.tmp.derive_dequant_ts(qp, log2_tb_w, log2_tb_h);
                for (int i = 0; i < ((1 << (log2_tb_w + log2_tb_h)) >> 4); i++) {
                    deq_prms.dequant_sb(&dst[16*i], deq_prms.scale, deq_prms.shift);
                }
            }
        } else {
            memcpy(dst, src, sizeof(int16_t) << (log2_tb_w + log2_tb_h));
        }
    }
}


static void
rcn_res_c(OVCTUDec *const ctu_dec, const struct TUInfo *tu_info,
          uint8_t x0, uint8_t y0,
          uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t cbf_mask, uint8_t lfnst_flag, CUFlags cu_flags)
{
    const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_funcs;

    if (cbf_mask & 0x2) {
        const OVBuffInfo *const ctu_buff = &ctu_dec->rcn_ctx.ctu_buff;
        OVSample *const dst_cb = &ctu_buff->cb[(x0) + (y0 * ctu_buff->stride_c)];
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        int16_t *const coeffs_cb = ctu_dec->residual_cb + tu_info->pos_offset;
        int16_t *tr_buff = ctu_dec->transform_buff;

        if (!(tu_info->tr_skip_mask & 0x2)) {
            const struct TBInfo *const tb_info_cb = &tu_info->tb_info[0];
            uint8_t qp = ctu_dec->dequant_cb.qp;

            rcn_residual_c(ctu_dec, tr_buff, coeffs_cb,
                           x0, y0, log2_tb_w, log2_tb_h,
                           tb_info_cb->last_pos, lfnst_flag, tu_info->lfnst_idx, tb_info_cb->sig_sb_map, qp);
        } else {
            const struct TBInfo *const tb_info = &tu_info->tb_info[0];
            int qp = ctu_dec->dequant_cb_skip.qp;

            rcn_transform_skip_tb_c(ctu_dec, tr_buff, coeffs_cb, tb_info,
                                    log2_tb_w, log2_tb_h, qp, cu_flags);

        }

        if (log2_tb_w + log2_tb_h > 2) {
            rcn_func->ict.ict[log2_tb_w][0](tr_buff, dst_cb, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        } else {
            rcn_func->ict.add[log2_tb_w](tr_buff, dst_cb, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        }

        if (!(cu_flags & flg_intra_bdpcm_chroma_flag)) {
            fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        }
    }

    if (cbf_mask & 0x1) {
        const OVBuffInfo *const ctu_buff = &ctu_dec->rcn_ctx.ctu_buff;
        OVSample *const dst_cr = &ctu_buff->cr[(x0) + (y0 * ctu_buff->stride_c)];
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        int16_t *const coeffs_cr = ctu_dec->residual_cr + tu_info->pos_offset;
        int16_t *tr_buff = ctu_dec->transform_buff;

        if (!(tu_info->tr_skip_mask & 0x1)) {
            const struct TBInfo *const tb_info_cr = &tu_info->tb_info[1];
            uint8_t qp = ctu_dec->dequant_cr.qp;
            rcn_residual_c(ctu_dec, tr_buff, coeffs_cr,
                           x0, y0, log2_tb_w, log2_tb_h,
                           tb_info_cr->last_pos, lfnst_flag, tu_info->lfnst_idx, tb_info_cr->sig_sb_map, qp);
        } else {
            const struct TBInfo *const tb_info = &tu_info->tb_info[1];
            int qp = ctu_dec->dequant_cr_skip.qp;

            rcn_transform_skip_tb_c(ctu_dec, ctu_dec->transform_buff, coeffs_cr, tb_info,
                                    log2_tb_w, log2_tb_h, qp, cu_flags);

        }

        if (log2_tb_w + log2_tb_h > 2) {
            rcn_func->ict.ict[log2_tb_w][0](tr_buff, dst_cr, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        } else {
            rcn_func->ict.add[log2_tb_w](tr_buff, dst_cr, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        }

        if (!(cu_flags & flg_intra_bdpcm_chroma_flag)) {
            fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        }
    }
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
}

static void
rcn_jcbcr(OVCTUDec *const ctu_dec, const struct TUInfo *const tu_info,
          uint8_t x0, uint8_t y0, uint8_t log2_tb_w, uint8_t log2_tb_h,
          uint8_t cbf_mask, uint8_t lfnst_flag, CUFlags cu_flags)
{
    const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_funcs;
    const OVBuffInfo *const ctu_buff = &ctu_dec->rcn_ctx.ctu_buff;
    OVSample *const dst_cb = &ctu_buff->cb[x0 + (y0 * ctu_buff->stride_c)];
    OVSample *const dst_cr = &ctu_buff->cr[x0 + (y0 * ctu_buff->stride_c)];
    int16_t *tr_buff = ctu_dec->transform_buff;

    if (!(tu_info->tr_skip_mask & 0x1)) {
        const struct TBInfo *const tb_info = &tu_info->tb_info[0];
        int16_t *const coeffs_jcbcr = ctu_dec->residual_cb + tu_info->pos_offset;
        uint8_t qp;

        if ((cbf_mask&0x3) == 3) {
            qp = ctu_dec->dequant_joint_cb_cr.qp;
        } else if (cbf_mask == 1) {
            qp = ctu_dec->dequant_cr.qp;
        } else {
            qp = ctu_dec->dequant_cb.qp;
        }

        rcn_residual_c(ctu_dec, tr_buff, coeffs_jcbcr,
                       x0, y0, log2_tb_w, log2_tb_h,
                       tb_info->last_pos, lfnst_flag, tu_info->lfnst_idx, tb_info->sig_sb_map, qp);
    } else {
        int16_t *const coeffs_jcbcr = ctu_dec->residual_cb + tu_info->pos_offset;
        int qp = (cbf_mask == 3) ? ctu_dec->dequant_jcbcr_skip.qp
                                 : (cbf_mask == 2) ? ctu_dec->dequant_cb_skip.qp
                                                   : ctu_dec->dequant_cr_skip.qp;

        const struct TBInfo *const tb_info = &tu_info->tb_info[0];

        rcn_transform_skip_tb_c(ctu_dec, tr_buff, coeffs_jcbcr, tb_info,
                                log2_tb_w, log2_tb_h, qp, cu_flags);

    }

    if (!(cu_flags & flg_intra_bdpcm_chroma_flag)) {
        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
        fill_bs_map(&ctu_dec->dbf_info.bs1_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    }

    fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    if ((cbf_mask & 0x3) == 0x3) {
        int qp_bd_offset = ctu_dec->qp_ctx.qp_bd_offset;
        uint8_t    qp = ctu_dec->dequant_joint_cb_cr.qp - qp_bd_offset;

        dbf_fill_qp_map(&ctu_dec->dbf_info.qp_map_cb, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1, qp);
        dbf_fill_qp_map(&ctu_dec->dbf_info.qp_map_cr, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1, qp);
    }

    /* FIXME better organisation based on cbf_mask */
    if (cbf_mask == 3) {
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        if (log2_tb_w + log2_tb_h == 2) scale = 1<<11;
        rcn_func->ict.ict[log2_tb_w][0](tr_buff, dst_cb, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        rcn_func->ict.ict[log2_tb_w][1](tr_buff, dst_cr, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
    } else if (cbf_mask == 2) {
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        if (log2_tb_w + log2_tb_h == 2) scale = 1<<11;
        rcn_func->ict.ict[log2_tb_w][0](tr_buff, dst_cb, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        rcn_func->ict.ict[log2_tb_w][2](tr_buff, dst_cr, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
    } else {
        int16_t scale  =  ctu_dec->lmcs_info.scale_c_flag ? ctu_dec->lmcs_info.lmcs_chroma_scale : 1<< 11;
        if (log2_tb_w + log2_tb_h == 2) scale = 1<<11;
        rcn_func->ict.ict[log2_tb_w][0](tr_buff, dst_cr, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
        rcn_func->ict.ict[log2_tb_w][2](tr_buff, dst_cb, ctu_buff->stride_c, log2_tb_w, log2_tb_h, scale);
    }

}

static void
rcn_isp_tu(OVCTUDec *const ctudec, const struct TBInfo *const tb_info, uint8_t log2_tb_w, uint8_t log2_tb_h, int16_t *coeffs_y, const struct ISPTUInfo *const tu_info, uint8_t type_v, uint8_t type_h, int8_t lfnst_intra_mode)
{
    DECLARE_ALIGNED(32, int16_t, tmp)[32*64];
    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[32*64];

    const struct TRFunctions *const TRFunc = &ctudec->rcn_funcs.tr;

    int16_t *src = coeffs_y;
    int16_t *dst = ctudec->transform_buff;

    int nb_col = derive_nb_cols(tb_info->sig_sb_map);
    int nb_row = derive_nb_rows(tb_info->sig_sb_map);

    int tb_h = 1 << log2_tb_h;
    int tb_w = 1 << log2_tb_w;

    int qp = ctudec->dequant_luma.qp;

    dequant_4x4_sb(ctudec, dequant_coeffs, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp, 1, 0, tu_info->lfnst_flag);

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (tu_info->lfnst_flag) {
        uint8_t lfnst_idx = tu_info->lfnst_idx;
        int16_t lfnst_sb[16];
        uint8_t log2_stride = OVMIN(5, log2_tb_w);
        uint8_t is_8x8 = log2_tb_w >= 3 && log2_tb_h >= 3;

        memcpy(lfnst_sb     , &dequant_coeffs[0               ], sizeof(int16_t) * 4);
        memcpy(lfnst_sb +  4, &dequant_coeffs[1 << log2_stride], sizeof(int16_t) * 4);
        memcpy(lfnst_sb +  8, &dequant_coeffs[2 << log2_stride], sizeof(int16_t) * 4);
        memcpy(lfnst_sb + 12, &dequant_coeffs[3 << log2_stride], sizeof(int16_t) * 4);

        process_lfnst_luma(ctudec, dequant_coeffs, lfnst_sb, OVMIN(5,log2_tb_w), OVMIN(5,log2_tb_h),
                           lfnst_idx, lfnst_intra_mode);

        nb_row = 4 << is_8x8;
        nb_col = 4 << is_8x8;

        /* lfnst forces IDCT II usage */
        type_v = type_h = DCT_II;
    }

    tb_w = OVMIN(32, tb_w);

    TRFunc->func[type_v][OVMIN(log2_tb_h, 6)](dequant_coeffs, tmp, tb_w, nb_row, nb_col, TR_SHIFT_V);
    TRFunc->func[type_h][OVMIN(log2_tb_w, 6)](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);
}

static void
rcn_2xX_tb(OVCTUDec *const ctudec, const struct TBInfo *const tb_info, uint8_t log2_tb_h,
           int16_t *coeffs_y, uint8_t type_v, uint8_t type_h)
{
    const struct TRFunctions *const TRFunc = &ctudec->rcn_funcs.tr;

    DECLARE_ALIGNED(32, int16_t, tmp)[2*64];
    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[2*64];

    int16_t *src = coeffs_y;
    int16_t *dst = ctudec->transform_buff;

    int nb_col = (derive_nb_cols(tb_info->sig_sb_map) >> 2) << 3;
    int nb_row = (derive_nb_rows(tb_info->sig_sb_map) >> 2) << 1;

    const uint8_t log2_tb_w = 1;

    int tb_h = 1 << log2_tb_h;
    int tb_w = 1 << log2_tb_w;

    int qp = ctudec->dequant_luma.qp;

    dequant_non_4x4_sb(ctudec, dequant_coeffs, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp);

    memset(tmp, 0, sizeof(int16_t) << (1 + log2_tb_h));

    TRFunc->func[type_v][OVMIN(log2_tb_h,6)](dequant_coeffs, tmp, tb_w, nb_row, nb_col, TR_SHIFT_V);
    TRFunc->func[type_h][OVMIN(log2_tb_w,6)](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);
}

static void
rcn_1xX_tb(OVCTUDec *const ctudec, const struct TBInfo *const tb_info, uint8_t log2_tb_h, int16_t *src, uint8_t type_v)
{
    const struct TRFunctions *const TRFunc = &ctudec->rcn_funcs.tr;
    const uint8_t log2_tb_w = 0;
    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[2*64];

    int tb_h = 1 << log2_tb_h;
    int tb_w = 1 << log2_tb_w;

    int qp = ctudec->dequant_luma.qp;

    dequant_non_4x4_sb(ctudec, dequant_coeffs, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp);

    TRFunc->func[type_v][OVMIN(log2_tb_h,6)](dequant_coeffs, ctudec->transform_buff, tb_w, tb_w, tb_info->last_pos & 0x1F, TR_SHIFT_H + 1);

}


static void
rcn_tu_isp_v(OVCTUDec *const ctudec,
             unsigned int x0, unsigned int y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t lfnst_intra_mode,
             const struct ISPTUInfo *const tu_info, uint8_t i,
             uint8_t type_v, uint8_t type_h, uint8_t offset)
{
    const struct TRFunctions *TRFunc = &ctudec->rcn_funcs.tr;
    const struct RCNFunctions *const rcn_func = &ctudec->rcn_funcs;
    const struct TBInfo *const tb_info = &tu_info->tb_info[i];
    const struct OVBuffInfo *const ctu_buff = &ctudec->rcn_ctx.ctu_buff;

    int16_t *coeffs_y = ctudec->residual_y + (i << (log2_tb_w + log2_tb_h));

    if (log2_tb_w) {

        if (log2_tb_w == 1) {
            rcn_2xX_tb(ctudec, tb_info, log2_tb_h, coeffs_y, type_v, type_h);
        } else {
            rcn_isp_tu(ctudec, tb_info, log2_tb_w, log2_tb_h, coeffs_y, tu_info,
                       type_v, type_h, lfnst_intra_mode);
        }

    } else {
        rcn_1xX_tb(ctudec, tb_info, log2_tb_h, coeffs_y, type_v);
    }

    OVSample *dst  = &ctu_buff->y[x0 + y0 * ctu_buff->stride];
    const int16_t *src  = ctudec->transform_buff;

    rcn_func->ict.add[log2_tb_w](src, dst, ctu_buff->stride, log2_tb_w, log2_tb_h, 0);
}

static void
rcn_Xx2_tb(OVCTUDec *const ctudec, const struct TBInfo *const tb_info, uint8_t log2_tb_w,
           int16_t *coeffs_y, uint8_t type_v, uint8_t type_h)
{
    const struct TRFunctions *const TRFunc = &ctudec->rcn_funcs.tr;

    DECLARE_ALIGNED(32, int16_t, tmp)[2*64];
    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[2*64];

    int16_t *src = coeffs_y;
    int16_t *dst = ctudec->transform_buff;

    int nb_col = (derive_nb_cols(tb_info->sig_sb_map) >> 2) << 1;
    int nb_row = (derive_nb_rows(tb_info->sig_sb_map) >> 2) << 3;

    const uint8_t log2_tb_h = 1;

    int tb_h = 1 << log2_tb_h;
    int tb_w = 1 << log2_tb_w;

    int qp = ctudec->dequant_luma.qp;

    dequant_non_4x4_sb(ctudec, dequant_coeffs, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp);

    memset(tmp, 0, sizeof(int16_t) << (1 + log2_tb_w));

    TRFunc->func[type_v][OVMIN(log2_tb_h,6)](dequant_coeffs, tmp, tb_w, nb_row, nb_col, TR_SHIFT_V);
    TRFunc->func[type_h][OVMIN(log2_tb_w,6)](tmp, dst, tb_h, tb_h, nb_row, TR_SHIFT_H);
}

static void
rcn_Xx1_tb(OVCTUDec *const ctudec, const struct TBInfo *const tb_info, uint8_t log2_tb_w, int16_t *src, uint8_t type_h)
{
    const struct TRFunctions *const TRFunc = &ctudec->rcn_funcs.tr;

    DECLARE_ALIGNED(32, int16_t, dequant_coeffs)[2*64];

    const uint8_t log2_tb_h = 0;

    int tb_h = 1 << log2_tb_h;
    int tb_w = 1 << log2_tb_w;

    int qp = ctudec->dequant_luma.qp;

    dequant_non_4x4_sb(ctudec, dequant_coeffs, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp);

    TRFunc->func[type_h][OVMIN(log2_tb_w,6)](dequant_coeffs, ctudec->transform_buff, tb_h, tb_h, tb_info->last_pos & 0x1F, TR_SHIFT_H + 1);

}

static void
rcn_tu_isp_h(OVCTUDec *const ctudec,
             unsigned int x0, unsigned int y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             uint8_t lfnst_intra_mode,
             const struct ISPTUInfo *const tu_info, uint8_t i,
             uint8_t type_v, uint8_t type_h, uint8_t offset)
{
    const struct TRFunctions *TRFunc = &ctudec->rcn_funcs.tr;
    const struct RCNFunctions *const rcn_func = &ctudec->rcn_funcs;
    const struct TBInfo *const tb_info = &tu_info->tb_info[i];
    const struct OVBuffInfo *const ctu_buff = &ctudec->rcn_ctx.ctu_buff;

    int16_t *coeffs_y = ctudec->residual_y + (i << (log2_tb_w + log2_tb_h));

    if (log2_tb_h) {

        if (log2_tb_h == 1) {
            rcn_Xx2_tb(ctudec, tb_info, log2_tb_w, coeffs_y, type_v, type_h);
        } else {
            rcn_isp_tu(ctudec, tb_info, log2_tb_w, log2_tb_h, coeffs_y, tu_info,
                       type_v, type_h, lfnst_intra_mode);
        }

    } else {
        rcn_Xx1_tb(ctudec, tb_info, log2_tb_w, coeffs_y, type_h);
    }

    OVSample *dst  = &ctu_buff->y[x0 + y0 * ctu_buff->stride];
    const int16_t *src  = ctudec->transform_buff;

    rcn_func->ict.add[log2_tb_w](src, dst, ctu_buff->stride, log2_tb_w, log2_tb_h, 0);
}


static void
recon_isp_subtree_v(OVCTUDec *const ctudec,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    uint8_t intra_mode,
                    const struct ISPTUInfo *const tu_info)
{
    uint8_t cbf_flags = tu_info->cbf_mask;
    uint8_t lfnst_flag = tu_info->lfnst_flag;

    int offset_x;
    int log2_pb_w = log2_cb_w - 2;
    int nb_pb;
    int pb_w;

    if (log2_cb_h < 4 && (log2_pb_w <= (4 - log2_cb_h))) {
        log2_pb_w = 4 - log2_cb_h;
    }

    nb_pb = (1 << log2_cb_w) >> log2_pb_w;
    pb_w =  1 << log2_pb_w;
    offset_x = 0;

    uint8_t type_h = ctudec->mts_enabled && log2_pb_w <= 4 && log2_pb_w > 1 ? DST_VII : DCT_II;
    uint8_t type_v = ctudec->mts_enabled && log2_cb_h <= 4 ? DST_VII : DCT_II;
    int8_t lfnst_intra_mode;
    const struct OVBuffInfo *const ctu_buff = &ctudec->rcn_ctx.ctu_buff;
    int i;

    if (lfnst_flag) {
        lfnst_intra_mode = drv_lfnst_mode_l(log2_cb_w, log2_cb_h, intra_mode);
    }

    for (i = 0; i < nb_pb; ++i) {
        uint8_t cbf = (cbf_flags >> (nb_pb - i - 1)) & 0x1;

        /* On vertical ISP the intra prediction can only be performed with a min width of 4
           thus requiring to perform prediction on 2 PBs for size 2xX and all PBs in the
           case of 1xX. Even when transforms are smaller */

        if (!(offset_x & 0x3)) {

            uint8_t log2_pred_w = log2_pb_w < 2 ? 2 : log2_pb_w;

            ctudec->rcn_funcs.intra_pred_isp(ctudec, &ctu_buff->y[0],
                                             ctu_buff->stride, intra_mode, x0, y0,
                                             log2_pred_w, log2_cb_h,
                                             log2_cb_w, log2_cb_h, offset_x, 0);

            fill_ctb_bound(&ctudec->dbf_info,      x0, y0, log2_pred_w, log2_cb_h);
            fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_pred_w, log2_cb_h);
        }

        if (cbf) {
            rcn_tu_isp_v(ctudec, x0, y0, log2_pb_w, log2_cb_h,
                         lfnst_intra_mode, tu_info, i,
                         type_v, type_h, offset_x);
        }

        x0       += pb_w;
        offset_x += pb_w;
    }
}

static void
recon_isp_subtree_h(OVCTUDec *const ctudec,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    uint8_t intra_mode,
                    const struct ISPTUInfo *const tu_info)
{
    const struct TRFunctions *TRFunc = &ctudec->rcn_funcs.tr;
    const struct RCNFunctions *const rcn_func = &ctudec->rcn_funcs;

    int log2_pb_h = log2_cb_h - 2;
    int nb_pb;
    uint8_t cbf_flags = tu_info->cbf_mask;
    uint8_t lfnst_flag = tu_info->lfnst_flag;
    int8_t lfnst_intra_mode;
    int pb_h, offset_y;

    if (lfnst_flag) {
        lfnst_intra_mode = drv_lfnst_mode_l(log2_cb_w, log2_cb_h, intra_mode);
    }
    // width < 16 imposes restrictions on split numbers
    if (log2_cb_w < 4 && (log2_pb_h <= (4 - log2_cb_w))) {
        log2_pb_h = 4 - log2_cb_w;
    }

    nb_pb = (1 << log2_cb_h) >> log2_pb_h;
    pb_h =  1 << log2_pb_h;
    offset_y = 0;

    uint8_t type_h = ctudec->mts_enabled && log2_cb_w <= 4 ? DST_VII : DCT_II;
    uint8_t type_v = ctudec->mts_enabled && log2_pb_h <= 4 && log2_pb_h > 1 ? DST_VII : DCT_II;
    const struct OVBuffInfo *const ctu_buff = &ctudec->rcn_ctx.ctu_buff;
    int i;

    for (i = 0; i < nb_pb; ++i) {

        uint8_t cbf = (cbf_flags >> (nb_pb - i - 1)) & 0x1;

        if (!(offset_y & 0x3) ) {
            fill_ctb_bound(&ctudec->dbf_info, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
            fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_cb_w, log2_pb_h >=2 ? log2_pb_h : 2);
        }

        ctudec->rcn_funcs.intra_pred_isp(ctudec, &ctu_buff->y[0],
                                         ctu_buff->stride, intra_mode, x0, y0,
                                         log2_cb_w, log2_pb_h, log2_cb_w, log2_cb_h, 0, offset_y);
        if (cbf) {
            rcn_tu_isp_h(ctudec, x0, y0, log2_cb_w, log2_pb_h,
                         lfnst_intra_mode, tu_info, i,
                         type_v, type_h, offset_y);
        }
        y0 += pb_h;
        offset_y += pb_h;
    }
}

static void
rcn_transform_skip_tb_l(OVCTUDec *const ctu_dec, int16_t *dst, int16_t *src,
                        const struct TBInfo *const tb_info,
                        uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t qp, CUFlags cu_flags)
{
    if (cu_flags & flg_intra_bdpcm_luma_flag) {
        uint8_t bdpcm_dir = !!(cu_flags & flg_intra_bdpcm_luma_dir);

        rcn_bdpcm_tb(ctu_dec, dst, src, tb_info, log2_tb_w, log2_tb_h,
                     bdpcm_dir, qp);

    } else {
        if (ctu_dec->sh_ts_disabled) {
            dequant_4x4_ts(ctu_dec, dst, src, tb_info->sig_sb_map, log2_tb_w, log2_tb_h, qp);
        } else {
            memcpy(dst, src, sizeof(int16_t) << (log2_tb_w + log2_tb_h));
        }
    }
}

static void
rcn_tu_st(OVCTUDec *const ctu_dec,
          uint8_t x0, uint8_t y0,
          uint8_t log2_tb_w, uint8_t log2_tb_h,
          CUFlags cu_flags, uint8_t cbf_mask,
          const struct TUInfo *const tu_info)
{
    uint8_t cbf_flag_l = cbf_mask & 0x10;
    uint8_t jcbcr_flag = cbf_mask & 0x8;
    uint8_t cbf_mask_c = cbf_mask & 0x3;

    if (cbf_flag_l) {
        const struct TBInfo *const tb_info = &tu_info->tb_info[2];
        const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_funcs;
        const struct OVBuffInfo *const ctu_buff = &ctu_dec->rcn_ctx.ctu_buff;
        int16_t *tr_buff = ctu_dec->transform_buff;

        if (!(tu_info->tr_skip_mask & 0x10)) {
            int lim_sb_s = ((((tb_info->last_pos >> 8)) >> 2) + (((tb_info->last_pos & 0xFF))>> 2) + 1) << 2;
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            uint8_t is_mip = !!(cu_flags & flg_mip_flag);
            uint8_t is_intra = !!(cu_flags & flg_pred_mode_flag);
            is_mip |= !is_intra;
            rcn_residual(ctu_dec, tr_buff, coeffs_y, x0, y0, log2_tb_w, log2_tb_h,
                         tu_info->cu_mts_flag, tu_info->cu_mts_idx,
                         !tb_info->last_pos, tu_info->lfnst_flag, is_mip, tu_info->lfnst_idx, tb_info->sig_sb_map);

        } else {
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            int qp = ctu_dec->dequant_luma_skip.qp;
            rcn_transform_skip_tb_l(ctu_dec, tr_buff, coeffs_y, tb_info,
                                    log2_tb_w, log2_tb_h, qp, cu_flags);
        }

        /* FIXME use transform add optimization */
        rcn_func->ict.add[log2_tb_w](tr_buff, &ctu_buff->y[x0 + y0 * ctu_buff->stride], ctu_buff->stride, log2_tb_w, log2_tb_h, 0);
            fill_bs_map(&ctu_dec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
        if (!(cu_flags & flg_intra_bdpcm_luma_flag)) {
            if (cu_flags & flg_pred_mode_flag) {
                fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_tb_w, log2_tb_h);
            }
        }
    }

    /* FIXME Avoid reprocessing CCLM from here by recontructing at the end of transform tree */
    if (cu_flags & flg_pred_mode_flag) {
        uint8_t x0_unit = (x0) >> LOG2_MIN_CU_S;
        uint8_t y0_unit = (y0) >> LOG2_MIN_CU_S;
        uint8_t nb_unit_w = (1 << log2_tb_w) >> LOG2_MIN_CU_S;
        uint8_t nb_unit_h = (1 << log2_tb_h) >> LOG2_MIN_CU_S;

        ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field_c, x0_unit,
                                    y0_unit, nb_unit_w, nb_unit_h);

        if (!(cu_flags & flg_intra_bdpcm_chroma_flag)) {
            fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0, y0, log2_tb_w, log2_tb_h);
        }

        ctu_dec->rcn_funcs.intra_pred_c(&ctu_dec->rcn_ctx, ctu_dec->intra_mode_c, x0 >> 1, y0 >> 1,
                                        log2_tb_w - 1, log2_tb_h - 1, cu_flags);
    }

    if (jcbcr_flag) {

        rcn_jcbcr(ctu_dec, tu_info, x0 >> 1, y0 >> 1, log2_tb_w - 1, log2_tb_h - 1, cbf_mask_c, 0, cu_flags);

    } else if (cbf_mask_c) {

        rcn_res_c(ctu_dec, tu_info, x0 >> 1, y0 >> 1, log2_tb_w - 1, log2_tb_h - 1, cbf_mask_c, 0, cu_flags);

    }

    fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
    fill_ctb_bound_c(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
}


static void
rcn_tu_l(OVCTUDec *const ctu_dec,
         uint8_t x0, uint8_t y0,
         uint8_t log2_tb_w, uint8_t log2_tb_h,
         CUFlags cu_flags, uint8_t cbf_mask,
         const struct TUInfo *const tu_info)
{
    const struct TBInfo *const tb_info = &tu_info->tb_info[2];
    if (cbf_mask) {
        const struct RCNFunctions *const rcn_func = &ctu_dec->rcn_funcs;
        const OVBuffInfo *const ctu_buff = &ctu_dec->rcn_ctx.ctu_buff;
        int16_t dst_stride = ctu_buff->stride;

        if (!(tu_info->tr_skip_mask & 0x10)) {
            int lim_sb_s = ((((tb_info->last_pos >> 8)) >> 2) + (((tb_info->last_pos & 0xFF))>> 2) + 1) << 2;
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            uint8_t is_mip = !!(cu_flags & flg_mip_flag);
            uint8_t is_intra = !!(cu_flags & flg_pred_mode_flag);
            is_mip |= !is_intra;
            rcn_residual(ctu_dec, ctu_dec->transform_buff, coeffs_y, x0, y0, log2_tb_w, log2_tb_h,
                         tu_info->cu_mts_flag, tu_info->cu_mts_idx,
                         !tb_info->last_pos, tu_info->lfnst_flag, is_mip, tu_info->lfnst_idx, tb_info->sig_sb_map);

        } else {
            int16_t *const coeffs_y = ctu_dec->residual_y + tu_info->pos_offset;
            int qp = ctu_dec->dequant_luma_skip.qp;
            rcn_transform_skip_tb_l(ctu_dec, ctu_dec->transform_buff, coeffs_y, tb_info,
                                    log2_tb_w, log2_tb_h, qp, cu_flags);
        }

        if (!(cu_flags & flg_intra_bdpcm_luma_flag)) {
            if (cu_flags & flg_pred_mode_flag) {
                fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_tb_w, log2_tb_h);
            }
        }
            fill_bs_map(&ctu_dec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);

        /* FIXME use transform add optimization */
        rcn_func->ict.add[log2_tb_w](ctu_dec->transform_buff, &ctu_buff->y[x0 + y0 * dst_stride],
                                     dst_stride, log2_tb_w, log2_tb_h, 0);
    }
    fill_ctb_bound(&ctu_dec->dbf_info, x0, y0, log2_tb_w, log2_tb_h);
}

static void
rcn_tu_c(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
         uint8_t log2_tb_w, uint8_t log2_tb_h,
         CUFlags cu_flags, uint8_t cbf_mask,
         const struct TUInfo *const tu_info)
{
    uint8_t jcbcr_flag = cbf_mask & 0x8;
    uint8_t cbf_mask_c = cbf_mask & 0x3;

    uint8_t x0_unit = (x0 << 1) >> LOG2_MIN_CU_S;
    uint8_t y0_unit = (y0 << 1) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (2 << log2_tb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (2 << log2_tb_h) >> LOG2_MIN_CU_S;

    ctu_field_set_rect_bitfield(&ctu_dec->rcn_ctx.progress_field_c, x0_unit,
                                y0_unit, nb_unit_w, nb_unit_h);

    ctu_dec->rcn_funcs.intra_pred_c(&ctu_dec->rcn_ctx, ctu_dec->intra_mode_c, x0, y0, log2_tb_w, log2_tb_h, cu_flags);

    fill_ctb_bound_c(&ctu_dec->dbf_info, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    if (!(cu_flags & flg_intra_bdpcm_chroma_flag)) {
        fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0 << 1, y0 << 1, log2_tb_w + 1, log2_tb_h + 1);
    }

    if (jcbcr_flag) {

        rcn_jcbcr(ctu_dec, tu_info, x0, y0, log2_tb_w, log2_tb_h, cbf_mask_c, tu_info->lfnst_flag, cu_flags);

    } else if (cbf_mask_c) {

        rcn_res_c(ctu_dec, tu_info, x0, y0, log2_tb_w, log2_tb_h, cbf_mask_c, tu_info->lfnst_flag, cu_flags);

    }
}

static void
rcn_intra_tu(OVCTUDec *const ctudec, uint8_t x0, uint8_t y0,
             uint8_t log2_tb_w, uint8_t log2_tb_h, CUFlags cu_flags)
{
    uint8_t mip_flag = cu_flags & flg_mip_flag;

    uint8_t intra_mode = ctudec->intra_mode;

    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_tb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_tb_h) >> LOG2_MIN_CU_S;

    if (mip_flag){
        const struct OVRCNCtx *rcn = &ctudec->rcn_ctx;
        ctudec->rcn_funcs.mip.rcn_intra_mip(rcn, x0, y0, log2_tb_w, log2_tb_h, ctudec->cu_opaque);

    } else {
        const OVBuffInfo *const ctu_buff = &ctudec->rcn_ctx.ctu_buff;
        uint8_t isp_flag = !!(cu_flags & flg_isp_flag);

        if (!isp_flag) {
            uint8_t mrl_flag = !!(cu_flags & flg_mrl_flag);
            if (!mrl_flag) {
                ctudec->rcn_funcs.intra_pred(&ctudec->rcn_ctx,
                                             &ctudec->rcn_ctx.ctu_buff,
                                             intra_mode, x0, y0,
                                             log2_tb_w, log2_tb_h, cu_flags);
            }

            if (mrl_flag){
                uint8_t mrl_idx = ctudec->cu_opaque;
                ctudec->rcn_funcs.intra_pred_mrl(ctudec, ctu_buff->y,
                                                 ctu_buff->stride, intra_mode, x0, y0,
                                                 log2_tb_w, log2_tb_h,
                                                 mrl_idx);
            }
        }
    }

    if (!(cu_flags & flg_intra_bdpcm_luma_flag)) {
        fill_bs_map(&ctudec->dbf_info.bs2_map, x0, y0, log2_tb_w, log2_tb_h);
    }

    ctu_field_set_rect_bitfield(&ctudec->rcn_ctx.progress_field,
                                x0_unit, y0_unit,
                                nb_unit_w, nb_unit_h);
}

static void
rcn_res_wrap(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
             uint8_t log2_tb_w, uint8_t log2_tb_h, CUFlags cu_flags,
             const struct TUInfo *const tu_info)
{
    uint8_t cbf_mask = tu_info->cbf_mask;
    if (ctu_dec->transform_unit == &transform_unit_st) {
        if (cu_flags & flg_pred_mode_flag) {
            rcn_intra_tu(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags);
        }
        rcn_tu_st(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, cbf_mask, tu_info);
    } else if (ctu_dec->transform_unit == &transform_unit_l) {
        if (cu_flags & flg_pred_mode_flag) {
            rcn_intra_tu(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags);
        }
        rcn_tu_l(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, cbf_mask, tu_info);
    } else if (ctu_dec->transform_unit == &transform_unit_c) {
        rcn_tu_c(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, cbf_mask, tu_info);
    }
}

static void
rcn_transform_tree(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
                   uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t log2_max_tb_s,
                   uint8_t tr_depth, CUFlags cu_flags, const struct TUInfo *const tu_info)
{
    uint8_t split_v = log2_tb_w > log2_max_tb_s;
    uint8_t split_h = log2_tb_h > log2_max_tb_s;
    uint8_t nb_subtrees = tr_depth ? 1 : (1 << (split_v + split_h));

    if (log2_tb_w > 6 && log2_tb_h < 7) {

        rcn_transform_tree(ctu_dec, x0, y0, 6, log2_tb_h,
                       log2_max_tb_s, tr_depth + 1, cu_flags,
                        &tu_info[0]);

        rcn_transform_tree(ctu_dec, x0 + 64, y0, 6, log2_tb_h,
                       log2_max_tb_s, tr_depth + 1, cu_flags,
                        &tu_info[8]);
        return;
    } else if (log2_tb_h > 6 && log2_tb_w < 7) {
        rcn_transform_tree(ctu_dec, x0, y0, log2_tb_w, 6,
                       log2_max_tb_s, tr_depth + 1, cu_flags,
                        &tu_info[0]);

        rcn_transform_tree(ctu_dec, x0, y0 + 64, log2_tb_w, 6,
                       log2_max_tb_s, tr_depth + 1, cu_flags,
                        &tu_info[8]);
        return;
    }

    if (split_v || split_h) {
        unsigned int tb_w1 = ((1 << log2_tb_w) >> split_v);
        unsigned int tb_h1 = ((1 << log2_tb_h) >> split_h);

        unsigned int log2_tb_w1 = log2_tb_w - split_v;
        unsigned int log2_tb_h1 = log2_tb_h - split_h;

        rcn_transform_tree(ctu_dec, x0, y0,
                           log2_tb_w1, log2_tb_h1,
                           log2_max_tb_s, tr_depth + 1,  cu_flags, &tu_info[0 * nb_subtrees]);
        if (split_v) {
            rcn_transform_tree(ctu_dec, x0 + tb_w1, y0,
                               log2_tb_w1, log2_tb_h1,
                               log2_max_tb_s,  tr_depth + 1, cu_flags, &tu_info[1 * nb_subtrees]);
        }

        if (split_h) {
            rcn_transform_tree(ctu_dec, x0, y0 + tb_h1,
                               log2_tb_w1, log2_tb_h1,
                               log2_max_tb_s, tr_depth + 1,  cu_flags, &tu_info[2 * nb_subtrees]);
        }

        if (split_h && split_v) {
            rcn_transform_tree(ctu_dec, x0 + tb_w1, y0 + tb_h1,
                               log2_tb_w1, log2_tb_h1,
                               log2_max_tb_s, tr_depth + 1,  cu_flags, &tu_info[3 * nb_subtrees]);
        }

    } else {
        rcn_res_wrap(ctu_dec, x0, y0, log2_tb_w, log2_tb_h, cu_flags, tu_info);
        if (ctu_dec->tmp_ciip) {
            fill_bs_map(&ctu_dec->dbf_info.bs2_map, x0, y0, log2_tb_w, log2_tb_h);
            fill_bs_map(&ctu_dec->dbf_info.bs2_map_c, x0, y0, log2_tb_w, log2_tb_h);
        }
    }
}

void
BD_DECL(rcn_init_transform_trees)(struct RCNFunctions *rcn_funcs)
{

    rcn_funcs->tmp.rcn_transform_tree = &rcn_transform_tree;

    rcn_funcs->tmp.rcn_tu_c = &rcn_tu_c;

    rcn_funcs->tmp.rcn_tu_st = &rcn_tu_st;

    rcn_funcs->tmp.recon_isp_subtree_h = &recon_isp_subtree_h;

    rcn_funcs->tmp.recon_isp_subtree_v = &recon_isp_subtree_v;
}
