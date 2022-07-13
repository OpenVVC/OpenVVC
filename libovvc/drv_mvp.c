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
#include <string.h>

#include "ovdefs.h"
#include "ovutils.h"

/* FIXME find a place for mv structures */
#include "ctudec.h"
#include "drv_utils.h"
#include "dec_structures.h"
#include "ovdpb.h"

#define POS_MASK(x, w) ((uint64_t) 1 << ((((x + 1)) + ((w)))))

#define TMVP_POS_MASK(y) ((uint64_t) 1 << ((y) + 1))

#define OFFSET_BUFF(x,y) (35 + x + (y) * 34)
#define MV_CMP(mv0,mv1) ((mv0.x == mv1.x) && (mv0.y == mv1.y))
#define MV_CMP2(mv0,mv1) ((mv0.x == mv1.x) && (mv0.y == mv1.y) && (mv0.ref_idx == mv1.ref_idx))

#define PB_POS_IN_BUF(x,y) (35 + (x) + ((y) * 34))
#define TMVP_BUFF_STRIDE 17
#define TMVP_POS_IN_BUF(x,y) ((x >> 1) + (((y >> 1)) * TMVP_BUFF_STRIDE))

#define TMVP_POS_IN_BUF2(x,y) ((x >> 1) + (((y >> 1)) * tmvp->pln0_stride))

#define MV_MANTISSA_BITCOUNT 6
#define MV_MANTISSA_UPPER_LIMIT ((1 << (MV_MANTISSA_BITCOUNT - 1)) - 1)
#define MV_MANTISSA_LIMIT (1 << (MV_MANTISSA_BITCOUNT - 1))

#define MV_BITS  18
#define MV_MAX   ((1 << (MV_BITS - 1)) - 1)
#define MV_MIN  (-(1 << (MV_BITS - 1)))

#define LOG2_MIN_CU_S 2

static int16_t
tmvp_compute_scale(int16_t dist_current, int16_t dist_colocated)
{
    int32_t scale;
    if (dist_current == dist_colocated || !dist_colocated)
        return 256;

    scale = dist_current * ((0x4000 + OVABS(dist_colocated >> 1)) / dist_colocated);
    scale += 32;
    scale >>= 6;

    scale = ov_clip_intp2(scale, 13);

    return (int16_t)scale;
}

OVMV
drv_change_precision_mv(OVMV mv, int src, int dst)
{
    int shift = dst - src;

    if (shift >= 0) {
        mv.x = (uint32_t)mv.x << shift;
        mv.y = (uint32_t)mv.y << shift;
    } else {
        shift *= -1;
        const int offset = 1 << (shift - 1);
        mv.x = mv.x >= 0 ? (mv.x + offset - 1) >> shift : (mv.x + offset) >> shift;
        mv.y = mv.y >= 0 ? (mv.y + offset - 1) >> shift : (mv.y + offset) >> shift;
    }

    return mv;
}

static inline struct MV
drv_change_precision_mv2(struct MV mv, int src, int dst)
{
    int shift = dst - src;

    if (shift >= 0) {
        mv.x = (uint32_t)mv.x << shift;
        mv.y = (uint32_t)mv.y << shift;
    } else {
        shift *= -1;
        const int offset = 1 << (shift - 1);
        mv.x = mv.x >= 0 ? (mv.x + offset - 1) >> shift : (mv.x + offset) >> shift;
        mv.y = mv.y >= 0 ? (mv.y + offset - 1) >> shift : (mv.y + offset) >> shift;
    }

    return mv;
}

OVMV
drv_round_to_precision_mv(OVMV mv, int src, int dst)
{
    OVMV mv1 = drv_change_precision_mv(mv, src, dst);
    return drv_change_precision_mv(mv1, dst, src);
}

static inline struct MV
drv_round_to_precision_mv2(struct MV mv, int src, int dst)
{
    struct MV mv1 = drv_change_precision_mv2(mv, src, dst);
    return drv_change_precision_mv2(mv1, dst, src);
}

static inline uint8_t
mi_cmp(const VVCMergeInfo *const a, const VVCMergeInfo *const b)
{
    if (a->inter_dir == b->inter_dir) {
        uint8_t is_neq = 0;
        if (a->inter_dir & 0x1) {
            is_neq += !MV_CMP(a->mv0, b->mv0);
            is_neq += a->mv0.ref_idx != b->mv0.ref_idx;
        }
        if (a->inter_dir & 0x2) {
            is_neq += !MV_CMP(a->mv1, b->mv1);
            is_neq += a->mv1.ref_idx != b->mv1.ref_idx;
        }

        return !is_neq;
    }
    return 0;
}

/* Compress and uncompress Motion vectors used by TMVP
 * FIXME there might be a more straight forward way of
 * doing this
 */
static inline int
tmvp_round(int32_t val)
{
  int sign  = val >> 31;

  if ((val ^ sign) - !!sign > 31) {
      int scale = floor_log2((val ^ sign) | MV_MANTISSA_UPPER_LIMIT) - (MV_MANTISSA_BITCOUNT - 1);
      int round = (1 << scale) >> 1;
      int n     = (val + round) >> scale;
      int exponent  = scale + ((n ^ sign) >> (MV_MANTISSA_BITCOUNT - 1));
      int mantissa  = (n & MV_MANTISSA_UPPER_LIMIT) | ((uint32_t)sign << (MV_MANTISSA_BITCOUNT - 1));
      return (uint32_t)(mantissa ^ MV_MANTISSA_LIMIT) << (exponent - !!exponent);
  } else {
      return val;
  }
}

static inline struct MV
tmvp_round_mv(struct MV mv)
{
    mv.x = tmvp_round(mv.x);
    mv.y = tmvp_round(mv.y);
    return mv;
}

static void
hmvp_add_cand_1(const struct HMVPLUT *const hmvp_lut,
                OVMV *const cand_list,
                int *const nb_cand, uint8_t inter_dir, uint8_t ref_idx, uint8_t opp_ref_idx)
{
    int max_nb_cand = OVMIN(4, hmvp_lut->nb_mv);
    int i;
    for (i = 1; i <= max_nb_cand && *nb_cand < 2;  ++i) {
        if(hmvp_lut->dir[i - 1] & inter_dir) {
            OVMV hmvp_cand = inter_dir & 0x1 ? hmvp_lut->hmv0[i - 1]
                                             : hmvp_lut->hmv1[i - 1];
            if (hmvp_cand.ref_idx == ref_idx) {
                cand_list[(*nb_cand)++] = hmvp_cand;
            }
        }

        if (*nb_cand == 2) {
            return;
        }

        if (hmvp_lut->dir[i - 1] & (3 - inter_dir)) {
            OVMV hmvp_cand = (3 - inter_dir) & 0x1 ? hmvp_lut->hmv0[i - 1]
                                                   : hmvp_lut->hmv1[i - 1];
            if (hmvp_cand.ref_idx == opp_ref_idx) {
                cand_list[(*nb_cand)++] = hmvp_cand;
            }
        }
    }
}

static uint8_t
hmvp_add_merge_cand(const struct HMVPLUT *const hmvp_lut,
                    OVMV *const cand_list,
                    OVMV *const cand_amvp,
                    uint8_t got_B1, uint8_t got_A1,
                    int *const nb_cand, uint8_t inter_dir, uint8_t mvp_idx,
                    uint8_t max_nb_mrg_cand_min1)
{
    int i;
    for (i = 1; i <= hmvp_lut->nb_mv;  ++i) {
        int mv_lut_idx = hmvp_lut->nb_mv - i;
        if (hmvp_lut->dir[mv_lut_idx] & inter_dir) {
            if (i > 2 || ((!got_B1 || !MV_CMP2(hmvp_lut->hmv0[mv_lut_idx], cand_amvp[0])) &&
                          (!got_A1 || !MV_CMP2(hmvp_lut->hmv0[mv_lut_idx], cand_amvp[1])))) {
                cand_list[(*nb_cand)++] = inter_dir & 0x1 ? hmvp_lut->hmv0[mv_lut_idx]
                                                          : hmvp_lut->hmv1[mv_lut_idx];
                if (*nb_cand == mvp_idx + 1) {
                    return 1;
                }

                if (*nb_cand == max_nb_mrg_cand_min1) {
                   return 0;
                }
            }
        }
    }
    return 0;
}

static uint8_t
hmvp_add_merge_cand_b(const struct HMVPLUT *const hmvp_lut,
                      VVCMergeInfo *const cand_list,
                      VVCMergeInfo *const cand_amvp, uint8_t got_B1, uint8_t got_A1,
                      int *const nb_cand, uint8_t mvp_idx, uint8_t max_nb_mrg_min1)
{
    int i;
    for (i = 1; i <= hmvp_lut->nb_mv;  ++i) {
        int mv_lut_idx = hmvp_lut->nb_mv - i;
        VVCMergeInfo lut_mi;
        lut_mi.inter_dir = hmvp_lut->dir[mv_lut_idx];
        lut_mi.mv0 = hmvp_lut->hmv0[mv_lut_idx];
        lut_mi.mv1 = hmvp_lut->hmv1[mv_lut_idx];
        if (i > 2 || ((!got_B1 || !mi_cmp(&lut_mi, &cand_amvp[0])) &&
                      (!got_A1 || !mi_cmp(&lut_mi, &cand_amvp[1])))) {

            cand_list[(*nb_cand)++] = lut_mi;

            if (*nb_cand == mvp_idx + 1) {
                return 1;
            }

            if (*nb_cand == max_nb_mrg_min1) {
                return 0;
            }
        }
    }
    return 0;
}

static void
hmvp_update_lut_b(struct HMVPLUT *const hmvp_lut, OVMV mv0, OVMV mv1, uint8_t inter_dir)
{
    int max_nb_cand = OVMIN(5, hmvp_lut->nb_mv);
    uint8_t duplicated_mv = 0;
    int i;
    uint8_t *const dir_lut = hmvp_lut->dir;

    for (i = 0; i < max_nb_cand;  ++i) {
        if (dir_lut[i] == inter_dir) {
            switch (inter_dir) {
                case 0x1:
                    duplicated_mv = MV_CMP(mv0, hmvp_lut->hmv0[i]) && mv0.ref_idx == hmvp_lut->hmv0[i].ref_idx;
                    break;
                case 0x2:
                    duplicated_mv = MV_CMP(mv1, hmvp_lut->hmv1[i]) && mv1.ref_idx == hmvp_lut->hmv1[i].ref_idx;
                    break;
                case 0x3:
                    duplicated_mv = MV_CMP(mv0, hmvp_lut->hmv0[i]) &&
                        MV_CMP(mv1, hmvp_lut->hmv1[i]) && mv0.ref_idx == hmvp_lut->hmv0[i].ref_idx && mv1.ref_idx == hmvp_lut->hmv1[i].ref_idx;
                    break;
            }
            if (duplicated_mv) {
                break;
            }
        } else duplicated_mv = 0;
    }

    if (duplicated_mv) {
       int j;
       for (j = i; j < max_nb_cand - 1; ++j) {
           hmvp_lut->hmv0[j] = hmvp_lut->hmv0[j + 1];
           hmvp_lut->hmv1[j] = hmvp_lut->hmv1[j + 1];
           hmvp_lut->dir[j] = hmvp_lut->dir[j + 1];
       }
       hmvp_lut->hmv0[j] = mv0;
       hmvp_lut->hmv1[j] = mv1;
       hmvp_lut->dir[j]= inter_dir;
    } else if (hmvp_lut->nb_mv == 5) {
       int j;
       for (j = 1; j < 5; ++j) {
           hmvp_lut->hmv0[j - 1] = hmvp_lut->hmv0[j];
           hmvp_lut->hmv1[j - 1] = hmvp_lut->hmv1[j];
           hmvp_lut->dir[j - 1] = hmvp_lut->dir[j];
       }
       hmvp_lut->hmv0[4] = mv0;
       hmvp_lut->hmv1[4] = mv1;
       hmvp_lut->dir[4] = inter_dir;
    } else {
        hmvp_lut->hmv0[hmvp_lut->nb_mv] = mv0;
        hmvp_lut->hmv1[hmvp_lut->nb_mv] = mv1;
        hmvp_lut->dir[hmvp_lut->nb_mv++] = inter_dir;
    }
}

void
tmvp_inter_synchronization(const OVPicture *ref_pic, int ctb_x, int ctb_y, int log2_ctu_s)
{
    const int pic_w = ref_pic->frame->width;
    const int pic_h = ref_pic->frame->height;
    
    /*Frame thread synchronization to ensure data is available
    */
    int nb_ctb_pic_w = (pic_w + ((1 << log2_ctu_s) - 1)) >> log2_ctu_s;
    int nb_ctb_pic_h = (pic_h + ((1 << log2_ctu_s) - 1)) >> log2_ctu_s;
    int br_ctu_x = OVMIN(ctb_x + 1, nb_ctb_pic_w-1);
    int br_ctu_y = OVMIN(ctb_y, nb_ctb_pic_h-1);
    uint16_t idx = atomic_load(ref_pic->idx_function);
    ref_pic->ovdpb_frame_synchro[idx](ref_pic, ctb_x, ctb_y, br_ctu_x, br_ctu_y);
}


static void
load_ctb_tmvp(OVCTUDec *const ctudec, int ctb_x, int ctb_y)
{

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct VVCTMVP *const tmvp_ctx = &inter_ctx->tmvp_ctx;

    const struct MVPlane *plane0 = tmvp_ctx->col_plane0;
    const struct MVPlane *plane1 = tmvp_ctx->col_plane1;

    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;

    uint8_t nb_unit_ctb = (1 << log2_ctb_s) >> LOG2_MIN_CU_S;
    uint16_t nb_ctb_w = ctudec->nb_ctb_pic_w;

    const OVPicture* ref_pic = tmvp_ctx->col_ref;
    if (ref_pic) {
        tmvp_inter_synchronization(ref_pic, ctb_x, ctb_y, log2_ctb_s);
    }

    uint8_t is_border_pic = nb_ctb_w - 1 == ctb_x;

    if (plane0 || plane1) {
        uint16_t ctb_addr_rs = ctb_x + ctb_y * nb_ctb_w;
        int32_t nb_tmvp_unit = (1 << log2_ctb_s) >> 3;
        int32_t pln_stride = nb_tmvp_unit * nb_ctb_w;
        int32_t ctb_offset = ctb_x * nb_tmvp_unit + (ctb_y * nb_tmvp_unit * pln_stride);

        if (is_border_pic) {
            memset(tmvp_ctx->dir_map_v0, 0, sizeof(uint64_t) * (nb_unit_ctb + 2));
            memset(tmvp_ctx->dir_map_v1, 0, sizeof(uint64_t) * (nb_unit_ctb + 2));
        }

        if (plane0 && plane0->dirs) {
            const uint64_t *src_map = plane0->dirs + ctb_addr_rs * nb_unit_ctb;
                  uint64_t *dst_map = tmvp_ctx->dir_map_v0 + 1;
            const struct TMVPMV *src_mvs = plane0->mvs + ctb_offset;
                  struct TMVPMV *dst_mvs = tmvp_ctx->mvs0;
            int i;

            tmvp_ctx->ctb_mv0 = src_mvs;
            tmvp_ctx->pln0_stride = pln_stride;
            memcpy(dst_map, src_map, sizeof(uint64_t) * (nb_unit_ctb + !is_border_pic));
            #if 0
            for (i = 0; i < nb_tmvp_unit; ++i) {
                memcpy(dst_mvs, src_mvs, sizeof(*dst_mvs) * (nb_tmvp_unit + !is_border_pic));
                dst_mvs += TMVP_BUFF_STRIDE;
                src_mvs += pln_stride;
            }
            #endif
        }

        if (plane1 && plane1->dirs) {
            const uint64_t *src_map = plane1->dirs + ctb_addr_rs * nb_unit_ctb;
                  uint64_t *dst_map = tmvp_ctx->dir_map_v1 + 1;
            const struct TMVPMV *src_mvs = plane1->mvs + ctb_offset;
                  struct TMVPMV *dst_mvs = tmvp_ctx->mvs1;
            int i;

            /*FIXME memory could be spared with smaller map size when possible */
            memcpy(dst_map, src_map, sizeof(uint64_t) * (nb_unit_ctb + !is_border_pic));
            tmvp_ctx->ctb_mv1 = src_mvs;
            tmvp_ctx->pln1_stride = pln_stride;
            tmvp_ctx->pln0_stride = pln_stride;
            #if 0
            for (i = 0; i < nb_tmvp_unit; ++i) {
                memcpy(dst_mvs, src_mvs, sizeof(*dst_mvs) * (nb_tmvp_unit + !is_border_pic));
                dst_mvs += TMVP_BUFF_STRIDE;
                src_mvs += pln_stride;
            }
            #endif
        }
    }

    inter_ctx->tmvp_avail |= 1;
}

static inline struct MV
tmvp_scale_mv(int scale, struct MV mv)
{
    mv.x = ov_clip((scale * mv.x + 128 - (scale * mv.x >= 0)) >> 8, MV_MIN, MV_MAX);
    mv.y = ov_clip((scale * mv.y + 128 - (scale * mv.y >= 0)) >> 8, MV_MIN, MV_MAX);
    return mv;
}

static int16_t
derive_tmvp_scale(int32_t dist_ref, int32_t dist_col)
{
    int32_t scale = 0;

    if (dist_ref == dist_col || !dist_col)
        return 256;

    /*FIXME POW2 clip */
    dist_ref = ov_clip(dist_ref, -128, 127);
    dist_col = ov_clip(dist_col, -128, 127);

    scale = dist_ref * ((0x4000 + OVABS(dist_col >> 1)) / dist_col);
    scale += 32;
    scale >>= 6;
    /* FIXME pow2_clip */
    scale = ov_clip(scale, -4096, 4095);

    return (int16_t)scale;
}

static uint8_t
derive_tmvp_status(const uint64_t *rpl0_vmap, const uint64_t *rpl1_vmap,
                   uint8_t c0_x, uint8_t c0_y,
                   uint8_t c1_x, uint8_t c1_y)
{
    uint64_t c0_vmsk = TMVP_POS_MASK(c0_y);
    uint64_t c1_vmsk = TMVP_POS_MASK(c1_y);

    uint8_t c0_col_idx = c0_x + 1;
    uint8_t c1_col_idx = c1_x + 1;

    uint64_t c0_rpl0 = rpl0_vmap[c0_col_idx] & c0_vmsk;
    uint64_t c0_rpl1 = rpl1_vmap[c0_col_idx] & c0_vmsk;
    uint64_t c1_rpl0 = rpl0_vmap[c1_col_idx] & c1_vmsk;
    uint64_t c1_rpl1 = rpl1_vmap[c1_col_idx] & c1_vmsk;

    uint8_t status = 0;

    status |= !!c0_rpl0;
    status |= !!c0_rpl1 << 1;
    status |= !!c1_rpl0 << 2;
    status |= !!c1_rpl1 << 3;

    return status;
}

static uint8_t
derive_tmvp_cand(const struct InterDRVCtx *const inter_ctx, const struct OVMVCtx *const mv_ctx,
                 OVMV *const cand, uint8_t pb_x, uint8_t pb_y, uint8_t nb_pb_w, uint8_t nb_pb_h,
                 int8_t ref_idx, int8_t opp_ref_idx, uint8_t prec_amvr)
{
    uint32_t msk_8x8 = ~1;

    const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;

    int c1_x = msk_8x8 & (pb_x + (nb_pb_w >> 1));
    int c1_y = msk_8x8 & (pb_y + (nb_pb_h >> 1));
    int c0_x = msk_8x8 & (pb_x + nb_pb_w);
    int c0_y = msk_8x8 & (pb_y + nb_pb_h);

    uint8_t status = derive_tmvp_status(tmvp->dir_map_v0, tmvp->dir_map_v1,
                                        c0_x, c0_y, c1_x, c1_y);

    if (status) {
        uint8_t is_c0 = status & 0x3;

        int pos_in_buff = is_c0 ? TMVP_POS_IN_BUF2(c0_x, c0_y) : TMVP_POS_IN_BUF2(c1_x, c1_y);

        uint8_t is_rpl0 = (mv_ctx == &inter_ctx->mv_ctx0);
        int16_t dist_ref = is_rpl0 ? inter_ctx->dist_ref_0[ref_idx] : inter_ctx->dist_ref_1[ref_idx];

        struct TMVPMV col_mv;
        int16_t scale;

        if (!(tmvp->col_ref_l0 || tmvp->ldc) || (tmvp->ldc && is_rpl0)) {

            if (status & 0x1) {
                col_mv = tmvp->ctb_mv0[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            } else if (status & 0x2) {
                col_mv = tmvp->ctb_mv1[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            } else if (status & 0x4) {
                col_mv = tmvp->ctb_mv0[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            } else if (status & 0x8) {
                col_mv = tmvp->ctb_mv1[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            }
        } else {
            if (status & 0x2) {
                col_mv = tmvp->ctb_mv1[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            } else if (status & 0x1) {
                col_mv = tmvp->ctb_mv0[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            } else if (status & 0x8) {
                col_mv = tmvp->ctb_mv1[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            } else if (status & 0x4) {
                col_mv = tmvp->ctb_mv0[pos_in_buff];
                scale = derive_tmvp_scale(dist_ref, col_mv.z);
            }
        }

        col_mv.mv = tmvp_round_mv(col_mv.mv);
        /* Discard candidate when only one is from long term ref */
        if ((dist_ref == 0) ^ (col_mv.z == 0)) return 0;

        col_mv.mv = tmvp_scale_mv(scale, col_mv.mv);
        col_mv.mv = drv_round_to_precision_mv2(col_mv.mv, MV_PRECISION_INTERNAL, prec_amvr);

        cand[0].x = col_mv.mv.x;
        cand[0].y = col_mv.mv.y;
        cand[0].ref_idx = ref_idx;
        cand[0].bcw_idx_plus1 = 0;
        cand[0].prec_amvr = 0;

        return 1;
    }
    return 0;

}

static OVMV
derive_mvp_candidates_1(struct InterDRVCtx *const inter_ctx,
                        const struct OVMVCtx *const mv_ctx,
                        uint8_t ref_idx,
                        uint8_t pb_x, uint8_t pb_y,
                        uint8_t nb_pb_w, uint8_t nb_pb_h,
                        uint8_t mvp_idx, uint8_t inter_dir,
                        const struct OVMVCtx *const mv_ctx_opp, uint8_t opp_ref_idx,
                        uint8_t prec_amvr, uint8_t is_small)
{
    const OVMV *const mv_buff     = mv_ctx->mvs;
    const OVMV *const mv_buff_opp = mv_ctx_opp->mvs;
    uint64_t lft_col = mv_ctx->map.vfield[pb_x];
    uint64_t abv_row = mv_ctx->map.hfield[pb_y];

    uint64_t lft_col0 = mv_ctx_opp->map.vfield[pb_x];
    uint64_t abv_row0 = mv_ctx_opp->map.hfield[pb_y];

    OVMV cand[2] = {0};
    int nb_cand = 0;

    /* Derive candidates availability based on CTU inter fields */
    uint8_t cand_bl = !!(lft_col & POS_MASK(pb_y, nb_pb_h));     /*AO*/
    uint8_t cand_l  = !!(lft_col & POS_MASK(pb_y, nb_pb_h - 1)); /*A1*/
    uint8_t cand_tr = !!(abv_row & POS_MASK(pb_x, nb_pb_w));     /*B0*/
    uint8_t cand_t  = !!(abv_row & POS_MASK(pb_x, nb_pb_w - 1)); /*B1*/
    uint8_t cand_tl = !!(abv_row & POS_MASK(pb_x - 1, 0));      /*B2*/

    uint8_t cand_bl0 = !!(lft_col0 & POS_MASK(pb_y, nb_pb_h));     /*AO*/
    uint8_t cand_l0  = !!(lft_col0 & POS_MASK(pb_y, nb_pb_h - 1)); /*A1*/
    uint8_t cand_tr0 = !!(abv_row0 & POS_MASK(pb_x, nb_pb_w));     /*B0*/
    uint8_t cand_t0  = !!(abv_row0 & POS_MASK(pb_x, nb_pb_w - 1)); /*B1*/
    uint8_t cand_tl0 = !!(abv_row0 & POS_MASK(pb_x - 1, 0));      /*B2*/
    uint8_t found = 0;

    if (cand_bl | cand_bl0) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h);
        if (cand_bl && mv_buff[pos_in_buff].ref_idx == ref_idx) {
            cand[nb_cand++] = mv_buff[pos_in_buff];
            found = 1;
        } else if (cand_bl0 && mv_buff_opp[pos_in_buff].ref_idx == opp_ref_idx) {
            cand[nb_cand++] = mv_buff_opp[pos_in_buff];
            found = 1;
        }
    }

    if (!found && (cand_l | cand_l0)) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h - 1);
        if (cand_l && mv_buff[pos_in_buff].ref_idx == ref_idx) {
            cand[nb_cand++] = mv_buff[pos_in_buff];
        } else if (cand_l0 && mv_buff_opp[pos_in_buff].ref_idx == opp_ref_idx) {
            cand[nb_cand++] = mv_buff_opp[pos_in_buff];
        }
    }

    found = 0;

    if (!found && (cand_tr | cand_tr0)) {
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w, pb_y - 1);
        if (cand_tr && mv_buff[pos_in_buff].ref_idx == ref_idx) {
            cand[nb_cand++] = mv_buff[pos_in_buff];
            found = 1;
        } else if (cand_tr0 && mv_buff_opp[pos_in_buff].ref_idx == opp_ref_idx) {
            cand[nb_cand++] = mv_buff_opp[pos_in_buff];
            found = 1;
        }
    }

    if (!found && (cand_t | cand_t0)) {
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w - 1, pb_y - 1);
        if (cand_t && mv_buff[pos_in_buff].ref_idx == ref_idx) {
            cand[nb_cand++] = mv_buff[pos_in_buff];
            found = 1;
        } else if (cand_t0 && mv_buff_opp[pos_in_buff].ref_idx == opp_ref_idx) {
            cand[nb_cand++] = mv_buff_opp[pos_in_buff];
            found = 1;
        }
    }

    if (!found && (cand_tl | cand_tl0)) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y - 1);
        if (cand_tl && mv_buff[pos_in_buff].ref_idx == ref_idx) {
            cand[nb_cand++] = mv_buff[pos_in_buff];
            found = 1;
        } else if (cand_tl0 && mv_buff_opp[pos_in_buff].ref_idx == opp_ref_idx) {
            cand[nb_cand++] = mv_buff_opp[pos_in_buff];
            found = 1;
        }
    }

    cand[0]  = drv_round_to_precision_mv(cand[0], MV_PRECISION_INTERNAL, prec_amvr);
    cand[1]  = drv_round_to_precision_mv(cand[1], MV_PRECISION_INTERNAL, prec_amvr);

    /* Remove on candidates if duplicated */
    if (nb_cand == 2) {
        if (MV_CMP(cand[0], cand[1])) {
            --nb_cand;
        }
    }

    if (inter_ctx->tmvp_enabled && nb_cand < 2 && !is_small) {
        if (!inter_ctx->tmvp_avail) {
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }

        nb_cand += derive_tmvp_cand(inter_ctx, mv_ctx, &cand[nb_cand], pb_x, pb_y, nb_pb_w, nb_pb_h,
                                    ref_idx, opp_ref_idx, prec_amvr);
    }

    if (nb_cand < 2) {
        const struct HMVPLUT *hmvp_lut = &inter_ctx->hmvp_lut;
        hmvp_add_cand_1(hmvp_lut, cand, &nb_cand, inter_dir, ref_idx, opp_ref_idx);
    }

    while (nb_cand < 2) {
        OVMV zmv = {0};
        zmv.ref_idx = ref_idx;
        cand[nb_cand++] = zmv;
    }

    cand[0]  = drv_round_to_precision_mv(cand[0], MV_PRECISION_INTERNAL, prec_amvr);
    cand[1]  = drv_round_to_precision_mv(cand[1], MV_PRECISION_INTERNAL, prec_amvr);

    return cand[mvp_idx];
}

static uint8_t
derive_tmvp_merge_cand(const struct InterDRVCtx *const inter_ctx,
                        uint8_t pb_x, uint8_t pb_y,
                        uint8_t nb_pb_w, uint8_t nb_pb_h,
                        OVMV *cand)
{
    const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;

    uint32_t msk_8x8 = ~1;

    int c1_x = msk_8x8 & (pb_x + (nb_pb_w >> 1));
    int c1_y = msk_8x8 & (pb_y + (nb_pb_h >> 1));
    int c0_x = msk_8x8 & (pb_x + nb_pb_w);
    int c0_y = msk_8x8 & (pb_y + nb_pb_h);

    uint8_t status = derive_tmvp_status(tmvp->dir_map_v0, tmvp->dir_map_v1,
                                        c0_x, c0_y, c1_x, c1_y);

    if (status) {
        int8_t dist_ref0 = inter_ctx->dist_ref_0[0];
        struct TMVPMV col_mv;
        int16_t scale;

        if (status & 0x1) {
            int pos_in_buff = TMVP_POS_IN_BUF2(c0_x, c0_y);
            col_mv = tmvp->ctb_mv0[pos_in_buff];
            scale = derive_tmvp_scale(dist_ref0, col_mv.z);
        } else if (status & 0x2) {
            int pos_in_buff = TMVP_POS_IN_BUF2(c0_x, c0_y);
            col_mv = tmvp->ctb_mv1[pos_in_buff];
            scale = derive_tmvp_scale(dist_ref0, col_mv.z);
        } else if (status & 0x4) {
            int pos_in_buff = TMVP_POS_IN_BUF2(c1_x, c1_y);
            col_mv = tmvp->ctb_mv0[pos_in_buff];
            scale = derive_tmvp_scale(dist_ref0, col_mv.z);
        } else if (status & 0x8) {
            int pos_in_buff = TMVP_POS_IN_BUF2(c1_x, c1_y);
            col_mv = tmvp->ctb_mv1[pos_in_buff];
            scale = derive_tmvp_scale(dist_ref0, col_mv.z);
        }
found:
        /* Discard candidate when only one is from long term ref */
        if ((dist_ref0 == 0) ^ (col_mv.z == 0)) return 0;
        col_mv.mv = tmvp_round_mv(col_mv.mv);
        col_mv.mv = tmvp_scale_mv(scale , col_mv.mv);

        cand->x = col_mv.mv.x;
        cand->y = col_mv.mv.y;
        cand->ref_idx = 0;
        cand->bcw_idx_plus1 = 0;
        cand->prec_amvr = 0;

        return 1;
    }
    return 0;
}

static OVMV
vvc_derive_merge_mvp(const struct InterDRVCtx *const inter_ctx,
                     const struct OVMVCtx *const mv_ctx,
                     uint8_t pb_x, uint8_t pb_y,
                     uint8_t nb_pb_w, uint8_t nb_pb_h,
                     uint8_t merge_idx, uint8_t max_nb_merge_cand,
                     uint8_t is_small)
{
    const OVMV *const mv_buff = mv_ctx->mvs;

    uint64_t lft_col = mv_ctx->map.vfield[pb_x];
    uint64_t abv_row = mv_ctx->map.hfield[pb_y];

    uint8_t cand_bl = !!(lft_col & POS_MASK(pb_y, nb_pb_h));     /*AO*/
    uint8_t cand_l  = !!(lft_col & POS_MASK(pb_y, nb_pb_h - 1)); /*A1*/
    uint8_t cand_tr = !!(abv_row & POS_MASK(pb_x, nb_pb_w));     /*B0*/
    uint8_t cand_t  = !!(abv_row & POS_MASK(pb_x, nb_pb_w - 1)); /*B1*/
    uint8_t cand_tl = !!(abv_row & POS_MASK(pb_x - 1, 0));      /*B2*/
    OVMV cand[6];
    OVMV cand_amvp[5];
    OVMV mv_z = {-1};

    int nb_cand = 0;

    cand_amvp[0] = mv_z;
    if (cand_t) { /* B1 */
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w - 1, pb_y - 1);
        OVMV mv_B1 = mv_buff[pos_in_buff];
        cand_amvp[0] = mv_buff[pos_in_buff];
        cand[nb_cand] = mv_B1;
        if (nb_cand++ == merge_idx)
            return mv_B1;
    }

    cand_amvp[1] = mv_z;
    if (cand_l) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h - 1);
        OVMV mv_A1 = mv_buff[pos_in_buff];
        cand_amvp[1] = mv_buff[pos_in_buff];
        if (!cand_t || !MV_CMP2(mv_A1, cand_amvp[0])) {
            cand[nb_cand] = mv_A1;
            if (nb_cand++ == merge_idx)
                return mv_A1;
        }
    }

    cand_amvp[2] = mv_z;
    if (cand_tr) {
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w, pb_y - 1);
        OVMV mv_B0 = mv_buff[pos_in_buff];
        cand_amvp[2] = mv_buff[pos_in_buff];
        if (!cand_t || !MV_CMP2(mv_B0, cand_amvp[0])) {
            cand[nb_cand] = mv_B0;
            if (nb_cand++ == merge_idx)
                return mv_B0;
        }
    }

    cand_amvp[3] = mv_z;
    if (cand_bl) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h);
        OVMV mv_A0 = mv_buff[pos_in_buff];
        cand_amvp[3] = mv_buff[pos_in_buff];
        if (!cand_l || !MV_CMP2(mv_A0, cand_amvp[1])) {
            cand[nb_cand] = mv_A0;
            if (nb_cand++ == merge_idx)
                return mv_A0;
        }
    }

    cand_amvp[4] = mv_z;
    if (nb_cand < 4) {
        if (cand_tl) {
            int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y - 1);
            OVMV mv_B2   = mv_buff[pos_in_buff];
            cand_amvp[4] = mv_buff[pos_in_buff];
            if ((!cand_l || !MV_CMP2(mv_B2, cand_amvp[1])) &&
                (!cand_t || !MV_CMP2(mv_B2, cand_amvp[0]))) {
                cand[nb_cand] = mv_B2;
                if (nb_cand++ == merge_idx)
                    return mv_B2;
            }
        }
    }

    if (inter_ctx->tmvp_enabled && !is_small) {

        if (!inter_ctx->tmvp_avail) {
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }

        nb_cand += derive_tmvp_merge_cand(inter_ctx, pb_x, pb_y, nb_pb_w, nb_pb_h, &cand[nb_cand]);
        if (nb_cand - 1 == merge_idx) {
            return cand[nb_cand - 1];
        }

    }

    if (nb_cand != max_nb_merge_cand - 1) {
        uint8_t found;
        const struct HMVPLUT *hmvp_lut = &inter_ctx->hmvp_lut;

        found = hmvp_add_merge_cand(hmvp_lut, cand, cand_amvp, cand_t, cand_l, &nb_cand, 1, merge_idx, max_nb_merge_cand - 1);
        if (found) {
            return cand[nb_cand - 1];
        }
    }


    if (nb_cand > 1 && nb_cand < max_nb_merge_cand) {

        OVMV avg_mv = cand[0];
        avg_mv.x += cand[1].x;
        avg_mv.y += cand[1].y;
        avg_mv.x += 1 - (avg_mv.x >= 0);
        avg_mv.y += 1 - (avg_mv.y >= 0);
        avg_mv.x >>= 1;
        avg_mv.y >>= 1;

        if (cand[0].prec_amvr != cand[1].prec_amvr) {
            avg_mv.prec_amvr = 0;
        }

        if (nb_cand == merge_idx)
            return avg_mv;

        nb_cand++;
    }

    int8_t diff_cand = merge_idx - nb_cand;
    int8_t num_min_ref = inter_ctx->nb_active_ref0;
    int8_t ref_idx = 0;
    if (diff_cand <= num_min_ref - 1) ref_idx = diff_cand;

    /*FIXME ref_idx*/
    while (nb_cand < max_nb_merge_cand) {
        OVMV zmv = {0};
        zmv.ref_idx = ref_idx;
        cand[nb_cand++] = zmv;
    }

    return cand[nb_cand - 1];
}

static uint8_t
derive_merge_enable_cand_list(uint8_t x0_u, uint8_t y0_u, uint8_t nb_units_w, uint8_t nb_units_h,
                              uint8_t log2_pmerge_lvl)
{
    uint8_t enable_msk = 0;
    int16_t x0 = x0_u << LOG2_MIN_CU_S;
    int16_t y0 = y0_u << LOG2_MIN_CU_S;
    int16_t cb_w = nb_units_w << LOG2_MIN_CU_S;
    int16_t cb_h = nb_units_h << LOG2_MIN_CU_S;

    uint8_t cand_a0 = ((y0 + cb_h)     >> log2_pmerge_lvl) != (y0 >> log2_pmerge_lvl);
    uint8_t cand_a1 = ((y0 + cb_h - 1) >> log2_pmerge_lvl) != (y0 >> log2_pmerge_lvl);

    uint8_t cand_b0 = ((x0 + cb_w)     >> log2_pmerge_lvl) != (x0 >> log2_pmerge_lvl);
    uint8_t cand_b1 = ((x0 + cb_w - 1) >> log2_pmerge_lvl) != (x0 >> log2_pmerge_lvl);
    uint8_t cand_b2 = ((x0 - 1)        >> log2_pmerge_lvl) != (x0 >> log2_pmerge_lvl);

    cand_a0 |= ((x0 - 1) >> log2_pmerge_lvl) != (x0 >> log2_pmerge_lvl);
    cand_a1 |= ((x0 - 1) >> log2_pmerge_lvl) != (x0 >> log2_pmerge_lvl);

    cand_b0 |= ((y0 - 1) >> log2_pmerge_lvl) != (y0 >> log2_pmerge_lvl);
    cand_b1 |= ((y0 - 1) >> log2_pmerge_lvl) != (y0 >> log2_pmerge_lvl);
    cand_b2 |= ((y0 - 1) >> log2_pmerge_lvl) != (y0 >> log2_pmerge_lvl);

    enable_msk |= cand_b2;
    enable_msk <<= 1;
    enable_msk |= cand_b1;
    enable_msk <<= 1;
    enable_msk |= cand_b0;
    enable_msk <<= 1;
    enable_msk |= cand_a1;
    enable_msk <<= 1;
    enable_msk |= cand_a0;

    return enable_msk;
}

static uint8_t
derive_tmvp_merge(const struct InterDRVCtx *const inter_ctx,
                  uint8_t pb_x, uint8_t pb_y,
                  uint8_t nb_pb_w, uint8_t nb_pb_h,
                  VVCMergeInfo *cand)
{
    const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;

    uint32_t msk_8x8 = ~1;

    int c1_x = msk_8x8 & (pb_x + (nb_pb_w >> 1));
    int c1_y = msk_8x8 & (pb_y + (nb_pb_h >> 1));
    int c0_x = msk_8x8 & (pb_x + nb_pb_w);
    int c0_y = msk_8x8 & (pb_y + nb_pb_h);

    uint8_t status = derive_tmvp_status(tmvp->dir_map_v0, tmvp->dir_map_v1,
                                        c0_x, c0_y, c1_x, c1_y);

    if (status) {
        uint8_t is_c0 = status & 0x3;
        int pos_in_buff = is_c0 ? TMVP_POS_IN_BUF2(c0_x, c0_y) : TMVP_POS_IN_BUF2(c1_x, c1_y);

        uint8_t cand_c0_cur = !tmvp->col_ref_l0 ? status & 0x1 : status & 0x2;
        uint8_t cand_c0_opp = !tmvp->col_ref_l0 ? status & 0x2 : status & 0x1 ;
        uint8_t cand_c1_cur = !tmvp->col_ref_l0 ? status & 0x4 : status & 0x8;
        uint8_t cand_c1_opp = !tmvp->col_ref_l0 ? status & 0x8 : status & 0x4 ;

        uint8_t cand_cur = is_c0 ? cand_c0_cur : cand_c1_cur;
        uint8_t cand_opp = is_c0 ? cand_c0_opp : cand_c1_opp;

        /* FIXME most of those could be swapped at slice init */
        const struct TMVPMV *mv_cur = !tmvp->col_ref_l0 ? tmvp->ctb_mv0 : tmvp->ctb_mv1;
        const struct TMVPMV *mv_opp = !tmvp->col_ref_l0 ? tmvp->ctb_mv1 : tmvp->ctb_mv0;

        int8_t dist_ref_cur = !tmvp->col_ref_l0 ? inter_ctx->dist_ref_0[0] : inter_ctx->dist_ref_1[0];
        int8_t dist_ref_opp = !tmvp->col_ref_l0 ? inter_ctx->dist_ref_1[0] : inter_ctx->dist_ref_0[0];

        OVMV *dst_cur = !tmvp->col_ref_l0 ? &cand->mv0 : &cand->mv1;
        OVMV *dst_opp = !tmvp->col_ref_l0 ? &cand->mv1 : &cand->mv0;
        int16_t scale_cur;
        int16_t scale_opp;

        struct TMVPMV tmp_cur;
        struct TMVPMV tmp_opp;
        struct TMVPMV col_mv_opp;
        struct TMVPMV col_mv;

        if (cand_cur) {
            col_mv  = mv_cur[pos_in_buff];
            if (cand_opp && tmvp->ldc) {
                col_mv_opp  = mv_opp[pos_in_buff];
            } else {
                col_mv_opp  = col_mv;
            }
            goto found;
        } else if (cand_opp) {
            col_mv_opp  = mv_opp[pos_in_buff];
            col_mv = col_mv_opp;
            goto found;
        }
found:
        scale_cur = derive_tmvp_scale(dist_ref_cur, col_mv.z);
        scale_opp = derive_tmvp_scale(dist_ref_opp, col_mv_opp.z);

        /* Discard candidate when only one is from long term ref */
        if ((dist_ref_cur == 0) ^ (col_mv.z == 0)) {
            if ((dist_ref_opp == 0) ^ (col_mv_opp.z == 0))
                return 0;
            col_mv = col_mv_opp;
            scale_cur = scale_opp;
        } else if ((dist_ref_opp == 0) ^ (col_mv_opp.z == 0)) {
            col_mv_opp = col_mv;
            scale_opp = scale_cur;
        }

        tmp_cur.mv = tmvp_round_mv(col_mv.mv);
        tmp_opp.mv = tmvp_round_mv(col_mv_opp.mv);

        tmp_cur.mv = tmvp_scale_mv(scale_cur, tmp_cur.mv);
        tmp_opp.mv = tmvp_scale_mv(scale_opp, tmp_opp.mv);

        dst_cur->x = tmp_cur.mv.x;
        dst_cur->y = tmp_cur.mv.y;
        dst_opp->x = tmp_opp.mv.x;
        dst_opp->y = tmp_opp.mv.y;

        cand->inter_dir = 3;

        cand->mv0.ref_idx = 0;
        cand->mv0.bcw_idx_plus1 = 0;
        cand->mv0.prec_amvr = 0;

        cand->mv1.ref_idx = 0;
        cand->mv1.bcw_idx_plus1 = 0;
        cand->mv1.prec_amvr = 0;

        return 1;
    }

    return 0;
}

static VVCMergeInfo
vvc_derive_merge_mvp_b(const struct InterDRVCtx *const inter_ctx,
                       uint8_t pb_x, uint8_t pb_y,
                       uint8_t nb_pb_w, uint8_t nb_pb_h,
                       uint8_t merge_idx, uint8_t max_nb_merge_cand,
                       uint8_t is_small)
{
    const struct OVMVCtx *mv_ctx0 = &inter_ctx->mv_ctx0;
    const struct OVMVCtx *mv_ctx1 = &inter_ctx->mv_ctx1;

    const OVMV *mv_buff0 = mv_ctx0->mvs;
    const OVMV *mv_buff1 = mv_ctx1->mvs;

    uint64_t lft_col0 = mv_ctx0->map.vfield[pb_x];
    uint64_t abv_row0 = mv_ctx0->map.hfield[pb_y];

    uint64_t lft_col1 = mv_ctx1->map.vfield[pb_x];
    uint64_t abv_row1 = mv_ctx1->map.hfield[pb_y];
    uint8_t log2_pmerge_lvl = inter_ctx->log2_parallel_merge_level;
    uint8_t enable_msk = derive_merge_enable_cand_list(pb_x, pb_y, nb_pb_w, nb_pb_h, log2_pmerge_lvl);

    /*FIXME use flags for inter_dir and availability*/
    uint8_t cand_bl0 = !!(lft_col0 & POS_MASK(pb_y, nb_pb_h))      && (enable_msk & (1 << 0));     /*AO*/
    uint8_t cand_bl1 = !!(lft_col1 & POS_MASK(pb_y, nb_pb_h))      && (enable_msk & (1 << 0));     /*AO*/
    uint8_t cand_l0  = !!(lft_col0 & POS_MASK(pb_y, nb_pb_h - 1))      && (enable_msk & (1 << 1)); /*A1*/
    uint8_t cand_l1  = !!(lft_col1 & POS_MASK(pb_y, nb_pb_h - 1))      && (enable_msk & (1 << 1)); /*A1*/

    uint8_t cand_tr0 = !!(abv_row0 & POS_MASK(pb_x, nb_pb_w))      && (enable_msk & (1 << 2));     /*B0*/
    uint8_t cand_tr1 = !!(abv_row1 & POS_MASK(pb_x, nb_pb_w))      && (enable_msk & (1 << 2));     /*B0*/
    uint8_t cand_t0  = !!(abv_row0 & POS_MASK(pb_x, nb_pb_w - 1))      && (enable_msk & (1 << 3)); /*B1*/
    uint8_t cand_t1  = !!(abv_row1 & POS_MASK(pb_x, nb_pb_w - 1))      && (enable_msk & (1 << 3)); /*B1*/
    uint8_t cand_tl0 = !!(abv_row0 & POS_MASK(pb_x - 1, 0))      && (enable_msk & (1 << 4));      /*B2*/
    uint8_t cand_tl1 = !!(abv_row1 & POS_MASK(pb_x - 1, 0))      && (enable_msk & (1 << 4));      /*B2*/

    VVCMergeInfo cand[6];
    VVCMergeInfo cand_amvp[5];

    VVCMergeInfo mv_z = { .mv0 = {0}, .mv1 = {0}, .inter_dir = 3};

    int nb_cand = 0;

    cand_amvp[0] = mv_z;
    if (cand_t0 | cand_t1) { /* B1 */
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w - 1, pb_y - 1);
        VVCMergeInfo mi_B1;
        mi_B1.mv0 = mv_buff0[pos_in_buff];
        mi_B1.mv1 = mv_buff1[pos_in_buff];
        mi_B1.inter_dir =  cand_t0 | (cand_t1 << 1);

        cand_amvp[0]  = mi_B1;

        cand[nb_cand] = mi_B1;
        if (nb_cand++ == merge_idx)
            return mi_B1;
    }

    cand_amvp[1] = mv_z;
    if (cand_l0 | cand_l1) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h - 1);
        VVCMergeInfo mv_A1;
        mv_A1.mv0 = mv_buff0[pos_in_buff];
        mv_A1.mv1 = mv_buff1[pos_in_buff];
        mv_A1.inter_dir = cand_l0 | (cand_l1 << 1);

        cand_amvp[1] = mv_A1;
        if (!(cand_t0 | cand_t1) || !mi_cmp(&mv_A1, &cand_amvp[0])) {
            cand[nb_cand] = mv_A1;
            if (nb_cand++ == merge_idx)
                return mv_A1;
        }
    }

    if (cand_tr0 | cand_tr1) {
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w, pb_y - 1);
        VVCMergeInfo mv_B0;
        mv_B0.mv0 = mv_buff0[pos_in_buff];
        mv_B0.mv1 = mv_buff1[pos_in_buff];
        mv_B0.inter_dir =  cand_tr0 | (cand_tr1 << 1);

        cand_amvp[2] = mv_B0;
        if (!(cand_t0 | cand_t1) || !mi_cmp(&mv_B0, &cand_amvp[0])) {
            cand[nb_cand] = mv_B0;
            if (nb_cand++ == merge_idx)
                return mv_B0;
        }
    }

    if (cand_bl0 | cand_bl1) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h);
        VVCMergeInfo mv_A0;
        mv_A0.mv0 = mv_buff0[pos_in_buff];
        mv_A0.mv1 = mv_buff1[pos_in_buff];
        mv_A0.inter_dir = cand_bl0 | (cand_bl1 << 1);

        cand_amvp[3] = mv_A0;
        if (!(cand_l0 | cand_l1) || !mi_cmp(&mv_A0, &cand_amvp[1])) {
            cand[nb_cand] = mv_A0;
            if (nb_cand++ == merge_idx)
                return mv_A0;
        }
    }

    if (nb_cand < 4) {
        if (cand_tl0 | cand_tl1) {
            int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y - 1);
            VVCMergeInfo mv_B2;
            mv_B2.mv0 = mv_buff0[pos_in_buff];
            mv_B2.mv1 = mv_buff1[pos_in_buff];
            mv_B2.inter_dir =  cand_tl0 | (cand_tl1 << 1);

            cand_amvp[4] = mv_B2;
            if ((!(cand_l0 | cand_l1) || !mi_cmp(&mv_B2, &cand_amvp[1])) &&
                (!(cand_t0 | cand_t1) || !mi_cmp(&mv_B2, &cand_amvp[0]))) {
                cand[nb_cand] = mv_B2;
                if (nb_cand++ == merge_idx)
                    return mv_B2;
            }
        }
    }

    if (inter_ctx->tmvp_enabled && !is_small) {
        if (!inter_ctx->tmvp_avail) {
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }

        nb_cand += derive_tmvp_merge(inter_ctx, pb_x, pb_y, nb_pb_w, nb_pb_h, &cand[nb_cand]);
        if (nb_cand - 1 == merge_idx) {
            return cand[nb_cand - 1];
        }
    }

    if (nb_cand != max_nb_merge_cand - 1) {
        const struct HMVPLUT *const hmvp_lut = &inter_ctx->hmvp_lut;
        uint8_t found;

        found = hmvp_add_merge_cand_b(hmvp_lut, cand, cand_amvp, cand_t0 | cand_t1, cand_l0 | cand_l1, &nb_cand, merge_idx, max_nb_merge_cand - 1);
        if (found) {
            return cand[nb_cand - 1];
        }
    }


    if (nb_cand > 1 && nb_cand < max_nb_merge_cand) {
        VVCMergeInfo avg_mv = cand[0];

        avg_mv.inter_dir = cand[0].inter_dir & cand[1].inter_dir;

        if (avg_mv.inter_dir & 0x1) {
            avg_mv.mv0.x += cand[1].mv0.x;
            avg_mv.mv0.y += cand[1].mv0.y;
            avg_mv.mv0.x += 1 - (avg_mv.mv0.x >= 0);
            avg_mv.mv0.y += 1 - (avg_mv.mv0.y >= 0);
            avg_mv.mv0.x >>= 1;
            avg_mv.mv0.y >>= 1;
        } else if (cand[1].inter_dir & 0x1) {
            avg_mv.mv0 = cand[1].mv0;
            avg_mv.inter_dir |= 1;
            avg_mv.mv0.ref_idx = cand[1].mv0.ref_idx;
        } else if (cand[0].inter_dir & 0x1) {
            avg_mv.inter_dir |= 1;
        }

        if (avg_mv.inter_dir & 0x2) {
            avg_mv.mv1.x += cand[1].mv1.x;
            avg_mv.mv1.y += cand[1].mv1.y;
            avg_mv.mv1.x += 1 - (avg_mv.mv1.x >= 0);
            avg_mv.mv1.y += 1 - (avg_mv.mv1.y >= 0);
            avg_mv.mv1.x >>= 1;
            avg_mv.mv1.y >>= 1;
        } else if (cand[1].inter_dir & 0x2) {
            avg_mv.mv1 = cand[1].mv1;
            avg_mv.inter_dir |= 2;
            avg_mv.mv1.ref_idx = cand[1].mv1.ref_idx;
        } else if (cand[0].inter_dir & 0x2) {
            avg_mv.inter_dir |= 2;
        }

        if (nb_cand == merge_idx){
            uint8_t prec_amvr0 = cand[0].inter_dir & 0x1 ? cand[0].mv0.prec_amvr : cand[0].mv1.prec_amvr;
            uint8_t prec_amvr1 = cand[1].inter_dir & 0x1 ? cand[1].mv0.prec_amvr : cand[1].mv1.prec_amvr;
            avg_mv.mv0.prec_amvr = (prec_amvr0 == prec_amvr1) ? prec_amvr0 : 0 ;
            avg_mv.mv1.prec_amvr = avg_mv.mv0.prec_amvr ;
            avg_mv.mv0.bcw_idx_plus1 = 0;
            avg_mv.mv1.bcw_idx_plus1 = 0;
            return avg_mv;
        }

        nb_cand++;
    }
    /*FIXME saturation or ridx based on nb_active_ref*/
    int8_t diff_cand = merge_idx - nb_cand;
    int8_t num_min_ref = OVMIN(inter_ctx->nb_active_ref0, inter_ctx->nb_active_ref1);
    int8_t ridx = 0;
    if (diff_cand <= num_min_ref - 1)
        ridx = diff_cand;

    mv_z.inter_dir = 3;
    mv_z.mv0.ref_idx = ridx;
    mv_z.mv0.bcw_idx_plus1 = 0;
    mv_z.mv0.prec_amvr = 0;
    mv_z.mv1.ref_idx = ridx;
    mv_z.mv1.bcw_idx_plus1 = 0;
    mv_z.mv1.prec_amvr = 0;

    return mv_z;
}

static void
fill_mvp_map(struct OVMVCtx *const mv_ctx, OVMV mv,
             int pb_x, int pb_y, int nb_pb_w, int nb_pb_h)
{
    int i, j;
    OVMV *const mv_start = &mv_ctx->mvs[PB_POS_IN_BUF(pb_x + 0, pb_y + 0)];
    OVMV *mv_l = mv_start;

    ctu_field_set_rect_bitfield(&mv_ctx->map, pb_x, pb_y, nb_pb_w, nb_pb_h);

    for (j = 0; j < nb_pb_h; ++j) {
        mv_l[0          ] = mv;
        mv_l[nb_pb_w - 1] = mv;
        mv_l += 34;
    }
    mv_l = mv_start + 34 * (nb_pb_h - 1);

    for (i = 0; i < nb_pb_w; ++i) {
        mv_start[i] = mv;
        mv_l[i]     = mv;
    }
}

static void
fill_tmvp_map(struct TMVPMV *const tmvp_mv, struct TMVPMV mv,
              int pb_x, int pb_y, int nb_pb_w, int nb_pb_h)
{
    int i, j;
    struct TMVPMV *dst_mv = tmvp_mv;
    uint8_t skip_first_x = pb_x & 0x1;
    uint8_t skip_first_y = pb_y & 0x1;

    /* Align MVs on 8x8 grid */
    struct TMVPMV *mv_line = &dst_mv[((pb_x + skip_first_x) >> 1) + ((pb_y + skip_first_y) >> 1) * 16];

    for (j = 0; j < (nb_pb_h + !skip_first_y >> 1); ++j) {
        struct TMVPMV *dst = mv_line;
        for (i = 0; i < (nb_pb_w + !skip_first_x >> 1); ++i) {
            *dst = mv;
            ++dst;
        }
        mv_line += 16;
    }
}

static void
update_gpm_mv_ctx_b(struct InterDRVCtx *const inter_ctx,
                    const OVMV mv0, const OVMV mv1,
                    uint8_t pb_x, uint8_t  pb_y,
                    uint8_t nb_pb_w, uint8_t nb_pb_h,
                    uint8_t inter_dir)
{
    /*FIXME Use specific DBF update function if DBF is disabled */
    /*FIXME Find a better way to retrieve dbf_info */
    if (inter_dir == 3) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv0 = {.mv.x=mv0.x, .mv.y=mv0.y, .z=inter_ctx->dist_ref_0[mv0.ref_idx]};
        struct TMVPMV tmvpmv1 = {.mv.x=mv1.x, .mv.y=mv1.y, .z=inter_ctx->dist_ref_1[mv1.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[0].mvs, tmvpmv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
        fill_tmvp_map(inter_ctx->tmvp_mv[1].mvs, tmvpmv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
        fill_mvp_map(mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x2) {
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv1 = {.mv.x=mv1.x, .mv.y=mv1.y, .z=inter_ctx->dist_ref_1[mv1.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[1].mvs, tmvpmv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct TMVPMV tmvpmv0 = {.mv.x=mv0.x, .mv.y=mv0.y, .z=inter_ctx->dist_ref_0[mv0.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[0].mvs, tmvpmv0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
    }
}

static void
update_mv_ctx_b(struct InterDRVCtx *const inter_ctx,
                const OVMV mv0, const OVMV mv1,
                uint8_t pb_x, uint8_t  pb_y,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t inter_dir)
{
    if (inter_dir == 3) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv0 = {.mv.x=mv0.x, .mv.y=mv0.y, .z=inter_ctx->dist_ref_0[mv0.ref_idx]};
        struct TMVPMV tmvpmv1 = {.mv.x=mv1.x, .mv.y=mv1.y, .z=inter_ctx->dist_ref_1[mv1.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[0].mvs, tmvpmv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
        fill_tmvp_map(inter_ctx->tmvp_mv[1].mvs, tmvpmv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
        fill_mvp_map(mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x2) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv1 = {.mv.x=mv1.x, .mv.y=mv1.y, .z=inter_ctx->dist_ref_1[mv1.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[1].mvs, tmvpmv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv0 = {.mv.x=mv0.x, .mv.y=mv0.y, .z=inter_ctx->dist_ref_0[mv0.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[0].mvs, tmvpmv0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
    }
}

static uint8_t
enable_hmvp_storage(uint8_t x0_u, uint8_t y0_u, uint8_t nb_unit_w, uint8_t nb_unit_h, uint8_t log2_pmerge_lvl)
{
    int16_t x0 = x0_u << LOG2_MIN_CU_S;
    int16_t y0 = y0_u << LOG2_MIN_CU_S;
    int16_t cb_w = nb_unit_w << LOG2_MIN_CU_S;
    int16_t cb_h = nb_unit_h << LOG2_MIN_CU_S;
    uint8_t enable_hmvp  = ((x0 + cb_w) >> log2_pmerge_lvl) > (x0 >> log2_pmerge_lvl);
            enable_hmvp &= ((y0 + cb_h) >> log2_pmerge_lvl) > (y0 >> log2_pmerge_lvl);
    
    return enable_hmvp;
    //return 1;
}

static void
update_mv_ctx(struct InterDRVCtx *const inter_ctx,
              const OVMV mv,
              uint8_t pb_x, uint8_t  pb_y,
              uint8_t nb_pb_w, uint8_t nb_pb_h,
              uint8_t inter_dir)
{
    uint8_t log2_pmerge_lvl = inter_ctx->log2_parallel_merge_level;
    uint8_t enable_hmvp = enable_hmvp_storage(pb_x, pb_y, nb_pb_w, nb_pb_h, log2_pmerge_lvl);
    if (inter_dir & 0x2) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv = {.mv.x=mv.x, .mv.y=mv.y, .z=inter_ctx->dist_ref_1[mv.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[1].mvs, tmvpmv, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx1, mv, pb_x, pb_y, nb_pb_w, nb_pb_h);

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct TMVPMV tmvpmv = {.mv.x=mv.x, .mv.y=mv.y, .z=inter_ctx->dist_ref_0[mv.ref_idx]};

        fill_tmvp_map(inter_ctx->tmvp_mv[0].mvs, tmvpmv, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx0, mv, pb_x, pb_y, nb_pb_w, nb_pb_h);
    }

    if (enable_hmvp)
        hmvp_update_lut_b(&inter_ctx->hmvp_lut, mv, mv, 0x1);
}

static void
update_gpm_mv_ctx(struct InterDRVCtx *const inter_ctx,
                const OVMV mv0, const OVMV mv1,
                VVCMergeInfo mv_info0, VVCMergeInfo mv_info1,
                uint8_t pb_x, uint8_t  pb_y,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t inter_dir0, uint8_t inter_dir1)
{   
    VVCMergeInfo mv_info;
    const struct VVCGPM *const gpm_info = &inter_ctx->gpm_ctx;

    uint8_t inter_dir = inter_dir0 | inter_dir1;

    /* FIXME can probably be simplified
     * check usage from loop
     */
    if (inter_dir == 0x1) {
        mv_info.mv0 = mv_info1.mv0;
    } else if (inter_dir == 0x2) {
        mv_info.mv1 = mv_info1.mv1;
    } else {
        if (inter_dir0 == 1 && inter_dir1 == 2) {
            mv_info.mv0     = mv_info0.mv0;
            mv_info.mv1     = mv_info1.mv1;
        } else if (inter_dir0 == 2 && inter_dir1 == 1) {
            mv_info.mv0     = mv_info1.mv0;
            mv_info.mv1     = mv_info0.mv1;
        }
    }

    mv_info.inter_dir = inter_dir;
    mv_info0.inter_dir = inter_dir0;
    mv_info1.inter_dir = inter_dir1;

    int split_dir = gpm_info->split_dir;

    int16_t angle = g_GeoParams[split_dir][0];

    int d_idx = g_GeoParams[split_dir][1];

    int x_dis = g_Dis[angle];
    int y_dis = g_Dis[(angle + (GEO_NUM_ANGLES >> 2)) % GEO_NUM_ANGLES];

    uint8_t flip = angle >= 13 && angle <= 27;

    /* FIXME use absolute coordinates instead */
    int offset_x = (-(int)nb_pb_w * 4) >> 1;
    int offset_y = (-(int)nb_pb_h * 4) >> 1;

    if (d_idx > 0) {
        if ((angle & 0xF) == 8 || ((angle & 0xF) && nb_pb_h >= nb_pb_w)) {
            offset_y += angle < 16 ? ((d_idx * nb_pb_h) >> 1) : -(int)((d_idx * nb_pb_h) >> 1);
        } else {
            offset_x += angle < 16 ? ((d_idx * nb_pb_w) >> 1) : -(int)((d_idx * nb_pb_w) >> 1);
        }
    }

    for (int y = 0; y < nb_pb_h; y++) {
        int lookup_y = (((4 * y + offset_y) * 2) + 5) * y_dis;

        for (int x = 0; x < nb_pb_w; x++) {
            int motion_idx = (((4 * x + offset_x) * 2) + 5) * x_dis + lookup_y;
            int tpm_mask = abs(motion_idx) < 32 ? 2 : (motion_idx <= 0 ? (1 - flip) : flip);

            VVCMergeInfo sbmv;

            if (tpm_mask == 2) {
                sbmv = mv_info;
            } else if (tpm_mask == 0) {
                sbmv = mv_info0;
            } else {
                sbmv = mv_info1;
            }
            update_gpm_mv_ctx_b(inter_ctx, sbmv.mv0, sbmv.mv1, pb_x + x, pb_y + y, 
                                1, 1, sbmv.inter_dir);
        }
    }
}


/* Derive motion vectors and update motion maps */
VVCMergeInfo
drv_mvp_b(struct InterDRVCtx *const inter_ctx,
          uint8_t x0, uint8_t y0,
          uint8_t log2_cb_w, uint8_t log2_cb_h,
          OVMV mvd0, OVMV mvd1, int prec_amvr,
          uint8_t mvp_idx0, uint8_t mvp_idx1, uint8_t bcw_idx,
          uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
          uint8_t is_small)
{
    OVMV mv0 = {0}, mv1 = {0};
    VVCMergeInfo mv_info;
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    uint8_t log2_pmerge_lvl = inter_ctx->log2_parallel_merge_level;
    uint8_t enable_hmvp = enable_hmvp_storage(x0_unit, y0_unit, nb_unit_w, nb_unit_h, log2_pmerge_lvl);
    /* FIXME can we combine mvp derivation for bi pred */
    if (inter_dir & 0x1) {
        uint8_t opp_ref_idx0 = inter_ctx->rpl0_opp[ref_idx0];
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        mv0 = derive_mvp_candidates_1(inter_ctx, mv_ctx0, ref_idx0,
                                      x0_unit, y0_unit, nb_unit_w, nb_unit_h,
                                      mvp_idx0, inter_dir & 0x1, mv_ctx1, opp_ref_idx0,
                                      prec_amvr, is_small);

        mvd0 = drv_change_precision_mv(mvd0, prec_amvr, MV_PRECISION_INTERNAL);

        mv0.x += mvd0.x;
        mv0.y += mvd0.y;

        mv0.ref_idx = ref_idx0;

        mv0.bcw_idx_plus1 = bcw_idx + 1;
        mv0.prec_amvr     = prec_amvr;
    }

    if (inter_dir & 0x2) {
        uint8_t opp_ref_idx1 = inter_ctx->rpl1_opp[ref_idx1];

        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

        mv1 = derive_mvp_candidates_1(inter_ctx, mv_ctx1, ref_idx1,
                                      x0_unit, y0_unit,
                                      nb_unit_w, nb_unit_h,
                                      mvp_idx1, inter_dir & 0x2, mv_ctx0, opp_ref_idx1,
                                      prec_amvr, is_small);

        mvd1 = drv_change_precision_mv(mvd1, prec_amvr, MV_PRECISION_INTERNAL);

        mv1.x += mvd1.x;
        mv1.y += mvd1.y;

        mv1.ref_idx = ref_idx1;

        mv1.bcw_idx_plus1 = bcw_idx + 1;
        mv1.prec_amvr = prec_amvr;
    }

    mv_info.inter_dir = inter_dir;
    mv_info.mv0 = mv0;
    mv_info.mv1 = mv1;

    /* Update for next pass */
    update_mv_ctx_b(inter_ctx, mv0, mv1, x0_unit, y0_unit, nb_unit_w,
                    nb_unit_h, inter_dir);

    if (enable_hmvp)
    hmvp_update_lut_b(&inter_ctx->hmvp_lut, mv0, mv1, inter_dir);

    return mv_info;
}


OVMV
drv_mmvd_merge_mvp(struct InterDRVCtx *const inter_ctx,
                   const struct OVMVCtx *const mv_ctx,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_cb_w, uint8_t log2_cb_h,
                   uint8_t merge_idx, uint8_t max_nb_cand)
{

    static const uint8_t ref_mvd_cands[8] = { 1, 2, 4, 8, 16, 32, 64, 128};
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    int smvd_mrg_idx = merge_idx / MMVD_MAX_REFINE_NUM;

    OVMV mv = vvc_derive_merge_mvp(inter_ctx, mv_ctx,
                                   x0_unit, y0_unit,
                                   nb_unit_w, nb_unit_h, smvd_mrg_idx,
                                   max_nb_cand, log2_cb_w + log2_cb_h <= 5);


    OVMV mvd0;

    int f_pos_group = merge_idx / (MMVD_BASE_MV_NUM * MMVD_MAX_REFINE_NUM);

    int idx = merge_idx;
    idx -= f_pos_group  * (MMVD_BASE_MV_NUM * MMVD_MAX_REFINE_NUM);
    idx -= smvd_mrg_idx * (MMVD_MAX_REFINE_NUM);

    int offset = (uint16_t)ref_mvd_cands[(idx >> 2)] << 2;

    offset <<= inter_ctx->mmvd_shift;

    int f_pos = idx - ((idx >> 2) << 2);

    if (f_pos == 0) {
        mvd0.x = offset;
        mvd0.y = 0;
    } else if (f_pos == 1) {
        mvd0.x = -offset;
        mvd0.y = 0;
    } else if (f_pos == 2) {
        mvd0.x = 0;
        mvd0.y = offset;
    } else {
        mvd0.x = 0;
        mvd0.y = -offset;
    }

    mv.x += mvd0.x;
    mv.y += mvd0.y;

    /* Force to RPL_0 */
    update_mv_ctx(inter_ctx, mv, x0_unit, y0_unit,
                    nb_unit_w, nb_unit_h, 0x1);



    return mv;
}


OVMV
drv_merge_mvp(struct InterDRVCtx *const inter_ctx,
              const struct OVMVCtx *const mv_ctx,
              uint8_t x0, uint8_t y0,
              uint8_t log2_cb_w, uint8_t log2_cb_h,
              uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    OVMV mv0 = vvc_derive_merge_mvp(inter_ctx, mv_ctx, x0_unit, y0_unit,
                                    nb_unit_w, nb_unit_h, merge_idx,
                                    max_nb_merge_cand,  log2_cb_w + log2_cb_h <= 5);

    update_mv_ctx(inter_ctx, mv0, x0_unit, y0_unit, nb_unit_w,
                  nb_unit_h, 1);
    return mv0;
}


OVMV
drv_mvp_mvd(struct InterDRVCtx *const inter_ctx,
            const struct OVMVCtx *const mv_ctx,
            OVMV mvd, int prec_amvr,
            uint8_t x0, uint8_t y0,
            uint8_t log2_cb_w, uint8_t log2_cb_h,
            uint8_t mvp_idx, uint8_t inter_dir,
            uint8_t ref_idx0, uint8_t ref_idx1)
{
    OVMV mv;

    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;
    uint8_t opp_ref_idx0 = inter_ctx->rpl0_opp[ref_idx0];
    uint8_t opp_ref_idx1 = inter_ctx->rpl1_opp[ref_idx1];

    const struct OVMVCtx *const mv_ctx_opp = &inter_ctx->mv_ctx0 == mv_ctx ? &inter_ctx->mv_ctx1 :
                                                                             &inter_ctx->mv_ctx0;

    uint8_t ref_idx = &inter_ctx->mv_ctx0 == mv_ctx  ? ref_idx0 : ref_idx1;
    uint8_t ref_idx_opp = &inter_ctx->mv_ctx0 == mv_ctx  ? opp_ref_idx0 : opp_ref_idx1;

    mv = derive_mvp_candidates_1(inter_ctx, mv_ctx, ref_idx,
                                 x0_unit, y0_unit, nb_unit_w, nb_unit_h,
                                 mvp_idx, 1, mv_ctx_opp, ref_idx_opp, prec_amvr, 0);

    mv.ref_idx = ref_idx;

    mvd = drv_change_precision_mv(mvd, prec_amvr, MV_PRECISION_INTERNAL);

    mv.x += mvd.x;
    mv.y += mvd.y;

    update_mv_ctx(inter_ctx, mv, x0_unit, y0_unit, nb_unit_w, nb_unit_h,
                  inter_dir);

   return mv;
}

VVCMergeInfo
drv_mmvd_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_pu_w, uint8_t log2_pu_h,
                     uint8_t merge_idx,
                     uint8_t max_nb_cand, uint8_t is_small)
{

    static const uint8_t ref_mvd_cands[8] = { 1, 2, 4, 8, 16, 32, 64, 128};
    /* FIXME better input to avoid div */
    int smvd_mrg_idx = merge_idx / MMVD_MAX_REFINE_NUM;
    uint8_t pb_x = x0 >> LOG2_MIN_CU_S;
    uint8_t pb_y = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_pb_w = (1 << log2_pu_w) >> LOG2_MIN_CU_S;
    uint8_t nb_pb_h = (1 << log2_pu_h) >> LOG2_MIN_CU_S;

    uint8_t log2_pmerge_lvl = inter_ctx->log2_parallel_merge_level;
    uint8_t enable_hmvp = enable_hmvp_storage(pb_x, pb_y, nb_pb_w, nb_pb_h, log2_pmerge_lvl);

    VVCMergeInfo mv_info = vvc_derive_merge_mvp_b(inter_ctx, pb_x, pb_y,
                                                  nb_pb_w, nb_pb_h, smvd_mrg_idx,
                                                  max_nb_cand, is_small);


    int idx = merge_idx;

    int f_pos_group = merge_idx / (MMVD_BASE_MV_NUM * MMVD_MAX_REFINE_NUM);

    idx -= f_pos_group * (MMVD_BASE_MV_NUM * MMVD_MAX_REFINE_NUM);
    idx -= smvd_mrg_idx * (MMVD_MAX_REFINE_NUM);

    int f_pos_step = idx >> 2;

    int f_pos = idx - (f_pos_step << 2);

    int offset = (uint16_t)ref_mvd_cands[f_pos_step] << 2;

    OVMV mvd0, mvd1;

    int ref0 = mv_info.mv0.ref_idx;
    int ref1 = mv_info.mv1.ref_idx;
    offset <<= inter_ctx->mmvd_shift;

    if (mv_info.inter_dir == 0x3){
        /* FIXME handle LT ref differently */

        if (f_pos == 0) {
            mvd0.x = offset;
            mvd0.y = 0;
        } else if (f_pos == 1) {
            mvd0.x = -offset;
            mvd0.y = 0;
        } else if (f_pos == 2) {
            mvd0.x = 0;
            mvd0.y = offset;
        } else {
            mvd0.x = 0;
            mvd0.y = -offset;
        }

        int32_t dist_ref0 = inter_ctx->dist_ref_0[ref0];
        int32_t dist_ref1 = inter_ctx->dist_ref_1[ref1];
        uint8_t is_lterm0 = !dist_ref0;
        uint8_t is_lterm1 = !dist_ref1;

        dist_ref0 = inter_ctx->poc - inter_ctx->rpl0[ref0]->poc;
        dist_ref1 = inter_ctx->poc - inter_ctx->rpl1[ref1]->poc;
        /* Same ref */
        if (dist_ref0 == dist_ref1){
            mvd1.x = mvd0.x;
            mvd1.y = mvd0.y;
        } else if (abs(dist_ref0) < abs(dist_ref1)){
            mvd1.x = mvd0.x;
            mvd1.y = mvd0.y;
            if (is_lterm0 || is_lterm1){
                if (dist_ref1 * dist_ref0 <= 0){
                    struct MV tmp = {.x = mvd1.x, .y = mvd1.y};
                    tmp = tmvp_scale_mv(-1, tmp);
                    mvd0.x = tmp.x;
                    mvd0.y = tmp.y;
                }
            } else {
                int scale = tmvp_compute_scale(dist_ref0, dist_ref1);
                struct MV tmp = {.x = mvd1.x, .y = mvd1.y};
                tmp = tmvp_scale_mv(scale, tmp);
                mvd0.x = tmp.x;
                mvd0.y = tmp.y;
            }
        } else {

            mvd1.x = mvd0.x;
            mvd1.y = mvd0.y;

            if (is_lterm0 || is_lterm1){
                if ((dist_ref1) * (dist_ref0) <= 0) {
                    struct MV tmp = {.x = mvd0.x, .y = mvd0.y};
                    tmp = tmvp_scale_mv(-1, tmp);
                    mvd1.x = tmp.x;
                    mvd1.y = tmp.y;
                }
            } else {
                int scale = tmvp_compute_scale(dist_ref1, dist_ref0);
                struct MV tmp = {.x = mvd0.x, .y = mvd0.y};
                tmp = tmvp_scale_mv(scale, tmp);
                mvd1.x = tmp.x;
                mvd1.y = tmp.y;
            }
        }
    } else if (mv_info.inter_dir == 0x1) {

        if (f_pos == 0) {
            mvd0.x = offset;
            mvd0.y = 0;
        } else if (f_pos == 1) {
            mvd0.x = -offset;
            mvd0.y = 0;
        } else if (f_pos == 2) {
            mvd0.x = 0;
            mvd0.y = offset;
        } else {
            mvd0.x = 0;
            mvd0.y = -offset;
        }

    } else if (mv_info.inter_dir == 0x2) {

        if (f_pos == 0) {
            mvd1.x = offset;
            mvd1.y = 0;
        } else if (f_pos == 1) {
            mvd1.x = -offset;
            mvd1.y = 0;
        } else if (f_pos == 2) {
            mvd1.x = 0;
            mvd1.y = offset;
        } else {
            mvd1.x = 0;
            mvd1.y = -offset;
        }

    }

    if (mv_info.inter_dir & 0x1) {
        mv_info.mv0.x += mvd0.x;
        mv_info.mv0.y += mvd0.y;
    }

    if (mv_info.inter_dir & 0x2) {
        mv_info.mv1.x += mvd1.x;
        mv_info.mv1.y += mvd1.y;
    }

    /* FIXME check before ? */
    if (is_small && mv_info.inter_dir == 0x3) {
        mv_info.inter_dir = 0x1;
    }

    update_mv_ctx_b(inter_ctx, mv_info.mv0, mv_info.mv1, pb_x, pb_y,
                    nb_pb_w, nb_pb_h, mv_info.inter_dir);


    if (enable_hmvp)
    hmvp_update_lut_b(&inter_ctx->hmvp_lut, mv_info.mv0, mv_info.mv1, mv_info.inter_dir);

    return mv_info;
}

void 
drv_gpm_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                    uint8_t x0, uint8_t y0,
                    uint8_t log2_pu_w, uint8_t log2_pu_h,
                    uint8_t max_nb_cand, uint8_t is_small)
{
    struct VVCGPM* gpm_ctx = &inter_ctx->gpm_ctx;
    VVCMergeInfo mv_info0, mv_info1;
    uint8_t pb_x = x0 >> LOG2_MIN_CU_S;
    uint8_t pb_y = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_pb_w = (1 << log2_pu_w) >> LOG2_MIN_CU_S;
    uint8_t nb_pb_h = (1 << log2_pu_h) >> LOG2_MIN_CU_S;

    mv_info0 = vvc_derive_merge_mvp_b(inter_ctx, pb_x, pb_y,
                                     nb_pb_w, nb_pb_h, gpm_ctx->merge_idx0,
                                     max_nb_cand, is_small);

    if (gpm_ctx->merge_idx0 != gpm_ctx->merge_idx1) { 
        mv_info1 = vvc_derive_merge_mvp_b(inter_ctx, pb_x, pb_y,
                                         nb_pb_w, nb_pb_h, gpm_ctx->merge_idx1,
                                         max_nb_cand, is_small);
    } else{
        mv_info1 = mv_info0;
    }

    mv_info0.mv0.bcw_idx_plus1 = 0;
    mv_info0.mv1.bcw_idx_plus1 = 0;
    mv_info1.mv0.bcw_idx_plus1 = 0;
    mv_info1.mv1.bcw_idx_plus1 = 0;

    mv_info0.mv0.prec_amvr = 0;
    mv_info0.mv1.prec_amvr = 0;
    mv_info1.mv0.prec_amvr = 0;
    mv_info1.mv1.prec_amvr = 0;

    uint8_t parity = gpm_ctx->merge_idx0 & 1;

    if (mv_info0.inter_dir & (0x01 + parity)) {
        gpm_ctx->inter_dir0 = 1 + parity;
        gpm_ctx->mv0 = parity ? mv_info0.mv1 : mv_info0.mv0; 
    } else if (mv_info0.inter_dir & (0x02 - parity)) {
        gpm_ctx->inter_dir0 = 2 - parity;
        gpm_ctx->mv0 = parity ? mv_info0.mv0 : mv_info0.mv1; 
    }   

    parity = gpm_ctx->merge_idx1 & 1;

    if(mv_info1.inter_dir & (0x01 + parity)) {
        gpm_ctx->inter_dir1 = 1 + parity;
        gpm_ctx->mv1 = parity ? mv_info1.mv1 : mv_info1.mv0; 
    } else if (mv_info1.inter_dir & (0x02 - parity)) {
        gpm_ctx->inter_dir1 = 2 - parity;
        gpm_ctx->mv1 = parity ? mv_info1.mv0 : mv_info1.mv1; 
    }  

    update_gpm_mv_ctx(inter_ctx, gpm_ctx->mv0, gpm_ctx->mv1, mv_info0, mv_info1, pb_x, pb_y,
                    nb_pb_w, nb_pb_h, gpm_ctx->inter_dir0, gpm_ctx->inter_dir1);
}

VVCMergeInfo
drv_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                uint8_t x0, uint8_t y0,
                uint8_t log2_cb_w, uint8_t log2_cb_h,
                uint8_t merge_idx,
                uint8_t max_nb_cand, uint8_t is_small)
{
    VVCMergeInfo mv_info;
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;
    uint8_t log2_pmerge_lvl = inter_ctx->log2_parallel_merge_level;
    uint8_t enable_hmvp = enable_hmvp_storage(x0_unit, y0_unit, nb_unit_w, nb_unit_h, log2_pmerge_lvl);

    mv_info = vvc_derive_merge_mvp_b(inter_ctx, x0_unit, y0_unit,
                                     nb_unit_w, nb_unit_h, merge_idx,
                                     max_nb_cand, is_small);

    if (is_small && mv_info.inter_dir == 3) {
        mv_info.inter_dir = 0x1;
    }

    update_mv_ctx_b(inter_ctx, mv_info.mv0, mv_info.mv1, x0_unit, y0_unit,
                    nb_unit_w, nb_unit_h, mv_info.inter_dir);

    if (enable_hmvp)
    hmvp_update_lut_b(&inter_ctx->hmvp_lut, mv_info.mv0, mv_info.mv1, mv_info.inter_dir);

    return mv_info;
}
