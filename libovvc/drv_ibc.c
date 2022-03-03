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

#include "ctudec.h"
#include "drv_utils.h"
#include "dec_structures.h"
#include "ovdpb.h"

#define POS_MASK(x, w) ((uint64_t) 1 << ((((x + 1)) + ((w)))))
#define LOG2_MIN_CU_S 2

/* FIXME simplify for IBC AMVR based on precision */
static IBCMV
drv_change_precision_mv(IBCMV mv, int src, int dst)
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

static IBCMV
drv_round_to_precision_mv(IBCMV mv, int src, int dst)
{
    IBCMV mv1 = drv_change_precision_mv(mv, src, dst);
    return drv_change_precision_mv(mv1, dst, src);
}

#define MV_CMP(mv0,mv1) ((mv0.x == mv1.x) && (mv0.y == mv1.y))
static inline uint8_t
ibc_mi_cmp(const IBCMV *const a, const IBCMV *const b)
{
    uint8_t is_equal = (a->x == b->x) & (a->y == b->y);
    return is_equal;
}

static void
ibc_update_hmvp_lut(struct IBCMVCtx *const ibc_ctx, const IBCMV mv)
{
    int max_nb_cand = OVMIN(5, ibc_ctx->nb_hmvp_cand);
    IBCMV *hmvp_lut = ibc_ctx->hmvp_lut;
    uint8_t duplicated_mv = 0;
    int i;

    for (i = 0; i < max_nb_cand; ++i) {
        duplicated_mv = ibc_mi_cmp(&hmvp_lut[i], &mv);
        if (duplicated_mv) {
            break;
        }
    }

    if (duplicated_mv) {
       int j;
       /* Erase duplicated MV in LUT and append new MV
        * at the end of LUT
        */
       for (j = i; j < max_nb_cand - 1; ++j) {
           hmvp_lut[j] = hmvp_lut[j + 1];
       }
       hmvp_lut[j] = mv;
    } else if (ibc_ctx->nb_hmvp_cand == 5) {
       /* Erase firts LUT element MV  and append new MV
        * at the end of LUT
        */
       int j;
       for (j = 1; j < 5; ++j) {
           hmvp_lut[j - 1] = hmvp_lut[j];
       }
       hmvp_lut[4] = mv;
    } else {
       /* Append new MV at the end of LUT and increase counter
        */
        hmvp_lut[ibc_ctx->nb_hmvp_cand++] = mv;
    }
}

static IBCMV
ibc_hmvp_mvp_cand(const struct IBCMVCtx *const ibc_ctx)
{
    IBCMV cand = { 0 };
    int nb_lut_cand = ibc_ctx->nb_hmvp_cand;
    if (nb_lut_cand) {
        cand = ibc_ctx->hmvp_lut[nb_lut_cand - 1];
    }
    return cand;
}

#define A1_MSK 0x1
#define B1_MSK 0x2

static IBCMV
ibc_derive_hmvp_merge_cand(const struct IBCMVCtx *const ibc_ctx,
                           const IBCMV const a1_b1[2],
                           uint8_t status_msk, int8_t nb_cand,
                           int8_t merge_idx)
{
    int i;
    IBCMV cand = { 0 };
    int8_t nb_lut_cand = ibc_ctx->nb_hmvp_cand;
    int8_t target_idx = merge_idx - nb_cand;
    if (target_idx < nb_lut_cand) {
        IBCMV first_cand = ibc_ctx->hmvp_lut[ibc_ctx->nb_hmvp_cand - 1];
        uint8_t already_cand = 0;
        if (status_msk & A1_MSK) already_cand |= ibc_mi_cmp(&a1_b1[0], &first_cand);
        if (status_msk & B1_MSK) already_cand |= ibc_mi_cmp(&a1_b1[1], &first_cand);
        if (already_cand) {
             /* skip first_cand */
             nb_lut_cand--;
        }

        int8_t lut_idx = nb_lut_cand - target_idx - 1;

        if (lut_idx >= 0) {
            cand = ibc_ctx->hmvp_lut[lut_idx];
        }
    }
    return cand;
}

static IBCMV
ibc_derive_mvp_mv(struct IBCMVCtx *const ibc_ctx,
                  uint8_t x0_unit, uint8_t y0_unit,
                  uint8_t nb_unit_w, uint8_t nb_unit_h,
                  uint8_t mvp_idx, uint8_t prec_amvr)
{
    uint64_t lft_col = ibc_ctx->ctu_map.vfield[x0_unit];
    uint64_t abv_row = ibc_ctx->ctu_map.hfield[y0_unit];

    IBCMV cand = {0};

    /* If small do not check for A1 and B1 candidates and look directly into
       HMVP LUT */
    if ((nb_unit_h | nb_unit_w) == 1) goto hmvp;

    uint8_t cand_a1 = !!(lft_col & POS_MASK(y0_unit, nb_unit_h - 1));
    uint8_t cand_b1 = !!(abv_row & POS_MASK(x0_unit, nb_unit_w - 1));

    if (cand_a1) {
        int cand_pos = y0_unit + nb_unit_h - 1;
        cand = ibc_ctx->lft_col[cand_pos];
        goto found;
    }

    if (cand_b1) {
        int cand_pos = x0_unit + nb_unit_w - 1;
        cand = ibc_ctx->abv_row[cand_pos];
        goto found;
    }

hmvp:
    cand = ibc_hmvp_mvp_cand(ibc_ctx);

found:
    cand = drv_round_to_precision_mv(cand, MV_PRECISION_INTERNAL, prec_amvr);
    return cand;
}

static IBCMV
ibc_derive_merge_mv(const struct IBCMVCtx *const ibc_ctx,
                    uint8_t x0_unit, uint8_t y0_unit,
                    uint8_t nb_unit_w, uint8_t nb_unit_h,
                    uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    IBCMV amvp_cand[2];
    IBCMV cand = { 0 };
    uint64_t lft_col = ibc_ctx->ctu_map.vfield[x0_unit];
    uint64_t abv_row = ibc_ctx->ctu_map.hfield[y0_unit];

    int nb_cand = 0;

    if ((nb_unit_h | nb_unit_w) == 1) goto hmvp;

    uint8_t cand_a1 = !!(lft_col & POS_MASK(y0_unit, nb_unit_h - 1));
    uint8_t cand_b1 = !!(abv_row & POS_MASK(x0_unit, nb_unit_w - 1));

    if (cand_a1) {
        int cand_pos = y0_unit + nb_unit_h - 1;
        amvp_cand[0] = ibc_ctx->lft_col[cand_pos];
        if (nb_cand++ == merge_idx) {
            cand = amvp_cand[0];
            goto found;
        }
    }

    if (cand_b1) {
        int cand_pos = x0_unit + nb_unit_w - 1;
        amvp_cand[1] = ibc_ctx->abv_row[cand_pos];
        if (!cand_a1 || !ibc_mi_cmp(&amvp_cand[0], &amvp_cand[1])) {
            if (nb_cand++ == merge_idx) {
                cand = amvp_cand[1];
                    goto found;
            }
        }
    }


hmvp:
    /* FIXME check if check on max_nb_merge_cand_min1 required */
    if (nb_cand != max_nb_merge_cand - 1) {
        uint8_t status_msk = cand_a1 | (cand_b1 << 1);
        /* FIXME small */
        cand = ibc_derive_hmvp_merge_cand(ibc_ctx, amvp_cand, status_msk, nb_cand, merge_idx);
    }

found:

    return cand;
}

static void
ibc_fill_mvp_map(struct IBCMVCtx *const ibc_ctx,
                 IBCMV mv,
                 int x0_unit, int y0_unit,
                 int nb_unit_w, int nb_unit_h)
{
    int i;
    IBCMV *abv_row = &ibc_ctx->abv_row[x0_unit + 0];
    IBCMV *lft_col = &ibc_ctx->lft_col[y0_unit + 0];

    ctu_field_set_rect_bitfield(&ibc_ctx->ctu_map, x0_unit, y0_unit, nb_unit_w, nb_unit_h);

    for (i = 0; i < nb_unit_w; ++i) {
        abv_row[i] = mv;
    }

    for (i = 0; i < nb_unit_h; ++i) {
        lft_col[i] = mv;
    }
}

static void
ibc_update_mv_ctx(struct IBCMVCtx *const ibc_ctx,
                  uint8_t x0_unit, uint8_t y0_unit,
                  uint8_t nb_unit_w, uint8_t nb_unit_h,
                  const IBCMV mv)
{
    ibc_fill_mvp_map(ibc_ctx, mv, x0_unit, y0_unit, nb_unit_w, nb_unit_h);

    ibc_update_hmvp_lut(ibc_ctx, mv);
}

IBCMV
drv_ibc_merge_mv(struct IBCMVCtx *const ibc_ctx,
                 uint8_t x0, uint8_t y0,
                 uint8_t log2_cb_w, uint8_t log2_cb_h,
                 uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    IBCMV mv = ibc_derive_merge_mv(ibc_ctx, x0_unit, y0_unit,
                                    nb_unit_w, nb_unit_h, merge_idx,
                                    max_nb_merge_cand);

    ibc_update_mv_ctx(ibc_ctx, x0_unit, y0_unit, nb_unit_w, nb_unit_h, mv);

    mv = drv_change_precision_mv(mv, MV_PRECISION_INTERNAL, MV_PRECISION_INT);

    return mv;
}

IBCMV
drv_ibc_mvp(struct IBCMVCtx *const ibc_ctx,
            uint8_t x0, uint8_t y0,
            uint8_t log2_cb_w, uint8_t log2_cb_h,
            IBCMV mvd, uint8_t mvp_idx, uint8_t prec_amvr)
{
    IBCMV mv;

    uint8_t x0_unit = x0 >> LOG2_MIN_CU_S;
    uint8_t y0_unit = y0 >> LOG2_MIN_CU_S;
    uint8_t nb_unit_w = (1 << log2_cb_w) >> LOG2_MIN_CU_S;
    uint8_t nb_unit_h = (1 << log2_cb_h) >> LOG2_MIN_CU_S;

    mv = ibc_derive_mvp_mv(ibc_ctx, x0_unit, y0_unit,
                           nb_unit_w, nb_unit_h, mvp_idx, prec_amvr);

    mvd = drv_change_precision_mv(mvd, prec_amvr, MV_PRECISION_INTERNAL);

    mv.x += mvd.x;
    mv.y += mvd.y;

    ibc_update_mv_ctx(ibc_ctx, x0_unit, y0_unit, nb_unit_w, nb_unit_h, mv);

    mv = drv_change_precision_mv(mv, MV_PRECISION_INTERNAL, MV_PRECISION_INT);

    return mv;
}
