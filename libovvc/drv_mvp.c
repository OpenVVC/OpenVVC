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

#define OFFSET_BUFF(x,y) (35 + x + (y) * 34)
#define MV_CMP(mv0,mv1) ((mv0.x == mv1.x) && (mv0.y == mv1.y))

#define PB_POS_IN_BUF(x,y) (35 + (x) + ((y) * 34))

#define MV_MANTISSA_BITCOUNT 6
#define MV_MANTISSA_UPPER_LIMIT ((1 << (MV_MANTISSA_BITCOUNT - 1)) - 1)
#define MV_MANTISSA_LIMIT (1 << (MV_MANTISSA_BITCOUNT - 1))

#define MV_BITS  18
#define MV_MAX   ((1 << (MV_BITS - 1)) - 1)
#define MV_MIN  (-(1 << (MV_BITS - 1)))

static inline OVMV
scale_mvd(OVMV mv)
{
    mv.x <<= 2;
    mv.y <<= 2;
    return mv;
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
static int
tmvp_round_mv(int32_t val)
{
  int sign  = val >> 31;

  if ((val ^ sign) - !!sign > 31) {
      int scale = floor_log2((val ^ sign) | MV_MANTISSA_UPPER_LIMIT) - (MV_MANTISSA_BITCOUNT - 1);
      int round = (1 << scale) >> 1;
      int n     = (val + round) >> scale;
      int exponent  = scale + ((n ^ sign) >> (MV_MANTISSA_BITCOUNT - 1));
      int mantissa  = (n & MV_MANTISSA_UPPER_LIMIT) | (sign << (MV_MANTISSA_BITCOUNT - 1));
      return (mantissa ^ MV_MANTISSA_LIMIT) << (exponent - !!exponent);
  } else {
      return val;
  }
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

static void
hmvp_add_cand(const struct HMVPLUT *const hmvp_lut,
              OVMV *const cand_list,
              int *const nb_cand, uint8_t inter_dir)
{
    int max_nb_cand = OVMIN(4, hmvp_lut->nb_mv);
    int i;
    for (i = 1; i <= max_nb_cand && *nb_cand < 2;  ++i) {
        if(hmvp_lut->dir[i - 1] & inter_dir) {
            cand_list[(*nb_cand)++] = inter_dir & 0x1 ? hmvp_lut->hmv0[i - 1]
                                                      : hmvp_lut->hmv1[i - 1];
        }
    }
}

static uint8_t
hmvp_add_merge_cand(const struct HMVPLUT *const hmvp_lut,
                    OVMV *const cand_list,
                    OVMV *const cand_amvp,
                    uint8_t got_B1, uint8_t got_A1,
                    int *const nb_cand, uint8_t inter_dir, uint8_t mvp_idx)
{
    #if 0
    int max_nb_cand = OVMIN(4, hmvp_lut->nb_mv);
    #endif
    int i;
    for (i = 1; i <= hmvp_lut->nb_mv;  ++i) {
        int mv_lut_idx = hmvp_lut->nb_mv - i;
        if (hmvp_lut->dir[mv_lut_idx] & inter_dir) {
            if (i > 2 || ((!got_B1 || !MV_CMP(hmvp_lut->hmv0[mv_lut_idx], cand_amvp[0])) &&
                (!got_A1 || !MV_CMP(hmvp_lut->hmv0[mv_lut_idx], cand_amvp[1])))) {
                cand_list[(*nb_cand)++] = inter_dir & 0x1 ? hmvp_lut->hmv0[mv_lut_idx]
                                                      : hmvp_lut->hmv1[mv_lut_idx];
                if (*nb_cand == mvp_idx + 1) {
                    return 1;
                }

                if (*nb_cand == 4) {
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
    #if 0
    int max_nb_cand = OVMIN(4, hmvp_lut->nb_mv);
    #endif
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
hmvp_update_lut(struct HMVPLUT *const hmvp_lut, OVMV mv)
{
    int max_nb_cand = OVMIN(5, hmvp_lut->nb_mv);
    int duplicated_mv = 0;
    int i;

    for (i = 0; i < max_nb_cand; ++i) {
        if (MV_CMP(mv, hmvp_lut->hmv0[i])) {
            duplicated_mv = 1;
            break;
        }
    }

    if (duplicated_mv) {
        int j;

        for (j = i; j < max_nb_cand - 1; ++j) {
            hmvp_lut->hmv0[j] = hmvp_lut->hmv0[j + 1];
        }

        hmvp_lut->hmv0[j] = mv;

    } else if (hmvp_lut->nb_mv == 5) {
        int j;

        for (j = 1; j < 5; ++j) {
            hmvp_lut->hmv0[j - 1] = hmvp_lut->hmv0[j];
        }

        hmvp_lut->hmv0[4] = mv;

    } else {
        hmvp_lut->hmv0[hmvp_lut->nb_mv++] = mv;
    }

    hmvp_lut->dir[0]= 0x1;
    hmvp_lut->dir[1]= 0x1;
    hmvp_lut->dir[2]= 0x1;
    hmvp_lut->dir[3]= 0x1;
    hmvp_lut->dir[4]= 0x1;
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


static void
load_ctb_tmvp(OVCTUDec *const ctudec, int ctb_x, int ctb_y)
{
    uint8_t log2_ctb_s = ctudec->part_ctx->log2_ctu_s;
    uint8_t log2_min_cb_s = ctudec->part_ctx->log2_min_cb_s;
    uint8_t nb_pb_ctb_w = (1 << log2_ctb_s) >> log2_min_cb_s;
    uint16_t nb_ctb_w = ctudec->nb_ctb_pic_w;
    uint16_t ctb_addr_rs = ctb_x + ctb_y * nb_ctb_w;
    uint8_t is_border_pic = nb_ctb_w - 1 == ctb_x;

    struct InterDRVCtx *const inter_ctx = &ctudec->drv_ctx.inter_ctx;
    struct VVCTMVP *const tmvp_ctx = &inter_ctx->tmvp_ctx;

    const struct MVPlane *plane0 = tmvp_ctx->col_plane0;
    const struct MVPlane *plane1 = tmvp_ctx->col_plane1;

    if (is_border_pic) {
        memset(tmvp_ctx->dir_map_v0, 0, sizeof(uint64_t) * 33);
        memset(tmvp_ctx->dir_map_v1, 0, sizeof(uint64_t) * 33);
    }

    if (plane0)
    if (plane0->dirs) {
        uint64_t *src_dirs = plane0->dirs + ctb_addr_rs * nb_pb_ctb_w;

        OVMV *src_mv = plane0->mvs + ctb_x * nb_pb_ctb_w + (ctb_y * nb_pb_ctb_w *nb_pb_ctb_w) * nb_ctb_w;
        #if 0
        struct OVMVCtx *mv_ctx = &tmvp_ctx->tmvp_mv.mv_ctx0;
        #else
        OVMV *mvs = tmvp_ctx->mvs0;
        #endif
        int i;

        memcpy(&tmvp_ctx->dir_map_v0[1], src_dirs, sizeof(uint64_t) * (nb_pb_ctb_w + !is_border_pic));
        for (i = 0; i < nb_pb_ctb_w; ++i) {
            memcpy(&mvs[1 + 34 * (i + 1)], src_mv, sizeof(*src_mv) * (nb_pb_ctb_w + !is_border_pic));
            src_mv += nb_pb_ctb_w * nb_ctb_w;
        }
    }

    if (plane1)
    if (plane1->dirs) {
        #if 0
        struct OVMVCtx *mv_ctx = &tmvp_ctx->tmvp_mv.mv_ctx1;
        #else
        OVMV *mvs = tmvp_ctx->mvs1;
        #endif
        uint64_t *src_dirs = plane1->dirs + ctb_addr_rs * nb_pb_ctb_w;
        int i;

        OVMV *src_mv = plane1->mvs + ctb_x * nb_pb_ctb_w + (ctb_y * nb_pb_ctb_w *nb_pb_ctb_w) * nb_ctb_w;

        /*FIXME memory could be spared with smaller map size when possible */
        memcpy(&tmvp_ctx->dir_map_v1[1], src_dirs, sizeof(uint64_t) * (nb_pb_ctb_w + !is_border_pic));
        for (i = 0; i < nb_pb_ctb_w; ++i) {
            memcpy(&mvs[1 + 34 * (i + 1)], src_mv, sizeof(*src_mv) * (nb_pb_ctb_w + !is_border_pic));
            src_mv += nb_pb_ctb_w * nb_ctb_w;
        }
    }

    inter_ctx->tmvp_avail |= 1;
}

static inline OVMV
tmvp_scale_mv(int scale, OVMV mv)
{
    mv.x = ov_clip((scale * mv.x + 128 - (scale * mv.x >= 0)) >> 8, MV_MIN, MV_MAX);
    mv.y = ov_clip((scale * mv.y + 128 - (scale * mv.y >= 0)) >> 8, MV_MIN, MV_MAX);
    return mv;
}

static OVMV
derive_mvp_candidates_1(struct InterDRVCtx *const inter_ctx,
                        const struct OVMVCtx *const mv_ctx,
                        uint8_t ref_idx,
                        uint8_t pb_x, uint8_t pb_y,
                        uint8_t nb_pb_w, uint8_t nb_pb_h,
                        uint8_t mvp_idx, uint8_t inter_dir,
                        const struct OVMVCtx *const mv_ctx_opp, uint8_t opp_ref_idx)
{
    const OVMV *const mv_buff = mv_ctx->mvs;
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

    cand[0].x += 2 - (cand[0].x >= 0);
    cand[0].y += 2 - (cand[0].y >= 0);

    cand[0].x >>= 2;
    cand[0].y >>= 2;

    cand[0].x <<= 2;
    cand[0].y <<= 2;

    cand[1].x += 2 - (cand[1].x >= 0);
    cand[1].y += 2 - (cand[1].y >= 0);

    cand[1].x >>= 2;
    cand[1].y >>= 2;

    cand[1].x <<= 2;
    cand[1].y <<= 2;


    /* Remove on candidates if duplicated */
    if (nb_cand == 2) {
        if (MV_CMP(cand[0], cand[1])) {
            --nb_cand;
        }
    }

    #if 0
    if (inter_ctx->tmvp_enabled && nb_cand < 2) {
    }
    #else
    if (inter_ctx->tmvp_enabled && nb_cand < 2) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        uint64_t c1_col;
        uint64_t c0_col;
        uint64_t c1_col1;
        uint64_t c0_col1;
        uint8_t cand_c0;
        uint8_t cand_c1;
        uint8_t cand_c01;
        uint8_t cand_c11;
        int c1_x = pb_x + (nb_pb_w >> 1);
        int c1_y = pb_y + (nb_pb_h >> 1);
        int c0_x = pb_x + nb_pb_w;
        int c0_y = pb_y + nb_pb_h;
        int scale0, scale1;
        #if 1
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }

        #endif

        /* Derive availability based on CTB inter fields */
        c0_col  = tmvp->dir_map_v0[c0_x + 1];
        c0_col1 = tmvp->dir_map_v1[c0_x + 1];
        c1_col  = tmvp->dir_map_v0[c1_x + 1];
        c1_col1 = tmvp->dir_map_v1[c1_x + 1];

        cand_c0  = !!(c0_col  & POS_MASK(pb_y, nb_pb_h));
        cand_c01 = !!(c0_col1 & POS_MASK(pb_y, nb_pb_h));
        cand_c1  = !!(c1_col  & POS_MASK(pb_y, nb_pb_h >> 1));
        cand_c11 = !!(c1_col1 & POS_MASK(pb_y, nb_pb_h >> 1));

        /*FIXME there might be an issue considering the order of RPL check 
         * for TMVP candidates
         */
        if (mv_ctx == &inter_ctx->mv_ctx0) {
            scale0 = tmvp->scale00;
            scale1 = tmvp->scale01;
        } else {
            scale0 = tmvp->scale10;
            scale1 = tmvp->scale11;
        }

        if (cand_c0) {
            /* Candidate 0 in collocated picture 0 */
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            #if 0
            OVMV c0 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c0 = tmvp->mvs0[pos_in_buff];
            #endif
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale0, c0);
            c0.x = ((c0.x + 2 - (c0.x >= 0)) >> 2) << 2;
            c0.y = ((c0.y + 2 - (c0.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c0;
        } else if (cand_c01) {
            /* Candidate 0 in collocated picture 1 */
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            #if 0
            OVMV c0 = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
            #else
            OVMV c0 = tmvp->mvs1[pos_in_buff];
            #endif
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale1, c0);
            c0.x = ((c0.x + 2 - (c0.x >= 0)) >> 2) << 2;
            c0.y = ((c0.y + 2 - (c0.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c0;
        } else if (cand_c1) {
            /* Candidate 1 in collocated picture 0 */
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            #if 0
            OVMV c1 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c1 = tmvp->mvs0[pos_in_buff];
            #endif
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale0, c1);
            c1.x = ((c1.x + 2 - (c1.x >= 0)) >> 2) << 2;
            c1.y = ((c1.y + 2 - (c1.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c1;
        } else if (cand_c11) {
            /* Candidate 1 in collocated picture 1 */
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            #if 0
            OVMV c1 = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
            #else
            OVMV c1 = tmvp->mvs1[pos_in_buff];
            #endif
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale1, c1);
            c1.x = ((c1.x + 2 - (c1.x >= 0)) >> 2) << 2;
            c1.y = ((c1.y + 2 - (c1.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c1;
        }
    }
    #endif

    if (nb_cand < 2) {
        const struct HMVPLUT *hmvp_lut = &inter_ctx->hmvp_lut;
        hmvp_add_cand_1(hmvp_lut, cand, &nb_cand, inter_dir, ref_idx, opp_ref_idx);
    }

    while (nb_cand < 2) {
        OVMV zmv = {0};
        zmv.ref_idx = ref_idx;
        cand[nb_cand++] = zmv;
    }

    cand[0].x += 2 - (cand[0].x >= 0);
    cand[0].y += 2 - (cand[0].y >= 0);
    cand[1].x += 2 - (cand[1].x >= 0);
    cand[1].y += 2 - (cand[1].y >= 0);
    cand[0].x >>= 2;
    cand[0].y >>= 2;
    cand[1].x >>= 2;
    cand[1].y >>= 2;
    cand[0].x <<= 2;
    cand[0].y <<= 2;
    cand[1].x <<= 2;
    cand[1].y <<= 2;

    return cand[mvp_idx];
}

static OVMV
derive_mvp_candidates(struct InterDRVCtx *const inter_ctx,
                      const struct OVMVCtx *const mv_ctx,
                      uint8_t pb_x, uint8_t pb_y,
                      uint8_t nb_pb_w, uint8_t nb_pb_h,
                      uint8_t mvp_idx, uint8_t inter_dir)
{
    const OVMV *const mv_buff = mv_ctx->mvs;
    uint64_t lft_col = mv_ctx->map.vfield[pb_x];
    uint64_t abv_row = mv_ctx->map.hfield[pb_y];

    OVMV cand[2] = {0};
    int nb_cand = 0;

    /* Derive candidates availability based on CTU inter fields */
    uint8_t cand_bl = !!(lft_col & POS_MASK(pb_y, nb_pb_h));     /*AO*/
    uint8_t cand_l  = !!(lft_col & POS_MASK(pb_y, nb_pb_h - 1)); /*A1*/
    uint8_t cand_tr = !!(abv_row & POS_MASK(pb_x, nb_pb_w));     /*B0*/
    uint8_t cand_t  = !!(abv_row & POS_MASK(pb_x, nb_pb_w - 1)); /*B1*/
    uint8_t cand_tl = !!(abv_row & POS_MASK(pb_x - 1, 0));      /*B2*/

    if (cand_bl) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    } else if (cand_l) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h - 1);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    }

    if (cand_tr) {
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w, pb_y - 1);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    } else if (cand_t) {
        int pos_in_buff = OFFSET_BUFF(pb_x + nb_pb_w - 1, pb_y - 1);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    } else if (cand_tl) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y - 1);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    }

    cand[0].x += 2 - (cand[0].x >= 0);
    cand[0].y += 2 - (cand[0].y >= 0);

    cand[0].x >>= 2;
    cand[0].y >>= 2;

    cand[0].x <<= 2;
    cand[0].y <<= 2;

    cand[1].x += 2 - (cand[1].x >= 0);
    cand[1].y += 2 - (cand[1].y >= 0);

    cand[1].x >>= 2;
    cand[1].y >>= 2;

    cand[1].x <<= 2;
    cand[1].y <<= 2;


    /* Remove on candidates if duplicated */
    if (nb_cand == 2) {
        if (MV_CMP(cand[0], cand[1])) {
            --nb_cand;
        }
    }

    #if 0
    if (inter_ctx->tmvp_enabled && nb_cand < 2) {
    }
    #else
    if (inter_ctx->tmvp_enabled && nb_cand < 2) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        uint64_t c1_col;
        uint64_t c0_col;
        uint64_t c1_col1;
        uint64_t c0_col1;
        uint8_t cand_c0;
        uint8_t cand_c1;
        uint8_t cand_c01;
        uint8_t cand_c11;
        int c1_x = pb_x + (nb_pb_w >> 1);
        int c1_y = pb_y + (nb_pb_h >> 1);
        int c0_x = pb_x + nb_pb_w;
        int c0_y = pb_y + nb_pb_h;
        int scale0, scale1;
        #if 1
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }

        #endif

        /* Derive availability based on CTB inter fields */
        c0_col  = tmvp->dir_map_v0[c0_x + 1];
        c0_col1 = tmvp->dir_map_v1[c0_x + 1];
        c1_col  = tmvp->dir_map_v0[c1_x + 1];
        c1_col1 = tmvp->dir_map_v1[c1_x + 1];

        cand_c0  = !!(c0_col  & POS_MASK(pb_y, nb_pb_h));
        cand_c01 = !!(c0_col1 & POS_MASK(pb_y, nb_pb_h));
        cand_c1  = !!(c1_col  & POS_MASK(pb_y, nb_pb_h >> 1));
        cand_c11 = !!(c1_col1 & POS_MASK(pb_y, nb_pb_h >> 1));

        /*FIXME there might be an issue considering the order of RPL check 
         * for TMVP candidates
         */
        if (mv_ctx == &inter_ctx->mv_ctx0) {
            scale0 = tmvp->scale00;
            scale1 = tmvp->scale01;
        } else {
            scale0 = tmvp->scale10;
            scale1 = tmvp->scale11;
        }

        if (cand_c0) {
            /* Candidate 0 in collocated picture 0 */
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            #if 0
            OVMV c0 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c0 = tmvp->mvs0[pos_in_buff];
            #endif
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale0, c0);
            c0.x = ((c0.x + 2 - (c0.x >= 0)) >> 2) << 2;
            c0.y = ((c0.y + 2 - (c0.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c0;
        } else if (cand_c01) {
            /* Candidate 0 in collocated picture 1 */
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            #if 0
            OVMV c0 = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
            #else
            OVMV c0 = tmvp->mvs1[pos_in_buff];
            #endif
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale1, c0);
            c0.x = ((c0.x + 2 - (c0.x >= 0)) >> 2) << 2;
            c0.y = ((c0.y + 2 - (c0.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c0;
        } else if (cand_c1) {
            /* Candidate 1 in collocated picture 0 */
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            #if 0
            OVMV c1 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c1 = tmvp->mvs0[pos_in_buff];
            #endif
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale0, c1);
            c1.x = ((c1.x + 2 - (c1.x >= 0)) >> 2) << 2;
            c1.y = ((c1.y + 2 - (c1.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c1;
        } else if (cand_c11) {
            /* Candidate 1 in collocated picture 1 */
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            #if 0
            OVMV c1 = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
            #else
            OVMV c1 = tmvp->mvs1[pos_in_buff];
            #endif
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale1, c1);
            c1.x = ((c1.x + 2 - (c1.x >= 0)) >> 2) << 2;
            c1.y = ((c1.y + 2 - (c1.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c1;
        }
    }
    #endif

    if (nb_cand < 2) {
        const struct HMVPLUT *hmvp_lut = &inter_ctx->hmvp_lut;
        hmvp_add_cand(hmvp_lut, cand, &nb_cand, inter_dir);
    }

    while (nb_cand < 2) {
        OVMV zmv = {0};
        cand[nb_cand++] = zmv;
    }

    cand[0].x += 2 - (cand[0].x >= 0);
    cand[0].y += 2 - (cand[0].y >= 0);
    cand[1].x += 2 - (cand[1].x >= 0);
    cand[1].y += 2 - (cand[1].y >= 0);
    cand[0].x >>= 2;
    cand[0].y >>= 2;
    cand[1].x >>= 2;
    cand[1].y >>= 2;
    cand[0].x <<= 2;
    cand[0].y <<= 2;
    cand[1].x <<= 2;
    cand[1].y <<= 2;

    return cand[mvp_idx];
}

static OVMV
vvc_derive_merge_mvp(const struct InterDRVCtx *const inter_ctx,
                     const struct OVMVCtx *const mv_ctx,
                     uint8_t pb_x, uint8_t pb_y,
                     uint8_t nb_pb_w, uint8_t nb_pb_h,
                     uint8_t merge_idx, uint8_t max_nb_merge_cand)
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
        if (!cand_t || !MV_CMP(mv_A1, cand_amvp[0])) {
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
        if (!cand_t || !MV_CMP(mv_B0, cand_amvp[0])) {
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
        if (!cand_l || !MV_CMP(mv_A0, cand_amvp[1])) {
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
            if ((!cand_l || !MV_CMP(mv_B2, cand_amvp[1])) &&
                (!cand_t || !MV_CMP(mv_B2, cand_amvp[0]))) {
                cand[nb_cand] = mv_B2;
                if (nb_cand++ == merge_idx)
                    return mv_B2;
            }
        }
    }


    #if 1
    /*FIXME TMVP disabled for 4x8 8x4 blocks*/
    if (inter_ctx->tmvp_enabled) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        uint64_t c1_col;
        uint64_t c0_col;
        uint64_t c1_col1;
        uint64_t c0_col1;
        uint8_t cand_c0;
        uint8_t cand_c1;
        uint8_t cand_c01;
        uint8_t cand_c11;
        int c1_x = pb_x + (nb_pb_w >> 1);
        int c1_y = pb_y + (nb_pb_h >> 1);
        int c0_x = pb_x + nb_pb_w;
        int c0_y = pb_y + nb_pb_h;
        /*FIXME determine whether or not RPL1 might be use when 
          collocated picture from P picture is a B picture */
        #if 1
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }
        #endif

        c0_col  = tmvp->dir_map_v0[c0_x + 1];
        c1_col  = tmvp->dir_map_v0[c1_x + 1];
        c0_col1 = tmvp->dir_map_v1[c0_x + 1];
        c1_col1 = tmvp->dir_map_v1[c1_x + 1];

        cand_c0  = !!(c0_col  & POS_MASK(pb_y, nb_pb_h));
        cand_c01 = !!(c0_col1 & POS_MASK(pb_y, nb_pb_h));
        cand_c1  = !!(c1_col  & POS_MASK(pb_y, nb_pb_h >> 1));
        cand_c11 = !!(c1_col1 & POS_MASK(pb_y, nb_pb_h >> 1));

        if (cand_c0) {
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            int scale = tmvp->scale00;
            #if 0
            OVMV c0 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c0 = tmvp->mvs0[pos_in_buff];
            #endif
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale, c0);
            cand[nb_cand] = c0;
            if (nb_cand++ == merge_idx)
                return c0;

        } else if (cand_c01) {
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            int scale = tmvp->scale01;
            #if 0
            OVMV c0 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c0 = tmvp->mvs1[pos_in_buff];
            #endif
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale, c0);
            cand[nb_cand] = c0;
            if (nb_cand++ == merge_idx)
                return c0;
        } else if (cand_c1) {
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            int scale = tmvp->scale00;
            #if 0
            OVMV c1 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c1 = tmvp->mvs0[pos_in_buff];
            #endif
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale , c1);
            cand[nb_cand] = c1;
            if (nb_cand++ == merge_idx)
                return c1;
        } else if (cand_c11) {
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            int scale = tmvp->scale01;
            #if 0
            OVMV c1 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            #else
            OVMV c1 = tmvp->mvs1[pos_in_buff];
            #endif
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale , c1);
            cand[nb_cand] = c1;
            if (nb_cand++ == merge_idx)
                return c1;
        }
    }
    #endif

    if (nb_cand != max_nb_merge_cand - 1) {
        uint8_t found;
        const struct HMVPLUT *hmvp_lut = &inter_ctx->hmvp_lut;

        found = hmvp_add_merge_cand(hmvp_lut, cand, cand_amvp, cand_t, cand_l, &nb_cand, 1, merge_idx);
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
        if (nb_cand == merge_idx)
            return avg_mv;
    }

    while (nb_cand < max_nb_merge_cand) {
        OVMV zmv = {0};
        cand[nb_cand++] = zmv;
    }

    return cand[nb_cand - 1];
}


static VVCMergeInfo
vvc_derive_merge_mvp_b(const struct InterDRVCtx *const inter_ctx,
                       uint8_t pb_x, uint8_t pb_y,
                       uint8_t nb_pb_w, uint8_t nb_pb_h,
                       uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    const struct OVMVCtx *mv_ctx0 = &inter_ctx->mv_ctx0;
    const struct OVMVCtx *mv_ctx1 = &inter_ctx->mv_ctx1;

    const OVMV *mv_buff0 = mv_ctx0->mvs;
    const OVMV *mv_buff1 = mv_ctx1->mvs;

    uint64_t lft_col0 = mv_ctx0->map.vfield[pb_x];
    uint64_t abv_row0 = mv_ctx0->map.hfield[pb_y];

    uint64_t lft_col1 = mv_ctx1->map.vfield[pb_x];
    uint64_t abv_row1 = mv_ctx1->map.hfield[pb_y];

    /*FIXME use flags for inter_dir and availability*/
    uint8_t cand_bl0 = !!(lft_col0 & POS_MASK(pb_y, nb_pb_h));     /*AO*/
    uint8_t cand_bl1 = !!(lft_col1 & POS_MASK(pb_y, nb_pb_h));     /*AO*/
    uint8_t cand_l0  = !!(lft_col0 & POS_MASK(pb_y, nb_pb_h - 1)); /*A1*/
    uint8_t cand_l1  = !!(lft_col1 & POS_MASK(pb_y, nb_pb_h - 1)); /*A1*/

    uint8_t cand_tr0 = !!(abv_row0 & POS_MASK(pb_x, nb_pb_w));     /*B0*/
    uint8_t cand_tr1 = !!(abv_row1 & POS_MASK(pb_x, nb_pb_w));     /*B0*/
    uint8_t cand_t0  = !!(abv_row0 & POS_MASK(pb_x, nb_pb_w - 1)); /*B1*/
    uint8_t cand_t1  = !!(abv_row1 & POS_MASK(pb_x, nb_pb_w - 1)); /*B1*/
    uint8_t cand_tl0 = !!(abv_row0 & POS_MASK(pb_x - 1, 0));      /*B2*/
    uint8_t cand_tl1 = !!(abv_row1 & POS_MASK(pb_x - 1, 0));      /*B2*/

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

    #if 1
    if (inter_ctx->tmvp_enabled) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        uint8_t cand_c0;
        uint8_t cand_c1;
        uint8_t cand_c01;
        uint8_t cand_c11;
        uint64_t c0_col;
        uint64_t c0_col1;
        uint64_t c1_col;
        uint64_t c1_col1;

        int c1_x = pb_x + (nb_pb_w >> 1);
        int c1_y = pb_y + (nb_pb_h >> 1);
        int c0_x = pb_x + nb_pb_w;
        int c0_y = pb_y + nb_pb_h;
        #if 1
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }
        #endif


        #if 0
        if (tmvp->col_ref == &ctudec->ref0) {
        #endif
            c0_col  = tmvp->dir_map_v0[c0_x + 1];
            c0_col1 = tmvp->dir_map_v1[c0_x + 1];
            c1_col  = tmvp->dir_map_v0[c1_x + 1];
            c1_col1 = tmvp->dir_map_v1[c1_x + 1];
            cand_c0  = !!(c0_col  & POS_MASK(pb_y, nb_pb_h));
            cand_c01 = !!(c0_col1 & POS_MASK(pb_y, nb_pb_h));
            cand_c1  = !!(c1_col  & POS_MASK(pb_y, nb_pb_h >> 1));
            cand_c11 = !!(c1_col1 & POS_MASK(pb_y, nb_pb_h >> 1));
        #if 0
        } else {
            uint64_t c0_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c0_x + 1];
            uint64_t c0_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c0_x + 1];
            uint64_t c1_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c1_x + 1];
            uint64_t c1_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c1_x + 1];
            cand_c0  = !!(c0_col  & POS_MASK(pb_y, nb_pb_h));
            cand_c01 = !!(c0_col1 & POS_MASK(pb_y, nb_pb_h));
            cand_c1  = !!(c1_col  & POS_MASK(pb_y, nb_pb_h >> 1));
            cand_c11 = !!(c1_col1 & POS_MASK(pb_y, nb_pb_h >> 1));
        }
        #endif

        /*FIXME check whether TMVP candidates RPL order is correct*/

        if (cand_c0 | cand_c01) {
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            cand[nb_cand].inter_dir = 3;
            if (cand_c0) {
                #if 0
                OVMV c0  = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
                #else
                OVMV c0  = tmvp->mvs0[pos_in_buff];
                #endif
                c0.x = tmvp_round_mv(c0.x);
                c0.y = tmvp_round_mv(c0.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale00, c0);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale10, c0);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            } else {
                #if 0
                OVMV c0  = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
                #else
                OVMV c0  = tmvp->mvs1[pos_in_buff];
                #endif
                c0.x = tmvp_round_mv(c0.x);
                c0.y = tmvp_round_mv(c0.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale01, c0);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale11, c0);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            }
        } else if (cand_c1 | cand_c11) {
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            cand[nb_cand].inter_dir = 3;
            if (cand_c1) {
                #if 0
                OVMV c1  = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
                #else
                OVMV c1  = tmvp->mvs0[pos_in_buff];
                #endif
                c1.x = tmvp_round_mv(c1.x);
                c1.y = tmvp_round_mv(c1.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale00, c1);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale10, c1);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            } else {
                #if 0
                OVMV c1  = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
                #else
                OVMV c1  = tmvp->mvs1[pos_in_buff];
                #endif
                c1.x = tmvp_round_mv(c1.x);
                c1.y = tmvp_round_mv(c1.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale01, c1);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale11, c1);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            }
        }
    }
    #endif

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
            avg_mv.mv0.x += cand[1].mv0.x + 1;
            avg_mv.mv0.y += cand[1].mv0.y + 1;
            avg_mv.mv0.x -= avg_mv.mv0.x >= 0;
            avg_mv.mv0.y -= avg_mv.mv0.y >= 0;
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
            avg_mv.mv1.x += cand[1].mv1.x + 1;
            avg_mv.mv1.y += cand[1].mv1.y + 1;
            avg_mv.mv1.x -= avg_mv.mv1.x >= 0;
            avg_mv.mv1.y -= avg_mv.mv1.y >= 0;
            avg_mv.mv1.x >>= 1;
            avg_mv.mv1.y >>= 1;
        } else if (cand[1].inter_dir & 0x2) {
            avg_mv.mv1 = cand[1].mv1;
            avg_mv.inter_dir |= 2;
            avg_mv.mv1.ref_idx = cand[1].mv1.ref_idx;
        } else if (cand[0].inter_dir & 0x2) {
            avg_mv.inter_dir |= 2;
        }

        if (nb_cand == merge_idx)
            return avg_mv;

        nb_cand++;
    }
    /*FIXME saturation or ridx based on nb_active_ref*/
    int8_t ridx = OVMIN(merge_idx - nb_cand, OVMIN(inter_ctx->nb_active_ref0, inter_ctx->nb_active_ref1));

    mv_z.inter_dir = 3;
    mv_z.mv0.ref_idx = ridx;
    mv_z.mv1.ref_idx = ridx;

    return mv_z;
}

static void
fill_mvp_map(struct OVMVCtx *const mv_ctx, OVMV mv,
             int pb_x, int pb_y, int nb_pb_w, int nb_pb_h)
{
    int i, j;

    ctu_field_set_rect_bitfield(&mv_ctx->map, pb_x, pb_y, nb_pb_w, nb_pb_h);

    for (j = 0; j < nb_pb_h; ++j) {
        for (i = 0; i < nb_pb_w; ++i) {
            memcpy(&mv_ctx->mvs[PB_POS_IN_BUF(pb_x + i, pb_y + j)], &mv, sizeof(OVMV));
        }
    }
}

/* FIXME DBF MV related */

#if 1
#define LF_MV_THRESHOLD 8
static void
fill_dbf_mv_map_b(struct DBFInfo *const dbf_info, struct OVMVCtx *const mv_ctx,
                  struct OVMVCtx *const mv_ctx1, OVMV mv,
                  int pb_x, int pb_y, int nb_pb_w, int nb_pb_h)
{
    int i, j;
    int log2_diff_min_cu = 1;
    int mask = (1 << (log2_diff_min_cu + 1)) - 1;
    int shift_v = 1 + (pb_y << log2_diff_min_cu);
    int shift_h = 2 + (pb_x << log2_diff_min_cu);

    uint64_t val = dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)];

    uint64_t tmp_mask_h = (uint64_t)mask << shift_h;
    uint64_t tmp_mask_v = (uint64_t)mask << shift_v;

    for (j = 0; j < nb_pb_w; ++j) {
        OVMV mv_above = mv_ctx->mvs[PB_POS_IN_BUF(pb_x + j, pb_y - 1)];
        int64_t above_avail = -((!!(mv_ctx->map.hfield[pb_y]  & POS_MASK(pb_x + j, 0))
                                & !(mv_ctx1->map.hfield[pb_y] & POS_MASK(pb_x + j, 0))));
        int64_t abv_th = -((abs(mv_above.x - mv.x) >= LF_MV_THRESHOLD) |
                           (abs(mv_above.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_h & abv_th & above_avail) | (tmp_mask_h & (-(!above_avail)));
        tmp_mask_h  <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)] |= val;

    val = dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)];

    for (i = 0; i < nb_pb_h; ++i) {
        OVMV mv_left = mv_ctx->mvs[PB_POS_IN_BUF(pb_x - 1, pb_y + i)];
        int64_t left_avail = -(!!(mv_ctx->map.vfield[pb_x]  & POS_MASK(pb_y + i, 0))
                              & !(mv_ctx1->map.vfield[pb_x] & POS_MASK(pb_y + i, 0)));
        int64_t abv_th = -((abs(mv_left.x - mv.x) >= LF_MV_THRESHOLD) |
                           (abs(mv_left.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_v & abv_th & left_avail) | (tmp_mask_v & (-(!left_avail)));
        tmp_mask_v <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)] |= val;
}

static void
fill_dbf_mv_map(struct DBFInfo *const dbf_info, struct OVMVCtx *const mv_ctx, OVMV mv,
                int pb_x, int pb_y, int nb_pb_w, int nb_pb_h)
{
    int i, j;
    int log2_diff_min_cu = 1;
    int mask = (1 << (log2_diff_min_cu + 1)) - 1;
    int shift_v = 1 + (pb_y << log2_diff_min_cu);
    int shift_h = 2 + (pb_x << log2_diff_min_cu);

    uint64_t val = dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)];

    uint64_t tmp_mask_h = (uint64_t)mask << shift_h;
    uint64_t tmp_mask_v = (uint64_t)mask << shift_v;
    for (j = 0; j < nb_pb_w; ++j) {
        OVMV mv_above = mv_ctx->mvs[PB_POS_IN_BUF(pb_x + j, pb_y - 1)];
        int64_t above_avail = -(!!(mv_ctx->map.hfield[pb_y] & POS_MASK(pb_x + j, 0)));
        int64_t abv_th = -((abs(mv_above.x - mv.x) >= LF_MV_THRESHOLD) |
                           (abs(mv_above.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_h & abv_th & above_avail) | (tmp_mask_h & (-(!above_avail)));
        tmp_mask_h  <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)] |= val;

    val = dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)];

    for (i = 0; i < nb_pb_h; ++i) {
        OVMV mv_left = mv_ctx->mvs[PB_POS_IN_BUF(pb_x - 1, pb_y + i)];
        int64_t left_avail = -(!!(mv_ctx->map.vfield[pb_x] & POS_MASK(pb_y + i, 0)));
        int64_t abv_th = -((abs(mv_left.x - mv.x) >= LF_MV_THRESHOLD) |
                           (abs(mv_left.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_v & abv_th & left_avail) | (tmp_mask_v & (-(!left_avail)));
        tmp_mask_v <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)] |= val;
}
#endif

static void
update_mv_ctx_b(struct InterDRVCtx *const inter_ctx,
                const OVMV mv0, const OVMV mv1,
                uint8_t pb_x, uint8_t  pb_y,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t inter_dir)
{
    #if 1
    /*FIXME Use specific DBF update function if DBF is disabled */
    /*FIXME Find a better way to retrieve dbf_info */
    struct DBFInfo *const dbf_info = &inter_ctx->tmvp_ctx.ctudec->dbf_info;
    #endif
    if (inter_dir == 3) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        #if 1
        fill_dbf_mv_map(dbf_info, mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        fill_dbf_mv_map(dbf_info, mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);
        #endif

    } else if (inter_dir & 0x2) {
        #if 1
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #endif
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx1, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        #if 1
        fill_dbf_mv_map_b(dbf_info, mv_ctx1, mv_ctx0, mv1, pb_x, pb_y, nb_pb_w, nb_pb_h);
        #endif

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #if 1
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        #endif

        fill_mvp_map(mv_ctx0, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        #if 1
        fill_dbf_mv_map_b(dbf_info, mv_ctx0, mv_ctx1, mv0, pb_x, pb_y, nb_pb_w, nb_pb_h);
        #endif

    }

    hmvp_update_lut_b(&inter_ctx->hmvp_lut, mv0, mv1, inter_dir);
}

static void
update_mv_ctx(struct InterDRVCtx *const inter_ctx,
              const OVMV mv,
              uint8_t pb_x, uint8_t  pb_y,
              uint8_t nb_pb_w, uint8_t nb_pb_h,
              uint8_t inter_dir)
{
    /*FIXME Use specific DBF update function if DBF is disabled */
    #if 1
    /*FIXME Use specific DBF update function if DBF is disabled */
    /*FIXME Find a better way to retrieve dbf_info */
    struct DBFInfo *const dbf_info = &inter_ctx->tmvp_ctx.ctudec->dbf_info;
    #endif
    if (inter_dir & 0x2) {
        #if 1
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #endif
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx1, mv, pb_x, pb_y, nb_pb_w, nb_pb_h);

        #if 1
        fill_dbf_mv_map_b(dbf_info, mv_ctx1, mv_ctx0, mv, pb_x, pb_y, nb_pb_w, nb_pb_h);
        #endif

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #if 1
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        #endif

        fill_mvp_map(mv_ctx0, mv, pb_x, pb_y, nb_pb_w, nb_pb_h);

        #if 1
        fill_dbf_mv_map_b(dbf_info, mv_ctx0, mv_ctx1, mv, pb_x, pb_y, nb_pb_w, nb_pb_h);
        #endif

    }

    hmvp_update_lut(&inter_ctx->hmvp_lut, mv);
}

/* Derive motion vectors and update motion maps */
VVCMergeInfo
drv_mvp_b(struct InterDRVCtx *const inter_ctx,
          uint8_t pb_x, uint8_t pb_y,
          uint8_t nb_pb_w, uint8_t nb_pb_h,
          OVMV mvd0, OVMV mvd1,
          uint8_t mvp_idx0, uint8_t mvp_idx1,
          uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1)
{
    OVMV mv0 = {0}, mv1 = {0};
    VVCMergeInfo mv_info;

    uint8_t opp_ref_idx0 = 0xFF;
    uint8_t opp_ref_idx1 = 0xFF;

    //uint8_t same_ref = inter_ctx->rpl0[ref_idx0] == inter_ctx->rpl1[ref_idx1];
    for (int i = 0; i < inter_ctx->nb_active_ref1; i ++) {
         if (inter_ctx->rpl0[ref_idx0] == inter_ctx->rpl1[i])
             opp_ref_idx0 = i;
    }

    for (int i = 0; i < inter_ctx->nb_active_ref0; i ++) {
         if (inter_ctx->rpl1[ref_idx1] == inter_ctx->rpl0[i])
             opp_ref_idx1 = i;
    }

    /* FIXME can we combine mvp derivation for bi pred */
    if (1) {
        if (inter_dir & 0x1) {
            struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
            struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

            mv0 = derive_mvp_candidates_1(inter_ctx, mv_ctx0, ref_idx0,
                                          pb_x, pb_y, nb_pb_w, nb_pb_h,
                                          mvp_idx0, inter_dir & 0x1, mv_ctx1, opp_ref_idx0);

            mvd0 = scale_mvd(mvd0);

            mv0.x += mvd0.x;
            mv0.y += mvd0.y;
            mv0.ref_idx = ref_idx0;
        }

        if (inter_dir & 0x2) {
            struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
            struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

            mv1 = derive_mvp_candidates_1(inter_ctx, mv_ctx1, ref_idx1,
                                          pb_x, pb_y,
                                          nb_pb_w, nb_pb_h,
                                          mvp_idx1, inter_dir & 0x2, mv_ctx0, opp_ref_idx1);

            mvd1 = scale_mvd(mvd1);

            mv1.x += mvd1.x;
            mv1.y += mvd1.y;
            mv1.ref_idx = ref_idx1;
        }
    } else {
        if (inter_dir & 0x1) {
            struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

            mv0 = derive_mvp_candidates(inter_ctx, mv_ctx0,
                                        pb_x, pb_y, nb_pb_w, nb_pb_h,
                                        mvp_idx0, inter_dir & 0x1);

            mvd0 = scale_mvd(mvd0);

            mv0.x += mvd0.x;
            mv0.y += mvd0.y;
        }

        if (inter_dir & 0x2) {
            struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

            mv1 = derive_mvp_candidates(inter_ctx, mv_ctx1,
                                        pb_x, pb_y,
                                        nb_pb_w, nb_pb_h,
                                        mvp_idx1, inter_dir & 0x2);

            mvd1 = scale_mvd(mvd1);

            mv1.x += mvd1.x;
            mv1.y += mvd1.y;
        }

    }

    mv_info.inter_dir = inter_dir;
    mv_info.mv0 = mv0;
    mv_info.mv1 = mv1;

    /* Update for next pass */
    update_mv_ctx_b(inter_ctx, mv0, mv1, pb_x, pb_y, nb_pb_w,
                    nb_pb_h, inter_dir);
    return mv_info;
}

OVMV
drv_merge_mvp(struct InterDRVCtx *const inter_ctx,
              const struct OVMVCtx *const mv_ctx,
              uint8_t pb_x, uint8_t pb_y,
              uint8_t nb_pb_w, uint8_t nb_pb_h,
              uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    OVMV mv0 = vvc_derive_merge_mvp(inter_ctx, mv_ctx, pb_x, pb_y,
                                    nb_pb_w, nb_pb_h, merge_idx,
                                    max_nb_merge_cand);

    update_mv_ctx(inter_ctx, mv0, pb_x, pb_y, nb_pb_w,
                  nb_pb_h, 1);
    return mv0;
}

OVMV
drv_mvp_mvd(struct InterDRVCtx *const inter_ctx,
            const struct OVMVCtx *const mv_ctx,
            OVMV mvd,
            uint8_t pb_x, uint8_t pb_y,
            uint8_t nb_pb_w, uint8_t nb_pb_h,
            uint8_t mvp_idx, uint8_t inter_dir,
            uint8_t ref_idx0, uint8_t ref_idx1)
{
    OVMV mv;
    //uint8_t same_ref = inter_ctx->rpl0[ref_idx0] == inter_ctx->rpl1[ref_idx1];
    uint8_t opp_ref_idx0 = 0xFF;
    uint8_t opp_ref_idx1 = 0xFF;

    //uint8_t same_ref = inter_ctx->rpl0[ref_idx0] == inter_ctx->rpl1[ref_idx1];
    for (int i = 0; i < inter_ctx->nb_active_ref1; i ++) {
         if (inter_ctx->rpl0[ref_idx0] == inter_ctx->rpl1[i])
             opp_ref_idx0 = i;
    }

    for (int i = 0; i < inter_ctx->nb_active_ref0; i ++) {
         if (inter_ctx->rpl1[ref_idx1] == inter_ctx->rpl0[i])
             opp_ref_idx1 = i;
    }

    if (1) {

        const struct OVMVCtx *const mv_ctx_opp = &inter_ctx->mv_ctx0 == mv_ctx ? &inter_ctx->mv_ctx1 :
                                                 &inter_ctx->mv_ctx0;

        uint8_t ref_idx = &inter_ctx->mv_ctx0 == mv_ctx  ? ref_idx0 : ref_idx1;
        uint8_t ref_idx_opp = &inter_ctx->mv_ctx0 == mv_ctx  ? opp_ref_idx0 : opp_ref_idx1;

        mv = derive_mvp_candidates_1(inter_ctx, mv_ctx, ref_idx,
                                     pb_x, pb_y, nb_pb_w, nb_pb_h,
                                     mvp_idx, 1, mv_ctx_opp, ref_idx_opp);
        mv.ref_idx = ref_idx;
    } else {
        mv = derive_mvp_candidates(inter_ctx, mv_ctx,
                                   pb_x, pb_y, nb_pb_w, nb_pb_h,
                                   mvp_idx, 1);
    }
    mvd = scale_mvd(mvd);

    mv.x += mvd.x;
    mv.y += mvd.y;

#if 1
    update_mv_ctx(inter_ctx, mv, pb_x, pb_y, nb_pb_w,
                  nb_pb_h, inter_dir);
#endif
   return mv;
}

VVCMergeInfo
drv_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                uint8_t pb_x, uint8_t pb_y,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t merge_idx,
                uint8_t max_nb_cand)
{
    VVCMergeInfo mv_info;

    mv_info = vvc_derive_merge_mvp_b(inter_ctx, pb_x, pb_y,
                                     nb_pb_w, nb_pb_h, merge_idx,
                                     max_nb_cand);

    update_mv_ctx_b(inter_ctx, mv_info.mv0, mv_info.mv1, pb_x, pb_y,
                    nb_pb_w, nb_pb_h, mv_info.inter_dir);

    return mv_info;
}
