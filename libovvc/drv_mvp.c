#include <stdint.h>
#include <string.h>

#include "ovdefs.h"
#include "ovutils.h"

/* FIXME find a place for mv structures */
#include "ctudec.h"
#include "drv_utils.h"
#include "dec_structures.h"

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
        if (a->inter_dir & 0x1)
            is_neq += !MV_CMP(a->mv0, b->mv0);
        if (a->inter_dir & 0x2)
            is_neq += !MV_CMP(a->mv1, b->mv1);

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
                      int *const nb_cand, uint8_t mvp_idx)
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

            if (*nb_cand == 4) {
                return 0;
            }
        }
    }
    return 0;
}


void
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

void
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
                    duplicated_mv = MV_CMP(mv0, hmvp_lut->hmv0[i]);
                    break;
                case 0x2:
                    duplicated_mv = MV_CMP(mv1, hmvp_lut->hmv1[i]);
                    break;
                case 0x3:
                    duplicated_mv = MV_CMP(mv0, hmvp_lut->hmv0[i]) &&
                        MV_CMP(mv1, hmvp_lut->hmv1[i]);
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


#if 0
static struct VVCTMVPStorage
derive_tmvp_ctx_storage(const struct VVCTMVPStorageinfo *const info, const uint8_t *const data)
{
    struct VVCTMVPStorage storage;
    storage.col0 = data;
    storage.col1 = data + info->tmvp_map_s;
    storage.mv0  = data + (info->tmvp_map_s << 1);
    storage.mv1  = data + (info->tmvp_map_s << 1) + info->tmvp_mvs_s;
    return storage;
}

static void
load_ctb_tmvp(const VVCLocalContext *const lc_ctx, int ctb_x, int ctb_y)
{
    VVCContext *const vvc_ctx = (VVCContext *) lc_ctx->parent;
    struct InterDRVCtx *const inter_ctx = &lc_ctx->inter_ctx;
    struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
    const struct VVCTMVPStorageinfo const *size_info = &tmvp->tmvp_size;
    int nb_pb = 1 << size_info->log2_nb_ctb_pb;
    const int ctb_offset = ctb_y * vvc_ctx->nb_ctu_w + ctb_x;

    const uint8_t *const ctb_tmvp = (uint8_t *)tmvp->data_ref + (ctb_offset * size_info->tmvp_s);

    struct VVCTMVPStorage storage = derive_tmvp_ctx_storage(size_info, ctb_tmvp);

    struct OVMVCtx *const mv_ctx0 = &tmvp->tmvp_mv.mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &tmvp->tmvp_mv.mv_ctx1;

    if ((vvc_ctx->threads_type & FF_THREAD_FRAME)) {
        ff_thread_await_progress(&tmvp->col_ref->tf, lc_ctx->ctb_y, 0);
    }

    OVMV *dst0 = mv_ctx0->mvs + 35;
    OVMV *dst1 = mv_ctx1->mvs + 35;

    memcpy(mv_ctx0->map.vfields + 1, storage.col0, size_info->tmvp_map_s);
    memcpy(mv_ctx1->map.vfields + 1, storage.col1, size_info->tmvp_map_s);

    for (int i = 0; i < (1 << size_info->log2_nb_ctb_pb); ++i) {
        memcpy(dst0, storage.mv0, (sizeof(OVMV) << size_info->log2_nb_ctb_pb));
        memcpy(dst1, storage.mv1, (sizeof(OVMV) << size_info->log2_nb_ctb_pb));
        storage.mv0 += 1 << size_info->log2_nb_ctb_pb;
        storage.mv1 += 1 << size_info->log2_nb_ctb_pb;
        dst0 += 34;
        dst1 += 34;
    }

    if (lc_ctx->ctb_x != vvc_ctx->nb_ctu_w - 1) {
        int i;
        uint8_t *const ctb_tmvp_r = ctb_tmvp + size_info->tmvp_s;
        storage = derive_tmvp_ctx_storage(size_info, ctb_tmvp_r);

        mv_ctx0->map.vfield[nb_pb + 1] = storage.col0[0];
        mv_ctx1->map.vfield[nb_pb + 1] = storage.col1[0];

        OVMV *mv_right0 = storage.mv0;
        OVMV *mv_right1 = storage.mv1;

        for (i = 0; i < nb_pb + (lc_ctx->ctb_y != vvc_ctx->nb_ctu_h - 1); ++i) {
            int pos_offset_l = PB_POS_IN_BUF(nb_pb, i);
            mv_ctx0->mvs[pos_offset_l] = *mv_right0;
            mv_ctx1->mvs[pos_offset_l] = *mv_right1;
            mv_right0 += 1 << size_info->log2_nb_ctb_pb;
            mv_right1 += 1 << size_info->log2_nb_ctb_pb;
        }

    } else {
        mv_ctx0->map.vfield[nb_pb + 1] = 0;
        mv_ctx1->map.vfield[nb_pb + 1] = 0;
    }
    inter_ctx->tmvp_avail |= 1;
}
#endif

static inline OVMV
tmvp_scale_mv(int scale, OVMV mv)
{
    mv.x = ov_clip((scale * mv.x + 128 - (scale * mv.x >= 0)) >> 8, MV_MIN, MV_MAX);
    mv.y = ov_clip((scale * mv.y + 128 - (scale * mv.y >= 0)) >> 8, MV_MIN, MV_MAX);
    return mv;
}

OVMV
derive_mvp_candidates(struct InterDRVCtx *const inter_ctx,
                      const struct OVMVCtx *const mv_ctx,
                      uint8_t pb_x, uint8_t pb_y,
                      uint8_t n_pb_w, uint8_t n_pb_h,
                      uint8_t mvp_idx, uint8_t inter_dir)
{
    const OVMV *const mv_buff = mv_ctx->mvs;
    uint64_t lft_col = mv_ctx->map.vfield[pb_x];
    uint64_t abv_row = mv_ctx->map.hfield[pb_y];

    OVMV cand[2];
    int nb_cand = 0;

    /* Derive candidates availability based on CTU inter fields */
    uint8_t cand_bl = !!(lft_col & POS_MASK(pb_y, n_pb_h));     /*AO*/
    uint8_t cand_l  = !!(lft_col & POS_MASK(pb_y, n_pb_h - 1)); /*A1*/
    uint8_t cand_tr = !!(abv_row & POS_MASK(pb_x, n_pb_w));     /*B0*/
    uint8_t cand_t  = !!(abv_row & POS_MASK(pb_x, n_pb_w - 1)); /*B1*/
    uint8_t cand_tl = !!(abv_row & POS_MASK(pb_x - 1, 0));      /*B2*/

    if (cand_bl) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + n_pb_h);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    } else if (cand_l) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + n_pb_h - 1);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    }

    if (cand_tr) {
        int pos_in_buff = OFFSET_BUFF(pb_x + n_pb_w, pb_y - 1);
        cand[nb_cand++] = mv_buff[pos_in_buff];
    } else if (cand_t) {
        int pos_in_buff = OFFSET_BUFF(pb_x + n_pb_w - 1, pb_y - 1);
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
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        uint64_t c1_col;
        uint64_t c0_col;
        uint64_t c1_col1;
        uint64_t c0_col1;
        uint8_t cand_c0;
        uint8_t cand_c1;
        uint8_t cand_c01;
        uint8_t cand_c11;
        int c1_x = pb_x + (n_pb_w >> 1);
        int c1_y = pb_y + (n_pb_h >> 1);
        int c0_x = pb_x + n_pb_w;
        int c0_y = pb_y + n_pb_h;
        int scale0, scale1;
        #if 0
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            load_ctb_tmvp(lc_ctx, lc_ctx->ctb_x, lc_ctx->ctb_y);
        }
        #endif

        /* Derive availability based on CTB inter fields */
        c0_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c0_x + 1];
        c0_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c0_x + 1];
        c1_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c1_x + 1];
        c1_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c1_x + 1];

        cand_c0  = !!(c0_col  & POS_MASK(pb_y, n_pb_h));
        cand_c01 = !!(c0_col1 & POS_MASK(pb_y, n_pb_h));
        cand_c1  = !!(c1_col  & POS_MASK(pb_y, n_pb_h >> 1));
        cand_c11 = !!(c1_col1 & POS_MASK(pb_y, n_pb_h >> 1));

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
            OVMV c0 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale0, c0);
            c0.x = ((c0.x + 2 - (c0.x >= 0)) >> 2) << 2;
            c0.y = ((c0.y + 2 - (c0.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c0;
        } else if (cand_c01) {
            /* Candidate 0 in collocated picture 1 */
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            OVMV c0 = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale1, c0);
            c0.x = ((c0.x + 2 - (c0.x >= 0)) >> 2) << 2;
            c0.y = ((c0.y + 2 - (c0.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c0;
        } else if (cand_c1) {
            /* Candidate 1 in collocated picture 0 */
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            OVMV c1 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            c1.x = tmvp_round_mv(c1.x);
            c1.y = tmvp_round_mv(c1.y);
            c1 = tmvp_scale_mv(scale0, c1);
            c1.x = ((c1.x + 2 - (c1.x >= 0)) >> 2) << 2;
            c1.y = ((c1.y + 2 - (c1.y >= 0)) >> 2) << 2;
            cand[nb_cand++] = c1;
        } else if (cand_c11) {
            /* Candidate 1 in collocated picture 1 */
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            OVMV c1 = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
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

OVMV
vvc_derive_merge_mvp(const struct InterDRVCtx *const inter_ctx,
                     const struct OVMVCtx *const mv_ctx,
                     uint8_t pb_x, uint8_t pb_y,
                     uint8_t n_pb_w, uint8_t n_pb_h,
                     uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    const OVMV *const mv_buff = mv_ctx->mvs;

    uint64_t lft_col = mv_ctx->map.vfield[pb_x];
    uint64_t abv_row = mv_ctx->map.hfield[pb_y];

    uint8_t cand_bl = !!(lft_col & POS_MASK(pb_y, n_pb_h));     /*AO*/
    uint8_t cand_l  = !!(lft_col & POS_MASK(pb_y, n_pb_h - 1)); /*A1*/
    uint8_t cand_tr = !!(abv_row & POS_MASK(pb_x, n_pb_w));     /*B0*/
    uint8_t cand_t  = !!(abv_row & POS_MASK(pb_x, n_pb_w - 1)); /*B1*/
    uint8_t cand_tl = !!(abv_row & POS_MASK(pb_x - 1, 0));      /*B2*/
    OVMV cand[6];
    OVMV cand_amvp[5];
    OVMV mv_z = {-1};

    int nb_cand = 0;

    cand_amvp[0] = mv_z;
    if (cand_t) { /* B1 */
        int pos_in_buff = OFFSET_BUFF(pb_x + n_pb_w - 1, pb_y - 1);
        OVMV mv_B1 = mv_buff[pos_in_buff];
        cand_amvp[0] = mv_buff[pos_in_buff];
        cand[nb_cand] = mv_B1;
        if (nb_cand++ == merge_idx)
            return mv_B1;
    }

    cand_amvp[1] = mv_z;
    if (cand_l) {
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + n_pb_h - 1);
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
        int pos_in_buff = OFFSET_BUFF(pb_x + n_pb_w, pb_y - 1);
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
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + n_pb_h);
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
            cand_amvp[4] = mv_buff[pos_in_buff];
            OVMV mv_B2 = mv_buff[pos_in_buff];
            if ((!cand_l || !MV_CMP(mv_B2, cand_amvp[1])) &&
                (!cand_t || !MV_CMP(mv_B2, cand_amvp[0]))) {
                cand[nb_cand] = mv_B2;
                if (nb_cand++ == merge_idx)
                    return mv_B2;
            }
        }
    }


    #if 0
    if (inter_ctx->tmvp_enabled) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        uint64_t c1_col;
        uint64_t c0_col;
        uint8_t cand_c0;
        uint8_t cand_c1;
        int c1_x = pb_x + (n_pb_w >> 1);
        int c1_y = pb_y + (n_pb_h >> 1);
        int c0_x = pb_x + n_pb_w;
        int c0_y = pb_y + n_pb_h;
        /*FIXME determine whether or not RPL1 might be use when 
          collocated picture from P picture is a B picture */
        #if 0
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            load_ctb_tmvp(lc_ctx, lc_ctx->ctb_x, lc_ctx->ctb_y);
        }
        #endif

        c0_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c0_x + 1];
        c1_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c1_x + 1];

        cand_c0  = !!(c0_col  & POS_MASK(pb_y, n_pb_h));
        cand_c1  = !!(c1_col  & POS_MASK(pb_y, n_pb_h >> 1));

        if (cand_c0) {
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            int scale = tmvp->scale00;
            OVMV c0 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
            c0.x = tmvp_round_mv(c0.x);
            c0.y = tmvp_round_mv(c0.y);
            c0 = tmvp_scale_mv(scale, c0);
            cand[nb_cand] = c0;
            if (nb_cand++ == merge_idx)
                return c0;

        } else if (cand_c1) {
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            int scale = tmvp->scale00;
            OVMV c1 = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
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


VVCMergeInfo
vvc_derive_merge_mvp_b(const struct InterDRVCtx *const inter_ctx,
                       uint8_t pb_x, uint8_t pb_y,
                       uint8_t n_pb_w, uint8_t n_pb_h,
                       uint8_t merge_idx, uint8_t max_nb_merge_cand)
{
    const OVMV *mv_buff0 = inter_ctx->mv_ctx0.mvs;
    const OVMV *mv_buff1 = inter_ctx->mv_ctx1.mvs;
    uint64_t lft_col0  = inter_ctx->mv_ctx0.map.vfield[pb_x];
    uint64_t abv_row0 = inter_ctx->mv_ctx0.map.hfield[pb_y];
    uint64_t lft_col1  = inter_ctx->mv_ctx1.map.vfield[pb_x];
    uint64_t abv_row1 = inter_ctx->mv_ctx1.map.hfield[pb_y];

    /*FIXME use flags for inter_dir and availability*/
    uint8_t cand_bl0 = !!(lft_col0 & POS_MASK(pb_y, n_pb_h));     /*AO*/
    uint8_t cand_bl1 = !!(lft_col1 & POS_MASK(pb_y, n_pb_h));     /*AO*/
    uint8_t cand_l0  = !!(lft_col0 & POS_MASK(pb_y, n_pb_h - 1)); /*A1*/
    uint8_t cand_l1  = !!(lft_col1 & POS_MASK(pb_y, n_pb_h - 1)); /*A1*/

    uint8_t cand_tr0 = !!(abv_row0 & POS_MASK(pb_x, n_pb_w));     /*B0*/
    uint8_t cand_tr1 = !!(abv_row1 & POS_MASK(pb_x, n_pb_w));     /*B0*/
    uint8_t cand_t0  = !!(abv_row0 & POS_MASK(pb_x, n_pb_w - 1)); /*B1*/
    uint8_t cand_t1  = !!(abv_row1 & POS_MASK(pb_x, n_pb_w - 1)); /*B1*/
    uint8_t cand_tl0 = !!(abv_row0 & POS_MASK(pb_x - 1, 0));      /*B2*/
    uint8_t cand_tl1 = !!(abv_row1 & POS_MASK(pb_x - 1, 0));      /*B2*/

    VVCMergeInfo cand[6];
    VVCMergeInfo cand_amvp[5];

    VVCMergeInfo mv_z = { .mv0 = {0}, .mv1 = {0}, .inter_dir = 3};

    int nb_cand = 0;

    cand_amvp[0] = mv_z;
    if (cand_t0 | cand_t1) { /* B1 */
        int pos_in_buff = OFFSET_BUFF(pb_x + n_pb_w - 1, pb_y - 1);
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
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + n_pb_h - 1);
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
        int pos_in_buff = OFFSET_BUFF(pb_x + n_pb_w, pb_y - 1);
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
        int pos_in_buff = OFFSET_BUFF(pb_x - 1, pb_y + n_pb_h);
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

    #if 0
    if (inter_ctx->tmvp_enabled) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        #if 0
        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            load_ctb_tmvp(lc_ctx, lc_ctx->ctb_x, lc_ctx->ctb_y);
        }
        #endif

        uint8_t cand_c0;
        uint8_t cand_c1;
        uint8_t cand_c01;
        uint8_t cand_c11;

        int c1_x = pb_x + (n_pb_w >> 1);
        int c1_y = pb_y + (n_pb_h >> 1);
        int c0_x = pb_x + n_pb_w;
        int c0_y = pb_y + n_pb_h;

        #if 0
        if (tmvp->col_ref == &lc_ctx->ref0) {
        #endif
            uint64_t c0_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c0_x + 1];
            uint64_t c0_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c0_x + 1];
            uint64_t c1_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c1_x + 1];
            uint64_t c1_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c1_x + 1];
            cand_c0  = !!(c0_col  & POS_MASK(pb_y, n_pb_h));
            cand_c01 = !!(c0_col1 & POS_MASK(pb_y, n_pb_h));
            cand_c1  = !!(c1_col  & POS_MASK(pb_y, n_pb_h >> 1));
            cand_c11 = !!(c1_col1 & POS_MASK(pb_y, n_pb_h >> 1));
        #if 0
        } else {
            uint64_t c0_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c0_x + 1];
            uint64_t c0_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c0_x + 1];
            uint64_t c1_col  = tmvp->tmvp_mv.mv_ctx0.map.vfield[c1_x + 1];
            uint64_t c1_col1 = tmvp->tmvp_mv.mv_ctx1.map.vfield[c1_x + 1];
            cand_c0  = !!(c0_col  & POS_MASK(pb_y, n_pb_h));
            cand_c01 = !!(c0_col1 & POS_MASK(pb_y, n_pb_h));
            cand_c1  = !!(c1_col  & POS_MASK(pb_y, n_pb_h >> 1));
            cand_c11 = !!(c1_col1 & POS_MASK(pb_y, n_pb_h >> 1));
        }
        #endif

        /*FIXME check whether TMVP candidates RPL order is correct*/

        if (cand_c0 | cand_c01) {
            cand[nb_cand].inter_dir = 3;
            int pos_in_buff = PB_POS_IN_BUF(c0_x, c0_y);
            if (cand_c0) {
                OVMV c0  = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
                c0.x = tmvp_round_mv(c0.x);
                c0.y = tmvp_round_mv(c0.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale00, c0);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale10, c0);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            } else {
                OVMV c0  = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
                c0.x = tmvp_round_mv(c0.x);
                c0.y = tmvp_round_mv(c0.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale01, c0);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale11, c0);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            }
        } else if (cand_c1 | cand_c11) {
            cand[nb_cand].inter_dir = 3;
            int pos_in_buff = PB_POS_IN_BUF(c1_x, c1_y);
            if (cand_c1) {
                OVMV c1  = tmvp->tmvp_mv.mv_ctx0.mvs[pos_in_buff];
                c1.x = tmvp_round_mv(c1.x);
                c1.y = tmvp_round_mv(c1.y);
                cand[nb_cand].mv0 = tmvp_scale_mv(tmvp->scale00, c1);
                cand[nb_cand].mv1 = tmvp_scale_mv(tmvp->scale10, c1);
                if (nb_cand++ == merge_idx)
                    return cand[merge_idx];
            } else {
                OVMV c1  = tmvp->tmvp_mv.mv_ctx1.mvs[pos_in_buff];
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

        found = hmvp_add_merge_cand_b(hmvp_lut, cand, cand_amvp, cand_t0 | cand_t1, cand_l0 | cand_l1, &nb_cand, merge_idx);
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
        } else if (cand[0].inter_dir & 0x2) {
            avg_mv.inter_dir |= 2;
        }

        if (nb_cand == merge_idx)
            return avg_mv;
    }

    return mv_z;
}

void
fill_mvp_map(struct OVMVCtx *const mv_ctx, OVMV mv,
             int pb_x, int pb_y, int n_pb_w, int n_pb_h)
{
    int i, j;

    ctu_field_set_rect_bitfield(&mv_ctx->map, pb_x, pb_y, n_pb_w, n_pb_h);

    for (j = 0; j < n_pb_h; ++j) {
        for (i = 0; i < n_pb_w; ++i) {
            memcpy(&mv_ctx->mvs[PB_POS_IN_BUF(pb_x + i, pb_y + j)], &mv, sizeof(OVMV));
        }
    }
}

/* FIXME
 * DBF MV related
 */

#if 0
#define LF_MV_THRESHOLD 8
static void
fill_dbf_mv_map_b(struct DBFInfo *const dbf_info, struct OVMVCtx *const mv_ctx, struct OVMVCtx *const mv_ctx1, OVMV mv,
                  int pb_x, int pb_y, int n_pb_w, int n_pb_h)
{
    int i, j;
    int log2_diff_min_cu = 1;
    int mask = (1 << (log2_diff_min_cu + 1)) - 1;
    int shift_v = 1 + (pb_y << log2_diff_min_cu);
    int shift_h = 2 + (pb_x << log2_diff_min_cu);

    uint64_t val = dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)];

    uint64_t tmp_mask_h = (uint64_t)mask << shift_h;
    uint64_t tmp_mask_v = (uint64_t)mask << shift_v;

    for (j = 0; j < n_pb_w; ++j) {
        OVMV mv_above = mv_ctx->mvs[PB_POS_IN_BUF(pb_x + j, pb_y - 1)];
        int64_t above_avail = -((!!(mv_ctx->map.hfield[pb_y]  & POS_MASK(pb_x + j, 0))
                                & !(mv_ctx1->map.hfield[pb_y] & POS_MASK(pb_x + j, 0))));
        int64_t abv_th = -((FFABS(mv_above.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_above.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_h & abv_th & above_avail) | (tmp_mask_h & (-(!above_avail)));
        tmp_mask_h  <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)] |= val;

    val = dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)];

    for (i = 0; i < n_pb_h; ++i) {
        OVMV mv_left = mv_ctx->mvs[PB_POS_IN_BUF(pb_x - 1, pb_y + i)];
        int64_t left_avail = -(!!(mv_ctx->map.vfield[pb_x]  & POS_MASK(pb_y + i, 0))
                              & !(mv_ctx1->map.vfield[pb_x] & POS_MASK(pb_y + i, 0)));
        int64_t abv_th = -((FFABS(mv_left.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_left.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_v & abv_th & left_avail) | (tmp_mask_v & (-(!left_avail)));
        tmp_mask_v <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)] |= val;
}

static void
fill_dbf_mv_map(struct DBFInfo *const dbf_info, struct OVMVCtx *const mv_ctx, OVMV mv,
                int pb_x, int pb_y, int n_pb_w, int n_pb_h)
{
    int i, j;
    int log2_diff_min_cu = 1;
    int mask = (1 << (log2_diff_min_cu + 1)) - 1;
    int shift_v = 1 + (pb_y << log2_diff_min_cu);
    int shift_h = 2 + (pb_x << log2_diff_min_cu);

    uint64_t val = dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)];

    uint64_t tmp_mask_h = (uint64_t)mask << shift_h;
    uint64_t tmp_mask_v = (uint64_t)mask << shift_v;
    for (j = 0; j < n_pb_w; ++j) {
        OVMV mv_above = mv_ctx->mvs[PB_POS_IN_BUF(pb_x + j, pb_y - 1)];
        int64_t above_avail = -(!!(mv_ctx->map.hfield[pb_y] & POS_MASK(pb_x + j, 0)));
        int64_t abv_th = -((FFABS(mv_above.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_above.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_h & abv_th & above_avail) | (tmp_mask_h & (-(!above_avail)));
        tmp_mask_h  <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.hor[(pb_y << log2_diff_min_cu)] |= val;

    val = dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)];

    for (i = 0; i < n_pb_h; ++i) {
        OVMV mv_left = mv_ctx->mvs[PB_POS_IN_BUF(pb_x - 1, pb_y + i)];
        int64_t left_avail = -(!!(mv_ctx->map.vfield[pb_x] & POS_MASK(pb_y + i, 0)));
        int64_t abv_th = -((FFABS(mv_left.x - mv.x) >= LF_MV_THRESHOLD) |
                           (FFABS(mv_left.y - mv.y) >= LF_MV_THRESHOLD));
        val |= (tmp_mask_v & abv_th & left_avail) | (tmp_mask_v & (-(!left_avail)));
        tmp_mask_v <<= (1 << log2_diff_min_cu);
    }
    dbf_info->bs1_map.ver[(pb_x << log2_diff_min_cu)] |= val;
}
#endif

void
update_mv_ctx_b(struct InterDRVCtx *const inter_ctx,
                const OVMV mv0, const OVMV mv1,
                uint8_t x_pu, uint8_t  y_pu,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t inter_dir)
{
    /*FIXME Use specific DBF update function if DBF is disabled */
    #if 0
    struct DBFInfo *const dbf_info = &lc_ctx->dbf_info;
    #endif
    if (inter_dir == 3) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_mvp_map(mv_ctx1, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);

        #if 0
        fill_dbf_mv_map(dbf_info, mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

        fill_dbf_mv_map(dbf_info, mv_ctx1, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);
        #endif

    } else if (inter_dir & 0x2) {
        #if 0
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #endif
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx1, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);

        #if 0
        fill_dbf_mv_map_b(dbf_info, mv_ctx1, mv_ctx0, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);
        #endif

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #if 0
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        #endif

        fill_mvp_map(mv_ctx0, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);

        #if 0
        fill_dbf_mv_map_b(dbf_info, mv_ctx0, mv_ctx1, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);
        #endif

    }

    hmvp_update_lut_b(&inter_ctx->hmvp_lut, mv0, mv1, inter_dir);
}

void
update_mv_ctx(struct InterDRVCtx *const inter_ctx,
              const OVMV mv,
              uint8_t x_pu, uint8_t  y_pu,
              uint8_t nb_pb_w, uint8_t nb_pb_h,
              uint8_t inter_dir)
{
    /*FIXME Use specific DBF update function if DBF is disabled */
    #if 0
    struct DBFInfo *const dbf_info = &lc_ctx->dbf_info;
    #endif
    if (inter_dir & 0x2) {
        #if 0
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #endif
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        fill_mvp_map(mv_ctx1, mv, x_pu, y_pu, nb_pb_w, nb_pb_h);

        #if 0
        fill_dbf_mv_map_b(dbf_info, mv_ctx1, mv_ctx0, mv1, x_pu, y_pu, nb_pb_w, nb_pb_h);
        #endif

    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        #if 0
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        #endif

        fill_mvp_map(mv_ctx0, mv, x_pu, y_pu, nb_pb_w, nb_pb_h);

        #if 0
        fill_dbf_mv_map_b(dbf_info, mv_ctx0, mv_ctx1, mv0, x_pu, y_pu, nb_pb_w, nb_pb_h);
        #endif

    }

    hmvp_update_lut(&inter_ctx->hmvp_lut, mv);
}

/* Derive motion vectors and update motion maps */
VVCMergeInfo
derive_mvp_b(struct InterDRVCtx *const inter_ctx,
             const OVPartInfo *const part_ctx,
             unsigned int x0, unsigned int y0,
             unsigned int log2_pb_w, unsigned int log2_pb_h,
             OVMV mvd0, OVMV mvd1,
             uint8_t mvp_idx0, uint8_t mvp_idx1,
             uint8_t inter_dir)
{
    /* FIXME replace part_ctx with something in inter CTX */
    uint8_t y_pu = y0 >> part_ctx->log2_min_cb_s;
    uint8_t x_pu = x0 >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_w = (1 << log2_pb_w) >> part_ctx->log2_min_cb_s;
    uint8_t nb_pb_h = (1 << log2_pb_h) >> part_ctx->log2_min_cb_s;
    OVMV mv0 = {0}, mv1 = {0};
    VVCMergeInfo mv_info;

    /* FIXME can we combine mvp derivation for bi pred */
    if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

        mv0 = derive_mvp_candidates(inter_ctx, mv_ctx0,
                                    x_pu, y_pu, nb_pb_w, nb_pb_h,
                                    mvp_idx0, inter_dir & 0x1);

        mvd0 = scale_mvd(mvd0);

        mv0.x += mvd0.x;
        mv0.y += mvd0.y;
    }

    if (inter_dir & 0x2) {
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

        mv1 = derive_mvp_candidates(inter_ctx, mv_ctx1,
                                    x_pu, y_pu,
                                    nb_pb_w, nb_pb_h,
                                    mvp_idx1, inter_dir & 0x2);

        mvd1 = scale_mvd(mvd1);

        mv1.x += mvd1.x;
        mv1.y += mvd1.y;
    }

    mv_info.inter_dir = inter_dir;
    mv_info.mv0 = mv0;
    mv_info.mv1 = mv1;

    /* Update for next pass */
    update_mv_ctx_b(inter_ctx, mv0, mv1, x_pu, y_pu, nb_pb_w,
                    nb_pb_h, inter_dir);
    return mv_info;
}

OVMV
derive_mvp_mvd(struct InterDRVCtx *const inter_ctx,
               const struct OVMVCtx *const mv_ctx,
               OVMV mvd,
               uint8_t pb_x, uint8_t pb_y,
               uint8_t nb_pb_w, uint8_t nb_pb_h,
               uint8_t mvp_idx, uint8_t inter_dir)
{
    OVMV mv;
    mv = derive_mvp_candidates(inter_ctx, mv_ctx,
                               pb_x, pb_y, nb_pb_w, nb_pb_h,
                               mvp_idx, 1);
    mvd = scale_mvd(mvd);

    mv.x += mvd.x;
    mv.y += mvd.y;

#if 1
    update_mv_ctx(inter_ctx, mv, pb_x, pb_y, nb_pb_w,
                  nb_pb_h, inter_dir);
#endif
   return mv;
}
