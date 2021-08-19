#include <string.h>

#include "ovutils.h"
#include "ctudec.h"
#include "drv_utils.h"
#include "drv.h"
#include "dec_structures.h"
#include "dbf_utils.h"
#include "ovdpb.h"
#include "rcn.h"

#define NB_TAP 6
#define NB_TAP_PLUS3 (NB_TAP + 3)
#define NB_TAP_PLUS5 (NB_TAP + 5)
#define NB_TAP_PLUS9 (NB_TAP + 9)

#define LOG2_MIN_CU_S 2

#define RND_AFF (4 << 11)

#define SB_SIZE 4
#define HALF_SB_SIZE (4 >> 1)

#define AFFINE_SHIFT 7

#define POS_MASK(x, w) ((uint64_t) 1 << ((((x + 1)) + ((w)))))

#define OFFSET_BUFF(x,y) (35 + x + (y) * 34)

#define PB_POS_IN_BUF(x,y) (35 + (x) + ((y) * 34))
#define TMVP_BUFF_STRIDE 17
#define TMVP_POS_IN_BUF(x,y) ((x >> 1) + (((y) >> 1) * TMVP_BUFF_STRIDE))

#define MAX_NB_AMVP_CAND 2

#define MV_BITS 18
#define CLIP_PERIOD  (1 << MV_BITS)
#define HALF_CLIP_PERIOD  (1 << (MV_BITS - 1))

#define MV_MAX   ((1 << (MV_BITS - 1)) - 1)
#define MV_MIN  (-(1 << (MV_BITS - 1)))

#define TMVP_POS_MASK(y) ((uint64_t) 1 << ((y) + 1))

#define MV_MANTISSA_BITCOUNT 6
#define MV_MANTISSA_UPPER_LIMIT ((1 << (MV_MANTISSA_BITCOUNT - 1)) - 1)
#define MV_MANTISSA_LIMIT (1 << (MV_MANTISSA_BITCOUNT - 1))


enum AffineType
{
    AFFINE_2CP = 0,
    AFFINE_3CP = 1
};

enum ControlPointIdx
{
   CP_LT = 0,
   CP_RT = 1,
   CP_LB = 2,
   CP_RB = 3,
};

enum ControlPointMask
{
   CP_LT_MASK = (1 << CP_LT),
   CP_RT_MASK = (1 << CP_RT),
   CP_LB_MASK = (1 << CP_LB),
   CP_RB_MASK = (1 << CP_RB),
};

enum ControlPointCandMask
{
   CP3_MASK0 = (1 << CP_LT) | (1 << CP_RT) | (1 << CP_LB),
   CP3_MASK1 = (1 << CP_LT) | (1 << CP_RT) | (1 << CP_RB),
   CP3_MASK2 = (1 << CP_LT) | (1 << CP_LB) | (1 << CP_RB),
   CP3_MASK3 = (1 << CP_RT) | (1 << CP_LB) | (1 << CP_RB),
   CP2_MASK0 = (1 << CP_LT) | (1 << CP_RT),
   CP2_MASK1 = (1 << CP_LT) | (1 << CP_LB)
};

/* FIXME only keep points */
struct AffineDeltaMV
{
    OVMV h;
    OVMV v;
};

struct Affine3CP
{
    OVMV cp[3];
};

struct Affine2CP
{
    OVMV cp[2];
};

struct AffineCPInfo
{
    uint8_t type;
    union {
        struct Affine3CP type3;
        struct Affine2CP type2;
    }cps;
};

enum CandName
{
   A0 = 0, 
   A1 = 1, 
   A2 = 2, 
   A3 = 3, 

   B0 = 4, 
   B1 = 5, 
   B2 = 6, 
   B3 = 7, 
};

enum CandListMask
{
   A0_MSK = 1 << A0, 
   A1_MSK = 1 << A1, 
   A2_MSK = 1 << A2, 
   A3_MSK = 1 << A3, 
                   
   B0_MSK = 1 << B0, 
   B1_MSK = 1 << B1, 
   B2_MSK = 1 << B2, 
   B3_MSK = 1 << B3
};

struct TMVPPos {
    uint8_t c0_x;
    uint8_t c0_y;
    uint8_t c1_x;
    uint8_t c1_y;
};

/* Generic Motion vector clipping and rounding helper functions
 */

static inline OVMV
clip_mv(OVMV src)
{
    OVMV dst;
    const int mv_min = -(1 << 17);
    const int mv_max = (1 << 17) - 1;

    dst.x = ov_clip(src.x, mv_min, mv_max);
    dst.y = ov_clip(src.y, mv_min, mv_max);

    return dst;
}

/*FIXME is different from clip_mv ? */
static inline OVMV
mv_clip_periodic(OVMV src)
{
    OVMV dst;

    dst.x = (src.x + CLIP_PERIOD) & (CLIP_PERIOD - 1);
    dst.y = (src.y + CLIP_PERIOD) & (CLIP_PERIOD - 1);

    dst.x = (dst.x >= HALF_CLIP_PERIOD) ? (dst.x - CLIP_PERIOD) : dst.x;
    dst.y = (dst.y >= HALF_CLIP_PERIOD) ? (dst.y - CLIP_PERIOD) : dst.y;

    return dst;
}

//TODOamvr
static inline OVMV
round_affine_mv(const OVMV mv, uint8_t prec_amvr)
{
    OVMV tmp;

    tmp  = drv_round_to_precision_mv(mv, MV_PRECISION_INTERNAL, prec_amvr);

    // tmp.x = (mv.x + 1 + (mv.x < 0)) >> 2;
    // tmp.y = (mv.y + 1 + (mv.y < 0)) >> 2;

    // tmp.x = tmp.x << 2;
    // tmp.y = tmp.y << 2;

    return tmp;
}

static inline OVMV
round_affine_mv2(OVMV mv)
{
  const int rnd = 1 << (AFFINE_SHIFT - 1);
  OVMV tmp;

  tmp.x = mv.x + rnd;
  tmp.y = mv.y + rnd;

  tmp.x -= mv.x >= 0;
  tmp.y -= mv.y >= 0;

  tmp.x >>= AFFINE_SHIFT;
  tmp.y >>= AFFINE_SHIFT;

  return tmp;
}

/* Generic Ref Pic List  helper functions
 */

enum RPLIndex 
{
    RPL_0 = 0,
    RPL_1 = 1,
};

static inline enum RPLIndex
opposit_rpl_idx(enum RPLIndex rpl_idx)
{
    return !rpl_idx;
}

/* Generic candidates derivation helper functions
 */

static uint8_t
is_above_cand(enum CandName cand_name)
{
    /* Note we use equal since A3 denotes same cand as B2*/
    return (cand_name >= A3);
}

static inline uint8_t
check_cand_available(uint64_t abv_row, uint64_t lft_col, uint8_t pb_x, uint8_t pb_y,
                     uint8_t nb_pb_w, uint8_t nb_pb_h)
{
    /* Note A3 and B2 correspond to the same candidate */
    uint64_t a0_pos_msk = POS_MASK(pb_y, nb_pb_h);
    uint64_t a1_pos_msk = POS_MASK(pb_y, nb_pb_h - 1);
    uint64_t a2_pos_msk = POS_MASK(pb_y, 0);
    uint64_t a3_pos_msk = POS_MASK(pb_y - 1, 0);

    uint64_t b0_pos_msk = POS_MASK(pb_x, nb_pb_w);
    uint64_t b1_pos_msk = POS_MASK(pb_x, nb_pb_w - 1);
    uint64_t b2_pos_msk = POS_MASK(pb_x - 1, 0);

    uint64_t b3_pos_msk = POS_MASK(pb_x, 0);

    uint8_t cand_a0 = !!(lft_col & a0_pos_msk);
    uint8_t cand_a1 = !!(lft_col & a1_pos_msk);
    uint8_t cand_a2 = !!(lft_col & a2_pos_msk);
    uint8_t cand_a3 = !!(lft_col & a3_pos_msk);

    uint8_t cand_b0 = !!(abv_row & b0_pos_msk);
    uint8_t cand_b1 = !!(abv_row & b1_pos_msk);
    uint8_t cand_b2 = !!(abv_row & b2_pos_msk);
    uint8_t cand_b3 = !!(abv_row & b3_pos_msk);

    uint8_t cand_list = cand_a0;

    cand_list |= cand_a1 << 1;
    cand_list |= cand_a2 << 2;
    cand_list |= cand_a3 << 3;

    cand_list |= cand_b0 << 4;
    cand_list |= cand_b1 << 5;
    cand_list |= cand_b2 << 6;
    cand_list |= cand_b3 << 7;

    return cand_list;
}

static inline enum CandName
cand_mask_to_idx(enum CandListMask cand_msk)
{
    return (enum CandName)ov_ctz(cand_msk);
}

static inline int16_t
derive_cand_position(struct PBInfo pb_info, enum CandName cand_name)
{
    int16_t pos = 0;

    int16_t pb_x = pb_info.x_pb;
    int16_t pb_y = pb_info.y_pb;
    int16_t nb_pb_w = pb_info.nb_pb_w;
    int16_t nb_pb_h = pb_info.nb_pb_h;

    switch (cand_name) {
        case A0:
            pos = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h);
        break;
        case A1:
            pos = OFFSET_BUFF(pb_x - 1, pb_y + nb_pb_h - 1);
        break;
        case A2:
            pos = OFFSET_BUFF(pb_x - 1, pb_y);
        break;
        case A3:
            pos = OFFSET_BUFF(pb_x - 1, pb_y - 1);
        break;
        case B0:
            pos = OFFSET_BUFF(pb_x + nb_pb_w, pb_y - 1);
        break;
        case B1:
            pos = OFFSET_BUFF(pb_x + nb_pb_w - 1, pb_y - 1);
        break;
        case B2:
            pos = OFFSET_BUFF(pb_x - 1, pb_y - 1);
        break;
        case B3:
            pos = OFFSET_BUFF(pb_x, pb_y - 1);
        break;
    }

    return pos;
}

/* Compress and uncompress Motion vectors used by TMVP
 * FIXME there might be a more straight forward way of
 * doing this
 */
static inline int32_t
tmvp_compress_uncompress_mv(int32_t val)
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
        memset(tmvp_ctx->dir_map_v0, 0, sizeof(uint64_t) * 34);
        memset(tmvp_ctx->dir_map_v1, 0, sizeof(uint64_t) * 34);
    }
        memset(tmvp_ctx->mvs0, 0, sizeof(tmvp_ctx->mvs0));
        memset(tmvp_ctx->mvs1, 0, sizeof(tmvp_ctx->mvs1));

    if (plane0 && plane0->dirs) {
        uint64_t *src_dirs = plane0->dirs + ctb_addr_rs * nb_pb_ctb_w;

        int32_t nb_tmvp_unit = nb_pb_ctb_w >> 1;
        int32_t pln_stride = nb_tmvp_unit * nb_ctb_w;
        int32_t ctb_offset = ctb_x * nb_tmvp_unit + (ctb_y * nb_tmvp_unit * pln_stride);
        OVMV *src_mv = plane0->mvs + ctb_offset;
        OVMV *mvs = tmvp_ctx->mvs0;
        int i;

        memcpy(&tmvp_ctx->dir_map_v0[1], src_dirs, sizeof(uint64_t) * (nb_pb_ctb_w + !is_border_pic));
        for (i = 0; i < nb_pb_ctb_w; i += 2) {
            memcpy(mvs, src_mv, sizeof(*mvs) * (nb_tmvp_unit + !is_border_pic));
            mvs += TMVP_BUFF_STRIDE;
            src_mv += pln_stride;
        }
    }

    if (plane1 && plane1->dirs) {
        OVMV *mvs = tmvp_ctx->mvs1;
        uint64_t *src_dirs = plane1->dirs + ctb_addr_rs * nb_pb_ctb_w;
        int i;

        int32_t nb_tmvp_unit = nb_pb_ctb_w >> 1;
        int32_t pln_stride = nb_tmvp_unit * nb_ctb_w;
        int32_t ctb_offset = ctb_x * nb_tmvp_unit + (ctb_y * nb_tmvp_unit * pln_stride);
        OVMV *src_mv = plane1->mvs + ctb_offset;

        /*FIXME memory could be spared with smaller map size when possible */
        memcpy(&tmvp_ctx->dir_map_v1[1], src_dirs, sizeof(uint64_t) * (nb_pb_ctb_w + !is_border_pic));
        for (i = 0; i < nb_pb_ctb_w; i += 2) {
            memcpy(mvs, src_mv, sizeof(*mvs) * (nb_tmvp_unit + !is_border_pic));
            mvs += TMVP_BUFF_STRIDE;
            src_mv += pln_stride;
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

struct TMVPPos
compute_tmpv_coord(struct PBInfo pb, uint8_t log2_min_cb_s)
{
    struct TMVPPos pos;
    uint8_t pos_8x8 = log2_min_cb_s == 2;

    pos.c0_x = (pb.x_pb + pb.nb_pb_w) & (~pos_8x8);
    pos.c0_y = (pb.y_pb + pb.nb_pb_h) & (~pos_8x8);

    /* center */
    pos.c1_x = (pb.x_pb + (pb.nb_pb_w >> 1)) & (~pos_8x8);
    pos.c1_y = (pb.y_pb + (pb.nb_pb_h >> 1)) & (~pos_8x8);

    return pos;
}

static inline uint8_t
check_tmvp_cand(const uint64_t *v_map0, const uint64_t *v_map1,
                struct TMVPPos pos)
{
    /* Derive availability based on CTB inter fields */
    uint64_t c0_col  = v_map0[pos.c0_x + 1];
    uint64_t c0_col1 = v_map1[pos.c0_x + 1];

    uint64_t c1_col  = v_map0[pos.c1_x + 1];
    uint64_t c1_col1 = v_map1[pos.c1_x + 1];

    uint8_t cand_c0  = !!(c0_col  & TMVP_POS_MASK(pos.c0_y));
    uint8_t cand_c01 = !!(c0_col1 & TMVP_POS_MASK(pos.c0_y));

    uint8_t cand_c1  = !!(c1_col  & TMVP_POS_MASK(pos.c1_y));
    uint8_t cand_c11 = !!(c1_col1 & TMVP_POS_MASK(pos.c1_y));

    uint8_t cand_msk = cand_c0; 
    cand_msk |= cand_c01 << 1;

    cand_msk |= cand_c1  << 2;
    cand_msk |= cand_c11 << 3;

    return cand_msk; 
}

static OVMV tmvp_rescale(OVMV mv, int16_t scale)
{
    OVMV dst;

    mv.x = tmvp_compress_uncompress_mv(mv.x);
    mv.y = tmvp_compress_uncompress_mv(mv.y);

    dst = tmvp_scale_mv(scale, mv);

    return dst;
}

static uint8_t
tmvp_from_l0(const struct InterDRVCtx *const inter_ctx, const struct VVCTMVP *const tmvp, struct TMVPPos pos,
             uint8_t rpl_idx, uint8_t ref_idx,
             uint8_t cand_msk, OVMV *const dst)
{
    int32_t dist_ref = rpl_idx == RPL_0 ? inter_ctx->dist_ref_0[ref_idx]
                                        : inter_ctx->dist_ref_1[ref_idx];

    uint8_t cand_c0  = cand_msk & 0x1;
    uint8_t cand_c01 = cand_msk & 0x2;
    uint8_t cand_c1  = cand_msk & 0x4;
    uint8_t cand_c11 = cand_msk & 0x8;

    OVMV mv;

    int32_t dist_col;
    int16_t scale;

    if (cand_c0 | cand_c01) {
        int16_t c0_pos = TMVP_POS_IN_BUF(pos.c0_x, pos.c0_y);

        const OVMV *mvs    = cand_c0 ? tmvp->mvs0
                                     : tmvp->mvs1;

        const int16_t *dist_cols = cand_c0 ? tmvp->dist_col_0
                                           : tmvp->dist_col_1;

        mv       = mvs[c0_pos];
        dist_col = dist_cols[mv.ref_idx];

        goto found;

    } else if (cand_c1 | cand_c11) {
        int16_t c1_pos = TMVP_POS_IN_BUF(pos.c1_x, pos.c1_y);

        const OVMV *mvs    = cand_c1 ? tmvp->mvs0
                                     : tmvp->mvs1;

        const int16_t *dist_cols = cand_c1 ? tmvp->dist_col_0
                                           : tmvp->dist_col_1;


        mv       = mvs[c1_pos];
        dist_col = dist_cols[mv.ref_idx];

        goto found;
    }

    return 0;

found :
    scale = derive_tmvp_scale(dist_ref, dist_col);

    mv = tmvp_rescale(mv, scale);

    mv.ref_idx = ref_idx;
    mv.bcw_idx_plus1 = 0;
    mv.prec_amvr = 0;

    *dst = mv;

    return 1;
}

/* FIXME We could invert TMVP context buff to avoid duplicating this function */
static uint8_t
tmvp_from_l1(const struct InterDRVCtx *const inter_ctx, const struct VVCTMVP *const tmvp, struct TMVPPos pos,
             uint8_t rpl_idx, uint8_t ref_idx,
             uint8_t cand_msk, OVMV *const dst)
{
    int32_t dist_ref = rpl_idx == RPL_0 ? inter_ctx->dist_ref_0[ref_idx]
                                        : inter_ctx->dist_ref_1[ref_idx];
    uint8_t cand_c0  = cand_msk & 0x1;
    uint8_t cand_c01 = cand_msk & 0x2;
    uint8_t cand_c1  = cand_msk & 0x4;
    uint8_t cand_c11 = cand_msk & 0x8;

    OVMV mv;

    int32_t dist_col;
    int16_t scale;

    if (cand_c0 | cand_c01) {
        int16_t c0_pos = TMVP_POS_IN_BUF(pos.c0_x, pos.c0_y);

        const OVMV *mvs    = cand_c01 ? tmvp->mvs1
                                      : tmvp->mvs0;

        const int16_t *dist_cols = cand_c01 ? tmvp->dist_col_1
                                            : tmvp->dist_col_0;

        mv       = mvs[c0_pos];
        dist_col = dist_cols[mv.ref_idx];

        goto found;

    } else if (cand_c1 | cand_c11) {
        int16_t c1_pos = TMVP_POS_IN_BUF(pos.c1_x, pos.c1_y);

        const OVMV *mvs    = cand_c11 ? tmvp->mvs1
                                      : tmvp->mvs0;

        const int16_t *dist_cols = cand_c11 ? tmvp->dist_col_1
                                            : tmvp->dist_col_0;

        mv       = mvs[c1_pos];
        dist_col = dist_cols[mv.ref_idx];

        goto found;
    }

    return 0;

found :
    scale = derive_tmvp_scale(dist_ref, dist_col);

    mv = tmvp_rescale(mv, scale);

    mv.ref_idx = ref_idx;
    mv.bcw_idx_plus1 = 0;
    mv.prec_amvr = 0;

    *dst = mv;

    return 1;
}

static uint8_t
merge_tmvp_from_ldc(const struct InterDRVCtx *const inter_ctx, const struct VVCTMVP *const tmvp, struct TMVPPos pos,
                    uint8_t rpl_idx, uint8_t ref_idx,
                    uint8_t cand_msk, OVMV *const dst)
{
    int32_t dist_ref = rpl_idx == RPL_0 ? inter_ctx->dist_ref_0[ref_idx]
                                        : inter_ctx->dist_ref_1[ref_idx];

    int32_t dist_ref_opp = rpl_idx == RPL_0 ? inter_ctx->dist_ref_1[ref_idx]
                                            : inter_ctx->dist_ref_0[ref_idx];
    uint8_t cand_c0  = cand_msk & 0x1;
    uint8_t cand_c01 = cand_msk & 0x2;

    OVMV mv;

    int32_t dist_col;
    int16_t scale;

    if (cand_c0 | cand_c01) {
        int16_t c0_pos = TMVP_POS_IN_BUF(pos.c0_x, pos.c0_y);
        uint8_t dir = 0;

        if (cand_c0 && cand_c01 && !tmvp->col_ref_l0) {
            const OVMV *mvs    = tmvp->mvs0;
            const int16_t *dist_cols = tmvp->dist_col_0;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            mv       = tmvp->mvs1[c0_pos];
            dist_col = tmvp->dist_col_1[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[1] = mv;
        } else if (cand_c0 && cand_c01) {
            const OVMV *mvs    = tmvp->mvs1;
            const int16_t *dist_cols = tmvp->dist_col_1;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            mv       = tmvp->mvs0[c0_pos];
            dist_col = tmvp->dist_col_0[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        } else if (cand_c0) {
            const OVMV *mvs    = tmvp->mvs0;
            const int16_t *dist_cols = tmvp->dist_col_0;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            mv       = mvs[c0_pos];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[1] = mv;

        } else if (cand_c01) {
            const OVMV *mvs    = tmvp->mvs1;
            const int16_t *dist_cols = tmvp->dist_col_1;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            mv       = mvs[c0_pos];

            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[1] = mv;
        }

        return dir;
    }

    return 0;
}

static uint8_t
merge_tmvp_from_l0(const struct InterDRVCtx *const inter_ctx, const struct VVCTMVP *const tmvp, struct TMVPPos pos,
                   uint8_t rpl_idx, uint8_t ref_idx,
                   uint8_t cand_msk, OVMV *const dst)
{
    int32_t dist_ref = rpl_idx == RPL_0 ? inter_ctx->dist_ref_0[ref_idx]
                                        : inter_ctx->dist_ref_1[ref_idx];

    int32_t dist_ref_opp = rpl_idx == RPL_0 ? inter_ctx->dist_ref_1[ref_idx]
                                            : inter_ctx->dist_ref_0[ref_idx];
    uint8_t cand_c0  = cand_msk & 0x1;
    uint8_t cand_c01 = cand_msk & 0x2;

    OVMV mv;

    int32_t dist_col;
    int16_t scale;

    if (cand_c0 | cand_c01) {
        int16_t c0_pos = TMVP_POS_IN_BUF(pos.c0_x, pos.c0_y);
        uint8_t dir = 0;

        if (cand_c0) {
            const OVMV *mvs    = tmvp->mvs0;
            const int16_t *dist_cols = tmvp->dist_col_0;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            mv       = mvs[c0_pos];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[1] = mv;

        } else if (cand_c01) {
            const OVMV *mvs    = tmvp->mvs1;
            const int16_t *dist_cols = tmvp->dist_col_1;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            mv       = mvs[c0_pos];

            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        }

        return dir;
    }

    return 0;
}

/* FIXME We could invert TMVP context buff to avoid duplicating this function */
static uint8_t
merge_tmvp_from_l1(const struct InterDRVCtx *const inter_ctx, const struct VVCTMVP *const tmvp, struct TMVPPos pos,
                   uint8_t rpl_idx, uint8_t ref_idx,
                   uint8_t cand_msk, OVMV *const dst)
{
    int32_t dist_ref = rpl_idx == RPL_0 ? inter_ctx->dist_ref_0[ref_idx]
                                        : inter_ctx->dist_ref_1[ref_idx];

    int32_t dist_ref_opp = rpl_idx == RPL_0 ? inter_ctx->dist_ref_1[ref_idx]
                                            : inter_ctx->dist_ref_0[ref_idx];
    uint8_t cand_c0  = cand_msk & 0x1;
    uint8_t cand_c01 = cand_msk & 0x2;

    OVMV mv;

    int32_t dist_col;
    int16_t scale;

    if (cand_c0 | cand_c01) {
        int16_t c0_pos = TMVP_POS_IN_BUF(pos.c0_x, pos.c0_y);
        uint8_t dir = 0;

        if (cand_c01) {
            const OVMV *mvs    = tmvp->mvs1;
            const int16_t *dist_cols = tmvp->dist_col_1;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            dir |= 0x2;

            mv       = mvs[c0_pos];

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        } else if (cand_c0) {
            const OVMV *mvs    = tmvp->mvs0;
            const int16_t *dist_cols = tmvp->dist_col_0;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[0] = mv;

            dir |= 0x1;

            mv       = mvs[c0_pos];

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        }

        return dir;
    }

    return 0;
}

/* Generic Affine control points and delta MVs derivation functions
 */

static struct AffineDeltaMV
derive_affine_delta_mvs(const struct AffineControlInfo *const cinfo,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t affine_type)
{
    struct AffineDeltaMV dmv;
    OVMV delta_mv_h, delta_mv_v;
    delta_mv_h = cinfo->lt;
    delta_mv_v = cinfo->lt;

    const uint8_t scale_h = AFFINE_SHIFT - log2_pb_w;

    delta_mv_h.x = (cinfo->rt.x - cinfo->lt.x) << scale_h;
    delta_mv_h.y = (cinfo->rt.y - cinfo->lt.y) << scale_h;

    if (affine_type == AFFINE_3CP) {
        const uint8_t scale_v = AFFINE_SHIFT - log2_pb_h;

        delta_mv_v.x = (cinfo->lb.x - cinfo->lt.x) << scale_v;
        delta_mv_v.y = (cinfo->lb.y - cinfo->lt.y) << scale_v;

    } else {
        delta_mv_v.x = -delta_mv_h.y;
        delta_mv_v.y =  delta_mv_h.x;
    }

    dmv.h = delta_mv_h;
    dmv.v = delta_mv_v;

    return dmv;
}

/* FIXME can we derive this without neighbour x/y info
         using delta_mv information directly */
static struct AffineControlInfo
derive_cp_from_cand(const struct AffineControlInfo *const ngh_cp,
                    struct PBInfo pb_info, struct PBInfo ngh_pb,
                    uint8_t affine_type, uint8_t ngh_affine_type,
                    enum CandName cand_name)
{
    int ngh_x0 = ngh_pb.x_pb << 2;
    int ngh_y0 = ngh_pb.y_pb << 2;
    uint8_t log2_ngh_w = ngh_pb.log2_w;
    uint8_t log2_ngh_h = ngh_pb.log2_h;

    uint8_t log2_pb_w = pb_info.log2_w;
    uint8_t log2_pb_h = pb_info.log2_h;
    int x0 = pb_info.x_pb << 2;
    int y0 = pb_info.y_pb << 2;

    struct AffineControlInfo dst_cp = {0};
    uint8_t ref_idx = ngh_cp->lt.ref_idx;
    uint8_t bcw_idx_plus1 = ngh_cp->lt.bcw_idx_plus1;
    uint8_t prec_amvr = ngh_cp->lt.prec_amvr;

    /* FIXME avoid checking for above candidate here */
    uint8_t is_abv_ctu = y0 == 0 && is_above_cand(cand_name);
    uint8_t is_lft_ctu = x0 == 0 && (!is_above_cand(cand_name) || cand_name == B2);
    uint8_t is_abv_rgt_ctu = cand_name == B0 && is_abv_ctu && x0 + (1 << log2_pb_w) == 1 << 7;

    /* FIXME use correct log2_ctu_s */
    int delta_pos_x = ((is_lft_ctu << 7) + x0 - ((is_abv_rgt_ctu << 7) + ngh_x0));
    int delta_pos_y = is_abv_ctu ? 0 : (y0 - ngh_y0);

    /* FIXME determine clip from merge or mvp cand derivation */
    struct AffineDeltaMV delta_mv = derive_affine_delta_mvs(ngh_cp, log2_ngh_w, log2_ngh_h,
                                                            ngh_affine_type);

    OVMV tmp = {0}, lt_mv = {0};

    if (is_abv_ctu || ngh_affine_type == AFFINE_2CP) {
        delta_mv.v.x = -delta_mv.h.y;
        delta_mv.v.y =  delta_mv.h.x;
    }

    /* Derive LT based on candidate position */
    lt_mv.x = ngh_cp->lt.x << AFFINE_SHIFT;
    lt_mv.y = ngh_cp->lt.y << AFFINE_SHIFT;

    lt_mv.x += delta_mv.h.x * delta_pos_x;
    lt_mv.y += delta_mv.h.y * delta_pos_x;

    /* Note that delta_pos_y is 0 if not from same CTU line 
     * delta_mv.v.x will be -delta_mv_h.y and
     * delta_mv.v.y will be  delta_mv_h.x
     */
    lt_mv.x += delta_mv.v.x * delta_pos_y;
    lt_mv.y += delta_mv.v.y * delta_pos_y;

    dst_cp.lt = round_affine_mv2(lt_mv);
    dst_cp.lt = clip_mv(dst_cp.lt);
    dst_cp.lt.ref_idx = ref_idx;
    dst_cp.lt.bcw_idx_plus1 = bcw_idx_plus1;
    dst_cp.lt.prec_amvr = prec_amvr;

    tmp.x = lt_mv.x + (delta_mv.h.x << log2_pb_w);
    tmp.y = lt_mv.y + (delta_mv.h.y << log2_pb_w);

    dst_cp.rt = round_affine_mv2(tmp);
    dst_cp.rt = clip_mv(dst_cp.rt);
    dst_cp.rt.ref_idx = ref_idx;
    dst_cp.rt.bcw_idx_plus1 = bcw_idx_plus1;
    dst_cp.rt.prec_amvr = prec_amvr;

    if (affine_type == AFFINE_3CP) {
        tmp.x = lt_mv.x + (delta_mv.v.x << log2_pb_h);
        tmp.y = lt_mv.y + (delta_mv.v.y << log2_pb_h);

        dst_cp.lb = round_affine_mv2(tmp);
        dst_cp.lb = clip_mv(dst_cp.lb);
        dst_cp.lb.ref_idx = ref_idx;
        dst_cp.lb.bcw_idx_plus1 = bcw_idx_plus1;
        dst_cp.lb.prec_amvr = prec_amvr;
    }

    return dst_cp;
}

/* Affine MVP related functions
 */

static uint8_t
derive_affine_mvp_cand(const struct AffineDRVInfo *const affine_ctx,
                       struct AffineControlInfo *const dst_cp_info,
                       struct PBInfo pb, enum CandName cand_name,
                       uint8_t inter_dir, uint8_t ref_idx, uint8_t ref_opp_idx,
                       uint8_t rpl_msk0, uint8_t rpl_msk1, uint8_t aff_msk,
                       uint8_t prec_amvr, uint8_t affine_type)
{

    /* Note this is OK since inter_dir cannot be 3 from MVP */
    uint8_t avail_affine  = (aff_msk     & (1 << cand_name));
    struct AffineControlInfo cp_info;

    if (avail_affine) {
        enum RPLIndex rpl_idx = inter_dir - 1;
        enum RPLIndex rpl_opp_idx = opposit_rpl_idx(rpl_idx);

        uint8_t rpl_msk = rpl_idx ? rpl_msk1 : rpl_msk0;
        uint8_t rpl_opp_msk = rpl_idx ? rpl_msk0 : rpl_msk1;

        uint8_t avail_rpl     = (rpl_msk     & (1 << cand_name));
        uint8_t avail_rpl_opp = (rpl_opp_msk & (1 << cand_name));

        const int16_t cand_pos = derive_cand_position(pb, cand_name);
        const struct AffineInfo *const affine_info = &affine_ctx->affine_info[cand_pos];
        const struct PBInfo ngh_pb = affine_info->pb;

        if (avail_rpl) {
            struct AffineControlInfo ngh_cp_info = affine_info->cps[rpl_idx];
            if (ngh_cp_info.lt.ref_idx == ref_idx) {

                cp_info = derive_cp_from_cand(&ngh_cp_info, pb, ngh_pb, affine_type,
                                              affine_info->type, cand_name);

                goto found;
            }
        }

        if (avail_rpl_opp) {
            struct AffineControlInfo ngh_cp_info = affine_info->cps[rpl_opp_idx];

            if (ngh_cp_info.lt.ref_idx == ref_opp_idx) {

                cp_info = derive_cp_from_cand(&ngh_cp_info, pb, ngh_pb, affine_type,
                                              affine_info->type, cand_name);

                /* override ref_idx since it is taken from opposit RPL */
                cp_info.lt.ref_idx = ref_idx;
                cp_info.rt.ref_idx = ref_idx;

                if (affine_type == AFFINE_3CP) {
                    cp_info.lb.ref_idx = ref_idx;
                }

                goto found;
            }
        }
    }

    return 0;

found:
    cp_info.lt = round_affine_mv(cp_info.lt, prec_amvr);
    cp_info.rt = round_affine_mv(cp_info.rt, prec_amvr);

    if (affine_type == AFFINE_3CP) {
        cp_info.lb = round_affine_mv(cp_info.lb, prec_amvr);
    }

    *dst_cp_info = cp_info;

    return 1;
}

static uint8_t
derive_mvp_cand(const struct InterDRVCtx *const inter_ctx,
                struct PBInfo pb, enum CandName cand_name,
                uint8_t inter_dir, uint8_t ref_idx, uint8_t ref_opp_idx,
                uint8_t rpl0_list, uint8_t rpl1_list,
                OVMV *dst_mv)
{
    const int16_t cand_pos = derive_cand_position(pb, cand_name);

    /* Note this is OK since inter_dir cannot be 3 from MVP */
    enum RPLIndex rpl_idx = inter_dir - 1;
    enum RPLIndex rpl_opp_idx = opposit_rpl_idx(rpl_idx);

    const OVMV *const mv_ctx     = rpl_idx ? inter_ctx->mv_ctx1.mvs : inter_ctx->mv_ctx0.mvs;
    const OVMV *const mv_ctx_opp = rpl_opp_idx ? inter_ctx->mv_ctx1.mvs : inter_ctx->mv_ctx0.mvs;

    uint8_t rpl_list     = rpl_idx ? rpl1_list : rpl0_list;
    uint8_t rpl_list_opp = rpl_idx ? rpl0_list : rpl1_list;
    OVMV mv_cand;

    if (rpl_list & (1 << cand_name)) {
        mv_cand = mv_ctx[cand_pos];

        if (mv_cand.ref_idx == ref_idx) {

            goto found;
        }
    }

    if (rpl_list_opp & (1 << cand_name)) {
        mv_cand = mv_ctx_opp[cand_pos];

        if (mv_cand.ref_idx == ref_opp_idx) {

            mv_cand.ref_idx = ref_idx;
            goto found;
        }
    }

    return 0;

found:

    *dst_mv = mv_cand;

    return 1;
}

static struct AffineControlInfo
drv_affine_mvp(struct InterDRVCtx *const inter_ctx,
               struct AffineDRVInfo *affine_ctx,
               uint8_t x_pb, uint8_t y_pb,
               uint8_t nb_pb_w, uint8_t nb_pb_h,
               uint8_t log2_cu_w, uint8_t log2_cu_h,
               uint8_t ref_idx, uint8_t ref_opp_idx, uint8_t mvp_idx,
               uint8_t inter_dir, uint8_t affine_type)
{
    /*FIXME do not search for third control point when type flag is zero */
    /*FIXME early termination ?*/
    uint8_t nb_cand = 0;
    uint8_t cand_aff_lft;
    uint8_t cand_aff_abv;

    uint8_t cand_lt;
    uint8_t cand_rt;
    uint8_t cand_lb;

    int cand_mask = 0;
    uint8_t prec_amvr = inter_ctx->prec_amvr;
    OVMV lt_mv_cand = {0};
    OVMV rt_mv_cand = {0};
    OVMV lb_mv_cand = {0};

    OVMV mv_aff[3];

    struct PBInfo pb_info = { 
        .x_pb = x_pb, 
        .y_pb = y_pb, 
        .nb_pb_w = nb_pb_w, 
        .nb_pb_h = nb_pb_h,
        .log2_w = log2_cu_w,
        .log2_h = log2_cu_h
    };

    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t aff_lft_col = affine_ctx->map.vfield[x_pb];
    uint64_t aff_abv_row = affine_ctx->map.hfield[y_pb];
    uint64_t rpl0_lft_col = mv_ctx0->map.vfield[x_pb];
    uint64_t rpl0_abv_row = mv_ctx0->map.hfield[y_pb];
    uint64_t rpl1_lft_col = mv_ctx1->map.vfield[x_pb];
    uint64_t rpl1_abv_row = mv_ctx1->map.hfield[y_pb];

    const uint8_t aff_cand_list = check_cand_available(aff_abv_row, aff_lft_col, x_pb, y_pb,
                                                       nb_pb_w, nb_pb_h);

    const uint8_t rpl0_cand = check_cand_available(rpl0_abv_row, rpl0_lft_col, x_pb, y_pb,
                                                   nb_pb_w, nb_pb_h);

    const uint8_t rpl1_cand = check_cand_available(rpl1_abv_row, rpl1_lft_col, x_pb, y_pb,
                                                   nb_pb_w, nb_pb_h);

    /* Control points from affine neighbours */
    /* Affine left cand */
    struct AffineControlInfo cp_info[2];


    /* FIXME check affine */
    cand_aff_lft = derive_affine_mvp_cand(affine_ctx, cp_info, pb_info, A0,
                                          inter_dir, ref_idx, ref_opp_idx,
                                          rpl0_cand, rpl1_cand, aff_cand_list,
                                          prec_amvr, affine_type);
    if (!cand_aff_lft) {
        cand_aff_lft = derive_affine_mvp_cand(affine_ctx, cp_info, pb_info, A1,
                                              inter_dir, ref_idx, ref_opp_idx,
                                              rpl0_cand, rpl1_cand, aff_cand_list,
                                              prec_amvr, affine_type);
    }

    /* Affine above cand */
    cand_aff_abv = derive_affine_mvp_cand(affine_ctx, &cp_info[cand_aff_lft], pb_info, B0,
                                          inter_dir, ref_idx, ref_opp_idx,
                                          rpl0_cand, rpl1_cand, aff_cand_list,
                                          prec_amvr, affine_type);
    if (!cand_aff_abv) {
        cand_aff_abv = derive_affine_mvp_cand(affine_ctx, &cp_info[cand_aff_lft], pb_info, B1,
                                              inter_dir, ref_idx, ref_opp_idx,
                                              rpl0_cand, rpl1_cand, aff_cand_list,
                                              prec_amvr, affine_type);
        if (!cand_aff_abv) {
            cand_aff_abv = derive_affine_mvp_cand(affine_ctx, &cp_info[cand_aff_lft], pb_info, B2,
                                                  inter_dir, ref_idx, ref_opp_idx,
                                                  rpl0_cand, rpl1_cand, aff_cand_list,
                                                  prec_amvr, affine_type);
        }
    }

    nb_cand = cand_aff_lft + cand_aff_abv;

    if (nb_cand >= MAX_NB_AMVP_CAND) {
        /* FIXME is rounding needed since done in cand derivation */
        for (int i = 0; i < nb_cand; i++) {
            cp_info[i].lt = round_affine_mv(cp_info[i].lt, prec_amvr);
            cp_info[i].rt = round_affine_mv(cp_info[i].rt, prec_amvr);
            cp_info[i].lb = round_affine_mv(cp_info[i].lb, prec_amvr);
        }
        return cp_info[mvp_idx];
    }

    /* Control points from MVs */
    /* Cand LT */
    cand_lt = derive_mvp_cand(inter_ctx, pb_info, B2, inter_dir, ref_idx, ref_opp_idx,
                              rpl0_cand, rpl1_cand, &lt_mv_cand);
    if (!cand_lt) {
        cand_lt = derive_mvp_cand(inter_ctx, pb_info, B3, inter_dir, ref_idx, ref_opp_idx,
                                  rpl0_cand, rpl1_cand, &lt_mv_cand);
        if (!cand_lt) {
            cand_lt = derive_mvp_cand(inter_ctx, pb_info, A2, inter_dir, ref_idx, ref_opp_idx,
                                      rpl0_cand, rpl1_cand, &lt_mv_cand);
        }
    }

    cand_mask |= cand_lt;

    /* Cand RT */
    cand_rt = derive_mvp_cand(inter_ctx, pb_info, B1, inter_dir, ref_idx, ref_opp_idx,
                              rpl0_cand, rpl1_cand, &rt_mv_cand);
    if (!cand_rt) {
        cand_rt = derive_mvp_cand(inter_ctx, pb_info, B0, inter_dir, ref_idx, ref_opp_idx,
                                  rpl0_cand, rpl1_cand, &rt_mv_cand);
    }

    cand_mask |= cand_rt << 1;

    /*Cand LB */
    cand_lb = derive_mvp_cand(inter_ctx, pb_info, A1, inter_dir, ref_idx, ref_opp_idx,
                              rpl0_cand, rpl1_cand, &lb_mv_cand);
    if (!cand_lb) {
        cand_lb = derive_mvp_cand(inter_ctx, pb_info, A0, inter_dir, ref_idx, ref_opp_idx,
                                  rpl0_cand, rpl1_cand, &lb_mv_cand);
    }

    cand_mask |= cand_lb << 2;

    mv_aff[CP_LT] = round_affine_mv(lt_mv_cand, prec_amvr);
    mv_aff[CP_RT] = round_affine_mv(rt_mv_cand, prec_amvr);
    mv_aff[CP_LB] = round_affine_mv(lb_mv_cand, prec_amvr);

    if (cand_mask == 0x7 || (cand_mask == 0x3 && affine_type == AFFINE_2CP)) {
        cp_info[nb_cand].lt = mv_aff[CP_LT];
        cp_info[nb_cand].rt = mv_aff[CP_RT];
        cp_info[nb_cand].lb = mv_aff[CP_LB];
        nb_cand++;
    }

    if (nb_cand < 2) {
        if (cand_mask & 0x4) {
            cp_info[nb_cand].lt = mv_aff[CP_LB];
            cp_info[nb_cand].rt = mv_aff[CP_LB];
            cp_info[nb_cand].lb = mv_aff[CP_LB];
            nb_cand++;
        }
    }

    if (nb_cand < 2) {
        if (cand_mask & 0x2) {
            cp_info[nb_cand].lt = mv_aff[CP_RT];
            cp_info[nb_cand].rt = mv_aff[CP_RT];
            cp_info[nb_cand].lb = mv_aff[CP_RT];
            nb_cand++;
        }
    }

    if (nb_cand < 2) {
        if (cand_mask & 0x1) {
            cp_info[nb_cand].lt = mv_aff[CP_LT];
            cp_info[nb_cand].rt = mv_aff[CP_LT];
            cp_info[nb_cand].lb = mv_aff[CP_LT];
            nb_cand++;
        }
    }

    /* TMVP candidate */
    if (nb_cand < 2 && inter_ctx->tmvp_enabled) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        const struct TMVPPos pos = compute_tmpv_coord(pb_info, 2);

        if (!inter_ctx->tmvp_avail) {
            /* FIXME thread synchro */
            /*FIXME dirty ref to ctudec */
            OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
            load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
        }

        uint8_t cand_msk = check_tmvp_cand(tmvp->dir_map_v0, tmvp->dir_map_v1, pos);

        if (cand_msk) {
            OVMV col_mv = {0};
            uint8_t avail = 0;
            uint8_t rpl_idx  = inter_dir - 1;
            uint8_t col_ref_l0 = tmvp->col_ref_l0;
            OVMV dst;
            if ((!col_ref_l0 && !tmvp->ldc) || (tmvp->ldc && rpl_idx == RPL_0)) {
                avail = tmvp_from_l0(inter_ctx, tmvp, pos, rpl_idx, ref_idx, cand_msk, &dst);
            } else {
                avail = tmvp_from_l1(inter_ctx, tmvp, pos, rpl_idx, ref_idx, cand_msk, &dst);
            }

            if (avail) {
                col_mv = round_affine_mv(dst, prec_amvr);
                col_mv.ref_idx = ref_idx;

                cp_info[nb_cand].lt = col_mv;
                cp_info[nb_cand].rt = col_mv;
                cp_info[nb_cand].lb = col_mv;

                nb_cand++;
            }
        }

    }

    if (nb_cand < 2) {
        for (int i = nb_cand; i < MAX_NB_AMVP_CAND; i++) {
            OVMV zmv = {0};
            cp_info[nb_cand].lt = zmv;
            cp_info[nb_cand].rt = zmv;
            cp_info[nb_cand].lb = zmv;
            nb_cand++;
        }
    }

    /* FIXME round needed? */
    /* Round control points */
    for (int i = 0; i < nb_cand; i++) {
        cp_info[i].lt = round_affine_mv(cp_info[i].lt, prec_amvr);
        cp_info[i].rt = round_affine_mv(cp_info[i].rt, prec_amvr);
        if (affine_type)
        cp_info[i].lb = round_affine_mv(cp_info[i].lb, prec_amvr);
    }
    return cp_info[mvp_idx];
}

/* SBTMVP related function
 */

#define LOG2_SBTMVP_S 3

#define TMVP_POS_MSK (~(0x7))

struct MergeInfo {
    OVMV mv0;
    OVMV mv1;
    uint8_t inter_dir;
};

struct OVPos
{
    int16_t x;
    int16_t y;
};

static inline OVMV
mv_internal_to_integer(OVMV mv)
{
    mv.x += 7 + (mv.x < 0);
    mv.y += 7 + (mv.y < 0);
    mv.x >>= 4;
    mv.y >>= 4;
    return mv;
}

/* FIXME truncated CTU */
static inline struct OVPos
clip_sb_pos_to_col_ctu(struct OVPos pos, int16_t ctu_w, int16_t ctu_h, uint8_t is_bnd)
{
  struct OVPos sb_pos;

  sb_pos.x = ov_clip(pos.x, 0, ctu_w + 3 - (is_bnd << 2));
  sb_pos.y = ov_clip(pos.y, 0, ctu_h - 1);

  sb_pos.x &= TMVP_POS_MSK;
  sb_pos.y &= TMVP_POS_MSK;

  return sb_pos;
}

static inline struct OVPos
derive_sbtmvp_cand_pos(uint8_t x0, uint8_t y0, uint8_t log2_pu_w, uint8_t log2_pu_h,
                       OVMV mv_offset, int16_t ctu_w, int16_t ctu_h, uint8_t is_bnd)
{
    struct OVPos center_pos;

    center_pos.x = x0 + ((1 << log2_pu_w) >> 1);
    center_pos.y = y0 + ((1 << log2_pu_h) >> 1);

    center_pos.x += mv_offset.x;
    center_pos.y += mv_offset.y;

    center_pos = clip_sb_pos_to_col_ctu(center_pos, ctu_w, ctu_h, is_bnd);

    return center_pos;
}

static OVMV
derive_sbtmvp_mv_offset(const struct InterDRVCtx *inter_ctx,
                        const struct PBInfo *const pb_info,
                        uint8_t cand_rpl0, uint8_t cand_rpl1)
{
    uint8_t avail_rpl0 = cand_rpl0 & A1_MSK;
    uint8_t avail_rpl1 = cand_rpl1 & A1_MSK;
    OVMV a1_mv;

    if (avail_rpl0 | avail_rpl1) {
        const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
        const uint16_t a1_pos = derive_cand_position(*pb_info, A1);
        /* FIXME derive this */
        const struct ColInfo *const col_info = &tmvp->col_info;

        if (avail_rpl0) {
            a1_mv = inter_ctx->mv_ctx0.mvs[a1_pos];
            uint8_t a1_ref0_is_col_pic = a1_mv.ref_idx == col_info->ref_idx_rpl0;
            if (a1_ref0_is_col_pic) {
                goto found;
            }
        }

        if (avail_rpl1) {
            a1_mv = inter_ctx->mv_ctx1.mvs[a1_pos];
            uint8_t a1_ref1_is_col_pic = a1_mv.ref_idx == col_info->ref_idx_rpl1;
            if (a1_ref1_is_col_pic) {
                goto found;
            }
        }
    }

    OVMV zmv = {0};

    return zmv;

found:
    a1_mv = mv_internal_to_integer(a1_mv);
    return a1_mv;
}

static uint8_t
sbtmvp_from_ldc(const struct InterDRVCtx *inter_ctx, const struct VVCTMVP *const tmvp, struct OVPos pos,
                uint8_t rpl_idx, uint8_t ref_idx,
                uint8_t cand_msk, OVMV *const dst)
{
    int32_t dist_ref = rpl_idx == RPL_0 ? inter_ctx->dist_ref_0[ref_idx]
                                        : inter_ctx->dist_ref_1[ref_idx];

    int32_t dist_ref_opp = rpl_idx == RPL_0 ? inter_ctx->dist_ref_1[ref_idx]
                                            : inter_ctx->dist_ref_0[ref_idx];
    uint8_t cand_c0  = cand_msk & 0x1;
    uint8_t cand_c01 = cand_msk & 0x2;

    OVMV mv;

    int32_t dist_col;
    int16_t scale;

    if (cand_c0 | cand_c01) {
        int16_t c0_pos = TMVP_POS_IN_BUF((pos.x >> 2), (pos.y >> 2));
        uint8_t dir = 0;

        if (cand_c0 && cand_c01 && !tmvp->col_ref_l0) {
            const OVMV *mvs    = tmvp->mvs0;
            const int16_t *dist_cols = tmvp->dist_col_0;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[0] = mv;

            mv       = tmvp->mvs1[c0_pos];
            dist_col = tmvp->dist_col_1[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        } else if (cand_c0 && cand_c01) {
            const OVMV *mvs    = tmvp->mvs1;
            const int16_t *dist_cols = tmvp->dist_col_1;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[0] = mv;

            mv       = tmvp->mvs0[c0_pos];
            dist_col = tmvp->dist_col_0[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        } else if (cand_c0) {
            const OVMV *mvs    = tmvp->mvs0;
            const int16_t *dist_cols = tmvp->dist_col_0;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];
            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[0] = mv;

            mv       = mvs[c0_pos];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            mv.ref_idx = ref_idx;

            dst[1] = mv;

        } else if (cand_c01) {
            const OVMV *mvs    = tmvp->mvs1;
            const int16_t *dist_cols = tmvp->dist_col_1;
            mv       = mvs[c0_pos];
            dist_col = dist_cols[mv.ref_idx];

            dir |= 0x2;

            scale = derive_tmvp_scale(dist_ref, dist_col);

            mv = tmvp_rescale(mv, scale);

            mv.ref_idx = ref_idx;
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            dst[0] = mv;

            mv       = mvs[c0_pos];

            dir |= 0x1;

            scale = derive_tmvp_scale(dist_ref_opp, dist_col);

            mv = tmvp_rescale(mv, scale);
            mv.bcw_idx_plus1 = 0;
            mv.prec_amvr = 0;

            mv.ref_idx = ref_idx;

            dst[1] = mv;
        }

        return dir;
    }

    return 0;
}

static uint8_t
sbtmvp_from_same_rpl(const struct InterDRVCtx *const inter_ctx, const struct VVCTMVP *const tmvp, struct OVPos pos,
                     uint8_t rpl_idx, uint8_t ref_idx,
                     uint8_t cand_msk, OVMV *const dst)
{

    int32_t dist_ref;
    const OVMV *mvs;
    const int16_t *dist_cols;
    uint8_t avail;
    OVMV mv;

    int32_t dist_col;
    int16_t scale;


    if (rpl_idx == RPL_0) {
        dist_ref  = inter_ctx->dist_ref_0[ref_idx];
        dist_cols = tmvp->dist_col_0;
        mvs       = tmvp->mvs0;

        avail  = cand_msk & 0x1;
    } else {
        dist_ref  = inter_ctx->dist_ref_1[ref_idx];
        dist_cols = tmvp->dist_col_1;
        mvs       = tmvp->mvs1;

        avail = cand_msk & 0x2;
    }

    if (avail) {
        int16_t c0_pos = TMVP_POS_IN_BUF((pos.x >> 2), (pos.y >> 2));

        mv       = mvs[c0_pos];
        dist_col = dist_cols[mv.ref_idx];

        goto found;
    }

    return 0;

found :
    scale = derive_tmvp_scale(dist_ref, dist_col);

    mv = tmvp_rescale(mv, scale);

    mv.ref_idx = ref_idx;
    mv.bcw_idx_plus1 = 0;
    mv.prec_amvr = 0;

    *dst = mv;

    return 1;
}

static inline uint8_t
check_sbtmvp_cand(const uint64_t *v_map0, const uint64_t *v_map1,
                  struct OVPos pos)
{
    /* Derive availability based on CTB inter fields */
    uint64_t col_rpl0 = v_map0[(pos.x >> 2) + 1];
    uint64_t col_rpl1 = v_map1[(pos.x >> 2) + 1];

    uint8_t cand_rpl0 = !!(col_rpl0 & TMVP_POS_MASK((pos.y >> 2)));
    uint8_t cand_rpl1 = !!(col_rpl1 & TMVP_POS_MASK((pos.y >> 2)));

    uint8_t cand_msk = cand_rpl0;
    cand_msk |= cand_rpl1 << 1;

    return cand_msk;
}

uint8_t
derive_sub_pu_merge_cand(const struct InterDRVCtx *inter_ctx,
                         uint8_t x0, uint8_t y0,
                         uint8_t log2_pu_w, uint8_t log2_pu_h,
                         struct MergeInfo *const mv_info,
                         OVMV *const a1_mv,
                         uint8_t cand_rpl0, uint8_t cand_rpl1)
{
    const struct PBInfo pb_info = {
        .x_pb = x0 >> 2,
        .y_pb = y0 >> 2,
        .nb_pb_w = (1 << log2_pu_w) >> 2,
        .nb_pb_h = (1 << log2_pu_h) >> 2,
        .log2_w = log2_pu_w,
        .log2_h = log2_pu_h
    };
    const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
    uint8_t is_bnd = tmvp->ctudec->ctb_x == tmvp->ctudec->nb_ctb_pic_w - 1;

    int16_t ctu_w = tmvp->ctu_w;
    int16_t ctu_h = tmvp->ctu_h;

    OVMV mv_offset = derive_sbtmvp_mv_offset(inter_ctx, &pb_info, cand_rpl0, cand_rpl1);

    struct OVPos center_pos = derive_sbtmvp_cand_pos(x0, y0, log2_pu_w, log2_pu_h, mv_offset,
                                                     ctu_w, ctu_h, is_bnd);

    uint8_t inter_dir = 0;

    if (!inter_ctx->tmvp_avail) {
        OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
        load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
    }

    uint8_t cand_msk = check_sbtmvp_cand(tmvp->dir_map_v0, tmvp->dir_map_v1, center_pos);

    if (cand_msk) {
        if (tmvp->ldc) {
            OVMV col_mv[2];
            inter_dir  = sbtmvp_from_ldc(inter_ctx, tmvp, center_pos,
                                         RPL_0, 0, cand_msk, col_mv);

            mv_info->mv0 = col_mv[0];
            mv_info->mv1 = col_mv[1];
        } else {
            OVMV col_mv[2];
            inter_dir  = sbtmvp_from_same_rpl(inter_ctx, tmvp, center_pos,
                                              RPL_0, 0, cand_msk, col_mv);
            inter_dir |= sbtmvp_from_same_rpl(inter_ctx, tmvp, center_pos,
                                              RPL_1, 0, cand_msk, &col_mv[1]) << 1;
            mv_info->mv0 = col_mv[0];
            mv_info->mv1 = col_mv[1];
        }
    }

    mv_info->inter_dir = inter_dir;
    *a1_mv = mv_offset;

    return !!inter_dir;
}

void
derive_sub_block_mvs(struct InterDRVCtx *inter_ctx,
                     const struct VVCTMVP *tmvp,
                     uint8_t x0, uint8_t y0,
                     uint8_t log2_pu_w, uint8_t log2_pu_h,
                     OVMV mv_offset, const struct MergeInfo *const main_mv)
{
    int nb_sb_w = OVMAX((1 << log2_pu_w) >> LOG2_SBTMVP_S, 1);
    int nb_sb_h = OVMAX((1 << log2_pu_h) >> LOG2_SBTMVP_S, 1);
    uint16_t ctu_w = tmvp->ctu_w;
    uint16_t ctu_h = tmvp->ctu_h;
    uint8_t is_bnd = tmvp->ctudec->ctb_x == tmvp->ctudec->nb_ctb_pic_w - 1;

    /* FIXME check if this clipping is needed */
    int sb_h = nb_sb_h == 1 ? 1 << log2_pu_h : 1 << LOG2_SBTMVP_S;
    int sb_w = nb_sb_w == 1 ? 1 << log2_pu_w : 1 << LOG2_SBTMVP_S;

    const uint8_t is_small = log2_pu_h + log2_pu_w <= 5;
    uint8_t x_pb = x0 >> 2;
    uint8_t y_pb = y0 >> 2;

    OVMV *mv_buff0 = &inter_ctx->mv_ctx0.mvs[35 + (x0 >> 2) + (y0 >> 2) * 34];
    OVMV *mv_buff1 = &inter_ctx->mv_ctx1.mvs[35 + (x0 >> 2) + (y0 >> 2) * 34];
    OVMV *tmvp_mv0 = &inter_ctx->tmvp_mv[0].mvs[((x0 + 4) >> 3) + (((y0 + 4) >> 3) << 4)];
    OVMV *tmvp_mv1 = &inter_ctx->tmvp_mv[1].mvs[((x0 + 4) >> 3) + (((y0 + 4) >> 3) << 4)];

    /* FIXME check start_x , start_y */
    int start_x = x0 + (sb_w >> 1) + mv_offset.x;
    int start_y = y0 + (sb_h >> 1) + mv_offset.y;

    int x, y;

    int i, j;

    for (y = start_y, i = 0; i < nb_sb_h; ++i, y += sb_h) {
        for (x = start_x, j = 0; j < nb_sb_w; ++j, x += sb_w) {
            struct OVPos col_pos = {
                .x = x,
                .y = y
            };

            uint8_t inter_dir = 0;

            OVMV mv0;
            OVMV mv1;

            /* FIXME avoid clipping + correct CTU dimension
             */
            col_pos = clip_sb_pos_to_col_ctu(col_pos, ctu_w, ctu_h, is_bnd);

            uint8_t cand_msk = check_sbtmvp_cand(tmvp->dir_map_v0, tmvp->dir_map_v1, col_pos);

            if (cand_msk) {
                if (tmvp->ldc) {
                    OVMV col_mv[2];
                    inter_dir  = sbtmvp_from_ldc(inter_ctx, tmvp, col_pos, RPL_0, 0, cand_msk, col_mv);

                    mv0 = col_mv[0];
                    mv1 = col_mv[1];
                } else {
                    OVMV col_mv[2];
                    inter_dir  = sbtmvp_from_same_rpl(inter_ctx, tmvp, col_pos, RPL_0, 0, cand_msk, col_mv);
                    inter_dir |= sbtmvp_from_same_rpl(inter_ctx, tmvp, col_pos, RPL_1, 0, cand_msk, &col_mv[1]) << 1;
                    mv0 = col_mv[0];
                    mv1 = col_mv[1];
                }
            }

            if (!inter_dir) {
                mv0 = main_mv->mv0;
                mv1 = main_mv->mv1;

                inter_dir = main_mv->inter_dir;
            }

            /* Disable B prediction on small blocks */
            if (is_small && inter_dir == 0x3) {
                inter_dir = 1;
            }

            if (inter_dir & 0x1) {
                ctu_field_set_rect_bitfield(&inter_ctx->mv_ctx0.map,
                                            x_pb + 2 * j, y_pb + 2 * i,
                                            2, 2);
                mv0.ref_idx = 0;
                mv0.bcw_idx_plus1 = 0;
                tmvp_mv0[j] = mv0;

                mv_buff0[j * 2]     = mv0;
                mv_buff0[j * 2 + 1] = mv0;
                mv_buff0[34 + j * 2]     = mv0;
                mv_buff0[34 + j * 2 + 1] = mv0;
            }

            if (inter_dir & 0x2) {
                ctu_field_set_rect_bitfield(&inter_ctx->mv_ctx1.map,
                                            x_pb + 2 * j, y_pb + 2 * i,
                                            2, 2);
                mv1.ref_idx = 0;
                mv1.bcw_idx_plus1 = 0;

                tmvp_mv1[j] = mv1;

                mv_buff1[j * 2]     = mv1;
                mv_buff1[j * 2 + 1] = mv1;
                mv_buff1[34 + j * 2]     = mv1;
                mv_buff1[34 + j * 2 + 1] = mv1;
            }

            /* FIXME Move somewhere else */
            rcn_mcp_b(tmvp->ctudec, tmvp->ctudec->rcn_ctx.ctu_buff,
                      inter_ctx, tmvp->ctudec->part_ctx,
                      mv0, mv1, x0 + 8 * j, y0 + 8 * i,
                      3, 3, inter_dir, 0, 0);

        }
        mv_buff0 += 34 * 2;
        mv_buff1 += 34 * 2;
        tmvp_mv0 += 16;
        tmvp_mv1 += 16;
    }
}

/* Affine merge related functions
 */

struct AffineMergeInfo
{
     struct AffineControlInfo cinfo[2];
     uint8_t inter_dir;
     uint8_t affine_type;
};

static inline uint8_t
avaiable_merge_affine_a0_a1(uint8_t cand_list_affine)
{
    if (cand_list_affine & A0_MSK) {
        return A0_MSK;
    }

    return (cand_list_affine & A1_MSK);
}

static inline int
available_merge_affine_b0_b1_b2(uint8_t cand_list_affine)
{
    if (cand_list_affine & B0_MSK) {
        return B0_MSK;
    }

    if (cand_list_affine & B1_MSK) {
        return B1_MSK;
    }

    return (cand_list_affine & B2_MSK);
}

static inline uint8_t
available_merge_b2_b3_a2(uint8_t cand_list_rpl0, uint8_t cand_list_rpl1)
{
    if ((cand_list_rpl0 | cand_list_rpl1) & B2_MSK) {
        return B2_MSK;
    }

    if ((cand_list_rpl0 | cand_list_rpl1) & B3_MSK) {
        return B3_MSK;
    }

    return ((cand_list_rpl0 | cand_list_rpl1) & A2_MSK);
}

/* TODO return for both list or only one */
static inline uint8_t
available_merge_b1_b0(uint8_t cand_list_rpl0, uint8_t cand_list_rpl1)
{
    if ((cand_list_rpl0 | cand_list_rpl1) & B1_MSK) {
        return B1_MSK;
    }

    return ((cand_list_rpl0 | cand_list_rpl1) & B0_MSK);
}

static inline uint8_t
available_merge_a1_a0(uint8_t cand_list_rpl0, uint8_t cand_list_rpl1)
{
    if ((cand_list_rpl0 | cand_list_rpl1) & A1_MSK) {
        return A1_MSK;
    }

    return ((cand_list_rpl0 | cand_list_rpl1) & A0_MSK);
}

static inline uint8_t
check_avail_cp(uint8_t avail_cp_map, enum ControlPointCandMask msk)
{
    return (msk & avail_cp_map) == msk;
}

struct ControlPointMVCand
{
     OVMV mv0[4];
     OVMV mv1[4];
     uint8_t dir;
};

static uint8_t
derive_affine_control_point_0(struct ControlPointMVCand mi, int model_idx,
                              uint8_t log2_cu_w, uint8_t log2_cu_h,
                              struct AffineMergeInfo *const aff_mrg_ctx)
{
    /* FIXME dir info should be derived from RPL cand list + ref_idx*/
    OVMV mv0[4];
    OVMV mv1[4];

    uint8_t dir = 0;

    switch (model_idx)
    {
        case 0:

            if (mi.mv0[CP_LT].ref_idx >= 0 && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_RT].ref_idx) {

                mv0[CP_LT] = mi.mv0[CP_LT];
                mv0[CP_RT] = mi.mv0[CP_RT];

                dir |= 0x1;
            }

            if (mi.mv1[CP_LT].ref_idx >= 0 && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_RT].ref_idx) {

                mv1[CP_LT] = mi.mv1[CP_LT];
                mv1[CP_RT] = mi.mv1[CP_RT];

                dir |= 0x2;
            }

            break;

        case 1:

            if (mi.mv0[CP_LT].ref_idx >= 0 && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_LB].ref_idx) {
                int shiftHtoW = AFFINE_SHIFT + log2_cu_w - log2_cu_h;
                OVMV tmp;

                mv0[CP_LT] = mi.mv0[CP_LT];
                mv0[CP_LB] = mi.mv0[CP_LB];

                tmp.x = (mv0[0].x << AFFINE_SHIFT) + ((mv0[2].y - mv0[0].y) << shiftHtoW);
                tmp.y = (mv0[0].y << AFFINE_SHIFT) - ((mv0[2].x - mv0[0].x) << shiftHtoW);

                mv0[1] = round_affine_mv2(tmp);

                mv0[1] = clip_mv(mv0[1]);

                dir |= 0x1;
            }

            if (mi.mv1[CP_LT].ref_idx >= 0 && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_LB].ref_idx) {
                int shift_hw = AFFINE_SHIFT + log2_cu_w - log2_cu_h;
                OVMV tmp;

                mv1[CP_LT] = mi.mv1[CP_LT];
                mv1[CP_LB] = mi.mv1[CP_LB];

                tmp.x = (mv1[0].x << AFFINE_SHIFT) + ((mv1[2].y - mv1[0].y) << shift_hw);
                tmp.y = (mv1[0].y << AFFINE_SHIFT) - ((mv1[2].x - mv1[0].x) << shift_hw);

                mv1[1] = round_affine_mv2(tmp);

                mv1[1] = clip_mv(mv1[1]);

                dir |= 0x2;
            }

            break;
    }

    if (dir == 0) {
        return 0;
    }

    aff_mrg_ctx[0].cinfo[0].lt = mv0[0];
    aff_mrg_ctx[0].cinfo[0].rt = mv0[1];
    aff_mrg_ctx[0].cinfo[0].lb = mv0[2];

    aff_mrg_ctx[0].cinfo[1].lt = mv1[0];
    aff_mrg_ctx[0].cinfo[1].rt = mv1[1];
    aff_mrg_ctx[0].cinfo[1].lb = mv1[2];

    aff_mrg_ctx[0].inter_dir   = dir;
    aff_mrg_ctx[0].affine_type = AFFINE_2CP;

    return 1;
}

static uint8_t
derive_affine_control_point_1(struct ControlPointMVCand mi, int model_idx,
                              struct AffineMergeInfo *const aff_mrg_ctx)
{
    OVMV mv0[4] = {0};
    OVMV mv1[4] = {0};

    uint8_t dir = 0;

    switch (model_idx)
    {
        case 0:

            if (mi.mv0[CP_LT].ref_idx >= 0 && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_RT].ref_idx
                                           && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_LB].ref_idx) {
                mv0[0] = mi.mv0[0];
                mv0[1] = mi.mv0[1];
                mv0[2] = mi.mv0[2];

                mv0[1].ref_idx = mi.mv0[0].ref_idx;
                mv0[2].ref_idx = mi.mv0[0].ref_idx;

                mv0[2].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[2].prec_amvr = mi.mv0[CP_LT].prec_amvr;
                mv0[0].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[0].prec_amvr = mi.mv0[CP_LT].prec_amvr;
                mv0[1].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[1].prec_amvr = mi.mv0[CP_LT].prec_amvr;

                dir |= 0x1;
            }

            if (mi.mv1[CP_LT].ref_idx >= 0 && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_RT].ref_idx
                                           && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_LB].ref_idx) {

                mv1[0] = mi.mv1[0];
                mv1[1] = mi.mv1[1];
                mv1[2] = mi.mv1[2];
                mv1[1].ref_idx = mi.mv1[0].ref_idx;
                mv1[2].ref_idx = mi.mv1[0].ref_idx;

                mv1[2].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[2].prec_amvr = mi.mv1[CP_LT].prec_amvr;
                mv1[0].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[0].prec_amvr = mi.mv1[CP_LT].prec_amvr;
                mv1[1].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[1].prec_amvr = mi.mv1[CP_LT].prec_amvr;

                dir |= 0x2;
            }

            break;

        case 1:
            if (mi.mv0[CP_LT].ref_idx >= 0 && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_RT].ref_idx
                                          && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_RB].ref_idx) {
                mv0[0] = mi.mv0[0];
                mv0[1] = mi.mv0[1];
                mv0[3] = mi.mv0[3];

                mv0[2].x = mv0[0].x - mv0[1].x + mv0[3].x;
                mv0[2].y = mv0[0].y - mv0[1].y + mv0[3].y;

                mv0[2] = clip_mv(mv0[2]);

                mv0[2].ref_idx = mi.mv0[CP_LT].ref_idx;
                mv0[2].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[2].prec_amvr = mi.mv0[CP_LT].prec_amvr;
                mv0[0].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[0].prec_amvr = mi.mv0[CP_LT].prec_amvr;
                mv0[1].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[1].prec_amvr = mi.mv0[CP_LT].prec_amvr;

                dir |= 0x1;
            }

            if (mi.mv1[CP_LT].ref_idx >= 0 && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_RT].ref_idx
                                          && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_RB].ref_idx) {

                mv1[0] = mi.mv1[0];
                mv1[1] = mi.mv1[1];
                mv1[3] = mi.mv1[3];

                mv1[2].x = mv1[0].x - mv1[1].x + mv1[3].x;
                mv1[2].y = mv1[0].y - mv1[1].y + mv1[3].y;

                mv1[2] = clip_mv(mv1[2]);

                mv1[2].ref_idx = mi.mv1[CP_LT].ref_idx;
                mv1[2].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[2].prec_amvr = mi.mv1[CP_LT].prec_amvr;
                mv1[0].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[0].prec_amvr = mi.mv1[CP_LT].prec_amvr;
                mv1[1].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[1].prec_amvr = mi.mv1[CP_LT].prec_amvr;

                dir |= 0x2;
            }

            break;

        case 2:
            if (mi.mv0[CP_LT].ref_idx >= 0 && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_LB].ref_idx
                                          && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_RB].ref_idx) {
                mv0[0] = mi.mv0[0];
                mv0[2] = mi.mv0[2];
                mv0[3] = mi.mv0[3];

                mv0[1].x = mv0[0].x - mv0[2].x + mv0[3].x;
                mv0[1].y = mv0[0].y - mv0[2].y + mv0[3].y;

                mv0[1] = clip_mv(mv0[1]);

                mv0[1].ref_idx = mi.mv0[CP_LT].ref_idx;
                mv0[1].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[1].prec_amvr = mi.mv0[CP_LT].prec_amvr;
                mv0[0].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[0].prec_amvr = mi.mv0[CP_LT].prec_amvr;
                mv0[2].bcw_idx_plus1 = mi.mv0[CP_LT].bcw_idx_plus1;
                mv0[2].prec_amvr = mi.mv0[CP_LT].prec_amvr;

                dir |= 0x1;
            }

            if (mi.mv1[CP_LT].ref_idx >= 0 && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_LB].ref_idx
                                          && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_RB].ref_idx) {

                mv1[0] = mi.mv1[0];
                mv1[2] = mi.mv1[2];
                mv1[3] = mi.mv1[3];

                mv1[1].x = mv1[0].x - mv1[2].x + mv1[3].x;
                mv1[1].y = mv1[0].y - mv1[2].y + mv1[3].y;

                mv1[1] = clip_mv(mv1[1]);

                mv1[1].ref_idx = mi.mv1[CP_LT].ref_idx;
                mv1[1].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[1].prec_amvr = mi.mv1[CP_LT].prec_amvr;
                mv1[0].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[0].prec_amvr = mi.mv1[CP_LT].prec_amvr;
                mv1[2].bcw_idx_plus1 = mi.mv1[CP_LT].bcw_idx_plus1;
                mv1[2].prec_amvr = mi.mv1[CP_LT].prec_amvr;

                dir |= 0x2;
            }

            break;

        case 3:
            if (mi.mv0[CP_RT].ref_idx >= 0 && mi.mv0[CP_RT].ref_idx == mi.mv0[CP_LB].ref_idx
                                          && mi.mv0[CP_RT].ref_idx == mi.mv0[CP_RB].ref_idx) {
                mv0[1] = mi.mv0[1];
                mv0[2] = mi.mv0[2];
                mv0[3] = mi.mv0[3];

                mv0[0].x = mv0[1].x + mv0[2].x - mv0[3].x;
                mv0[0].y = mv0[1].y + mv0[2].y - mv0[3].y;

                mv0[0] = clip_mv(mv0[0]);

                mv0[0].ref_idx = mi.mv0[CP_RT].ref_idx;
                mv0[0].bcw_idx_plus1 = mi.mv0[CP_RT].bcw_idx_plus1;
                mv0[0].prec_amvr = mi.mv0[CP_RT].prec_amvr;
                mv0[1].bcw_idx_plus1 = mi.mv0[CP_RT].bcw_idx_plus1;
                mv0[1].prec_amvr = mi.mv0[CP_RT].prec_amvr;
                mv0[2].bcw_idx_plus1 = mi.mv0[CP_RT].bcw_idx_plus1;
                mv0[2].prec_amvr = mi.mv0[CP_RT].prec_amvr;

                dir |= 0x1;
            }

            if (mi.mv1[CP_RT].ref_idx >= 0 && mi.mv1[CP_RT].ref_idx == mi.mv1[CP_LB].ref_idx
                                          && mi.mv1[CP_RT].ref_idx == mi.mv1[CP_RB].ref_idx) {

                mv1[1] = mi.mv1[1];
                mv1[2] = mi.mv1[2];
                mv1[3] = mi.mv1[3];

                mv1[0].x = mv1[1].x + mv1[2].x - mv1[3].x;
                mv1[0].y = mv1[1].y + mv1[2].y - mv1[3].y;

                mv1[0] = clip_mv(mv1[0]);

                mv1[0].ref_idx = mi.mv1[CP_RT].ref_idx;
                mv1[0].bcw_idx_plus1 = mi.mv1[CP_RT].bcw_idx_plus1;
                mv1[0].prec_amvr = mi.mv1[CP_RT].prec_amvr;
                mv1[1].bcw_idx_plus1 = mi.mv1[CP_RT].bcw_idx_plus1;
                mv1[1].prec_amvr = mi.mv1[CP_RT].prec_amvr;
                mv1[2].bcw_idx_plus1 = mi.mv1[CP_RT].bcw_idx_plus1;
                mv1[2].prec_amvr = mi.mv1[CP_RT].prec_amvr;

                dir |= 0x2;
            }
            break;
    }

    /* Not found */
    if (dir == 0) {
        return 0;
    }

    aff_mrg_ctx[0].cinfo[0].lt = mv0[0];
    aff_mrg_ctx[0].cinfo[0].rt = mv0[1];
    aff_mrg_ctx[0].cinfo[0].lb = mv0[2];

    aff_mrg_ctx[0].cinfo[1].lt = mv1[0];
    aff_mrg_ctx[0].cinfo[1].rt = mv1[1];
    aff_mrg_ctx[0].cinfo[1].lb = mv1[2];

    aff_mrg_ctx[0].inter_dir   = dir;
    aff_mrg_ctx[0].affine_type = AFFINE_3CP;

    return 1;
}

void
derive_affine_merge_mv(struct InterDRVCtx *const inter_ctx,
                       struct AffineDRVInfo *affine_ctx,
                       struct AffineMergeInfo *const aff_mrg_ctx,
                       uint8_t x_pb, uint8_t y_pb,
                       uint8_t nb_pb_w, uint8_t nb_pb_h,
                       uint8_t log2_cu_w, uint8_t log2_cu_h,
                       const int mrg_idx)
{
    uint8_t nb_cand = 0;

    struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    /* FIXME missing tool SBTMVP */

    uint64_t aff_lft_col = affine_ctx->map.vfield[x_pb];
    uint64_t aff_abv_row = affine_ctx->map.hfield[y_pb];
    uint64_t rpl0_lft_col = mv_ctx0->map.vfield[x_pb];
    uint64_t rpl0_abv_row = mv_ctx0->map.hfield[y_pb];
    uint64_t rpl1_lft_col = mv_ctx1->map.vfield[x_pb];
    uint64_t rpl1_abv_row = mv_ctx1->map.hfield[y_pb];

    struct PBInfo pb_info = { 
        .x_pb = x_pb, 
        .y_pb = y_pb, 
        .nb_pb_w = nb_pb_w,
        .nb_pb_h = nb_pb_h,
        .log2_w = log2_cu_w,
        .log2_h = log2_cu_h
    };

    const uint8_t aff_cand_list = check_cand_available(aff_abv_row, aff_lft_col, x_pb, y_pb,
                                                       nb_pb_w, nb_pb_h);

    const uint8_t rpl0_cand = check_cand_available(rpl0_abv_row, rpl0_lft_col, x_pb, y_pb,
                                                   nb_pb_w, nb_pb_h);

    const uint8_t rpl1_cand = check_cand_available(rpl1_abv_row, rpl1_lft_col, x_pb, y_pb,
                                                   nb_pb_w, nb_pb_h);

    /* derive affine control points from affine neighbours */
    uint8_t cand_lft = avaiable_merge_affine_a0_a1(aff_cand_list);
    uint8_t cand_abv = available_merge_affine_b0_b1_b2(aff_cand_list);

    if (cand_lft) {
        enum CandName cand_name = cand_mask_to_idx(cand_lft);
        const int16_t cand_pos = derive_cand_position(pb_info, cand_name);
        const struct AffineInfo *const affine_info = &affine_ctx->affine_info[cand_pos];
        const enum AffineType affine_type = affine_info->type;

        struct AffineControlInfo cp_info0 = {0};
        struct AffineControlInfo cp_info1 = {0};

        uint8_t dir = (!!(rpl0_cand & cand_lft)) | ((!!(rpl1_cand & cand_lft)) << 1);

        if (dir & 0x1) {
            struct AffineControlInfo cand_cp_info0 = affine_info->cps[RPL_0];
            struct PBInfo ngh_pb = affine_info->pb;
            cp_info0 = derive_cp_from_cand(&cand_cp_info0, pb_info, ngh_pb, affine_type,
                                           affine_info->type, cand_name);
        }

        /* Note we do not check for B slice since dir should already be 1 if P slice */
        if (dir & 0x2) {
            struct AffineControlInfo cand_cp_info1 = affine_info->cps[RPL_1];
            struct PBInfo ngh_pb = affine_info->pb;
            cp_info1 = derive_cp_from_cand(&cand_cp_info1, pb_info, ngh_pb, affine_type,
                                           affine_info->type, cand_name);
        }

        aff_mrg_ctx[0].cinfo[0] = cp_info0;
        aff_mrg_ctx[0].cinfo[1] = cp_info1;

        aff_mrg_ctx[0].inter_dir   = dir;
        aff_mrg_ctx[0].affine_type = affine_type;

        if (nb_cand == mrg_idx) {
            return;
        }

        nb_cand++;
    }


    if (cand_abv) {
        enum CandName cand_name = cand_mask_to_idx(cand_abv);
        const int16_t cand_pos = derive_cand_position(pb_info, cand_name);
        const struct AffineInfo *const affine_info = &affine_ctx->affine_info[cand_pos];
        const enum AffineType affine_type = affine_info->type;

        struct AffineControlInfo cp_info0 = {0};
        struct AffineControlInfo cp_info1 = {0};

        uint8_t dir = (!!(rpl0_cand & cand_abv)) | ((!!(rpl1_cand & cand_abv)) << 1);

        if (dir & 0x1) {
            struct AffineControlInfo cand_cp_info0 = affine_info->cps[RPL_0];
            struct PBInfo ngh_pb = affine_info->pb;
            cp_info0 = derive_cp_from_cand(&cand_cp_info0, pb_info, ngh_pb, affine_type,
                                           affine_info->type, cand_name);
        }

        /* Note we do not check for B slice since dir should be 1 if P slice */
        if (dir & 0x2) {
            struct AffineControlInfo cand_cp_info1 = affine_info->cps[RPL_1];
            struct PBInfo ngh_pb = affine_info->pb;
            cp_info1 = derive_cp_from_cand(&cand_cp_info1, pb_info, ngh_pb, affine_type,
                                           affine_info->type, cand_name);
        }

        aff_mrg_ctx[0].cinfo[0] = cp_info0;
        aff_mrg_ctx[0].cinfo[1] = cp_info1;

        aff_mrg_ctx[0].inter_dir   = dir;
        aff_mrg_ctx[0].affine_type = affine_type;

        if (nb_cand == mrg_idx) {
            return;
        }

        nb_cand++;
    }

    /* derive affine control points from neighbours */
    {
        uint8_t avail_cp_map = 0;
        uint8_t cand_msk;
        uint8_t dir[4];

        struct ControlPointMVCand mi;

        if ((cand_msk = available_merge_b2_b3_a2(rpl0_cand, rpl1_cand))) {
            enum CandName cand_id = cand_mask_to_idx(cand_msk);
            uint16_t pos = derive_cand_position(pb_info, cand_id);
            avail_cp_map |= 0x1;
            dir[0]  = !!(cand_msk & rpl0_cand);
            dir[0] |= (!!(cand_msk & rpl1_cand)) << 1;

            if (dir[0] & 0x1) {
                mi.mv0[0] = mv_ctx0->mvs[pos];
            } else {
                /* FIXME this is only so  we check for ref_idx < 0 in CPInfo derivation */
                mi.mv0[0].ref_idx = -1;
            }

            if (dir[0] & 0x2) {
                mi.mv1[0] = mv_ctx1->mvs[pos];
            } else {
                mi.mv1[0].ref_idx = -1;
            }

        }

        if ((cand_msk = available_merge_b1_b0(rpl0_cand, rpl1_cand))) {
            enum CandName cand_id = cand_mask_to_idx(cand_msk);
            uint16_t pos = derive_cand_position(pb_info, cand_id);
            avail_cp_map |= 0x2;
            dir[1]  = !!(cand_msk & rpl0_cand);
            dir[1] |= (!!(cand_msk & rpl1_cand)) << 1;

            if (dir[1] & 0x1) {
                mi.mv0[1] = mv_ctx0->mvs[pos];
            } else {
                mi.mv0[1].ref_idx = -1;
            }


            if (dir[1] & 0x2) {
                mi.mv1[1] = mv_ctx1->mvs[pos];
            } else {
                mi.mv1[1].ref_idx = -1;
            }

        }

        if ((cand_msk = available_merge_a1_a0(rpl0_cand, rpl1_cand))) {
            enum CandName cand_id = cand_mask_to_idx(cand_msk);
            uint16_t pos = derive_cand_position(pb_info, cand_id);
            avail_cp_map |= 0x4;
            dir[2]  = !!(cand_msk & rpl0_cand);
            dir[2] |= (!!(cand_msk & rpl1_cand)) << 1;

            if (dir[2] & 0x1) {
                mi.mv0[2] = mv_ctx0->mvs[pos];
            } else {
                mi.mv0[2].ref_idx = -1;
            }

            if (dir[2] & 0x2) {
                mi.mv1[2] = mv_ctx1->mvs[pos];
            } else {
                mi.mv1[2].ref_idx = -1;
            }

        }

        /* FIXME test if affine type enabled so we skip TMVP when not needed ? */
        if (inter_ctx->tmvp_enabled) {
            const struct VVCTMVP *const tmvp = &inter_ctx->tmvp_ctx;
            const struct TMVPPos pos = compute_tmpv_coord(pb_info, 2);

            if (!inter_ctx->tmvp_avail) {
                OVCTUDec *ctudec = inter_ctx->tmvp_ctx.ctudec;
                load_ctb_tmvp(ctudec, ctudec->ctb_x, ctudec->ctb_y);
            }

            uint8_t cand_msk = check_tmvp_cand(tmvp->dir_map_v0, tmvp->dir_map_v1, pos);

            /* We only check for C0 candidate */
            cand_msk &= 0x3;

            if (cand_msk) {
                OVMV c0_mv[2] = {0};
                uint8_t col_ref_l0 = tmvp->col_ref_l0;
                uint8_t rpl_idx  = RPL_0;
                uint8_t avail_dir = 0;
                dir[3] = 0;

                /* FIXME inter_dir */
                if (tmvp->ldc) {
                    avail_dir = merge_tmvp_from_ldc(inter_ctx, tmvp, pos,
                                                    rpl_idx, 0, cand_msk, c0_mv);
                } else if (!col_ref_l0) {
                    avail_dir = merge_tmvp_from_l0(inter_ctx, tmvp, pos,
                                                   rpl_idx, 0, cand_msk, c0_mv);
                } else {
                    avail_dir = merge_tmvp_from_l1(inter_ctx, tmvp, pos,
                                                   rpl_idx, 0, cand_msk, c0_mv);
                }


                if (avail_dir & 0x1) {
                    c0_mv[0].ref_idx = 0;
                    mi.mv0[3] = c0_mv[0];
                    avail_cp_map |= 0x8;
                    dir[3] |= 0x1;
                } else {
                    mi.mv0[3].ref_idx = -1;
                }

                if (avail_dir & 0x2) {
                    c0_mv[1].ref_idx = 0;
                    mi.mv1[3] = c0_mv[1];
                    avail_cp_map |= 0x8;
                    dir[3] |= 0x2;
                } else {
                    mi.mv1[3].ref_idx = -1;
                }
            }
        }

        /* Set control points */
        /* FIXME plug this */
        uint8_t sps_affine_type_flag = 1;
        if (sps_affine_type_flag) {
            if (check_avail_cp(avail_cp_map, CP3_MASK0)) {
                nb_cand += derive_affine_control_point_1(mi, 0, aff_mrg_ctx);

                if (nb_cand - 1 == mrg_idx) {
                    return;
                }
            }

            /* Note check those only if TMVP is enabled */
            if (check_avail_cp(avail_cp_map, CP3_MASK1)) {
                nb_cand += derive_affine_control_point_1(mi, 1, aff_mrg_ctx);

                if (nb_cand - 1 == mrg_idx) {
                    return;
                }
            }

            if (check_avail_cp(avail_cp_map, CP3_MASK2)) {
                nb_cand += derive_affine_control_point_1(mi, 2, aff_mrg_ctx);

                if (nb_cand - 1 == mrg_idx) {
                    return;
                }
            }

            if (check_avail_cp(avail_cp_map, CP3_MASK3)) {
                nb_cand += derive_affine_control_point_1(mi, 3, aff_mrg_ctx);

                if (nb_cand - 1 == mrg_idx) {
                    return;
                }
            }
        }

        if (check_avail_cp(avail_cp_map, CP2_MASK0)) {
            nb_cand += derive_affine_control_point_0(mi, 0, log2_cu_w, log2_cu_h,
                                                     aff_mrg_ctx);

            if (nb_cand - 1 == mrg_idx) {
                return;
            }
        }

        if (check_avail_cp(avail_cp_map, CP2_MASK1)) {
            nb_cand += derive_affine_control_point_0(mi, 1, log2_cu_w, log2_cu_h,
                                                     aff_mrg_ctx);

            if (nb_cand - 1 == mrg_idx) {
                return;
            }
        }
    }

    /* FIXME check if we can return zmv directly */
    /* FIXME check ref_idx increment as in classic merge */
    while (nb_cand <= mrg_idx) {
        struct AffineControlInfo z_cinfo = {0};

        aff_mrg_ctx[0].cinfo[0] = z_cinfo;
        aff_mrg_ctx[0].cinfo[1] = z_cinfo;

        aff_mrg_ctx[0].inter_dir = 3;

        aff_mrg_ctx[0].affine_type = AFFINE_2CP;

        nb_cand++;

        if (nb_cand - 1 == mrg_idx) {
            return;
        }
    }
}

/* Affine MV derivation related functions
 */

static uint8_t
broadcast_mv(struct AffineDeltaMV delta_mv, uint8_t inter_dir)
{
    int a = delta_mv.h.x * 4;
    int b = delta_mv.h.y * 4;
    int c = delta_mv.v.x * 4;
    int d = delta_mv.v.y * 4;

    if (inter_dir == 0x3) {
        int blk_w = OVMAX(OVMAX(0, a + RND_AFF), OVMAX(c, a + c + RND_AFF))
                  - OVMIN(OVMIN(0, a + RND_AFF), OVMIN(c, a + c + RND_AFF));
        int blk_h = OVMAX(OVMAX(0, b), OVMAX(d + RND_AFF, b + d + RND_AFF))
                  - OVMIN(OVMIN(0, b), OVMIN(d + RND_AFF, b + d + RND_AFF));

        blk_w = (blk_w >> 11) + NB_TAP_PLUS3;
        blk_h = (blk_h >> 11) + NB_TAP_PLUS3;

        if (blk_w * blk_h > NB_TAP_PLUS9 * NB_TAP_PLUS9) {
            return 1;
        }

    } else {
        int blk_w = OVMAX(0, a + RND_AFF) - OVMIN(0, a + RND_AFF);
        int blk_h = OVMAX(0, b) - OVMIN(0, b);

        blk_w = (blk_w >> 11) + NB_TAP_PLUS3;
        blk_h = (blk_h >> 11) + NB_TAP_PLUS3;

        if (blk_w * blk_h > NB_TAP_PLUS9 * NB_TAP_PLUS5) {
            return 1;
        }

        blk_w = OVMAX(0, c) - OVMIN(0, c);
        blk_h = OVMAX(0, d + RND_AFF) - OVMIN(0, d + RND_AFF);

        blk_h = (blk_h >> 11) + NB_TAP_PLUS3;
        blk_w = (blk_w >> 11) + NB_TAP_PLUS3;

        if (blk_w * blk_h > NB_TAP_PLUS5 * NB_TAP_PLUS9) {
            return 1;
        }
    }
    return 0;
}

void
compute_subblock_mvs(const struct AffineControlInfo *const cinfo,
                     const struct AffineDeltaMV delta_mv,
                     OVMV *mv_buff,
                     uint8_t log2_cu_w, uint8_t log2_cu_h,
                     uint8_t mv_broad)
{
    /* Compute delta_mv from control points */
    /* TODO call before and give as an argument */

    uint8_t nb_sb_w = (1 << log2_cu_w) >> LOG2_MIN_CU_S;
    uint8_t nb_sb_h = (1 << log2_cu_h) >> LOG2_MIN_CU_S;
    uint8_t ref_idx = cinfo->lt.ref_idx;

    if (!mv_broad) {
        int i, j;
        OVMV accu_mv_v;
        OVMV mv_dst;

        accu_mv_v.x = (cinfo->lt.x << AFFINE_SHIFT) + delta_mv.h.x * HALF_SB_SIZE
                                                    + delta_mv.v.x * HALF_SB_SIZE;
        accu_mv_v.y = (cinfo->lt.y << AFFINE_SHIFT) + delta_mv.h.y * HALF_SB_SIZE
                                                    + delta_mv.v.y * HALF_SB_SIZE;

        for (i = 0; i < nb_sb_h; ++i) {
            OVMV accu_mv_h;

            accu_mv_h.x = accu_mv_v.x;
            accu_mv_h.y = accu_mv_v.y;

            for (j = 0; j < nb_sb_w; ++j) {

                mv_dst.x = accu_mv_h.x;
                mv_dst.y = accu_mv_h.y;

                mv_dst = round_affine_mv2(mv_dst);
                mv_dst = clip_mv(mv_dst);

                mv_dst.ref_idx = ref_idx;
                mv_dst.bcw_idx_plus1 = cinfo->lt.bcw_idx_plus1;
                mv_dst.prec_amvr = cinfo->lt.prec_amvr;

                mv_buff[j] = mv_dst;

                accu_mv_h.x += SB_SIZE * delta_mv.h.x;
                accu_mv_h.y += SB_SIZE * delta_mv.h.y;
            }

            accu_mv_v.x += SB_SIZE * delta_mv.v.x;
            accu_mv_v.y += SB_SIZE * delta_mv.v.y;
            mv_buff += 34;
        }

    } else {

        /* Broadcast MV computed from PU center */
        int i, j;
        OVMV center_mv;

        center_mv.x = cinfo->lt.x << AFFINE_SHIFT;
        center_mv.y = cinfo->lt.y << AFFINE_SHIFT;

        center_mv.x += (delta_mv.h.x << log2_cu_w) >> 1;
        center_mv.y += (delta_mv.h.y << log2_cu_w) >> 1;

        center_mv.x += (delta_mv.v.x << log2_cu_h) >> 1;
        center_mv.y += (delta_mv.v.y << log2_cu_h) >> 1;

        center_mv = round_affine_mv2(center_mv);

        center_mv = clip_mv(center_mv);

        center_mv.ref_idx = ref_idx;
        center_mv.bcw_idx_plus1 = cinfo->lt.bcw_idx_plus1;
        center_mv.prec_amvr = cinfo->lt.prec_amvr;

        for (i = 0; i < nb_sb_h; ++i) {
            for (j = 0; j < nb_sb_w; ++j) {
                mv_buff[j] = center_mv;
            }
            mv_buff += 34;
        }
    }
}

#define LF_MV_THRESHOLD 8
static inline uint8_t
mv_threshold_check(OVMV a, OVMV b)
{
    uint32_t abs_delta_x = abs(a.x - b.x);
    uint32_t abs_delta_y = abs(a.y - b.y);

    uint8_t chk = (abs_delta_x >= LF_MV_THRESHOLD) || (abs_delta_y >= LF_MV_THRESHOLD);

    return chk;
}

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

static void
dbf_mv_check_p(const struct InterDRVCtx *const inter_ctx,
               struct DBFInfo *const dbf_info,
               struct OVMVCtx *const mv_ctx0, struct OVMVCtx *const mv_ctx1,
               int x0_unit, int y0_unit,
               int nb_unit_w, int nb_unit_h)
{
    uint64_t unit_msk_w = (uint64_t)((uint64_t)1 << nb_unit_w) - 1llu;
    uint64_t unit_msk_h = (uint64_t)((uint64_t)1 << nb_unit_h) - 1llu;

    uint64_t abv0_msk = (mv_ctx0->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_w;
    uint64_t abv1_msk = (mv_ctx1->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_w;

    uint64_t lft0_msk = (mv_ctx0->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_h;
    uint64_t lft1_msk = (mv_ctx1->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_h;

    uint64_t bs1_map_h = (dbf_info->bs1_map.hor[y0_unit] >> (2 + x0_unit)) & unit_msk_w;
    uint64_t bs1_map_v = (dbf_info->bs1_map.ver[x0_unit] >> y0_unit)       & unit_msk_h;
    const int16_t *dist_ref_q = mv_ctx0 == &inter_ctx->mv_ctx0 ? inter_ctx->dist_ref_0
                                                               : inter_ctx->dist_ref_1;

    /* Avoid checking already set bs1 or bs2.
     * Note if no MV in map the other PU is intra so boundary strength is already 2
     * This way we only check condition on MVs when required
     */
    uint64_t chk_abv = (abv0_msk ^ abv1_msk) & (~bs1_map_h);
    uint64_t chk_lft = (lft0_msk ^ lft1_msk) & (~bs1_map_v);

    /* Init checked boundary strength map part to 0 it will be disabled if all MVs
     * are in the same ref * and the delta MVs are inferior to
     * integer MV precision.
     */
    uint64_t dst_map_h = (~chk_abv) & unit_msk_w;
    uint64_t dst_map_v = (~chk_lft) & unit_msk_h;

    if (chk_abv) {
        const OVMV *mv_abv0 = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv_abv1 = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv0     = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_abv);
            abv0_msk >>= nb_skipped_blk;
            uint8_t is_l0 = abv0_msk & 0x1;
            const int16_t *dist_ref_p = is_l0 ? ((mv_ctx0 == &inter_ctx->mv_ctx0) ? inter_ctx->dist_ref_0 : inter_ctx->dist_ref_1) : ((mv_ctx0 == &inter_ctx->mv_ctx1) ? inter_ctx->dist_ref_0 : inter_ctx->dist_ref_1);

            mv_abv0   += nb_skipped_blk;
            mv_abv1   += nb_skipped_blk;
            mv0       += nb_skipped_blk;

            pos_shift += nb_skipped_blk;

            uint64_t abv_th = check_dbf_enabled_p(dist_ref_p, dist_ref_q, is_l0 ? *mv_abv0 : *mv_abv1, *mv0);

            dst_map_h |= abv_th << pos_shift;

            mv_abv0++;
            mv_abv1++;
            mv0++;

            pos_shift++;

            abv0_msk >>= 1;
            chk_abv >>= nb_skipped_blk + 1;

        } while (chk_abv);
    }

    if (chk_lft) {
        const OVMV *mv_lft0 = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv_lft1 = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv0     = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];

        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_lft);

            lft0_msk >>= nb_skipped_blk;

            uint8_t is_l0 = lft0_msk & 0x1;
            const int16_t *dist_ref_p = is_l0 ? ((mv_ctx0 == &inter_ctx->mv_ctx0) ? inter_ctx->dist_ref_0 : inter_ctx->dist_ref_1) : ((mv_ctx0 == &inter_ctx->mv_ctx1) ? inter_ctx->dist_ref_0 : inter_ctx->dist_ref_1);

            mv_lft0   += 34 * nb_skipped_blk;
            mv_lft1   += 34 * nb_skipped_blk;
            mv0       += 34 * nb_skipped_blk;

            pos_shift += nb_skipped_blk;

            uint64_t lft_th = check_dbf_enabled_p(dist_ref_p, dist_ref_q, is_l0 ? *mv_lft0 : *mv_lft1, *mv0);

            dst_map_v |= lft_th << pos_shift;

            mv_lft0 += 34;
            mv_lft1 += 34;
            mv0     += 34;

            pos_shift++;

            lft0_msk >>= 1;
            chk_lft  >>= nb_skipped_blk + 1;

        } while (chk_lft);
    }
    dbf_info->bs1_map.hor[y0_unit] |= (bs1_map_h | dst_map_h) << (x0_unit + 2);
    dbf_info->bs1_map.ver[x0_unit] |= dst_map_v << y0_unit;
}

static void 
dbf_sbtmvp_set_hedges(const struct InterDRVCtx *const inter_ctx,
                      struct DBFInfo *const dbf_info,
                      int x0_unit, int y0_unit,
                      int nb_unit_w, int nb_unit_h)
{
    const struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    const struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t unit_msk_h = (uint64_t)((uint64_t)1 << nb_unit_w) - 1llu;

    uint64_t abv0_p_msk = (mv_ctx0->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_h;
    uint64_t abv1_p_msk = (mv_ctx1->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_h;

    uint64_t abv0_q_msk = (mv_ctx0->map.hfield[y0_unit + 1] >> (x0_unit + 1)) & unit_msk_h;
    uint64_t abv1_q_msk = (mv_ctx1->map.hfield[y0_unit + 1] >> (x0_unit + 1)) & unit_msk_h;

    uint64_t bs1_map_h = (dbf_info->bs1_map.hor[y0_unit] >> (x0_unit + 2)) & unit_msk_h;

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

    chk_b  &= (~bs1_map_h) & unit_msk_h;
    chk_p0 &= (~bs1_map_h) & unit_msk_h;
    chk_p1 &= (~bs1_map_h) & unit_msk_h;

    uint64_t dst_map_h = (~(chk_p0 | chk_p1 | chk_b)) & unit_msk_h;

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

    dbf_info->bs1_map.hor[y0_unit] |= dst_map_h << (x0_unit + 2);
}

static void 
dbf_sbtmvp_set_vedges(const struct InterDRVCtx *const inter_ctx,
                      struct DBFInfo *const dbf_info,
                      int x0_unit, int y0_unit,
                      int nb_unit_w, int nb_unit_h)
{
    const struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
    const struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;

    uint64_t unit_msk_v = (uint64_t)((uint64_t)1 << nb_unit_h) - 1llu;

    uint64_t lft0_p_msk = (mv_ctx0->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_v;
    uint64_t lft1_p_msk = (mv_ctx1->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_v;

    uint64_t lft0_q_msk = (mv_ctx0->map.vfield[x0_unit + 1] >> (y0_unit + 1)) & unit_msk_v;
    uint64_t lft1_q_msk = (mv_ctx1->map.vfield[x0_unit + 1] >> (y0_unit + 1)) & unit_msk_v;

    uint64_t bs1_map_v = (dbf_info->bs1_map.ver[x0_unit] >> y0_unit)      & unit_msk_v;

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

    chk_b  &= (~bs1_map_v) & unit_msk_v;
    chk_p0 &= (~bs1_map_v) & unit_msk_v;
    chk_p1 &= (~bs1_map_v) & unit_msk_v;

    uint64_t dst_map_v = (~(chk_p0 | chk_p1 | chk_b)) & unit_msk_v;

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

    dbf_info->bs1_map.ver[x0_unit] |= dst_map_v << y0_unit;
}

static void
dbf_mv_check_sbtmvp(const struct InterDRVCtx *const inter_ctx,
                    struct DBFInfo *const dbf_info,
                    int x0_unit, int y0_unit,
                    int nb_unit_w, int nb_unit_h)
{
    int i;
    for (i = 0; i < nb_unit_w; i += 2) {
        dbf_sbtmvp_set_vedges(inter_ctx, dbf_info, x0_unit + i, y0_unit, nb_unit_w, nb_unit_h);
    }

    for (i = 0; i < nb_unit_h; i += 2) {
        dbf_sbtmvp_set_hedges(inter_ctx, dbf_info, x0_unit, y0_unit + i, nb_unit_w, nb_unit_h);
    }
}

static void
dbf_mv_check_b(const struct InterDRVCtx *const inter_ctx,
               struct DBFInfo *const dbf_info,
               struct OVMVCtx *const mv_ctx0, struct OVMVCtx *const mv_ctx1,
               int x0_unit, int y0_unit,
               int nb_unit_w, int nb_unit_h)
{
    uint64_t unit_msk_w = (uint64_t)((uint64_t)1 << nb_unit_w) - 1llu;
    uint64_t unit_msk_h = (uint64_t)((uint64_t)1 << nb_unit_h) - 1llu;

    uint64_t abv0_msk = (mv_ctx0->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_w;
    uint64_t lft0_msk = (mv_ctx0->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_h;

    uint64_t abv1_msk = (mv_ctx1->map.hfield[y0_unit] >> (x0_unit + 1)) & unit_msk_w;
    uint64_t lft1_msk = (mv_ctx1->map.vfield[x0_unit] >> (y0_unit + 1)) & unit_msk_h;

    uint64_t bs1_map_h = (dbf_info->bs1_map.hor[y0_unit] >> (2 + x0_unit)) & unit_msk_w;
    uint64_t bs1_map_v = (dbf_info->bs1_map.ver[x0_unit] >> y0_unit)       & unit_msk_h;

    /* Note if no MV in map the other PU is intra so boundary strength is already 2
     * This way we only check condition on MVs when required
     */
    uint64_t chk_abv = (abv0_msk & abv1_msk) & (~bs1_map_h);
    uint64_t chk_lft = (lft0_msk & lft1_msk) & (~bs1_map_v);

    /* Init checked boundary strength map part to 0 it will be disabled if all MVs
     * are in the same ref * and the delta MVs are inferior to
     * integer MV precision.
     */
    uint64_t dst_map_h = (~chk_abv) & unit_msk_w;
    uint64_t dst_map_v = (~chk_lft) & unit_msk_h;

    if (chk_abv) {
        const OVMV *mv_abv0 = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv_abv1 = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit - 1)];
        const OVMV *mv0     = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        const OVMV *mv1     = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_abv);

            mv_abv0   += nb_skipped_blk;
            mv_abv1   += nb_skipped_blk;
            mv0       += nb_skipped_blk;
            mv1       += nb_skipped_blk;
            pos_shift += nb_skipped_blk;

            uint64_t abv_th = check_dbf_enabled(inter_ctx, *mv_abv0, *mv_abv1, *mv0, *mv1);

            dst_map_h |= abv_th << pos_shift;

            mv_abv0++;
            mv_abv1++;
            mv0++;
            mv1++;
            pos_shift++;

            chk_abv >>= nb_skipped_blk + 1;

        } while (chk_abv);
    }

    if (chk_lft) {
        const OVMV *mv_lft0 = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv_lft1 = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit - 1, y0_unit)];
        const OVMV *mv0     = &mv_ctx0->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];
        const OVMV *mv1     = &mv_ctx1->mvs[PB_POS_IN_BUF(x0_unit, y0_unit)];

        uint8_t pos_shift = 0;
        do {
            uint8_t nb_skipped_blk = ov_ctz64(chk_lft);
            mv_lft0   += 34 * nb_skipped_blk;
            mv_lft1   += 34 * nb_skipped_blk;
            mv0       += 34 * nb_skipped_blk;
            mv1       += 34 * nb_skipped_blk;
            pos_shift += nb_skipped_blk;

            uint64_t lft_th = check_dbf_enabled(inter_ctx, *mv_lft0, *mv_lft1, *mv0, *mv1);

            dst_map_v |= lft_th << pos_shift;

            mv_lft0 += 34;
            mv_lft1 += 34;
            mv0     += 34;
            mv1     += 34;
            pos_shift++;

            chk_lft >>= nb_skipped_blk + 1;

        } while (chk_lft);
    }
    dbf_info->bs1_map.hor[y0_unit] |= dst_map_h << (x0_unit + 2);
    dbf_info->bs1_map.ver[x0_unit] |= dst_map_v << y0_unit;
}

static inline OVMV
round_dmv(OVMV mv)
{
  const int rnd = (1 << (5 - 1)) - 1;
  OVMV tmp;

  #if 1
  tmp.x = (mv.x) + rnd;
  tmp.y = (mv.y) + rnd;

  tmp.x += mv.x >= 0;
  tmp.y += mv.y >= 0;
  #else
  tmp.x = mv.x << 2;
  tmp.y = mv.y << 2;
  #endif

  tmp.x >>= 5;
  tmp.y >>= 5;

  return tmp;
}
static inline uint8_t
check_dbf_mv(const OVMV mv0)
{
    /* FIXME check if need of crossed with vertical delta */
    #if 0
    OVMV mv = {.x= mv0.x*SB_SIZE, .y= mv0.y*SB_SIZE};
    mv = round_affine_mv2(mv);
    //mv = clip_mv(mv);
    #else
    OVMV mv = {.x= mv0.x, .y= mv0.y};
    mv = round_dmv(mv);
    #endif
    uint8_t int_diff  = abs(mv.x) >= LF_MV_THRESHOLD;
            int_diff |= abs(mv.y) >= LF_MV_THRESHOLD;
    return int_diff;
}

static inline uint8_t
check_dbf_mv_b(const OVMV mv_0, const OVMV mv_1)
{
    /* FIXME check if need of crossed with vertical delta */
    #if 0
    OVMV mv0 = {.x= mv_0.x*SB_SIZE + 31, .y= mv_0.y*SB_SIZE + 31};
    OVMV mv1 = {.x= mv_1.x*SB_SIZE + 31, .y= mv_1.y*SB_SIZE + 31};

    mv0 = round_affine_mv2(mv0);
    mv1 = round_affine_mv2(mv1);

    //mv0 = clip_mv(mv0);
    //mv1 = clip_mv(mv1);
    #else
    OVMV mv0 = {.x= mv_0.x, .y= mv_0.y};
    OVMV mv1 = {.x= mv_1.x, .y= mv_1.y};
    mv0 = round_dmv(mv0);
    mv1 = round_dmv(mv1);
    #endif

    uint8_t int_diff  = abs(mv0.x) >= LF_MV_THRESHOLD;
            int_diff |= abs(mv0.y) >= LF_MV_THRESHOLD;

            int_diff |= abs(mv1.x) >= LF_MV_THRESHOLD;
            int_diff |= abs(mv1.y) >= LF_MV_THRESHOLD;

    return int_diff;
}

static void
dbf_update_internal_p(const struct InterDRVCtx *inter_ctx,
                      struct DBFInfo *const dbf_info,
                      const struct AffineDeltaMV dmv0,
                      uint8_t x0_u, uint8_t  y0_u,
                      uint8_t nb_unit_w, uint8_t nb_unit_h,
                      uint8_t mv_broad, const struct OVMVCtx *mv_ctx,
                      const int16_t *dist_ref0)
{
    /* vertical edge */
    /* FIXME in case of Affine MVP there might be a shortcut using
     * Affine MV deltas however MV rounding makes it hard to come with a clean
     * solution so we use MVs instead
     */
    if (mv_broad) {
        int i;
        for (i = 2; i < nb_unit_w; i += 2) {
#if 0
            if (check_dbf_mv(dmv0.h)) {
                dbf_info->bs1_map.ver[x0_u + i] |= mask_ver;
            }
#else
            if (mv_broad) {
                const OVMV *mv = &mv_ctx->mvs[35 + i + x0_u + y0_u * 34];
                int j;
                uint64_t tmp_msk = 0;
                for (j = 0; j < nb_unit_h; j++) {
                    uint64_t val =  check_dbf_enabled_p(dist_ref0, dist_ref0, *mv, mv[-1]);
                    tmp_msk |= val << (j + y0_u);
                    mv += 34;
                }
                dbf_info->bs1_map.ver[x0_u + i] |= tmp_msk;
            }
#endif
        }

        for (i = 2; i < nb_unit_h; i += 2) {
#if 0
            if (check_dbf_mv(dmv0.v)) {
                dbf_info->bs1_map.hor[y0_u + i] |= mask_hor;
            }
#else
            const OVMV *mv = &mv_ctx->mvs[35 + x0_u + (y0_u + i) * 34];
            int j;
            uint64_t tmp_msk = 0;
            for (j = 0; j < nb_unit_w; j++) {
                uint64_t val =  check_dbf_enabled_p(dist_ref0, dist_ref0, *mv, mv[-34]);
                tmp_msk |= val << (j + 2 + x0_u);
                mv++;
            }
            dbf_info->bs1_map.hor[y0_u + i] |= tmp_msk;
#endif
        }
    }
}

static void
dbf_set_sb_edges(struct DBFInfo *const dbf_info,
                 uint8_t x0_u, uint8_t  y0_u,
                 uint8_t nb_unit_w, uint8_t nb_unit_h)
{
    int i;
    uint64_t msk_v = (uint64_t)((uint64_t)1 << nb_unit_h) - 1;
    uint64_t msk_h = (uint64_t)((uint64_t)1 << nb_unit_w) - 1;

    msk_v <<= y0_u;
    msk_h <<= (2 + x0_u);

    for (i = 2; i < nb_unit_w; i += 2) {
        dbf_info->aff_edg_ver[8 + x0_u + i] |= msk_v;
    }

    for (i = 2; i < nb_unit_h; i += 2) {
        dbf_info->aff_edg_hor[8 + y0_u + i] |= msk_h;
    }
}

static void
dbf_update_internal_b(const struct InterDRVCtx *inter_ctx,
                      struct DBFInfo *const dbf_info,
                      const struct AffineDeltaMV dmv0, const struct AffineDeltaMV dmv1,
                      uint8_t x0_u, uint8_t  y0_u,
                      uint8_t nb_unit_w, uint8_t nb_unit_h, uint8_t mv_broad,
                      const struct OVMVCtx *mv_ctx0, const struct OVMVCtx *mv_ctx1)
{

    if (mv_broad) {
        #if 0
        uint64_t mask_ver = (uint64_t)((uint64_t)1 << nb_unit_h) - 1;
        uint64_t mask_hor = (uint64_t)((uint64_t)1 << nb_unit_w) - 1;
        mask_ver <<= y0_u;
        mask_hor <<= (2 + x0_u);
        #endif
        int i;
        for (i = 2; i < nb_unit_w; i += 2) {
#if 0
            if (check_dbf_mv_b(dmv0.h, dmv1.h)) {
                dbf_info->bs1_map.ver[x0_u + i] |= mask_ver;
            }
#else
            const OVMV *mv0 = &mv_ctx0->mvs[35 + i + x0_u + y0_u * 34];
            const OVMV *mv1 = &mv_ctx1->mvs[35 + i + x0_u + y0_u * 34];
            int j;
            uint64_t tmp_msk = 0;
            for (j = 0; j < nb_unit_h; j++) {
#if 0
                uint64_t val  = mv_threshold_check(*mv0, mv0[-1]);
                val |= mv_threshold_check(*mv1, mv1[-1]);
                uint64_t val2  = mv_threshold_check(*mv0, mv1[-1]);
                val2 |= mv_threshold_check(*mv1, mv0[-1]);
#else
                uint64_t val =  check_dbf_enabled(inter_ctx, *mv0, *mv1, mv0[-1], mv1[-1]);
#endif
                tmp_msk |= val << (j + y0_u);
                mv0 += 34;
                mv1 += 34;
            }
            dbf_info->bs1_map.ver[x0_u + i] |= tmp_msk;
#endif
        }

        for (i = 2; i < nb_unit_h; i += 2) {
#if 0
            if (check_dbf_mv_b(dmv0.v, dmv1.v)) {
                dbf_info->bs1_map.hor[y0_u + i] |= mask_hor;
            }
#else
            const OVMV *mv0 = &mv_ctx0->mvs[35 + x0_u + (y0_u + i) * 34];
            const OVMV *mv1 = &mv_ctx1->mvs[35 + x0_u + (y0_u + i) * 34];
            int j;
            uint64_t tmp_msk = 0;
            for (j = 0; j < nb_unit_w; j++) {
                uint64_t val =  check_dbf_enabled(inter_ctx, *mv0, *mv1, mv0[-34], mv1[-34]);
                tmp_msk |= val << (j + 2 + x0_u);
                mv0++;
                mv1++;
            }
            dbf_info->bs1_map.hor[y0_u + i] |= tmp_msk;
#endif
        }
    }
}

static uint8_t
update_mv_ctx_b2(struct InterDRVCtx *const inter_ctx,
                uint8_t pb_x, uint8_t  pb_y,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t log2_cu_w, uint8_t log2_cu_h,
                uint8_t inter_dir)
{
    struct DBFInfo *const dbf_info = &inter_ctx->tmvp_ctx.ctudec->dbf_info;

    dbf_set_sb_edges(dbf_info, pb_x, pb_y, nb_pb_w, nb_pb_h);
    dbf_fill_aff_map(&dbf_info->affine_map, pb_x, pb_y, nb_pb_w, nb_pb_h);
    dbf_mv_check_sbtmvp(inter_ctx, dbf_info, pb_x, pb_y, nb_pb_w, nb_pb_h);

    #if 0
    if (inter_dir == 0x3) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        const struct AffineDeltaMV dmv_0 = {0}; 
        const struct AffineDeltaMV dmv_1 = {0}; 

        const uint8_t mv_broad_0 = 0;
        const uint8_t mv_broad_1 = 0;



        return ((!mv_broad_0) | (!mv_broad_1 << 1));
    } else if (inter_dir & 0x2) {
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;

        const struct AffineDeltaMV dmv_1 = {0}; 

        const uint8_t mv_broad_1 = 0;

        dbf_mv_check_p(inter_ctx, dbf_info, mv_ctx1, mv_ctx0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        dbf_update_internal_p(inter_ctx, dbf_info, dmv_1, pb_x, pb_y, nb_pb_w, nb_pb_h, !mv_broad_1, mv_ctx1, inter_ctx->dist_ref_1);

        return (!mv_broad_1 << 1);
    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        const struct AffineDeltaMV dmv_0 = {0}; 

        const uint8_t mv_broad_0 = 0;

        dbf_mv_check_p(inter_ctx, dbf_info, mv_ctx0, mv_ctx1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        dbf_update_internal_p(inter_ctx, dbf_info, dmv_0, pb_x, pb_y, nb_pb_w, nb_pb_h, !mv_broad_0, mv_ctx0, inter_ctx->dist_ref_0);

        return (!mv_broad_0);
    }
    #endif
    return 0;
}

/* FIXME avoid duplicated dmv derivation */
static uint8_t
update_mv_ctx_b(struct InterDRVCtx *const inter_ctx,
                struct AffineMergeInfo mv_info,
                uint8_t pb_x, uint8_t  pb_y,
                uint8_t nb_pb_w, uint8_t nb_pb_h,
                uint8_t log2_cu_w, uint8_t log2_cu_h,
                uint8_t inter_dir)
{
    /*FIXME check for specific DBF application for sub blokc MVs */
    struct AffineDRVInfo *aff_info = &inter_ctx->affine_ctx;
    const struct AffineControlInfo *const cinfo = mv_info.cinfo;
    struct DBFInfo *const dbf_info = &inter_ctx->tmvp_ctx.ctudec->dbf_info;
    uint8_t affine_type = mv_info.affine_type;

    uint16_t pos = PB_POS_IN_BUF(pb_x, pb_y);

    /* FIXME check required */
    ctu_field_set_rect_bitfield(&aff_info->map, pb_x, pb_y, nb_pb_w, nb_pb_h);

    dbf_set_sb_edges(dbf_info, pb_x, pb_y, nb_pb_w, nb_pb_h);
    dbf_fill_aff_map(&dbf_info->affine_map, pb_x, pb_y, nb_pb_w, nb_pb_h);

    if (inter_dir == 0x3) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        const struct AffineDeltaMV dmv_0 = derive_affine_delta_mvs(&cinfo[0],
                                                                   log2_cu_w, log2_cu_h,
                                                                   affine_type);

        const struct AffineDeltaMV dmv_1 = derive_affine_delta_mvs(&cinfo[1],
                                                                   log2_cu_w, log2_cu_h,
                                                                   affine_type);

        const uint8_t mv_broad_0 = broadcast_mv(dmv_0, 0x3);
        const uint8_t mv_broad_1 = broadcast_mv(dmv_1, 0x3);

        ctu_field_set_rect_bitfield(&mv_ctx0->map, pb_x, pb_y, nb_pb_w, nb_pb_h);
        ctu_field_set_rect_bitfield(&mv_ctx1->map, pb_x, pb_y, nb_pb_w, nb_pb_h);

        compute_subblock_mvs(&cinfo[0], dmv_0, &mv_ctx0->mvs[pos],
                             log2_cu_w, log2_cu_h, mv_broad_0);
        compute_subblock_mvs(&cinfo[1], dmv_1, &mv_ctx1->mvs[pos],
                             log2_cu_w, log2_cu_h, mv_broad_1);

        dbf_update_internal_b(inter_ctx, dbf_info, dmv_0, dmv_1, pb_x, pb_y, nb_pb_w, nb_pb_h,
                              !(mv_broad_0 && mv_broad_1), mv_ctx0, mv_ctx1);

        dbf_mv_check_b(inter_ctx, dbf_info, mv_ctx0, mv_ctx1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        return ((!mv_broad_0) | (!mv_broad_1 << 1));
    } else if (inter_dir & 0x2) {
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        const struct AffineDeltaMV dmv_1 = derive_affine_delta_mvs(&cinfo[1],
                                                                   log2_cu_w, log2_cu_h,
                                                                   affine_type);

        const uint8_t mv_broad_1 = broadcast_mv(dmv_1, 0x2);

        ctu_field_set_rect_bitfield(&mv_ctx1->map, pb_x, pb_y, nb_pb_w, nb_pb_h);

        compute_subblock_mvs(&cinfo[1], dmv_1, &mv_ctx1->mvs[pos],
                             log2_cu_w, log2_cu_h, mv_broad_1);

        dbf_mv_check_p(inter_ctx, dbf_info, mv_ctx1, mv_ctx0, pb_x, pb_y, nb_pb_w, nb_pb_h);

        dbf_update_internal_p(inter_ctx, dbf_info, dmv_1, pb_x, pb_y, nb_pb_w, nb_pb_h,
                              !mv_broad_1, mv_ctx1, inter_ctx->dist_ref_1);

        return (!mv_broad_1 << 1);
    } else if (inter_dir & 0x1) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        const struct AffineDeltaMV dmv_0 = derive_affine_delta_mvs(&cinfo[0],
                                                                   log2_cu_w, log2_cu_h,
                                                                   affine_type);

        const uint8_t mv_broad_0 = broadcast_mv(dmv_0, 0x1);

        ctu_field_set_rect_bitfield(&mv_ctx0->map, pb_x, pb_y, nb_pb_w, nb_pb_h);

        compute_subblock_mvs(&cinfo[0], dmv_0, &mv_ctx0->mvs[pos],
                             log2_cu_w, log2_cu_h, mv_broad_0);

        dbf_mv_check_p(inter_ctx, dbf_info, mv_ctx0, mv_ctx1, pb_x, pb_y, nb_pb_w, nb_pb_h);

        dbf_update_internal_p(inter_ctx, dbf_info, dmv_0, pb_x, pb_y, nb_pb_w, nb_pb_h,
                              !mv_broad_0, mv_ctx0, inter_ctx->dist_ref_0);

        return (!mv_broad_0);
    }
    return 0;
}

void
store_affine_info(struct AffineDRVInfo *const affine_ctx, struct AffineInfo aff_info, uint8_t x_pb, uint8_t y_pb, uint8_t nb_pb_w, uint8_t nb_pb_h)
{
    int i, j;
    uint16_t pos = PB_POS_IN_BUF(x_pb, y_pb);
    struct AffineInfo *aff_buff = &affine_ctx->affine_info[pos];

    for (i = 0; i < nb_pb_h; ++i) {
        for (j = 0; j < nb_pb_w; ++j) {
            aff_buff[j] = aff_info;
        }
        aff_buff += 34;
    }
}
#define SB_W 4
#define SB_H 4

#define GRAD_SHIFT 6

/* Link to GRAD SHIFT */
#define PROF_DMV_MAX ((1 << 5) - 1)

#define PROF_MV_SHIFT 8
#define PROF_MV_RND (1 << (PROF_MV_SHIFT - 1))

#define BITDEPTH 10
#define PROF_SMP_SHIFT (14 - BITDEPTH)
#define PROF_SMP_RND (1 << (14 - 1))
#define PROF_SMP_OFFSET (1 << (PROF_SMP_SHIFT - 1)) + PROF_SMP_RND

#define PROF_BUFF_PADD_H 1
#define PROF_BUFF_PADD_W 1

struct OVDMV {
    int32_t x;
    int32_t y;
};

struct OVDMV
round_prof_dmv_scale(struct OVDMV dmv)
{
    dmv.x += PROF_MV_RND - (dmv.x >= 0);
    dmv.y += PROF_MV_RND - (dmv.y >= 0);
    dmv.x >>= PROF_MV_SHIFT;
    dmv.y >>= PROF_MV_SHIFT;
    return dmv;
}

static void
compute_prof_dmv_scale(struct AffineDeltaMV delta_mv,
                       int32_t dmv_scale_h[16], int32_t dmv_scale_v[16])
{
    int x, y, i;
    OVMV quad_dmv_h;
    OVMV quad_dmv_v;
    int32_t *dmv_scale_h_p = dmv_scale_h;
    int32_t *dmv_scale_v_p = dmv_scale_v;

    quad_dmv_h.x = delta_mv.h.x << 2;
    quad_dmv_h.y = delta_mv.h.y << 2;
    quad_dmv_v.x = delta_mv.v.x << 2;
    quad_dmv_v.y = delta_mv.v.y << 2;

    dmv_scale_h_p[0] = ((delta_mv.h.x + delta_mv.v.x) << 1) - ((quad_dmv_h.x + quad_dmv_v.x) << 1);
    dmv_scale_v_p[0] = ((delta_mv.h.y + delta_mv.v.y) << 1) - ((quad_dmv_h.y + quad_dmv_v.y) << 1);

    for (x = 1; x < SB_W; x++) {
        dmv_scale_h_p[x] = dmv_scale_h_p[x - 1] + quad_dmv_h.x;
        dmv_scale_v_p[x] = dmv_scale_v_p[x - 1] + quad_dmv_h.y;
    }

    dmv_scale_h_p += SB_W;
    dmv_scale_v_p += SB_W;

    for (y = 1; y < SB_H; y++) {
        for (x = 0; x < SB_W; x++) {
            dmv_scale_h_p[x] = dmv_scale_h_p[x - SB_W] + quad_dmv_v.x;
            dmv_scale_v_p[x] = dmv_scale_v_p[x - SB_W] + quad_dmv_v.y;
        }
        dmv_scale_h_p += SB_W;
        dmv_scale_v_p += SB_W;
    }

    for (i = 0; i < SB_W * SB_H; i++) {

        struct OVDMV dmv = {
            .x = dmv_scale_h[i],
            .y = dmv_scale_v[i]
        };

        dmv = round_prof_dmv_scale(dmv);

        dmv_scale_h[i] = ov_clip(dmv.x, -PROF_DMV_MAX, PROF_DMV_MAX);
        dmv_scale_v[i] = ov_clip(dmv.y, -PROF_DMV_MAX, PROF_DMV_MAX);
    }
}

static void
rcn_affine_mcp_b_l(OVCTUDec *const ctudec,
                   struct InterDRVCtx *const inter_ctx,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_cu_w, uint8_t log2_cu_h,
                   uint8_t inter_dir)
{
    int i, j;
    uint8_t nb_sb_w = (1 << log2_cu_w) >> LOG2_MIN_CU_S;
    uint8_t nb_sb_h = (1 << log2_cu_h) >> LOG2_MIN_CU_S;

    uint16_t pos = PB_POS_IN_BUF(x0 >> 2, y0 >> 2);

    const struct OVMV *mv_buff0 = &inter_ctx->mv_ctx0.mvs[pos];
    const struct OVMV *mv_buff1 = &inter_ctx->mv_ctx1.mvs[pos];
    uint8_t ref_idx0 = mv_buff0->ref_idx;
    uint8_t ref_idx1 = mv_buff1->ref_idx;


    OVMV *tmvp_mv0 = &inter_ctx->tmvp_mv[0].mvs[0];
    OVMV *tmvp_mv1 = &inter_ctx->tmvp_mv[1].mvs[0];

    for (i = 0; i < nb_sb_h; ++i) {
        for (j = 0; j < nb_sb_w; ++j) {
            OVMV mv0 = mv_buff0[j];
            OVMV mv1 = mv_buff1[j];

            #if 0
            if (!((j + start_x) & 0x1) && !((i + start_y) & 0x1)) {
            #else

               if (!(((x0+4*j)>>2)& 0x1) && !((((y0+4*i)>>2)& 0x1))) {
            #endif
               tmvp_mv0[((x0+4*j)>>3) + ((y0+4*i)>>3) *16] = mv0;
               tmvp_mv1[((x0+4*j)>>3) + ((y0+4*i)>>3) *16] = mv1;
            }

            rcn_mcp_b_l(ctudec, ctudec->rcn_ctx.ctu_buff, inter_ctx, ctudec->part_ctx,
                        mv0, mv1, x0 + 4*j, y0 + 4*i,
                        2, 2, inter_dir, ref_idx0, ref_idx1);
        }

        #if 0
        if (((i + start_y) & 0x1)) {
            tmvp_mv0 += 16;
            tmvp_mv1 += 16;
        }
        #endif

        mv_buff0 += 34;
        mv_buff1 += 34;
    }
}

struct PROFInfo {
    int32_t dmv_scale_h_0[16];
    int32_t dmv_scale_v_0[16];
    int32_t dmv_scale_h_1[16];
    int32_t dmv_scale_v_1[16];
};

static void
rcn_affine_prof_mcp_b_l(OVCTUDec *const ctudec,
                        struct InterDRVCtx *const inter_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_cu_w, uint8_t log2_cu_h,
                        uint8_t inter_dir, uint8_t prof_dir,
                        const struct AffineDeltaMV *const dmv_0,
                        const struct AffineDeltaMV *const dmv_1)
{
    int i, j;
    uint8_t nb_sb_w = (1 << log2_cu_w) >> LOG2_MIN_CU_S;
    uint8_t nb_sb_h = (1 << log2_cu_h) >> LOG2_MIN_CU_S;

    uint16_t pos = PB_POS_IN_BUF(x0 >> 2, y0 >> 2);

    const struct OVMV *mv_buff0 = &inter_ctx->mv_ctx0.mvs[pos];
    const struct OVMV *mv_buff1 = &inter_ctx->mv_ctx1.mvs[pos];
    uint8_t ref_idx0 = mv_buff0->ref_idx;
    uint8_t ref_idx1 = mv_buff1->ref_idx;
    struct PROFInfo prof_info;

    if (prof_dir & 0x1) {
        compute_prof_dmv_scale(*dmv_0, prof_info.dmv_scale_h_0, prof_info.dmv_scale_v_0);
    }

    if (prof_dir & 0x2) {
        compute_prof_dmv_scale(*dmv_1, prof_info.dmv_scale_h_1, prof_info.dmv_scale_v_1);
    }

    OVMV *tmvp_mv0 = &inter_ctx->tmvp_mv[0].mvs[0];
    OVMV *tmvp_mv1 = &inter_ctx->tmvp_mv[1].mvs[0];

    // inter_ctx->prec_amvr = inter_dir & 0x1 ? mv_buff0[0].prec_amvr : mv1.prec_amvr;
    for (i = 0; i < nb_sb_h; ++i) {
        for (j = 0; j < nb_sb_w; ++j) {
            OVMV mv0 = mv_buff0[j];
            OVMV mv1 = mv_buff1[j];

            if (!(((x0+4*j)>>2)& 0x1) && !((((y0+4*i)>>2)& 0x1))) {
               tmvp_mv0[((x0+4*j)>>3) + ((y0+4*i)>>3) *16] = mv0;
               tmvp_mv1[((x0+4*j)>>3) + ((y0+4*i)>>3) *16] = mv1;
            }

            //TODOrebase: take into account amvr and bcw
            rcn_prof_mcp_b_l(ctudec, ctudec->rcn_ctx.ctu_buff, inter_ctx, ctudec->part_ctx,
                             mv0, mv1, x0 + 4*j, y0 + 4*i,
                             2, 2, inter_dir, ref_idx0, ref_idx1,
                             prof_dir, &prof_info);
            // rcn_mcp_b_l(ctudec, ctudec->rcn_ctx.ctu_buff, inter_ctx, ctudec->part_ctx,
            //             mv0, mv1, x0 + 4*j, y0 + 4*i,
            //             2, 2, inter_dir, ref_idx0, ref_idx1);
        }
        #if 0
        if (((i + start_y) & 0x1)) {
            tmvp_mv0 += 16;
            tmvp_mv1 += 16;
        }
        #endif

        mv_buff0 += 34;
        mv_buff1 += 34;
    }
}


void
rcn_affine_mcp_b_c(OVCTUDec *const ctudec,
                   struct InterDRVCtx *const inter_ctx,
                   uint8_t x0, uint8_t y0,
                   uint8_t log2_cu_w, uint8_t log2_cu_h,
                   uint8_t inter_dir)
{
    uint8_t nb_sb_w = (1 << log2_cu_w) >> LOG2_MIN_CU_S;
    uint8_t nb_sb_h = (1 << log2_cu_h) >> LOG2_MIN_CU_S;

    uint16_t pos = PB_POS_IN_BUF(x0 >> 2, y0 >> 2);

    const struct OVMV *mv_buff0 = &inter_ctx->mv_ctx0.mvs[pos];
    const struct OVMV *mv_buff1 = &inter_ctx->mv_ctx1.mvs[pos];
    uint8_t ref_idx0 = mv_buff0->ref_idx;
    uint8_t ref_idx1 = mv_buff1->ref_idx;

    int i, j;

    for (i = 0; i < nb_sb_h; i += 2) {
        for (j = 0; j < nb_sb_w; j += 2) {
            OVMV mv0 = mv_buff0[j];
            OVMV mv1 = mv_buff1[j];

            mv0.x += mv_buff0[j + 35].x;
            mv0.y += mv_buff0[j + 35].y;

            mv0.x += mv0.x < 0;
            mv0.y += mv0.y < 0;
            mv0.x >>= 1;
            mv0.y >>= 1;

            mv1.x += mv_buff1[j + 35].x;
            mv1.y += mv_buff1[j + 35].y;

            mv1.x += mv1.x < 0;
            mv1.y += mv1.y < 0;
            mv1.x >>= 1;
            mv1.y >>= 1;
            
            rcn_mcp_b_c(ctudec, ctudec->rcn_ctx.ctu_buff, inter_ctx, ctudec->part_ctx,
                        mv0, mv1, x0+ 4*j, y0 + 4*i,
                        3, 3, inter_dir, ref_idx0, ref_idx1);
        }

        mv_buff0 += 34 * 2;
        mv_buff1 += 34 * 2;
    }
}


/* TODO P slices functions */
void
drv_affine_merge_mvp()
{
}

void
drv_affine_mvp_mvd()
{
}

static inline uint8_t
mv_cmp(const OVMV a, const OVMV b)
{
     uint8_t is_eq;
     is_eq = a.x == b.x;
     is_eq &= a.y == b.y;
     return is_eq;
}

uint8_t check_affine_prof(const struct AffineMergeInfo *const affine_info, uint8_t rpl_idx)
{
    uint8_t prof_enabled = 1;
    const struct AffineControlInfo *const cpinfo = &affine_info->cinfo[rpl_idx];

    if (affine_info->affine_type == AFFINE_3CP) {
        prof_enabled &= !(mv_cmp(cpinfo->lt, cpinfo->rt) && mv_cmp(cpinfo->lt, cpinfo->lb));
    } else {
        prof_enabled &= !mv_cmp(cpinfo->lt, cpinfo->rt);
    }

    return prof_enabled;
}

/* Derive motion vectors and update motion maps */
void
drv_affine_mvp_b(struct InterDRVCtx *const inter_ctx,
                 uint8_t x0, uint8_t y0,
                 uint8_t log2_cu_w, uint8_t log2_cu_h,
                 struct AffineControlInfo * cp_mvd0,
                 struct AffineControlInfo * cp_mvd1,
                 uint8_t mvp_idx0, uint8_t mvp_idx1, uint8_t bcw_idx,
                 uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                 uint8_t affine_type)
{
    struct AffineDRVInfo *affine_ctx = &inter_ctx->affine_ctx;
    struct AffineMergeInfo mv_info ={0};
    uint8_t prec_amvr = inter_ctx->prec_amvr;
    
    uint8_t x_pb = x0 >> 2;
    uint8_t y_pb = y0 >> 2;

    uint8_t nb_pb_w = (1 << log2_cu_w) >> 2;
    uint8_t nb_pb_h = (1 << log2_cu_h) >> 2;

    uint8_t opp_ref_idx0 = 0xFF;
    uint8_t opp_ref_idx1 = 0xFF;
    uint8_t prof_dir = inter_ctx->prof_enabled ? 0x3 : 0;

    for (int i = 0; i < inter_ctx->nb_active_ref1; i ++) {
         if (inter_ctx->rpl0[ref_idx0] == inter_ctx->rpl1[i])
             opp_ref_idx0 = i;
    }

    for (int i = 0; i < inter_ctx->nb_active_ref0; i ++) {
         if (inter_ctx->rpl1[ref_idx1] == inter_ctx->rpl0[i])
             opp_ref_idx1 = i;
    }

    /* FIXME can we combine mvp derivation for bi pred */
    if (inter_dir & 0x1) {
        struct AffineControlInfo *const cp_info = &mv_info.cinfo[0];

        *cp_info = drv_affine_mvp(inter_ctx, affine_ctx, x_pb, y_pb,
                                  nb_pb_w, nb_pb_h, log2_cu_w, log2_cu_h,
                                  ref_idx0, opp_ref_idx0, mvp_idx0,
                                  inter_dir & 0x1, affine_type);

        // cp_info->lt  = drv_round_to_precision_mv(cp_info->lt, MV_PRECISION_INTERNAL, prec_amvr);
        // cp_info->rt  = drv_round_to_precision_mv(cp_info->lt, MV_PRECISION_INTERNAL, prec_amvr);
        cp_mvd0->lt = drv_change_precision_mv(cp_mvd0->lt, prec_amvr, MV_PRECISION_INTERNAL);
        cp_mvd0->rt = drv_change_precision_mv(cp_mvd0->rt, prec_amvr, MV_PRECISION_INTERNAL);

        cp_info->lt.x +=  cp_mvd0->lt.x; // << 2;
        cp_info->lt.y +=  cp_mvd0->lt.y; // << 2;

        cp_info->rt.x +=  cp_mvd0->lt.x; // << 2;
        cp_info->rt.y +=  cp_mvd0->lt.y; // << 2;
        cp_info->rt.x +=  cp_mvd0->rt.x; // << 2;
        cp_info->rt.y +=  cp_mvd0->rt.y; // << 2;

        cp_info->lt = mv_clip_periodic(cp_info->lt);
        cp_info->rt = mv_clip_periodic(cp_info->rt);
        cp_info->lt.ref_idx = ref_idx0;
        cp_info->rt.ref_idx = ref_idx0;
        cp_info->lt.bcw_idx_plus1 = bcw_idx + 1;
        cp_info->rt.bcw_idx_plus1 = bcw_idx + 1;
        cp_info->lt.prec_amvr = prec_amvr ;
        cp_info->rt.prec_amvr = prec_amvr ;
        if (affine_type == AFFINE_3CP) {
            // cp_info->lb  = drv_round_to_precision_mv(cp_info->lb, MV_PRECISION_INTERNAL, prec_amvr);
            cp_mvd0->lb = drv_change_precision_mv(cp_mvd0->lb, prec_amvr, MV_PRECISION_INTERNAL);
            cp_info->lb.x +=  cp_mvd0->lt.x; // << 2;
            cp_info->lb.y +=  cp_mvd0->lt.y; // << 2;
            cp_info->lb.x +=  cp_mvd0->lb.x; // << 2;
            cp_info->lb.y +=  cp_mvd0->lb.y; // << 2;
            cp_info->lb = mv_clip_periodic(cp_info->lb);
            cp_info->lb.ref_idx = ref_idx0;
            cp_info->lb.bcw_idx_plus1 = bcw_idx + 1;
            cp_info->lb.prec_amvr = prec_amvr;
        }
    }

    if (inter_dir & 0x2) {
        struct AffineControlInfo *const cp_info = &mv_info.cinfo[1];

        *cp_info = drv_affine_mvp(inter_ctx, affine_ctx, x_pb, y_pb,
                                  nb_pb_w, nb_pb_h, log2_cu_w, log2_cu_h,
                                  ref_idx1, opp_ref_idx1, mvp_idx1,
                                  inter_dir & 0x2, affine_type);

        // cp_info->lt  = drv_round_to_precision_mv(cp_info->lt, MV_PRECISION_INTERNAL, prec_amvr);
        // cp_info->rt  = drv_round_to_precision_mv(cp_info->lt, MV_PRECISION_INTERNAL, prec_amvr);
        cp_mvd1->lt = drv_change_precision_mv(cp_mvd1->lt, prec_amvr, MV_PRECISION_INTERNAL);
        cp_mvd1->rt = drv_change_precision_mv(cp_mvd1->rt, prec_amvr, MV_PRECISION_INTERNAL);

        cp_info->lt.x +=  cp_mvd1->lt.x; // << 2;
        cp_info->lt.y +=  cp_mvd1->lt.y; // << 2;

        cp_info->rt.x +=  cp_mvd1->lt.x; // << 2;
        cp_info->rt.y +=  cp_mvd1->lt.y; // << 2;
        cp_info->rt.x +=  cp_mvd1->rt.x; // << 2;
        cp_info->rt.y +=  cp_mvd1->rt.y; // << 2;

        cp_info->lt = mv_clip_periodic(cp_info->lt);
        cp_info->rt = mv_clip_periodic(cp_info->rt);
        cp_info->lt.ref_idx = ref_idx1;
        cp_info->rt.ref_idx = ref_idx1;
        cp_info->lt.bcw_idx_plus1 = bcw_idx + 1;
        cp_info->rt.bcw_idx_plus1 = bcw_idx + 1;
        cp_info->lt.prec_amvr = prec_amvr ;
        cp_info->rt.prec_amvr = prec_amvr ;

        if (affine_type == AFFINE_3CP) {
            // cp_info->lb  = drv_round_to_precision_mv(cp_info->lb, MV_PRECISION_INTERNAL, prec_amvr);
            cp_mvd1->lb = drv_change_precision_mv(cp_mvd1->lb, prec_amvr, MV_PRECISION_INTERNAL);

            cp_info->lb.x +=  cp_mvd1->lt.x; // << 2;
            cp_info->lb.y +=  cp_mvd1->lt.y; // << 2;
            cp_info->lb.x +=  cp_mvd1->lb.x; // << 2;
            cp_info->lb.y +=  cp_mvd1->lb.y; // << 2;
            cp_info->lb = mv_clip_periodic(cp_info->lb);
            cp_info->lb.ref_idx = ref_idx1;
            cp_info->lb.bcw_idx_plus1 = bcw_idx + 1;
            cp_info->lb.prec_amvr = prec_amvr;
        }
    }

    mv_info.inter_dir = inter_dir;
    mv_info.affine_type = affine_type;
    inter_ctx->prec_amvr = 0;

    const struct AffineControlInfo *const cinfo = mv_info.cinfo;
    const struct AffineDeltaMV dmv_0 = derive_affine_delta_mvs(&cinfo[0],
                                                               log2_cu_w, log2_cu_h,
                                                               affine_type);

    const struct AffineDeltaMV dmv_1 = derive_affine_delta_mvs(&cinfo[1],
                                                               log2_cu_w, log2_cu_h,
                                                               affine_type);

    /* Update for next pass */
    prof_dir &= update_mv_ctx_b(inter_ctx, mv_info, x_pb, y_pb, nb_pb_w,
                                nb_pb_h, log2_cu_w, log2_cu_h, inter_dir);

    if (prof_dir) {
        uint8_t prof_0 = check_affine_prof(&mv_info, RPL_0);
        uint8_t prof_1 = check_affine_prof(&mv_info, RPL_1);

        prof_dir &= (prof_0) | (prof_1 << 1);

        prof_dir &= inter_dir;
    }

    if (!prof_dir) {
        rcn_affine_mcp_b_l(inter_ctx->tmvp_ctx.ctudec, inter_ctx, x0, y0,
                           log2_cu_w, log2_cu_h,
                           inter_dir);
    } else {
        rcn_affine_prof_mcp_b_l(inter_ctx->tmvp_ctx.ctudec, inter_ctx, x0, y0,
                                log2_cu_w, log2_cu_h,
                                inter_dir, prof_dir, &dmv_0, &dmv_1);
    }

    rcn_affine_mcp_b_c(inter_ctx->tmvp_ctx.ctudec, inter_ctx, x0, y0,
                       log2_cu_w, log2_cu_h,
                       inter_dir);

    struct PBInfo pb = {
        .x_pb = x_pb,
        .y_pb = y_pb,
        .nb_pb_w = nb_pb_w,
        .nb_pb_h = nb_pb_h,
        .log2_w = log2_cu_w,
        .log2_h = log2_cu_h
    };

    struct AffineInfo aff_info = {
        .cps[0] = mv_info.cinfo[0],
        .cps[1] = mv_info.cinfo[1],
        .pb = pb,
        .type = affine_type
    };

    store_affine_info(affine_ctx, aff_info, x_pb, y_pb, nb_pb_w, nb_pb_h);
}

void
drv_affine_merge_mvp_b(struct InterDRVCtx *const inter_ctx,
                       uint8_t x0, uint8_t y0,
                       uint8_t log2_cu_w, uint8_t log2_cu_h,
                       uint8_t merge_idx)
{
    struct AffineDRVInfo *affine_ctx = &inter_ctx->affine_ctx;
    struct AffineMergeInfo mv_info;

    uint8_t x_pb = x0 >> 2;
    uint8_t y_pb = y0 >> 2;

    uint8_t nb_pb_w = (1 << log2_cu_w) >> 2;
    uint8_t nb_pb_h = (1 << log2_cu_h) >> 2;
    uint8_t sbtmvp_enabled = inter_ctx->sbtmvp_enabled;
    uint8_t is_sbtmvp = 0;

    if (sbtmvp_enabled) {
        struct OVMVCtx *const mv_ctx0 = &inter_ctx->mv_ctx0;
        struct OVMVCtx *const mv_ctx1 = &inter_ctx->mv_ctx1;
        struct MergeInfo mv_info;
        struct AffineMergeInfo mv_info2;
        OVMV mv_offset = {0};

        /* FIXME avoid duplication with merge derivation */
        uint64_t rpl0_lft_col = mv_ctx0->map.vfield[x_pb];
        uint64_t rpl0_abv_row = mv_ctx0->map.hfield[y_pb];
        uint64_t rpl1_lft_col = mv_ctx1->map.vfield[x_pb];
        uint64_t rpl1_abv_row = mv_ctx1->map.hfield[y_pb];

        const uint8_t rpl0_cand = check_cand_available(rpl0_abv_row, rpl0_lft_col, x_pb, y_pb,
                                                       nb_pb_w, nb_pb_h);

        const uint8_t rpl1_cand = check_cand_available(rpl1_abv_row, rpl1_lft_col, x_pb, y_pb,
                                                       nb_pb_w, nb_pb_h);

        uint8_t sb_cand = derive_sub_pu_merge_cand(inter_ctx, x0, y0,
                                                   log2_cu_w, log2_cu_h,
                                                   &mv_info, &mv_offset,
                                                   rpl0_cand, rpl1_cand);

        if (sb_cand && merge_idx == 0) {
            derive_sub_block_mvs(inter_ctx, &inter_ctx->tmvp_ctx,
                                 x0, y0,
                                 log2_cu_w, log2_cu_h,
                                 mv_offset, &mv_info);

            #if 1
            update_mv_ctx_b2(inter_ctx, x_pb, y_pb,
                             nb_pb_w, nb_pb_h, log2_cu_w, log2_cu_h,
                             0x3);
                             #endif
            is_sbtmvp = 1;
        }

        merge_idx -= sb_cand;
    }

    if (!is_sbtmvp) {
        derive_affine_merge_mv(inter_ctx, affine_ctx, &mv_info,
                               x_pb, y_pb, nb_pb_w, nb_pb_h,
                               log2_cu_w, log2_cu_h,
                               merge_idx);
    }

    mv_info.cinfo[0].lt.prec_amvr = 0;
    mv_info.cinfo[0].rt.prec_amvr = 0;
    mv_info.cinfo[0].lb.prec_amvr = 0;
    mv_info.cinfo[1].lt.prec_amvr = 0;
    mv_info.cinfo[1].rt.prec_amvr = 0;
    mv_info.cinfo[1].lb.prec_amvr = 0;
    /* FIXME can we have small blocks bidir requiring inter_dir
     * override
     */
    if (!is_sbtmvp) {
        uint8_t prof_dir = inter_ctx->prof_enabled ? 0x3 : 0;

        uint8_t affine_type = mv_info.affine_type;
        const struct AffineControlInfo *const cinfo = mv_info.cinfo;

        const struct AffineDeltaMV dmv_0 = derive_affine_delta_mvs(&cinfo[0],
                                                                   log2_cu_w, log2_cu_h,
                                                                   affine_type);

        const struct AffineDeltaMV dmv_1 = derive_affine_delta_mvs(&cinfo[1],
                                                                   log2_cu_w, log2_cu_h,
                                                                   affine_type);
        prof_dir &= update_mv_ctx_b(inter_ctx, mv_info, x_pb, y_pb,
                                    nb_pb_w, nb_pb_h, log2_cu_w, log2_cu_h,
                                    mv_info.inter_dir);

        if (prof_dir) {
            uint8_t prof_0 = check_affine_prof(&mv_info, RPL_0);
            uint8_t prof_1 = check_affine_prof(&mv_info, RPL_1);

            prof_dir &= (prof_0) | (prof_1 << 1);

            prof_dir &= mv_info.inter_dir;
        }

        if (!prof_dir) {
            rcn_affine_mcp_b_l(inter_ctx->tmvp_ctx.ctudec, inter_ctx, x0, y0,
                               log2_cu_w, log2_cu_h,
                               mv_info.inter_dir);
        } else {
            rcn_affine_prof_mcp_b_l(inter_ctx->tmvp_ctx.ctudec, inter_ctx, x0, y0,
                                    log2_cu_w, log2_cu_h,
                                    mv_info.inter_dir, prof_dir, &dmv_0, &dmv_1);
        }

        rcn_affine_mcp_b_c(inter_ctx->tmvp_ctx.ctudec, inter_ctx, x0, y0,
                           log2_cu_w, log2_cu_h,
                           mv_info.inter_dir);

    }

    struct PBInfo pb = {
        .x_pb = x_pb,
        .y_pb = y_pb,
        .nb_pb_w = nb_pb_w,
        .nb_pb_h = nb_pb_h,
        .log2_w = log2_cu_w,
        .log2_h = log2_cu_h
    };

    struct AffineInfo aff_info = {
        .cps[0] = mv_info.cinfo[0],
        .cps[1] = mv_info.cinfo[1],
        .pb = pb,
        .type = mv_info.affine_type
    };

    store_affine_info(affine_ctx, aff_info, x_pb, y_pb, nb_pb_w, nb_pb_h);

}
