#include "ovutils.h"
#include "ctudec.h"

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

#define MAX_NB_AMVP_CAND 2

#define MV_BITS 18
#define CLIP_PERIOD  (1 << MV_BITS)
#define HALF_CLIP_PERIOD  (1 << (MV_BITS - 1))

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

static inline OVMV
round_affine_mv(const OVMV mv)
{
    OVMV tmp;

    tmp.x = (mv.x + 1 + (mv.x < 0)) >> 2;
    tmp.y = (mv.y + 1 + (mv.y < 0)) >> 2;

    tmp.x = tmp.x << 2;
    tmp.y = tmp.y << 2;

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
    int16_t pos;

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

/* Generic Affine control points and delta MVs derivation functions
 */

static struct AffineDeltaMV
derive_affine_delta_mvs(const struct AffineControlInfo *const cinfo,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t affine_type)
{
    struct AffineDeltaMV dmv;
    OVMV delta_mv_h, delta_mv_v;

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

    struct AffineControlInfo dst_cp;
    uint8_t ref_idx = ngh_cp->lt.ref_idx;

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

    OVMV tmp, lt_mv;

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

    tmp.x = lt_mv.x + (delta_mv.h.x << log2_pb_w);
    tmp.y = lt_mv.y + (delta_mv.h.y << log2_pb_w);

    dst_cp.rt = round_affine_mv2(tmp);
    dst_cp.rt = clip_mv(dst_cp.rt);
    dst_cp.rt.ref_idx = ref_idx;

    if (affine_type == AFFINE_3CP) {
        tmp.x = lt_mv.x + (delta_mv.v.x << log2_pb_h);
        tmp.y = lt_mv.y + (delta_mv.v.y << log2_pb_h);

        dst_cp.lb = round_affine_mv2(tmp);
        dst_cp.lb = clip_mv(dst_cp.lb);
        dst_cp.lb.ref_idx = ref_idx;
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
                       uint8_t affine_type)
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
    cp_info.lt = round_affine_mv(cp_info.lt);
    cp_info.rt = round_affine_mv(cp_info.rt);

    if (affine_type == AFFINE_3CP) {
        cp_info.lb = round_affine_mv(cp_info.lb);
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
        OVMV mv_cand = mv_ctx[cand_pos];

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

struct AffineControlInfo
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

    OVMV lt_mv_cand;
    OVMV rt_mv_cand;
    OVMV lb_mv_cand;

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
                                          affine_type);
    if (!cand_aff_lft) {
        cand_aff_lft = derive_affine_mvp_cand(affine_ctx, cp_info, pb_info, A1,
                                              inter_dir, ref_idx, ref_opp_idx,
                                              rpl0_cand, rpl1_cand, aff_cand_list,
                                              affine_type);
    }

    /* Affine above cand */
    cand_aff_abv = derive_affine_mvp_cand(affine_ctx, &cp_info[cand_aff_lft], pb_info, B0,
                                          inter_dir, ref_idx, ref_opp_idx,
                                          rpl0_cand, rpl1_cand, aff_cand_list,
                                          affine_type);
    if (!cand_aff_abv) {
        cand_aff_abv = derive_affine_mvp_cand(affine_ctx, &cp_info[cand_aff_lft], pb_info, B1,
                                              inter_dir, ref_idx, ref_opp_idx,
                                              rpl0_cand, rpl1_cand, aff_cand_list,
                                              affine_type);
        if (!cand_aff_abv) {
            cand_aff_abv = derive_affine_mvp_cand(affine_ctx, &cp_info[cand_aff_lft], pb_info, B2,
                                                  inter_dir, ref_idx, ref_opp_idx,
                                                  rpl0_cand, rpl1_cand, aff_cand_list,
                                                  affine_type);
        }
    }

    nb_cand = cand_aff_lft + cand_aff_abv;

    if (nb_cand >= MAX_NB_AMVP_CAND) {
        /* FIXME is rounding needed since done in cand derivation */
        for (int i = 0; i < nb_cand; i++) {
            cp_info[i].lt = round_affine_mv(cp_info[i].lt);
            cp_info[i].rt = round_affine_mv(cp_info[i].rt);
            cp_info[i].lb = round_affine_mv(cp_info[i].lb);
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

    mv_aff[CP_LT] = round_affine_mv(lt_mv_cand);
    mv_aff[CP_RT] = round_affine_mv(rt_mv_cand);
    mv_aff[CP_LB] = round_affine_mv(lb_mv_cand);

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
    uint8_t ph_tmvp_enabled_flag = 0;
    if (nb_cand < 2 && ph_tmvp_enabled_flag) {
        /* TODO retrieve TMVP cand C0 or C1 from current rpl */
        OVMV col_mv = {0};

        col_mv = round_affine_mv(col_mv);

        cp_info[nb_cand].lt = col_mv;
        cp_info[nb_cand].rt = col_mv;
        cp_info[nb_cand].lb = col_mv;

        nb_cand++;
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
        cp_info[i].lt = round_affine_mv(cp_info[i].lt);
        cp_info[i].rt = round_affine_mv(cp_info[i].rt);
        cp_info[i].lb = round_affine_mv(cp_info[i].lb);
    }
    return cp_info[mvp_idx];
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
    OVMV mv0[4];
    OVMV mv1[4];

    uint8_t dir = 0;

    switch (model_idx)
    {
        case 0:

            if (mi.mv0[CP_LT].ref_idx >= 0 && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_RT].ref_idx
                                           && mi.mv0[CP_LT].ref_idx == mi.mv0[CP_LB].ref_idx) {
                mv0[0] = mi.mv0[0];
                mv0[1] = mi.mv0[1];
                mv0[2] = mi.mv0[2];

                dir |= 0x1;
            }

            if (mi.mv1[CP_LT].ref_idx >= 0 && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_RT].ref_idx
                                           && mi.mv1[CP_LT].ref_idx == mi.mv1[CP_LB].ref_idx) {

                mv1[0] = mi.mv1[0];
                mv1[1] = mi.mv1[1];
                mv1[2] = mi.mv1[2];

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

        struct AffineControlInfo cp_info0;
        struct AffineControlInfo cp_info1;

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

        struct AffineControlInfo cp_info0;
        struct AffineControlInfo cp_info1;

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

        if (cand_msk = available_merge_b2_b3_a2(rpl0_cand, rpl1_cand)) {
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

        if (cand_msk = available_merge_b1_b0(rpl0_cand, rpl1_cand)) {
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

        if (cand_msk = available_merge_a1_a0(rpl0_cand, rpl1_cand)) {
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
            /* TODO retrieve TMVP cand C0 from both rpl */
            OVMV c0_mv = {0};
            uint8_t avail_dir;

            if (avail_dir & 0x1) {
                mi.mv0[3] = c0_mv;
                /* FIXME ref_idx scale*/
                avail_cp_map |= 0x8;
                dir[3] = 0x1;
            }

            if (avail_dir & 0x2) {
                mi.mv1[3] = c0_mv;
                /* FIXME ref_idx scale*/
                avail_cp_map |= 0x8;
                dir[3] |= 0x2;
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
                     OVMV *mv_buff,
                     uint8_t log2_cu_w, uint8_t log2_cu_h,
                     uint8_t inter_dir, uint8_t affine_type)
{
    /* Compute delta_mv from control points */
    /* TODO call before and give as an argument */
    const struct AffineDeltaMV delta_mv = derive_affine_delta_mvs(cinfo, log2_cu_w, log2_cu_h,
                                                                  affine_type);

    const uint8_t mv_broad = broadcast_mv(delta_mv, inter_dir);

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

                /* TODO Store sub blocks MVs for MCP */
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

        for (i = 0; i < nb_sb_h; ++i) {
            for (j = 0; j < nb_sb_w; ++j) {
                mv_buff[j] = center_mv;
            }
            mv_buff += 34;
        }
    }
}
        }
    }

    /* TODO Store CPInfo for later derivation */
    /* FIXME do this somewhere else ?*/
}
