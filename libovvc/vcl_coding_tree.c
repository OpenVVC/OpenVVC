#include <string.h>

#include "vcl.h"
#include "cabac_internal.h"
#include "dec_structures.h"
#include "ctudec.h"
#include "rcn_lmcs.h"

enum CUMode {
    OV_NA = 0xFF,
    OV_INTER = 1,
    OV_INTRA = 2,
    OV_INTER_SKIP = 3,
    OV_MIP = 4,
};

#if 0
static int coding_quadtree(OVCTUDec *const ctu_dec,
                           const OVPartInfo *const part_ctx,
                           unsigned int x0, unsigned int y0,
                           unsigned int log2_cb_s, unsigned int qt_depth);

static int coding_quadtree_implicit(OVCTUDec *const ctu_dec,
                                    const OVPartInfo *const part_ctx,
                                    unsigned int x0, unsigned int y0,
                                    unsigned int log2_cb_s, unsigned int qt_depth,
                                    unsigned int rem_w, unsigned int rem_h);

static int dual_tree(OVCTUDec *const ctu_dec,
                     const OVPartInfo *const part_ctx,
                     unsigned int x0, unsigned int y0,
                     unsigned int log2_cb_s, unsigned int qt_depth);

static int dual_tree_implicit(OVCTUDec *const ctu_dec,
                              const OVPartInfo *const part_ctx,
                              unsigned int x0, unsigned int y0,
                              unsigned int log2_cb_s, unsigned int qt_depth,
                              unsigned int rem_w,
                              unsigned int rem_h);

static void separate_tree_mtt(OVCTUDec *const ctu_dec,
                              unsigned int x0, unsigned int y0,
                              unsigned int log2_cb_w, unsigned int log2_cb_h,
                              unsigned int mtt_depth,
                              uint8_t implicit_mt_depth);
#endif
static uint8_t separate_trees_qt(OVCTUDec *ctudec, const OVPartInfo *const part_ctx,
                                 uint8_t x0, uint8_t y0,
                                 uint8_t log2_cb_w, uint8_t log2_cb_h, uint8_t split_cu_v);

/* FIXME cast to uint8_t and check result */
static int bt_split(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    unsigned int mtt_depth,
                    uint8_t implicit_mt_depth,
                    uint8_t split_cu_v);

/* FIXME cast to uint8_t and check result */
static int tt_split(OVCTUDec *const ctu_dec,
                    const OVPartInfo *const part_ctx,
                    unsigned int x0, unsigned int y0,
                    unsigned int log2_cb_w, unsigned int log2_cb_h,
                    unsigned int mtt_depth,
                    uint8_t implicit_mt_depth,
                    uint8_t split_cu_v);

static int multi_type_tree(OVCTUDec *const ctu_dec,
                           const OVPartInfo *const part_ctx,
                           uint8_t x0, uint8_t y0,
                           uint8_t log2_cb_w, uint8_t log2_cb_h,
                           uint8_t mtt_depth,
                           uint8_t middle_tt_flag, uint8_t implicit_mt_depth);

static int binary_tree_implicit_h(OVCTUDec *const ctu_dec,
                                  const OVPartInfo *const part_ctx,
                                  unsigned int x0, unsigned int y0,
                                  unsigned int log2_cb_w, unsigned int log2_cb_h,
                                  unsigned int mtt_depth, unsigned int rem_h);

static int binary_tree_implicit_v(OVCTUDec *const ctu_dec,
                                  const OVPartInfo *const part_ctx,
                                  unsigned int x0, unsigned int y0,
                                  unsigned int log2_cb_w, unsigned int log2_cb_h,
                                  unsigned int mtt_depth, unsigned int rem_w);


static uint8_t
ovcabac_read_ae_split_cu_flag(OVCABACCtx *const cabac_ctx,
                              uint8_t log2_cu_w_abv, uint8_t log2_cu_h_lft,
                              uint8_t log2_cu_w, uint8_t log2_cu_h,
                              uint8_t nb_split_cand)
{
    /* Note when not available log2_cu_h/w_lft/abv are set to 0xFF
     * so that they cannot be < to log2_cu_w/h
     */
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx_offset = (log2_cu_h_lft < log2_cu_h) + (log2_cu_w_abv < log2_cu_w) +
                     ((nb_split_cand >> 1) * 3);

    return ovcabac_ae_read(cabac_ctx, &cabac_state[SPLIT_FLAG_CTX_OFFSET + ctx_offset]);
}

static uint8_t
ovcabac_read_ae_split_qt_flag(OVCABACCtx *const cabac_ctx,
                         int8_t qt_depth_abv, int8_t qt_depth_lft,
                         int8_t qt_depth)
{
    /* Note when not available qt_depth_lft are set to 0 so that
     * qt_depth_lft cannot be > to current qt_depth
     */
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    uint8_t ctx_offset = (qt_depth_lft > qt_depth) + (qt_depth_abv > qt_depth) +
                         ((qt_depth < 2) ? 0 : 3);

    return ovcabac_ae_read(cabac_ctx, &cabac_state[SPLIT_QT_FLAG_CTX_OFFSET + ctx_offset]);
}

static uint8_t
ovcabac_read_ae_mtt_split_cu_vertical_flag(OVCABACCtx *const cabac_ctx,
                                      unsigned int log2_cu_w_abv,
                                      unsigned int log2_cu_h_lft,
                                      unsigned int log2_cu_w, unsigned int log2_cu_h,
                                      uint8_t nb_split_cand_v, uint8_t nb_split_cand_h)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx;

    /* FIXME try to avoid branch in this function */
    if (nb_split_cand_v == nb_split_cand_h) {
        int w_ratio = (1 << log2_cu_w) >> log2_cu_w_abv;
        int h_ratio = (1 << log2_cu_h) >> log2_cu_h_lft;

        if ((w_ratio == h_ratio) || (log2_cu_w_abv == 0xFF)
                                 || (log2_cu_h_lft == 0xFF)) {
            ctx = 0;
        } else if (w_ratio < h_ratio) {
            ctx = 1;
        } else {
            ctx = 2;
        }
    } else if (nb_split_cand_v < nb_split_cand_h) {
        ctx = 3;
    } else {
        ctx = 4;
    }

    return ovcabac_ae_read(cabac_ctx, &cabac_state[SPLIT_HV_FLAG_CTX_OFFSET + ctx]);
}

static uint8_t
ovcabac_read_ae_mtt_split_cu_binary_flag(OVCABACCtx *const cabac_ctx,
                                         uint8_t mtt_depth,
                                         uint8_t mtt_split_cu_vertical_flag)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    int ctx = (mtt_split_cu_vertical_flag << 1) | (mtt_depth <= 1);

    return ovcabac_ae_read(cabac_ctx, &cabac_state[SPLIT12_FLAG_CTX_OFFSET + ctx]);
}

static uint8_t
ovcabac_read_ae_mode_constraint(OVCABACCtx *const cabac_ctx, uint8_t intra_ngh)
{
    uint64_t *const cabac_state = cabac_ctx->ctx_table;
    return ovcabac_ae_read(cabac_ctx, &cabac_state[MODE_CONS_FLAG_CTX_OFFSET + intra_ngh]);
}

static void
store_qt_depth(OVCTUDec *const ctu_dec,
               const OVPartInfo *const part_ctx, uint8_t log2_cb_s,
               uint8_t x0, uint8_t y0, uint8_t qt_depth)
{
    struct PartMap *const part_map = ctu_dec->active_part_map;
    int log2_nb_pb_s = log2_cb_s - part_ctx->log2_min_cb_s;
    int x_cb = x0 >> part_ctx->log2_min_cb_s;
    int y_cb = y0 >> part_ctx->log2_min_cb_s;
    memset(&part_map->qt_depth_map_x[x_cb], qt_depth, sizeof(uint8_t) << log2_nb_pb_s);
    memset(&part_map->qt_depth_map_y[y_cb], qt_depth, sizeof(uint8_t) << log2_nb_pb_s);

    /* In single tree and when the conditions for separable tree have not been reached
     * we need to also store qt_depth for chroma so that the chroma partitioner
     * called in separable tree uses correct qt_depths
     */
    if (!ctu_dec->share && ctu_dec->coding_tree != &dual_tree &&
        ctu_dec->coding_tree_implicit != &dual_tree_implicit) {
        struct PartMap *const part_map_c = &ctu_dec->part_map_c;
        memset(&part_map_c->qt_depth_map_x[x_cb], qt_depth, sizeof(uint8_t) << log2_nb_pb_s);
        memset(&part_map_c->qt_depth_map_y[y_cb], qt_depth, sizeof(uint8_t) << log2_nb_pb_s);
    }
}


int
coding_quadtree(OVCTUDec *const ctu_dec,
               const OVPartInfo *const part_ctx,
               unsigned int x0, unsigned int y0,
               unsigned int log2_cb_s, unsigned int qt_depth)
{
    uint8_t split_cu_flag = 0;

    int x_cb = x0 >> part_ctx->log2_min_cb_s;
    int y_cb = y0 >> part_ctx->log2_min_cb_s;

    uint8_t allow_qt = log2_cb_s  > part_ctx->log2_min_qt_s;
    uint8_t allow_tt = log2_cb_s <= part_ctx->log2_max_tt_s;
    uint8_t allow_bt = log2_cb_s <= part_ctx->log2_max_bt_s;

    allow_tt &= ((log2_cb_s - 1) > part_ctx->log2_min_cb_s);
    allow_bt &= ( log2_cb_s      > part_ctx->log2_min_cb_s);
    /* Disable splits if less than 16 samples */
    allow_tt &= log2_cb_s > 2;
    allow_bt &= log2_cb_s > 2;

    if (ctu_dec->share == 2 && log2_cb_s + log2_cb_s == 5) {
        allow_bt = 0;
    }

    if (ctu_dec->share == 2 && log2_cb_s + log2_cb_s == 6) {
        allow_tt = 0;
    }

    uint8_t compute_chr_scale = (log2_cb_s == 6 && ctu_dec->lmcs_info.lmcs_enabled_flag) ;
    if (compute_chr_scale){
        rcn_lmcs_compute_chroma_scale(ctu_dec, x0, y0);
    }

    if (allow_qt | allow_bt | allow_tt) {
        OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
        const struct PartMap *const part_map = ctu_dec->active_part_map;
        uint8_t log2_cb_w_abv = part_map->log2_cu_w_map_x[x_cb];
        uint8_t log2_cb_h_lft = part_map->log2_cu_h_map_y[y_cb];
        int nb_split_cand = (allow_tt << 1) + (allow_bt << 1) + (allow_qt << 1) - 1;

        split_cu_flag = ovcabac_read_ae_split_cu_flag(cabac_ctx,
                                                 log2_cb_w_abv, log2_cb_h_lft,
                                                 log2_cb_s, log2_cb_s,
                                                 nb_split_cand);
        if (split_cu_flag) {
            uint8_t split_qt_flag = allow_qt;
            if (allow_qt && (allow_bt | allow_tt)) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                uint8_t depth_abv = part_map->qt_depth_map_x[x_cb];
                uint8_t depth_lft = part_map->qt_depth_map_y[y_cb];

                split_qt_flag = ovcabac_read_ae_split_qt_flag(cabac_ctx,
                                                              depth_abv, depth_lft,
                                                              qt_depth);
            }
            /* Note we can not merge both checks because split_qt_flag can be
             * implicit if both binary and ternary tree are disabled
             */
            if (split_qt_flag) {
                unsigned int x1 = x0 + (1 << (log2_cb_s - 1));
                unsigned int y1 = y0 + (1 << (log2_cb_s - 1));
                /* FIXME Find a better way to handle separable tree */
                #if 0
                uint8_t sep_tree = !ctu_dec->share && log2_cb_s == 3 &&
                                   ctu_dec->coding_tree_implicit != &dual_tree_implicit
                                   && ctu_dec->coding_tree != &dual_tree;
                #else
                uint8_t sep_tree = separate_trees_qt(ctu_dec, part_ctx, x0, y0, log2_cb_s, log2_cb_s, 1);
                #endif

                if (!ctu_dec->share) {
                    ctu_dec->share = sep_tree;
                }

                /* FIXME dirty hack for CCLM management */
                if (log2_cb_s == 5) {
                    ctu_dec->enable_cclm = 1;
                }

                if (sep_tree == 1) {
                    /*FIXME use specific function to launch chroma tree */
                    const OVPartInfo * part_ctx = ctu_dec->part_ctx;
                    void (*coding_unit_bkup)    = ctu_dec->coding_unit;
                    void (*transform_unit_bkup) = ctu_dec->transform_unit;

                    ctu_dec->coding_unit    = &coding_unit_intra;
                    ctu_dec->transform_unit = &transform_unit_l;
                    ctu_dec->active_part_map = &ctu_dec->part_map;
                    ctu_dec->tmp_disable_cclm     = 0;

                    coding_quadtree(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x1, y0, log2_cb_s - 1, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x0, y1, log2_cb_s - 1, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x1, y1, log2_cb_s - 1, qt_depth + 1);

                    part_ctx = ctu_dec->part_ctx_c;

                    ctu_dec->coding_unit   = &coding_unit_intra_c;
                    ctu_dec->transform_unit= &transform_unit_c;
                    ctu_dec->active_part_map = &ctu_dec->part_map_c;
                    ctu_dec->enable_cclm = 0;

                    coding_quadtree(ctu_dec, part_ctx, x0 >> 1, y0 >> 1, log2_cb_s - 2, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x1 >> 1, y0 >> 1, log2_cb_s - 2, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x0 >> 1, y1 >> 1, log2_cb_s - 2, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x1 >> 1, y1 >> 1, log2_cb_s - 2, qt_depth + 1);

                    ctu_dec->coding_unit = coding_unit_bkup;
                    ctu_dec->transform_unit = transform_unit_bkup;
                    ctu_dec->active_part_map = &ctu_dec->part_map;
                } else {
                    coding_quadtree(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x1, y0, log2_cb_s - 1, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x0, y1, log2_cb_s - 1, qt_depth + 1);
                    coding_quadtree(ctu_dec, part_ctx, x1, y1, log2_cb_s - 1, qt_depth + 1);
                }

                if (sep_tree) {
                    ctu_dec->share = 0;
                }

                return 1;
            }

            /* we enter the multi type tree we will not use QT anymore
             * the only allowed splits are bt or tt
             */
            store_qt_depth(ctu_dec, part_ctx, log2_cb_s, x0, y0, qt_depth);

            /* FIXME CCLM henadling */
            /* 64X64 luma node in dual tree will BT or TT split => disable CCLM */
            if (ctu_dec->coding_unit == coding_unit_intra && log2_cb_s == 6){
                ctu_dec->tmp_disable_cclm = 1;
            }

            multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_s, log2_cb_s,
                            0, 0, 0);
            return 1;
        }
    }

    /* FIXME CCLM henadling */
    if (log2_cb_s <= 5) {
        /* Last split before 32x32 was at least QT or 32x 32 did not split
         * => permit CCLM
         */
        ctu_dec->enable_cclm = 1;
    }

    coding_unit(ctu_dec, part_ctx, x0, y0, log2_cb_s, log2_cb_s);

    /* FIXME determine if we could store qt_depth before calling 
     * coding_unit function
     */
    store_qt_depth(ctu_dec, part_ctx, log2_cb_s, x0, y0, qt_depth);

    return 1;
}

int
coding_quadtree_implicit(OVCTUDec *const ctu_dec,
                        const OVPartInfo *const part_ctx,
                        unsigned int x0, unsigned int y0,
                        unsigned int log2_cb_s, unsigned int qt_depth,
                        unsigned int rem_w, unsigned int rem_h)
{
    #if 0
    uint8_t qt_split_cu_flag = 0;
    #endif
    uint8_t implicit_qt_split = 0;
    #if 0
    uint8_t implicit_bt_split = 0;
    #endif
    uint8_t force_implicit_qt = 0;

    unsigned int x1 = x0 + (1 << log2_cb_s);
    unsigned int y1 = y0 + (1 << log2_cb_s);

    uint8_t allow_qt = log2_cb_s >  part_ctx->log2_min_qt_s;
    uint8_t allow_bt = log2_cb_s <= part_ctx->log2_max_bt_s && log2_cb_s <= 6;

    uint8_t compute_chr_scale = (log2_cb_s == 6 && ctu_dec->lmcs_info.lmcs_enabled_flag) ;
    if (compute_chr_scale){
        rcn_lmcs_compute_chroma_scale(ctu_dec, x0, y0);
    }

    /* FIXME check everything is correct */
    if ((x1 > rem_w || y1 > rem_h)) {
        /* if remaining width or remaining height and min qt size is not reached
         * implicit split is quad tree split.
         */
        implicit_qt_split = 1;
            #if 0
        if (log2_cb_s <= part_ctx->log2_max_bt_s) {
            /* if bt size permits it implicit split is binary
             * according to required split direction
             * Priority is set to horizontal split if qt is not allowed and both
             * horizontal and vertical splits are required
             */
            implicit_bt_split = (x1 > rem_w) ? 1 : 2;
        }
            #endif
        force_implicit_qt = x1 > rem_w && y1 > rem_h;
    }

    if (implicit_qt_split) {
        uint8_t split_qt_flag = allow_qt;

        if (allow_bt && !force_implicit_qt && allow_qt) {

            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
            const struct PartMap *const part_map = ctu_dec->active_part_map;
            int x_pb = x0 >> part_ctx->log2_min_cb_s;
            int y_pb = y0 >> part_ctx->log2_min_cb_s;
            uint8_t depth_abv = part_map->qt_depth_map_x[x_pb];
            uint8_t depth_lft = part_map->qt_depth_map_y[y_pb];

            split_qt_flag = ovcabac_read_ae_split_qt_flag(cabac_ctx,
                                                          depth_abv, depth_lft,
                                                          qt_depth);

        }

        if (split_qt_flag) {
            x1 = x0 + (1 << (log2_cb_s - 1));
            y1 = y0 + (1 << (log2_cb_s - 1));

            /* FIXME dirty hack for CCLM management */
            if (log2_cb_s == 5) {
                ctu_dec->enable_cclm = 1;
            }

            if (x1 <= rem_w && y1 <= rem_h) {
                /* does not require implicit split anymore */
                coding_quadtree(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1,
                                qt_depth + 1);
            } else {
                /* we still need to check either one or the other direction for
                 * implicit split
                 */
                coding_quadtree_implicit(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1,
                                         qt_depth + 1, rem_w, rem_h);
            }

            /* if there is still some width to process continue  with splitting */
            if (x1 < rem_w)
                coding_quadtree_implicit(ctu_dec, part_ctx, x1, y0, log2_cb_s - 1,
                                         qt_depth + 1, rem_w, rem_h);
            if (y1 < rem_h)
                coding_quadtree_implicit(ctu_dec, part_ctx, x0, y1, log2_cb_s - 1,
                                         qt_depth + 1, rem_w, rem_h);

            if (x1 < rem_w && y1 < rem_h)
                coding_quadtree_implicit(ctu_dec, part_ctx, x1, y1, log2_cb_s - 1,
                                         qt_depth + 1, rem_w, rem_h);

            return 1;
        } else {

            store_qt_depth(ctu_dec, part_ctx, log2_cb_s, x0, y0, qt_depth);

            if (x1 > rem_w) {
                binary_tree_implicit_v(ctu_dec, part_ctx, x0, y0,
                                       log2_cb_s, log2_cb_s,
                                       0, rem_w);
            } else {
                binary_tree_implicit_h(ctu_dec, part_ctx, x0, y0,
                                       log2_cb_s, log2_cb_s,
                                       0, rem_h);
            }

            return 1;
        }
    } else {

        coding_quadtree(ctu_dec, part_ctx, x0, y0, log2_cb_s, qt_depth);

        return 1;
    }
}



int
dual_tree(OVCTUDec *const ctu_dec,
          const OVPartInfo *const part_ctx,
          unsigned int x0, unsigned int y0,
          unsigned int log2_cb_s, unsigned int qt_depth)
{
    if (log2_cb_s > 6) {
        unsigned int x1 = x0 + (1 << (log2_cb_s - 1));
        unsigned int y1 = y0 + (1 << (log2_cb_s - 1));

        dual_tree(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1, qt_depth + 1);
        dual_tree(ctu_dec, part_ctx, x1, y0, log2_cb_s - 1, qt_depth + 1);
        dual_tree(ctu_dec, part_ctx, x0, y1, log2_cb_s - 1, qt_depth + 1);
        dual_tree(ctu_dec, part_ctx, x1, y1, log2_cb_s - 1, qt_depth + 1);

        return 1;
    } else {
        const OVPartInfo * part_ctx = ctu_dec->part_ctx;

        ctu_dec->coding_unit    = &coding_unit_intra;
        ctu_dec->transform_unit = &transform_unit_l;
        ctu_dec->active_part_map = &ctu_dec->part_map;
        ctu_dec->tmp_disable_cclm     = 0;
        coding_quadtree(ctu_dec, part_ctx, x0, y0, log2_cb_s, qt_depth);

        part_ctx = ctu_dec->part_ctx_c;

        #if 0
        /* FIXME  Update LMCS for chroma */
        ctu_dec->compute_lmcs(ctu_dec, x0, y0);
        #endif

        ctu_dec->coding_unit   = &coding_unit_intra_c;
        ctu_dec->transform_unit= &transform_unit_c;
        ctu_dec->active_part_map = &ctu_dec->part_map_c;
        ctu_dec->enable_cclm = 0;

        coding_quadtree(ctu_dec, part_ctx, x0 >> 1, y0 >> 1, log2_cb_s - 1,
                        qt_depth);
    }
    return 1;
}

int
dual_tree_implicit(OVCTUDec *const ctu_dec,
                   const OVPartInfo *const part_ctx,
                   unsigned int x0, unsigned int y0,
                   unsigned int log2_cb_s, unsigned int qt_depth,
                   unsigned int rem_w,
                   unsigned int rem_h)
{

    if (log2_cb_s > 6) {
        unsigned int x1 = x0 + (1 << (log2_cb_s - 1));
        unsigned int y1 = y0 + (1 << (log2_cb_s - 1));

        if (x1 <= rem_w && y1 <= rem_h) {

            dual_tree(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1, qt_depth + 1);

        } else {

            dual_tree_implicit(ctu_dec, part_ctx, x0, y0, log2_cb_s - 1, qt_depth + 1,
                               rem_w, rem_h);
        }

        if (x1 < rem_w){
            dual_tree_implicit(ctu_dec, part_ctx, x1, y0, log2_cb_s - 1, qt_depth + 1,
                               rem_w, rem_h);
        }

        if (y1 < rem_h){
            dual_tree_implicit(ctu_dec, part_ctx, x0, y1, log2_cb_s - 1, qt_depth + 1,
                               rem_w, rem_h);
        }

        if (x1 < rem_w && y1 < rem_h){
            dual_tree_implicit(ctu_dec, part_ctx, x1, y1, log2_cb_s - 1, qt_depth + 1,
                               rem_w, rem_h);
        }

        return 1;

    } else {
        const OVPartInfo * part_ctx = ctu_dec->part_ctx;

        ctu_dec->coding_unit     = &coding_unit_intra;
        ctu_dec->transform_unit  = &transform_unit_l;
        ctu_dec->active_part_map = &ctu_dec->part_map;
        ctu_dec->tmp_disable_cclm     = 0;

        coding_quadtree_implicit(ctu_dec, part_ctx, x0, y0, log2_cb_s, qt_depth,
                                 rem_w, rem_h);

        #if 0
        ctu_dec->compute_lmcs(ctu_dec, x0, y0);
        #endif

        part_ctx = ctu_dec->part_ctx_c;

        ctu_dec->coding_unit     = &coding_unit_intra_c;
        ctu_dec->transform_unit  = &transform_unit_c;
        ctu_dec->active_part_map = &ctu_dec->part_map_c;
        ctu_dec->enable_cclm = 0;

        coding_quadtree_implicit(ctu_dec, part_ctx, x0 >> 1, y0 >> 1, log2_cb_s - 1, qt_depth,
                                 rem_w >> 1, rem_h >> 1);
    }

    return 1;
}

static void
separate_tree_mtt(OVCTUDec *const ctu_dec,
                   unsigned int x0, unsigned int y0,
                   unsigned int log2_cb_w, unsigned int log2_cb_h,
                   unsigned int mtt_depth,
                   uint8_t implicit_mt_depth)
{
    void (*coding_unit_bkup)    = ctu_dec->coding_unit;
    void (*transform_unit_bkup) = ctu_dec->transform_unit;

    ctu_dec->coding_unit     = &coding_unit_intra_c;
    ctu_dec->transform_unit  = &transform_unit_c;
    ctu_dec->active_part_map = &ctu_dec->part_map_c;

    multi_type_tree(ctu_dec, ctu_dec->part_ctx_c, x0 >> 1, y0 >> 1,
                    log2_cb_w - 1, log2_cb_h - 1,
                    mtt_depth, 0, implicit_mt_depth);

    ctu_dec->active_part_map = &ctu_dec->part_map;
    ctu_dec->coding_unit     = coding_unit_bkup;
    ctu_dec->transform_unit  = transform_unit_bkup;
    ctu_dec->share = 0;
}

/* FIXME cast to uint8_t and check result */
static int
bt_split(OVCTUDec *const ctu_dec,
         const OVPartInfo *const part_ctx,
         unsigned int x0, unsigned int y0,
         unsigned int log2_cb_w, unsigned int log2_cb_h,
         unsigned int mtt_depth,
         uint8_t implicit_mt_depth,
         uint8_t split_cu_v)
{
    const unsigned split_mask_h = -(!split_cu_v);
    const unsigned split_mask_v = ~split_mask_h;

    unsigned int log2_cb_w1 = log2_cb_w - (split_mask_v & 1);
    unsigned int log2_cb_h1 = log2_cb_h - (split_mask_h & 1);

    unsigned int x1 = x0 + ((1 << log2_cb_w1) & split_mask_v);
    unsigned int y1 = y0 + ((1 << log2_cb_h1) & split_mask_h);

    /* FIXME CCLM Handling */
    uint8_t cclm_state_bckp = ctu_dec->enable_cclm;

    multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w1, log2_cb_h1,
                    mtt_depth + 1, 0, implicit_mt_depth);

    /* FIXME CCLM Handling */
    /* restore first split is bt h */
    ctu_dec->enable_cclm = cclm_state_bckp;

    multi_type_tree(ctu_dec, part_ctx, x1, y1, log2_cb_w1, log2_cb_h1,
                    mtt_depth + 1, 0, implicit_mt_depth);

    return 1;
}

/* FIXME cast to uint8_t and check result */
static int
tt_split(OVCTUDec *const ctu_dec,
         const OVPartInfo *const part_ctx,
         unsigned int x0, unsigned int y0,
         unsigned int log2_cb_w, unsigned int log2_cb_h,
         unsigned int mtt_depth,
         uint8_t implicit_mt_depth,
         uint8_t split_cu_v)
{
    const unsigned split_mask_h = -(!split_cu_v);
    const unsigned split_mask_v = ~split_mask_h;

    unsigned int log2_cb_w1 = log2_cb_w - (split_mask_v & 2);
    unsigned int log2_cb_h1 = log2_cb_h - (split_mask_h & 2);
    unsigned int x1 = x0 + ((1 << log2_cb_w1) & split_mask_v);
    unsigned int y1 = y0 + ((1 << log2_cb_h1) & split_mask_h);

    unsigned int log2_cb_w2 = log2_cb_w - (split_mask_v & 1);
    unsigned int log2_cb_h2 = log2_cb_h - (split_mask_h & 1);
    unsigned int x2 = x1 + ((1 << log2_cb_w2) & split_mask_v);
    unsigned int y2 = y1 + ((1 << log2_cb_h2) & split_mask_h);

    multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w1, log2_cb_h1,
                    mtt_depth + 1, 0, implicit_mt_depth);

    multi_type_tree(ctu_dec, part_ctx, x1, y1, log2_cb_w2, log2_cb_h2,
                    mtt_depth + 1, 1 << (!split_cu_v), implicit_mt_depth);

    multi_type_tree(ctu_dec, part_ctx, x2, y2, log2_cb_w1, log2_cb_h1,
                    mtt_depth + 1, 0, implicit_mt_depth);

    return 1;
}

static uint8_t
separate_trees_tt(OVCTUDec *ctudec, const OVPartInfo *const part_ctx, uint8_t x0, uint8_t y0,
                  uint8_t log2_cb_w, uint8_t log2_cb_h, uint8_t split_cu_v)
{
    uint8_t log2_cb_s = log2_cb_w + log2_cb_h;
    uint16_t luma_area = (1 << log2_cb_s) >> 2;
    uint16_t chroma_area = luma_area >> 2;

    if (ctudec->share || ctudec->coding_tree == &dual_tree) {
        return 0;
    }

    if (chroma_area >= 16 && !(split_cu_v && log2_cb_w == 4)) {
        return 0;
    } else if (luma_area < 32 || ctudec->coding_unit == &coding_unit_intra_st) {
        return 1;
    } else {
       /* signal */
        OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
        int x_cb = x0 >> part_ctx->log2_min_cb_s;
        int y_cb = y0 >> part_ctx->log2_min_cb_s;
        uint8_t cu_type_abv = ctudec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctudec->part_map.cu_mode_y[y_cb];
        return 2 >> ovcabac_read_ae_mode_constraint(cabac_ctx, cu_type_abv == OV_INTRA || cu_type_abv == OV_MIP || cu_type_lft == OV_INTRA || cu_type_lft == OV_MIP);
    }
}

static uint8_t
separate_trees_qt(OVCTUDec *ctudec, const OVPartInfo *const part_ctx, uint8_t x0, uint8_t y0,
                  uint8_t log2_cb_w, uint8_t log2_cb_h, uint8_t split_cu_v)
{
    uint8_t log2_cb_s = log2_cb_w + log2_cb_h;
    uint16_t luma_area = (1 << log2_cb_s) >> 2;
    uint16_t chroma_area = luma_area >> 2;
    uint8_t sep_tree = 0;

    if (ctudec->share || ctudec->coding_tree == &dual_tree) {
        return 0;
    }

    if (chroma_area >= 16) {
        return 0;
    } else if (luma_area < 32 || ctudec->coding_unit == &coding_unit_intra_st) {
        /* infer */
        return 1;
    } else {
        /* signal */
        OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
        int x_cb = x0 >> part_ctx->log2_min_cb_s;
        int y_cb = y0 >> part_ctx->log2_min_cb_s;
        uint8_t cu_type_abv = ctudec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctudec->part_map.cu_mode_y[y_cb];
        return 2 >> ovcabac_read_ae_mode_constraint(cabac_ctx, cu_type_abv == OV_INTRA || cu_type_abv == OV_MIP || cu_type_lft == OV_INTRA || cu_type_lft == OV_MIP);
    }
}

static uint8_t
separate_trees_bt(OVCTUDec *ctudec, const OVPartInfo *const part_ctx, uint8_t x0, uint8_t y0,
                  uint8_t log2_cb_w, uint8_t log2_cb_h, uint8_t split_cu_v)
{
    uint8_t log2_cb_s = log2_cb_w + log2_cb_h;
    uint16_t luma_area = (1 << log2_cb_s) >> 1;
    uint16_t chroma_area = luma_area >> 2;
    uint8_t sep_tree = 0;

    if (ctudec->share || ctudec->coding_tree == &dual_tree) {
        return 0;
    }

    if (chroma_area >= 16 && !(split_cu_v && log2_cb_w == 3)) {
        return 0;
    } else if (luma_area < 32 || ctudec->coding_unit == &coding_unit_intra_st) {
        /* infer */
        return 1;
    } else {
        /* signal */
        OVCABACCtx *const cabac_ctx = ctudec->cabac_ctx;
        int x_cb = x0 >> part_ctx->log2_min_cb_s;
        int y_cb = y0 >> part_ctx->log2_min_cb_s;
        uint8_t cu_type_abv = ctudec->part_map.cu_mode_x[x_cb];
        uint8_t cu_type_lft = ctudec->part_map.cu_mode_y[y_cb];
        return 2 >> ovcabac_read_ae_mode_constraint(cabac_ctx, cu_type_abv == OV_INTRA || cu_type_abv == OV_MIP || cu_type_lft == OV_INTRA || cu_type_lft == OV_MIP);
    }
}

static int
multi_type_tree(OVCTUDec *const ctu_dec,
               const OVPartInfo *const part_ctx,
               uint8_t x0, uint8_t y0,
               uint8_t log2_cb_w, uint8_t log2_cb_h,
               uint8_t mtt_depth,
               uint8_t middle_tt_flag, uint8_t implicit_mt_depth)
{
    uint8_t can_split = mtt_depth - implicit_mt_depth < part_ctx->max_mtt_depth;

    if (ctu_dec->share == 1 && ctu_dec->active_part_map == &ctu_dec->part_map_c) {
        can_split = 0;
    }

    uint8_t allow_tt_v = 0, allow_tt_h = 0, allow_bt_h = 0, allow_bt_v = 0;

    if (can_split) {
        /* Ternary Tree needs both width and height to be less than max_tt_size
         * to be permitted
         */
        uint8_t allow_tt = log2_cb_w <= part_ctx->log2_max_tt_s &&
                           log2_cb_h <= part_ctx->log2_max_tt_s;

        allow_tt &= log2_cb_w <= 6 && log2_cb_h <= 6;
        allow_tt &= !(ctu_dec->share == 2 && log2_cb_w + log2_cb_h == 6);

        allow_tt_v = allow_tt && ((log2_cb_w - 1) > part_ctx->log2_min_cb_s);
        allow_tt_h = allow_tt && ((log2_cb_h - 1) > part_ctx->log2_min_cb_s);

        allow_tt_v &= !(ctu_dec->coding_tree == dual_tree && ctu_dec->coding_unit == coding_unit_intra_c && log2_cb_w == 3);

        /* BT split is disabled in the middle of TT partition of same direction
         */
        allow_bt_v = (log2_cb_w >  part_ctx->log2_min_cb_s) &&
                     (log2_cb_w <= part_ctx->log2_max_bt_s) &&
                     !(middle_tt_flag == 1);

        allow_bt_h = (log2_cb_h >  part_ctx->log2_min_cb_s) &&
                     (log2_cb_h <= part_ctx->log2_max_bt_s) &&
                     !(middle_tt_flag == 2);

        allow_bt_v &= !(ctu_dec->coding_tree == dual_tree && ctu_dec->coding_unit == coding_unit_intra_c && log2_cb_w == 2);

        /* Disable splits leading to CU with less than 16 samples */
        if ((log2_cb_h + log2_cb_w) <= 4) {
            allow_bt_v = allow_bt_h = 0;
        }

        if ((log2_cb_h + log2_cb_w - 1) <= 4) {
            allow_tt_v = allow_tt_h = 0;
        }

        if (ctu_dec->share == 2 && log2_cb_w + log2_cb_h == 5) {
            allow_bt_v = allow_bt_h = 0;
        }

        if (log2_cb_h > 6 && log2_cb_w <= 6) {
            allow_bt_v = 0;
        }

        if (log2_cb_w > 6 && log2_cb_h <= 6) {
            allow_bt_h = 0;
        }

        can_split = allow_bt_v | allow_bt_h | allow_tt_v | allow_tt_h;
    }

    if (can_split) {
        int x_cb = x0 >> part_ctx->log2_min_cb_s;
        int y_cb = y0 >> part_ctx->log2_min_cb_s;


        /* We force first split when mtt_depth equals to 0 since split_cu_flag
         * was already read in quadtree function
         */
        uint8_t split_cu_flag = !mtt_depth;
        if (!split_cu_flag) {
            OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;

            const struct PartMap *const part_map = ctu_dec->active_part_map;
            uint8_t log2_cb_w_abv = part_map->log2_cu_w_map_x[x_cb];
            uint8_t log2_cb_h_lft = part_map->log2_cu_h_map_y[y_cb];
            int nb_split_cand = allow_bt_v + allow_bt_h + allow_tt_v + allow_tt_h - 1;

            split_cu_flag = ovcabac_read_ae_split_cu_flag(cabac_ctx,
                                                          log2_cb_w_abv, log2_cb_h_lft,
                                                          log2_cb_w, log2_cb_h,
                                                          nb_split_cand);
        }

        if (split_cu_flag) {
            uint8_t can_split_v = allow_tt_v | allow_bt_v;
            uint8_t can_split_h = allow_tt_h | allow_bt_h;
            uint8_t split_cu_v = can_split_v;
            uint8_t split_cu_bt = allow_bt_h | allow_bt_v;

            if (can_split_v & can_split_h) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                const struct PartMap *const part_map = ctu_dec->active_part_map;
                uint8_t log2_cb_w_abv  = part_map->log2_cu_w_map_x[x_cb];
                uint8_t log2_cb_h_lft  = part_map->log2_cu_h_map_y[y_cb];
                split_cu_v = ovcabac_read_ae_mtt_split_cu_vertical_flag(cabac_ctx,
                                                                        log2_cb_w_abv,
                                                                        log2_cb_h_lft,
                                                                        log2_cb_w, log2_cb_h,
                                                                        allow_bt_v + allow_tt_v,
                                                                        allow_bt_h + allow_tt_h);
            }

            if (split_cu_v & (allow_tt_v & allow_bt_v) ||
               (!split_cu_v) & (allow_tt_h & allow_bt_h)) {
                OVCABACCtx *const cabac_ctx = ctu_dec->cabac_ctx;
                split_cu_bt = ovcabac_read_ae_mtt_split_cu_binary_flag(cabac_ctx,
                                                                  mtt_depth, split_cu_v);
            } else {
                split_cu_bt = ((!split_cu_v) & allow_bt_h) | (split_cu_v & allow_bt_v);
            }

            if (split_cu_bt) {
                unsigned int log2_nb_s = log2_cb_w + log2_cb_h;
                /* FIXME Separable tree */
#if 0
                uint8_t sep_tree = !ctu_dec->share && ((log2_nb_s == 6) || (log2_nb_s == 5) ||
                                    (split_cu_v && log2_cb_w == 3)) &&
                                    ctu_dec->coding_tree_implicit != &dual_tree_implicit
                                    && ctu_dec->coding_tree != &dual_tree;
#else
                uint8_t sep_tree = separate_trees_bt(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, split_cu_v);
#endif

                /* FIXME Separable tree */
                if (!ctu_dec->share && sep_tree) {
                    ctu_dec->share = sep_tree;
                }

                /* FIXME CCLM */
                if (!split_cu_v && (mtt_depth == 0) && ctu_dec->enable_cclm != 1) {
                    /* first split on 32x32 node is bt h */
                    ctu_dec->enable_cclm = 2;
                }

                if (split_cu_v && (mtt_depth == 1) && (ctu_dec->enable_cclm == 2)) {
                    /* second split after bt h is bt_v */
                    ctu_dec->enable_cclm = 1;
                }

                /* FIXME Separable tree */
                void (*transform_unit_bkup) = ctu_dec->transform_unit;
                void (*coding_unit_bkup) = ctu_dec->coding_unit;
                if (sep_tree == 1) {
                    ctu_dec->transform_unit = &transform_unit_l;
                    ctu_dec->coding_unit = &coding_unit_intra;
                }

                bt_split(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h,
                         mtt_depth, implicit_mt_depth, split_cu_v);

                /* FIXME Separable tree */
                ctu_dec->transform_unit = transform_unit_bkup;
                ctu_dec->coding_unit = coding_unit_bkup;

                /* FIXME Separable tree */
                if (sep_tree == 1) {
                    separate_tree_mtt(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                                      mtt_depth, implicit_mt_depth);
                }

                if  (sep_tree)  {
                    ctu_dec->share = 0;
                }

                return 1;

            } else {
                unsigned int log2_nb_s = log2_cb_w + log2_cb_h;
                /* FIXME separable tree */
#if 0
                uint8_t sep_tree = !ctu_dec->share && ((log2_nb_s == 7) || (log2_nb_s == 6) ||
                                    (split_cu_v && log2_cb_w == 4)) &&
                                    ctu_dec->coding_tree_implicit != &dual_tree_implicit
                                    && ctu_dec->coding_tree != &dual_tree;
#else
                uint8_t sep_tree = separate_trees_tt(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h, split_cu_v);
#endif

                if (!ctu_dec->share && sep_tree) {
                    ctu_dec->share = sep_tree;
                }

                void (*transform_unit_bkup) = ctu_dec->transform_unit;
                void (*coding_unit_bkup) = ctu_dec->coding_unit;
                if (sep_tree == 1) {
                    ctu_dec->transform_unit = &transform_unit_l;
                    ctu_dec->coding_unit = &coding_unit_intra;
                }
                /* end of FIXME */

                tt_split(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h,
                         mtt_depth, implicit_mt_depth, split_cu_v);

                /* FIXME separable tree */
                ctu_dec->transform_unit = transform_unit_bkup; 
                ctu_dec->coding_unit = coding_unit_bkup;

                if (sep_tree == 1) {
                    separate_tree_mtt(ctu_dec, x0, y0, log2_cb_w, log2_cb_h,
                                      mtt_depth, implicit_mt_depth);
                }

                if  (sep_tree)  {
                    ctu_dec->share = 0;
                }
                /* end of FIXME */
                return 1;
            }
        }
    }

    /* FIXME find an other way of testing these conditions */
    if ((mtt_depth == 1) && ctu_dec->enable_cclm == 2) {
        /* first split is bt_h and no further split */
        ctu_dec->enable_cclm = 1;
    }

    coding_unit(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h);

    return 0;
}

static int
binary_tree_implicit_h(OVCTUDec *const ctu_dec,
                       const OVPartInfo *const part_ctx,
                       unsigned int x0, unsigned int y0,
                       unsigned int log2_cb_w, unsigned int log2_cb_h,
                       unsigned int mtt_depth, unsigned int rem_h)
{

    /* Split is implicit if the end of current cb is outside of dimension */
    uint8_t implicit_split = y0 + (1 << log2_cb_h) <= rem_h ? 0 : 1;

    if (implicit_split) {

        unsigned int log2_cb_h1 = log2_cb_h - 1;
        unsigned int y1 = y0 + (1 << log2_cb_h1);

        if ((mtt_depth == 0) && ctu_dec->enable_cclm != 1)
            ctu_dec->enable_cclm = 2;/*first split is bt_h */

        int enable_cclm = ctu_dec->enable_cclm;

        /* Split is implicit if the end of current cb is inside the avalable area */
        /* we can call the multi type tree for the first split part */
        if (y1 <= rem_h) {

            multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h1,
                            mtt_depth + 1, 0, mtt_depth + 1);
        } else { /* We need to further split the current partition */

            binary_tree_implicit_h(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h1,
                                   mtt_depth + 1, rem_h);
        }

        ctu_dec->enable_cclm = enable_cclm;

        /* if there exists a second binary partition, we need to check whether */
        /* or not it requires more implicit splits */
        if (y1 < rem_h)
            binary_tree_implicit_h(ctu_dec, part_ctx, x0, y1, log2_cb_w, log2_cb_h1,
                                   mtt_depth + 1, rem_h);
        return 1;
    } else {
        multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h,
                        mtt_depth, 0, mtt_depth);
    }
    return 0;
}

static int
binary_tree_implicit_v(OVCTUDec *const ctu_dec,
                       const OVPartInfo *const part_ctx,
                       unsigned int x0, unsigned int y0,
                       unsigned int log2_cb_w, unsigned int log2_cb_h,
                       unsigned int mtt_depth, unsigned int rem_w)
{
    uint8_t implicit_split = x0 + (1 << log2_cb_w) <= rem_w ? 0 : 1;

    if (implicit_split) {
        unsigned int log2_cb_w1 = (log2_cb_w - 1);
        unsigned int x1 = x0 + (1 << log2_cb_w1);

        if (x1 <= rem_w) {
            multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w1, log2_cb_h,
                            mtt_depth + 1, 0, mtt_depth + 1);
        } else {
            binary_tree_implicit_v(ctu_dec, part_ctx, x0, y0, log2_cb_w1, log2_cb_h,
                                   mtt_depth + 1, rem_w);
        }

        if (x1 < rem_w)
            binary_tree_implicit_v(ctu_dec, part_ctx, x1, y0, log2_cb_w1, log2_cb_h,
                                   mtt_depth + 1, rem_w);
        return 1;
    } else {
        multi_type_tree(ctu_dec, part_ctx, x0, y0, log2_cb_w, log2_cb_h,
                        mtt_depth, 0, mtt_depth);
    }
    return 0;
}

