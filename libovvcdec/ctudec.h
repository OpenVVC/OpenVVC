#ifndef CTU_DEC_H
#define CTU_DEC_H

#include "ovdefs.h"

#include "ovvcdec.h"

    /* FIXME COMPAT old
     * Old structures to be removed
     * Those structures were imported as is from previous version
     * because they were mandatory for a quickly operational decoder
     * and are still used by the decoder however they will be removed
     * in the future
     */

#define FLG_STORE(name, dst) \
   dst |= -(!!name) & flg_##name; 

#define DECL_FLG(name, pos)\
    flg_##name = (1llu << (pos))

enum VVCCUFlag{
     DECL_FLG(cu_skip_flag,0),
     DECL_FLG(pred_mode_flag,1),
     DECL_FLG(mip_flag,2),
     DECL_FLG(isp_flag,3),
     DECL_FLG(mrl_flag,4),
     DECL_FLG(mpm_flag,5),
     DECL_FLG(cclm_flag,6),
     DECL_FLG(mpm_flag_c,7),
     DECL_FLG(merge_flag,2),
     DECL_FLG(inter_dir,3),
};

typedef struct VVCDeQuantCtx{
    uint8_t qp;
}VVCDeQuantCtx;

typedef struct VVCQPCTX{
    int8_t current_qp;
    int8_t min_qp_prime_ts;
    int8_t cb_offset;
    int8_t cr_offset;
    int8_t jcbcr_offset;
    const int8_t *chroma_qp_map_cb;
    const int8_t *chroma_qp_map_cr;
    const int8_t *chroma_qp_map_jcbcr;
}VVCQPCTX;

struct VVCCU{

/*  intra_flags
 *  _ _ _ _ _ _ _ _
 * |_|_|_|_|_|_|_|_|
 *  | | | | | | | |
 *  | | | | | | | cu_skip_flag
 *  | | | | | | pred_mode_flag
 *  | | | | | mip_flag
 *  | | | | isp_flag
 *  | | | mrl_flag
 *  | | mpm_flag
 *  | cclm_flag
 *  mpm_flag_c
 */

/*  inter_flag
 *  _ _ _ _ _ _ _ _
 * |_|_|_|_|_|_|_|_|
 *  | | | | | | | |
 *  | | | | | | | cu_skip_flag
 *  | | | | | | pred_mode_flag
 *  | | | | | merge_flag
 *  | | | | inter_dir
 *  | | |
 *  | |
 *  |
 */

    uint8_t cu_flags;
    uint8_t cu_mode_idx;
    uint8_t cu_opaque;
    uint8_t cu_mode_info;

};

/* Definitions of CTUDecoder structures */
struct InterDRVCtx
{
   struct RPLInfo *rpl_0;
   struct RPLInfo *rpl_1;
   /* Information to be used by motion vector predicition (MVP)
    * Contains Motion Vector and TMVP if needed
    */
   struct MVPCtx *mvp_context;

   /* HMVP lut info */
   #if 0
   struct HMVP hmvp_lut;
   #endif
};

struct IntraDRVCtx;

struct FiltersDRVCtx
{
    struct OVDBFInfo *dbf_ctx;
    struct OVALFInfo *alf_ctx;
    struct OVSAOInfo *sao_ctx;
    struct OVLMCSInfo *lmcs_ctx;
};

/* Storage for transform coefficient
 * TODO this is intended to replace 
 * it in ctu dec
 */
struct TrCoeffData{
    int16_t residual_y[64*64];
    int16_t residual_cb[64*64];
    int16_t residual_cr[64*64];
    int16_t lfnst_subblock[16*2];

};

struct OVCTUDec {
    /**
     * Pointer to the parent frame decoder
     * Note this should be only be used only to retrieve
     * parent decoder error handling
     */
    const struct OVDec *parent_dec;
    /* TODO decide whether or not we should add thread info
     * here
     */

    /* Associated cabac_context
     * Contains context tables, cabac status and position in the
     * current entry
     */
    struct OVCABACCtx *cabac_ctx;

    /* List of flag/modes to match activated tools */
    struct {
        uint8_t flag;
        uint8_t status;
    } tools_status;

    /* Derivations context according to activated tools
    */
    struct OVDrvCtx {
        /* Pointers to Information used by specific tools */
        struct InterDRVCtx *inter_ctx;
        struct IntraDRVCtx *intra_ctx;
        struct DeltaQPDRVCtx *delta_qp_ctx;
        struct FiltersDRVCtx *loop_filters_ctx;
    } drv_ctx;

    /* Reconstruction context */
    struct OVRCNCtx {
        /* Pointers to the first sample data of CTU in the current
         * picture
         */
        struct {
            uint16_t *data_y;
            uint16_t *data_cb;
            uint16_t *data_cr;
        } ctu_pos;

        /* Pointers to CTU reconstruction buffers to be used
         * to reconstruct current CTU.
         * These buffers will be written to the destination picture
         * before filtering operation
         */
        struct {
            uint16_t *data_y;
            uint16_t *data_cb;
            uint16_t *data_cr;
        } rcn_ctu_buff;

        /* Side Buffer to be used by reconstruction functions
         * when needed
         */
        struct {
            uint16_t *data;
        } tmp_buff;

        /* Bit fields corresponding to the decoding progress in
         * current CTU, and its borders those are used for example
         * in order to derive references samples for intra prediction
         */
        struct CTUBitField{
            uint64_t hfield[33];
            uint64_t vfield[33];
        } progress_field;
    } rcn_ctx;

    /* CTU neighbours availability flags
     * An aggregation of flag used to tell the decoder if
     * the CTU neighbours are supposed to be known from
     * the current CTU decoder.
     */
    uint8_t ctu_ngh_flags;

    /* FIXME COMPAT old passed this line
     * Old structures to be removed
     * Those structures were imported as is from previous version
     * because they were mandatory for a quickly operational decoder
     * and are still used by the decoder however they will be removed
     * in the future
     */
    /* FIXME
     * Check usage of those sort between flags and status
     * and move them to tools_status structure
     * (as an example transform_skip_flag can be replaced
     * by a max_log2_trskip_s to zere, mts_implicit and
     * mts_enabled could be replace by a mts status etc.)
     */
    uint8_t transform_skip_enabled;
    uint8_t max_log2_transform_skip_size;

    uint8_t enable_sdh;
    uint8_t jcbcr_enabled;
    uint8_t isp_enabled;
    uint8_t mts_implicit;
    uint8_t mts_enabled;
    uint8_t delta_qp_enabled;
    uint8_t enable_lfnst;
    uint8_t enabled_mip;
    uint8_t lm_chroma_enabled;
    uint8_t tmp_disable_cclm;
    uint8_t enable_cclm;
    uint8_t enable_mrl;
    uint8_t share;
    uint8_t max_num_merge_candidates;
    uint8_t dbf_disable;

    /**
     * Depths of left and up neighbours during in the decision tree
     * needed to derive cabac contexts for split decision
     * and mode context derivation
     * TODO rename as CABAC Maps
     */
    struct PartMap{
        uint8_t *qt_depth_map_x;
        uint8_t *log2_cu_w_map_x;
        uint8_t *cu_mode_x;
        uint8_t qt_depth_map_y[32];
        uint8_t log2_cu_h_map_y[32];
        uint8_t cu_mode_y[32];
    } part_map;

    struct PartMap part_map_c;

    /* Pointer to active part map to be used in current
     * tree to derive cabac context
     */
    struct PartMap *active_part_map;

    const struct OVPartInfo *part_ctx;

    /**
     * Chroma partition context for dual tree
     */
    const struct OVPartInfo *part_ctx_c;

    /**
     * Pointer to recursive coding tree structure based on the slice type and
     * qtbt qtbt_dual_intra_flag
     */
    int (*coding_tree)(struct OVCTUDec *const lc_ctx,
                       const OVPartInfo *const part_ctx,
                       unsigned int x0, unsigned int y0,
                       unsigned int log2_cb_size, unsigned int qt_depth);

    int (*coding_tree_implicit)(struct OVCTUDec *const lc_ctx,
                                const OVPartInfo *const part_ctx,
                                unsigned int x0, unsigned int y0,
                                unsigned int log2_cb_size, unsigned int qt_depth,
                                unsigned int remaining_width,
                                unsigned int remaining_height);

    VVCCU (*coding_unit)(struct OVCTUDec *const lc_ctx,
                         const OVPartInfo *const part_ctx,
                         uint8_t x0, uint8_t y0,
                         uint8_t log2_cb_w, uint8_t log2_cb_h);

    int (*transform_unit)(struct OVCTUDec *const lc_ctx,
                          unsigned int x0, unsigned int y0,
                          unsigned int log2_tb_w,
                          unsigned int log2_tb_h,
                          uint8_t cbf_ctx, uint8_t cu_flags);

    int (*prediction_unit)(struct OVCTUDec *const lc_ctx,
                           const OVPartInfo *const part_ctx,
                           uint8_t x0, uint8_t y0,
                           uint8_t log2_pb_w, uint8_t log2_pb_h,
                           uint8_t cu_skip_flag);

    int (*residual_coding_isp_h)(struct OVCTUDec *const lc_ctx, uint16_t *const dst,
                                 unsigned int log2_tb_w, unsigned int log2_tb_h,
                                 uint16_t last_pos);

    int (*residual_coding_isp_v)(struct OVCTUDec *const lc_ctx, uint16_t *const dst,
                                 unsigned int log2_tb_w, unsigned int log2_tb_h,
                                 uint16_t last_pos);

    int (*residual_coding_ts)(struct OVCTUDec *const lc_ctx,
                              unsigned int log2_tb_w,
                              unsigned int log2_tb_h);

    uint64_t (*residual_coding)(struct OVCTUDec *const lc_ctx, uint16_t *const dst,
                                unsigned int log2_tb_w, unsigned int log2_tb_h,
                                uint16_t last_pos);

    int (*residual_coding_chroma)(struct OVCTUDec *const lc_ctx, uint16_t *const dst,
                                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                                  uint16_t last_pos);

    int16_t residual_y[64*64];
    int16_t residual_cb[64*64];
    int16_t residual_cr[64*64];
    int16_t lfnst_subblock[16*2];
    int16_t transform_buff[64*64];


    uint8_t slice_qp;
    VVCQPCTX qp_ctx;
    VVCDeQuantCtx dequant_luma;
    VVCDeQuantCtx dequant_luma_skip;
    VVCDeQuantCtx dequant_cb;
    VVCDeQuantCtx dequant_cr;
    VVCDeQuantCtx dequant_joint_cb_cr;

    const VVCDeQuantCtx *dequant_chroma;
};

int ovdec_decode_ctu(OVVCDec *dec, OVCTUDec *ctu_dec);

int ctudec_init(OVCTUDec **ctudec_p);
int ctudec_uninit(OVCTUDec *ctudec_p);

/* TODO add hook functors for PU / TU init and update */

#endif
