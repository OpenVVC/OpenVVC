#ifndef CTU_DEC_H
#define CTU_DEC_H

#include "ovdefs.h"

#include "ovdec.h"
#include "rcn_structures.h"

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

enum CTUNGHFlags
{
     CTU_UP_FLG    = 1 << 0,
     CTU_LFT_FLG   = 1 << 1,
     CTU_UPLFT_FLG = 1 << 2,
     CTU_UPRGT_FLG = 1 << 3,
};

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

/* FIXME move inter struct definitions somewhere else */
struct CTUBitField{
    uint64_t hfield[33];
    uint64_t vfield[33];
};

struct OVMV
{
    int32_t x;
    int32_t y;
};

typedef struct VVCMergeInfo
{
    OVMV mv0;
    OVMV mv1;
    uint8_t inter_dir;
} VVCMergeInfo;

struct VVCTMVPStorage
{
    uint64_t *col0;
    uint64_t *col1;
    OVMV *mv0;
    OVMV *mv1;
};

struct OVMVCtx
{
    struct CTUBitField map;
    OVMV mvs[34*34];
};

struct HMVPLUT
{
    OVMV hmv0[5];
    OVMV hmv1[5];
    uint8_t dir[5];
    uint8_t nb_mv;
};
struct MVPlane;

struct InterDRVCtx
{
    /* References Pictures Lists */
    OVPicture *rpl0[16];
    OVPicture *rpl1[16];

    /* CTU Local Map Motion Vectors */
    struct OVMVCtx mv_ctx0;
    struct OVMVCtx mv_ctx1;

    /* History based Motion Vector Predicition
     * Look-Up table containing the five last
     * Motion Vectors used.
     */
    struct HMVPLUT hmvp_lut;

    /* Temporal Motion Vector Prediction Related
     * information
     */
    uint8_t tmvp_avail;
    uint8_t tmvp_enabled;

    struct VVCTMVP
    {
        /* FIXME tmp info */
        struct OVCTUDec *ctudec;
        /* FIXME we used Ref Pictures
         * but we only need motion vectors
         * information or ref idx
         */
        /* MV plane storage for current picture */
        struct MVPlane *plane0;
        struct MVPlane *plane1;

        #if 0
        OVPicture *col_ref;
        #endif

        /* MV plane storage for collocated reference picture */
        struct MVPlane *col_plane0;
        struct MVPlane *col_plane1;

        /* Scale info computed at slice start
         * based on the distance between collocated
         * and current picture in POC
         */
        int16_t scale00;
        int16_t scale01;
        int16_t scale10;
        int16_t scale11;

        /* FIXME do not start from 1 */
        /* Vertical bit map of motions available in collocated MV CTU */
        uint64_t dir_map_v0[33];
        uint64_t dir_map_v1[33];

        OVMV mvs0[34*34];
        OVMV mvs1[34*34];

        #if 0
        struct VVCCTBTMVP
        {
            struct OVMVCtx mv_ctx0;
            struct OVMVCtx mv_ctx1;
        }tmvp_mv;
        #endif
    } tmvp_ctx;
};


/* Structure reserved for Intra Modes dervivation
 */
struct IntraDRVInfo
{
    /*
     * FIXME Reset value should be 0 as PLANAR is used
     * when not available
     */

    /* Pointer to above mode line at CTU start */
    /* FIXME
     * This line is not required since mode
     * is considered zero
     */
    uint8_t *luma_mode_x;

    /* Storage for CTU left intra modes
     * reset to PLANAR at each CTU line start
     * FIXME Reduce size for according to partition
     * limits
     */
    uint8_t luma_mode_y[32];

    /* Luma Intra Modes
     * Used in Dual / Separable Tree to derive Chroma Intra Modes
     * FIXME Reduce size according to partition limits
     *       Or keep size constant
     *       Move to Trees Specific info ?
     */
    uint8_t luma_modes[32*32];

};

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

/* FIXME
 * We use enum here since its easier when
 * debugging however a const value
 * would probably be a better option
 */

enum RCNSizes
{
   /* Stride Used in CTU buffers MAX_CTU_S
    *     + 64 samples right used for intra
    *     +  4 samples for intra Multi Ref Lines
    *     + 12 samples for memory alignement purposes
    */
   RCN_CTB_STRIDE  = (128 + 16 + 64),

   /* A padding of 4 upper lines and 16 left
    * columns from buffer start to be used for
    * border copy for intra prediction
    */
   RCN_CTB_PADDING = (RCN_CTB_STRIDE * 4 + 16),

   /* Size of CTB Buffer in samples */
   RCN_CTB_SIZE    = (RCN_CTB_STRIDE * RCN_CTB_STRIDE),
};

struct CTURCNData
{
    uint16_t y_buff[RCN_CTB_SIZE];
    uint16_t cb_buff[RCN_CTB_SIZE];
    uint16_t cr_buff[RCN_CTB_SIZE];

    /* To be used for temporary storage
     * when we needed
     */
    uint16_t tmp_buff[RCN_CTB_SIZE];
};

struct OVCTUDec
{
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
        struct IntraDRVInfo intra_info;
        struct InterDRVCtx  inter_ctx;
        struct DeltaQPDRVCtx *delta_qp_ctx;
        struct FiltersDRVCtx *loop_filters_ctx;
        int8_t qp_map_x[32];
        int8_t qp_map_y[32];
    } drv_ctx;

    /* Reconstruction context */
    struct OVRCNCtx {

        /* FIXME
         * decide where we should store / alloc init this
         * since this should be used by ctudec
         */
        struct CTURCNData data;

        
        /* Pointers to the first sample data of CTU in the current
         * picture
         */
        struct OVBuffInfo{
            uint16_t *y;
            uint16_t *cb;
            uint16_t *cr;
            uint32_t stride;
            uint32_t stride_c;
        } frame_buff;

        /* Pointers to CTU reconstruction buffers to be used
         * to reconstruct current CTU.
         * These buffers will be written to the destination picture
         * before filtering operation
         */
        struct OVBuffInfo ctu_buff;

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
         /* FIXME move to drv */
        struct CTUBitField progress_field;

        struct CTUBitField progress_field_c;

        /* A structure containing various functions pointers
         * to block reconstruction function
         */
        struct RCNFunctions rcn_funcs;
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
    /* FIXME no need for cu_mode in chroma part */
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

    int (*residual_coding_isp_h)(struct OVCTUDec *const lc_ctx, int16_t *const dst,
                                 unsigned int log2_tb_w, unsigned int log2_tb_h,
                                 uint16_t last_pos);

    int (*residual_coding_isp_v)(struct OVCTUDec *const lc_ctx, int16_t *const dst,
                                 unsigned int log2_tb_w, unsigned int log2_tb_h,
                                 uint16_t last_pos);

    int (*residual_coding_ts)(struct OVCTUDec *const lc_ctx,
                              unsigned int log2_tb_w,
                              unsigned int log2_tb_h);

    uint64_t (*residual_coding)(struct OVCTUDec *const lc_ctx, int16_t *const dst,
                                unsigned int log2_tb_w, unsigned int log2_tb_h,
                                uint16_t last_pos);

    int (*residual_coding_chroma)(struct OVCTUDec *const lc_ctx, int16_t *const dst,
                                  unsigned int log2_tb_w, unsigned int log2_tb_h,
                                  uint16_t last_pos);

    /* FIXME
     *    -Reduce residual buff to 32x32
     *    -Might need to store 4 TUs in
     *   case of transform tree or ISP tree
     *    -Move to Transform Unit structure
     */
    int16_t residual_y[64*64];
    int16_t residual_cb[64*64];
    int16_t residual_cr[64*64];
    int16_t lfnst_subblock[16*2];
    int16_t transform_buff[64*64];


    uint8_t slice_qp;
    /* FIXME
     * harmonize qp info and dequant structures
     */
    VVCQPCTX qp_ctx;
    VVCDeQuantCtx dequant_luma;
    VVCDeQuantCtx dequant_luma_skip;
    VVCDeQuantCtx dequant_cb;
    VVCDeQuantCtx dequant_cr;
    VVCDeQuantCtx dequant_joint_cb_cr;

    const VVCDeQuantCtx *dequant_chroma;
    uint16_t ctb_x;
    uint16_t ctb_y;
    uint16_t nb_ctb_pic_w;
};

int ovdec_decode_ctu(OVVCDec *dec, OVCTUDec *ctu_dec);

int ctudec_init(OVCTUDec **ctudec_p);
int ctudec_uninit(OVCTUDec *ctudec_p);

/* TODO add hook functors for PU / TU init and update */

#endif
