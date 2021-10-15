#ifndef CTU_DEC_H
#define CTU_DEC_H

#include "ovdefs.h"
#include "ovframe.h"
#include "ovmem.h"
#include "ovdec.h"

#include "nvcl_structures.h"
#include "rcn_structures.h"
#include "rcn_alf.h"
#include "dec_structures.h"

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

enum VVCCUFlag
{
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

struct TUInfo;

/* FIXME: Move here for SSE */
struct LMParams
{
   int shift;
   int a;
   int b;
};

/* FIXME: Move here for SSE */
struct CCLMParams
{
   struct LMParams cb;
   struct LMParams cr;
};

struct DBFMap
{
    uint64_t ver[33]; /* Column map of vertical edges */
    uint64_t hor[33]; /* Row map of horizontal edges */
};

struct DBFQPMap
{
    uint8_t hor[34*33];
};

struct DBFInfo
{
    /*FIXME replace with DBFMap structure */
    /* FIXME this overlap with edge maps */
    uint64_t ctb_bound_ver[16 + 33]; /* Column map of vertical edges */
    uint64_t ctb_bound_hor[16 + 33]; /* Row map of horizontal edges */
    uint64_t ctb_bound_ver_c[16 + 33]; /* Column map of vertical edges */
    uint64_t ctb_bound_hor_c[16 + 33]; /* Row map of horizontal edges */

    /* FIXME reduce since internal */
    uint64_t aff_edg_ver[16 + 33]; /* Column map of affine internal vertical edges */
    uint64_t aff_edg_hor[16 + 33]; /* Row map of affine internal horizontal edges */


    /* FIXME this overlap with edge maps */
    struct DBFMap bs2_map;
    struct DBFMap bs2_map_c;

    struct DBFMap bs1_map;
    struct DBFMap bs1_map_cb;
    struct DBFMap bs1_map_cr;
    struct DBFMap affine_map;

    /* FIXME this overlap with edge maps */

    int16_t beta_offset;
    int16_t tc_offset;

    /*FIXME reduce those maps deriving qp in dbf_function 
     * Those tables are not mandatory if delta qp is disabled
     */
    struct DBFQPMap qp_map_y;
    struct DBFQPMap qp_map_cb;
    struct DBFQPMap qp_map_cr;
    uint8_t disable_v;
    uint8_t disable_h;
};


struct SAOInfo
{
    /* FIXME move flags to ctudec */
    uint8_t sao_luma_flag;
    uint8_t sao_chroma_flag;
    /* FIXME remove */
    uint8_t chroma_format_idc;

    /*array of SAO parameters structure for each ctu */ 
    SAOParamsCtu *sao_params;
};

struct ALFInfo
{
    /* FIXME move flags to ctu dec */
    uint8_t alf_luma_enabled_flag;
    uint8_t alf_cb_enabled_flag;
    uint8_t alf_cr_enabled_flag;
    uint8_t cc_alf_cb_enabled_flag;
    uint8_t cc_alf_cr_enabled_flag;

    uint8_t num_alf_aps_ids_luma;
    uint8_t left_ctb_alf_flag;

    //TODO: use width in ctu of image (or tile) instead of max value 32.
    uint8_t ctb_cc_alf_flag_line[2][256];

    const struct OVALFData* aps_alf_data[8];
    const struct OVALFData* aps_alf_data_c;
    const struct OVALFData* aps_cc_alf_data_cb;
    const struct OVALFData* aps_cc_alf_data_cr;

    uint8_t left_ctb_cc_alf_flag[2];
    //TODO: use width in ctu of image (or tile) instead of max value 32.
    uint8_t ctb_alf_flag_line[256];
    uint8_t* ctb_cc_alf_filter_idx[2];

    //TODO:  better here or in slicedec ?
     /* arrays of ALF parameters structure for each ctu*/
    ALFParamsCtu *ctb_alf_params;

    //ALF reconstruction structure
    RCNALF rcn_alf;
};

struct LMCSInfo
{
    /* FIXME use new LMCS structures */
    uint8_t  lmcs_enabled_flag;
    uint8_t  scale_c_flag;
    uint16_t lmcs_output_pivot[PIC_CODE_CW_BINS];
    uint16_t lmcs_chroma_scale;
    int16_t  lmcs_chroma_scaling_offset;
    uint8_t min_idx;
    uint8_t max_idx;
    uint16_t* lmcs_lut_inv_luma;
    uint16_t* lmcs_lut_fwd_luma;
    const OVLMCSData *data;
};


typedef struct VVCDeQuantCtx{
    uint8_t qp;
}VVCDeQuantCtx;

typedef struct VVCQPCTX
{
    /* FIXME move current qp outside */
    int8_t current_qp;
    int8_t min_qp_prime_ts; /* FIXME usefull ?*/
    int8_t cb_offset;
    int8_t cr_offset;
    int8_t jcbcr_offset;
    const int8_t *chroma_qp_map_cb;
    const int8_t *chroma_qp_map_cr;
    const int8_t *chroma_qp_map_jcbcr;
} VVCQPCTX;

struct VVCCU
{

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
    uint8_t cu_mode_idx_c;
    uint8_t cu_opaque;
    uint8_t cu_mode_info;

};

/* FIXME move inter struct definitions somewhere else */
struct CTUBitField {
    uint64_t hfield[33];
    uint64_t vfield[33];
};

struct OVMV
{
    /*FIXME move ref_idx, bcw_idx and prec_amvr outside of struct */
    int32_t x;
    int32_t y;
    int8_t ref_idx;
    uint8_t bcw_idx_plus1;
    uint8_t prec_amvr;
};

typedef struct VVCMergeInfo
{
    OVMV mv0;
    OVMV mv1;
    uint8_t inter_dir;
} VVCMergeInfo;

struct OVMVCtx
{
    struct CTUBitField map;
    OVMV mvs[34*34];
};

struct AffineControlInfo
{
    OVMV lt;
    OVMV rt;
    OVMV lb;
};

struct PBInfo
{
    uint8_t x_pb;
    uint8_t y_pb;
    uint8_t log2_w;
    uint8_t log2_h;
    uint8_t nb_pb_w;
    uint8_t nb_pb_h;
};

struct AffineInfo {
    struct AffineControlInfo cps[2];
    struct PBInfo pb;
    uint8_t type;
};

struct AffineDRVInfo {
    struct CTUBitField map;
    struct AffineInfo affine_info[34*34];
};

struct HMVPLUT
{
    OVMV hmv0[5];
    OVMV hmv1[5];
    uint8_t dir[5];
    uint8_t nb_mv;
};

struct MVPlane;

struct MVCTU {
    OVMV mvs[16*16];
};

struct InterDRVCtx
{
    /* References Pictures Lists */
    OVPicture *rpl0[16];
    OVPicture *rpl1[16];

    uint8_t nb_active_ref0;
    uint8_t nb_active_ref1;

    struct RPLInfo* rpl_info0;
    struct RPLInfo* rpl_info1;

    int16_t dist_ref_0[16];
    int16_t dist_ref_1[16];

    /* CTU Local Map Motion Vectors */
    struct OVMVCtx mv_ctx0;
    struct OVMVCtx mv_ctx1;
    struct MVCTU tmvp_mv[2];

    /* History based Motion Vector Predicition
     * Look-Up table containing the five last
     * Motion Vectors used.
     */
    struct HMVPLUT hmvp_lut;

    /* Symmetric MVD (SMVD) Related
     * information
     */
    uint8_t bi_dir_pred_flag;
    int ref_smvd_idx0;
    int ref_smvd_idx1;

    /* Combined Inter Intra Pred(CIIP) Related
     * information
     */
    uint8_t ciip_flag;

    /* Merge mode with MVD(MMVD) Related
     * information
     */
    uint8_t mmvd_flag;
    uint8_t mmvd_shift;

    /* Bi-prediction CU-level Weight (BCW) Related
     * information
     */
    uint8_t bcw_flag;

    /* Geometric Partition Mode(GPM) Related
     * information
     */
    uint8_t gpm_flag;
    int max_gpm_cand;

    /* Adaptive MV Refinement(AMVR) Related
     * information
     */
    uint8_t amvr_flag;
    uint8_t affine_amvr_flag;
    uint8_t affine_6params_enabled;
    uint8_t prec_amvr;

    struct VVCGPM
    {
        int split_dir;

        int merge_idx0;
        int merge_idx1;
        int inter_dir0;
        int inter_dir1;
        OVMV mv0;
        OVMV mv1;
    } gpm_ctx;

    uint8_t mvd1_zero_flag;

    uint8_t tmvp_enabled;
    uint8_t sbtmvp_enabled;

    uint8_t prof_enabled;

    uint8_t tmvp_avail;
    struct VVCTMVP
    {
        /* FIXME tmp info */
        struct OVCTUDec *ctudec;

        /* MV plane storage for current picture */
        const struct MVPlane *plane0;
        const struct MVPlane *plane1;

        struct ColInfo {
            int8_t ref_idx_rpl0;
            int8_t ref_idx_rpl1;
        }col_info;

        const OVPicture *col_ref;

        /* MV plane storage for collocated reference picture */
        const struct MVPlane *col_plane0;
        const struct MVPlane *col_plane1;
        uint8_t col_ref_l0;
        uint8_t ldc;

        /* Scale info computed at slice start
         * based on the distance between collocated
         * and current picture in POC
         */
        int16_t scale00;
        int16_t scale01;
        int16_t scale10;
        int16_t scale11;

        int16_t dist_col_0[16];
        int16_t dist_col_1[16];

        /* FIXME do not start from 1 */
        /* Vertical bit map of motions available in collocated MV CTU */
        uint64_t dir_map_v0[34];
        uint64_t dir_map_v1[34];

        OVMV mvs0[16*17];
        OVMV mvs1[16*17];

        uint8_t ctu_w;
        uint8_t ctu_h;

    } tmvp_ctx;

    struct AffineDRVInfo affine_ctx;
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
struct TrCoeffData
{
    int16_t residual_y[128*128];
    int16_t residual_cb[128*128];
    int16_t residual_cr[128*128];
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
    DECLARE_ALIGNED(32, uint16_t, y_buff)[RCN_CTB_SIZE];
    DECLARE_ALIGNED(32, uint16_t, cb_buff)[RCN_CTB_SIZE];
    DECLARE_ALIGNED(32, uint16_t, cr_buff)[RCN_CTB_SIZE];

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
    struct OVRCNCtx
    {
        /*FIXME tmp*/
        struct OVCTUDec *ctudec;

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

        /*Pointers to intra line reconstruction buffers*/
        struct OVBuffInfo intra_line_buff;

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


    struct DBFInfo dbf_info;
    
    struct SAOInfo sao_info;

    struct ALFInfo alf_info;

    struct LMCSInfo lmcs_info;


    struct OVFilterBuffers{
        int16_t* filter_region[3];
        int16_t  filter_region_h[3];
        int16_t  filter_region_w[3];
        int16_t  filter_region_stride[3];
        int16_t  filter_region_offset[3];

        int16_t* saved_rows_sao[3];
        int16_t* saved_rows_alf[3];
        int16_t* saved_cols[3];
        int16_t  saved_rows_stride[3];

        uint8_t  margin;

        //TODO: change alf/sao to use ctudec buffer instead of frame buffer.
        struct Frame* pic_frame; 
    }filter_buffers;

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
    uint8_t mts_implicit;
    uint8_t mts_enabled;
    uint8_t mts_explicit_intra;
    uint8_t mts_explicit_inter;
    uint8_t enable_lfnst;

    uint8_t delta_qp_enabled;

    uint8_t isp_enabled;
    uint8_t enabled_mip;
    uint8_t enable_mrl;

    uint8_t lm_chroma_enabled;
    uint8_t enable_cclm;
    uint8_t max_num_merge_candidates;
    uint8_t sbt_enabled;

    uint8_t affine_enabled;
    uint8_t affine_status;
    uint8_t affine_nb_merge_cand;

    uint8_t bdof_enabled;
    uint8_t dmvr_enabled;

    uint8_t dbf_disable;

    /* Separate chroma tree */
    uint8_t share;

    uint8_t tmp_disable_cclm;

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
                          uint8_t cbf_ctx, uint8_t cu_flags,
                          uint8_t tr_depth, struct TUInfo *const tu_info);

    int (*prediction_unit)(struct OVCTUDec *const lc_ctx,
                           const OVPartInfo *const part_ctx,
                           uint8_t x0, uint8_t y0,
                           uint8_t log2_pb_w, uint8_t log2_pb_h,
                           uint8_t skip_flag, uint8_t cu_merge_flag);

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
    DECLARE_ALIGNED(32, int16_t, residual_y)[128*128];
    DECLARE_ALIGNED(32, int16_t, residual_cb)[128*128];
    DECLARE_ALIGNED(32, int16_t, residual_cr)[128*128];
    int16_t lfnst_subblock[16*2];
    DECLARE_ALIGNED(32, int16_t, transform_buff)[64*64];

    int8_t slice_qp;
    /* FIXME
     * harmonize qp info and dequant structures
     */
    VVCQPCTX qp_ctx;
    VVCDeQuantCtx dequant_luma;
    VVCDeQuantCtx dequant_luma_skip;
    VVCDeQuantCtx dequant_cb;
    VVCDeQuantCtx dequant_cr;
    VVCDeQuantCtx dequant_joint_cb_cr;
    VVCDeQuantCtx dequant_cb_skip;
    VVCDeQuantCtx dequant_cr_skip;
    VVCDeQuantCtx dequant_jcbcr_skip;

    const VVCDeQuantCtx *dequant_chroma;
    const VVCDeQuantCtx *dequant_skip;
    uint16_t ctb_x;
    uint16_t ctb_y;
    uint16_t nb_ctb_pic_w;
    uint16_t prev_nb_ctu_w_rect_entry;
    
    //image height and width in luma samples
    uint16_t pic_h;
    uint16_t pic_w;
    uint16_t cur_poc;
    uint8_t intra_mode_c;
    /* FIXME to be removed */
    uint8_t tmp_ciip;
    uint8_t tmp_red;
};

int ovdec_decode_ctu(OVVCDec *dec, OVCTUDec *ctu_dec);

void ctudec_alloc_filter_buffers(OVCTUDec *const ctudec, int nb_ctu_w, int margin);
void ctudec_extend_filter_region(OVCTUDec *const ctudec, int16_t** saved_rows, int x_l, int x_pic_l, int y_pic_l, uint8_t is_border_rect);
void ctudec_save_last_rows(OVCTUDec *const ctudec, int16_t** saved_rows, int x_l, int x_pic_l, int y_pic_l, uint8_t is_border_rect);
void ctudec_save_last_cols(OVCTUDec *const ctudec, int x_pic_l, int y_pic_l, uint8_t is_border_rect);
void ctudec_free_filter_buffers(OVCTUDec *const ctudec);

void ctudec_alloc_intra_line_buff(OVCTUDec *const ctudec, int nb_ctu_w);
void ctudec_free_intra_line_buff(OVCTUDec *const ctudec);

int ctudec_init(OVCTUDec **ctudec_p);
int ctudec_uninit(OVCTUDec *ctudec_p);

/* TODO add hook functors for PU / TU init and update */

#endif
