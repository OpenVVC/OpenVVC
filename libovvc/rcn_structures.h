#ifndef RCN_STRUCTURES_H
#define RCN_STRUCTURES_H
#include <stddef.h>
#include <stdint.h>

#include <stdint.h>
#include <stddef.h>

struct ALFClassifier;
struct Area;
struct CCLMParams;
struct SAOParamsCtu;
struct LMCSInfo;

enum DCTType
{
    DST_VII = 0,
    DCT_VIII = 1,
    DCT_II = 2,
    NB_TR_TYPES = 3
};

#define NB_TR_SIZES 7

typedef void (*MCUniDirFunc)(uint16_t *_dst, ptrdiff_t _dststride,
                             const uint16_t *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBiDir0Func)(int16_t *_dst,
                             const uint16_t *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBiDir1Func)(uint16_t *_dst, ptrdiff_t _dststride,
                             const uint16_t *_src0, ptrdiff_t _srcstride,
                             const int16_t *_src1,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCUniDirWFunc)(uint8_t* dst, ptrdiff_t dststride, uint8_t* _src,
                              ptrdiff_t _srcstride, int height, int denom,
                              int wx, int ox, intptr_t mx, intptr_t my,
                              int width);
typedef void (*MCBiDirWFunc)(uint8_t* dst, ptrdiff_t dststride, uint8_t* _src,
                             ptrdiff_t _srcstride, int16_t* src2,
                             ptrdiff_t src2stride, int height, int denom,
                             int wx0, int wx1, intptr_t mx,
                             intptr_t my, int width);

typedef void (*LMsubsampleFunc)(const uint16_t *lm_src, uint16_t *dst_cb, uint16_t *dst_cr,
                               ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                               const struct CCLMParams *const lm_params,
                               int pb_w, int pb_h, uint8_t lft_avail);

typedef void (*CCLMFunc)( const uint16_t* const src_luma, uint16_t* const dst_cb,
                          uint16_t* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
                          int up_available, int left_available, LMsubsampleFunc const *compute_subsample);

typedef void (*MDLMFunc)(const uint16_t* const src_luma, uint16_t* const dst_cb,
                         uint16_t* const dst_cr, uint64_t intra_map_rows,
                         int log2_pb_w, int log2_pb_h, int x0, int y0,
                         uint8_t left_available, uint8_t up_available, LMsubsampleFunc const *compute_subsample);

typedef void (*ResidualAddScaleFunc)(const int16_t *src, uint16_t *dst,
                                     int log2_tb_w, int log2_tb_h,
                                     int scale);

typedef void (*TrFunc)(const int16_t *src, int16_t *dst,
                 ptrdiff_t src_stride,
                 int num_lines, int num_columns, int shift);

typedef void (*DCFunc)(const uint16_t* const src_above,
                 const uint16_t* const src_left, uint16_t* const dst,
                 ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h);

typedef void (*PlanarFunc)(const uint16_t* const src_above,
                     const uint16_t* const src_left, uint16_t* const dst,
                     ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h);

typedef void (*LFNSTFunc)(const int16_t* const src, int16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h);

typedef void (*MIPUpSample)(uint16_t *const dst, const int16_t *const src,
                            const uint16_t *ref,
                            int log2_upsampled_size_src, int log2_opposite_size,
                            int src_step, int src_stride,
                            int dst_step, int dst_stride,
                            int ref_step, int log2_scale);

typedef void (*MIPMatMult)(int16_t * src, uint16_t *dst, const int stride,
                           const uint8_t *matrix, int16_t offset, const int rnd,
                           uint8_t log2_src, uint8_t log2_red_w, uint8_t log2_red_h);

typedef void (*ALFClassifBlkFunc)(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr,
                                  int16_t *const src, const int stride, const struct Area blk,
                                  const int shift, const int ctu_height, int virbnd_pos);

typedef void (*ALFFilterBlkFunc)(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr, int16_t *const dst, int16_t *const src, const int dstStride, const int srcStride,
                        struct Area blk_dst, const int16_t *filter_set, const int16_t *clip_set,
                        const int ctu_height, int virbnd_pos);


typedef void (*ALFChromaFilterBlkFunc)(int16_t *const dst, const int16_t *const src,
                                       const int dstStride, const int srcStride,
                                       struct Area blk_dst,
                                       const int16_t *const filter_set, const int16_t *const clip_set,
                                       const int ctu_height, int virbnd_pos);

typedef void (*CCALFFilterBlkFunc)(int16_t * chroma_dst, int16_t * luma_src, const int chr_stride, const int luma_stride,
                        const struct Area blk_dst, const uint8_t c_id, const int16_t *filt_coeff,
                        const int vbCTUHeight, int vbPos);

typedef void (*SAOBandFilterFunc)(uint8_t* _dst, uint8_t* _src,
                                  ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                  struct SAOParamsCtu* sao, int width,
                                  int height, int c_idx);

typedef void (*SAOEdgeFilterFunc)(uint8_t* _dst, uint8_t* _src,
                                  ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                  struct SAOParamsCtu* sao, int width,
                                  int height, int c_idx);

typedef void (*LMCSReshapeFunc)(uint16_t *_dst, ptrdiff_t stride_dst, uint16_t* lmcs_lut_luma, int width, int height);

typedef uint64_t (*DMVRSADFunc)(const int16_t *ref0, const int16_t *ref1, int16_t dmvr_stride, int16_t pb_w, int16_t pb_h);

typedef uint8_t (*DMVRComputeSADsFunc)(const int16_t *ref0, const int16_t *ref1, uint64_t *sad_array, int sb_w, int sb_h);

typedef void (*PROFGradFunction)(const uint16_t* src, int src_stride, int sb_w, int sb_h, int grad_stride, int16_t* grad_x, int16_t* grad_y);

typedef void (*PROFFunction)(uint16_t* dst, int dst_stride, const uint16_t* src, int src_stride,
         const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
         const int32_t* dmv_scale_h, const int32_t* dmv_scale_v, uint8_t bidir);

typedef void (*BDOFSBFunction)(const int16_t* src0, int src0_stride,
                        const int16_t* src1, int src1_stride,
                        int16_t *dst, int dst_stride,
                        const int16_t *gradX0, const int16_t *gradX1,
                        const int16_t *gradY0, const int16_t *gradY1, int grad_stride,
                        int wgt_x, int wgt_y);

typedef void (*CIIPWeightedFuntion)(uint16_t* dst, int dststride, const uint16_t* src_intra,
                                    const uint16_t* src_inter, int srcstride, int width, int height, int wt);

typedef void (*DFFilterFunction)(int16_t *src, const int stride, const int tc);

/**
 * The Context put together all functions used by strategies.
 */

struct MCFunctions{

    MCUniDirFunc unidir[4][8];

    MCBiDir0Func bidir0[4][8];
    MCBiDir1Func bidir1[4][8];

    MCUniDirWFunc unidir_w[4][8];
    MCBiDirWFunc bidir_w[4][8];

    MCUniDirFunc bilinear[4][8];
};

struct CCLMFunctions
{
    CCLMFunc cclm;
    MDLMFunc mdlm_left;
    MDLMFunc mdlm_top;
    LMsubsampleFunc compute_subsample;
};

struct TRFunctions
{
   TrFunc func[NB_TR_TYPES][NB_TR_SIZES];
   void (*dc)(int16_t* const dst, int log2_tb_w, int log2_tb_h, int dc_val);
};

struct DCFunctions
{
  DCFunc func;
  DCFunc pdpc;
};

struct PlanarFunctions
{
  PlanarFunc func;
  PlanarFunc pdpc[2];
};

struct ICTFunctions
{
  ResidualAddScaleFunc add[7];
  ResidualAddScaleFunc ict[6][3];
};

struct LFNSTFunctions
{
  LFNSTFunc func[2][2];
};

struct MIPFunctions
{
  MIPUpSample upsample_h[2][3];
  MIPUpSample upsample_v[2][3];
  MIPMatMult matmult;
};

struct ALFFunctions{
  ALFClassifBlkFunc classif;
  ALFFilterBlkFunc luma[2];
  ALFChromaFilterBlkFunc chroma[2];
  CCALFFilterBlkFunc ccalf[2];
};

struct SAOFunctions{
    SAOBandFilterFunc band;
    SAOEdgeFilterFunc edge[2];
};

struct DMVRFunctions{
    DMVRSADFunc sad[2];
    DMVRComputeSADsFunc computeSB[2];
};

struct PROFFunctions{
    PROFGradFunction grad;
    PROFFunction rcn;
};

struct BDOFFunctions{
    PROFGradFunction grad;
    BDOFSBFunction subblock;
};

struct CIIPFunctions{
    CIIPWeightedFuntion weighted;
};

struct DFFunctions{
    DFFilterFunction filter_h[11];
    DFFilterFunction filter_v[11];
};

struct RCNFunctions
{
    /* Motion Compensation Luma */
    struct MCFunctions mc_l;

    /* Motion Compensation Chroma */
    struct MCFunctions mc_c;

    struct CCLMFunctions cclm;

    struct ICTFunctions ict;

    /* Transform Functions */
    struct TRFunctions tr;

    /* LFNST Functions */
    struct LFNSTFunctions lfnst;

    /* DC Functions */
    struct DCFunctions dc;

    /* Planar Functions */
    struct PlanarFunctions planar;

    /* MIP Functions */
    struct MIPFunctions mip;

    /* ALF Functions */
    struct ALFFunctions alf;

    /* SAO Functions */
    struct SAOFunctions sao;

    /* LMCS Functions */
    LMCSReshapeFunc lmcs_reshape;

    /* DMVR Functions */
    struct DMVRFunctions dmvr;

    /* PROF Functions */
    struct PROFFunctions prof;

    /* BDOF Functions */
    struct BDOFFunctions bdof;

    /* CIIP Functions */
    struct CIIPFunctions ciip;

    /* DF Functions */
    struct DFFunctions df;
};


#endif
