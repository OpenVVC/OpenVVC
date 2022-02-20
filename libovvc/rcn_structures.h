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

#ifndef RCN_STRUCTURES_H
#define RCN_STRUCTURES_H
#include <stddef.h>
#include <stdint.h>

#include "rcn_intra_angular.h"
#include "bitdepth.h"

struct ALFClassifier;
struct Area;
struct CCLMParams;
struct SAOParamsCtu;
struct LMCSLUTs;

struct CTUBitField;
struct RCNFunctions;
struct OVLMCSData;
struct LMCSInfo;
struct ISPTUInfo;
struct TUInfo;


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

struct CTUBitField
{
    uint64_t hfield[33];
    uint64_t vfield[33];
};

enum DCTType
{
    DST_VII = 0,
    DCT_VIII = 1,
    DCT_II = 2,
    NB_TR_TYPES = 3
};

#define NB_TR_SIZES 7

typedef void (*MCUniDirFunc)(OVSample *_dst, ptrdiff_t _dststride,
                             const OVSample *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBilinear)(uint16_t *_dst, ptrdiff_t _dststride,
                           const OVSample *_src, ptrdiff_t _srcstride,
                           int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBiDir0Func)(int16_t *_dst,
                             const OVSample *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBiDir1Func)(OVSample *_dst, ptrdiff_t _dststride,
                             const OVSample *_src0, ptrdiff_t _srcstride,
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

typedef void (*MCRPRHor)(uint16_t *_dst, ptrdiff_t _dststride,
                             const OVSample *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width, uint8_t filter_idx);

typedef void (*MCRPRVerUni)(OVSample *_dst, ptrdiff_t _dststride,
                             const uint16_t *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width, uint8_t filter_idx);

typedef void (*MCRPRVerBi)(uint16_t *_dst, ptrdiff_t _dststride,
                             const uint16_t *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width, uint8_t filter_idx);

typedef void (*MCRPRSum)(OVSample* _dst, ptrdiff_t _dststride,
                          const uint16_t* _src0, ptrdiff_t _src0stride,
                          const uint16_t* _src1, ptrdiff_t _src1stride,
                          int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCRPRWeighted)(OVSample* _dst, ptrdiff_t _dststride,
                      const uint16_t* _src, ptrdiff_t _srcstride,
                      const uint16_t* _src2, ptrdiff_t _src2stride,
                      int height, int denom, int wx0, int wx1,
                      intptr_t mx, intptr_t my, int width);

typedef void (*LMsubsampleFunc)(const OVSample *lm_src, OVSample *dst_cb, OVSample *dst_cr,
                                ptrdiff_t lm_src_stride, ptrdiff_t dst_stride_c,
                                const struct CCLMParams *const lm_params,
                                int pb_w, int pb_h, uint8_t lft_avail);

typedef void (*CCLMFunc)( const OVSample* const src_luma, OVSample* const dst_cb,
                          OVSample* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
                          int up_available, int left_available, LMsubsampleFunc const compute_subsample);

typedef void (*MDLMFunc)(const OVSample* const src_luma, OVSample* const dst_cb,
                         OVSample* const dst_cr, uint64_t intra_map_rows,
                         int log2_pb_w, int log2_pb_h, int x0, int y0,
                         uint8_t left_available, uint8_t up_available, LMsubsampleFunc const compute_subsample);

typedef void (*ResidualAddScaleFunc)(const int16_t *src, OVSample *dst,
                                     int log2_tb_w, int log2_tb_h,
                                     int scale);

typedef void (*TrFunc)(const int16_t *src, int16_t *dst,
                 ptrdiff_t src_stride,
                 int num_lines, int num_columns, int shift);

typedef void (*DCFunc)(const OVSample* const src_above,
                 const OVSample* const src_left, OVSample* const dst,
                 ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h);

typedef void (*PlanarFunc)(const OVSample* const src_above,
                     const OVSample* const src_left, OVSample* const dst,
                     ptrdiff_t dst_stride, int log2_pb_w, int log2_pb_h);

typedef void (*LFNSTFunc)(const int16_t* const src, int16_t* const dst,
                     const int8_t* const lfnst_matrix, int log2_tb_w,
                     int log2_tb_h);

typedef void (*MIPUpSample)(OVSample *const dst, const OVSample *const src,
                            const OVSample *ref,
                            int log2_upsampled_size_src, int log2_opposite_size,
                            int src_step, int src_stride,
                            int dst_step, int dst_stride,
                            int ref_step, int log2_scale);

typedef void (*MIPMatMult)(const int16_t *src, OVSample *dst,
                           const uint8_t *matrix, int16_t offset, int rnd,
                           uint8_t log2_src, uint8_t log2_red_w, uint8_t log2_red_h);

typedef void (*ALFClassifBlkFunc)(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr,
                                  OVSample *const src, const int stride, const struct Area blk,
                                  const int shift, const int ctu_height, int virbnd_pos);

typedef void (*ALFFilterBlkFunc)(uint8_t * class_idx_arr, uint8_t * transpose_idx_arr, OVSample *const dst, OVSample *const src, const int dstStride, const int srcStride,
                        struct Area blk_dst, const int16_t *filter_set, const int16_t *clip_set,
                        const int ctu_height, int virbnd_pos);


typedef void (*ALFChromaFilterBlkFunc)(OVSample *const dst, const OVSample *const src,
                                       const int dstStride, const int srcStride,
                                       struct Area blk_dst,
                                       const int16_t *const filter_set, const int16_t *const clip_set,
                                       const int ctu_height, int virbnd_pos);

typedef void (*CCALFFilterBlkFunc)(OVSample * chroma_dst, OVSample * luma_src, const int chr_stride, const int luma_stride,
                        const struct Area blk_dst, const uint8_t c_id, const int16_t *filt_coeff,
                        const int vbCTUHeight, int vbPos);

typedef void (*SAOBandFilterFunc)(OVSample* _dst, OVSample* _src,
                                  ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                  struct SAOParamsCtu* sao, int width,
                                  int height, int c_idx);

typedef void (*SAOEdgeFilterFunc)(OVSample* _dst, OVSample* _src,
                                  ptrdiff_t _stride_dst, ptrdiff_t _stride_src,
                                  struct SAOParamsCtu* sao, int width,
                                  int height, int c_idx);

typedef void (*LMCSReshapeFunc)(OVSample *_dst, ptrdiff_t stride_dst, const struct LMCSLUTs *const luts, int width, int height);

typedef uint64_t (*DMVRSADFunc)(const uint16_t *ref0, const uint16_t *ref1, int16_t dmvr_stride, int16_t pb_w, int16_t pb_h);

typedef uint8_t (*DMVRComputeSADsFunc)(const uint16_t *ref0, const uint16_t *ref1, uint64_t *sad_array, int sb_w, int sb_h);

typedef void (*PROFGradFunction)(const int16_t* src, int src_stride, int sb_w, int sb_h, int grad_stride, int16_t* grad_x, int16_t* grad_y);

typedef void (*PROFFunction)(OVSample* dst, int dst_stride, const int16_t* src, int src_stride,
                             const int16_t* grad_x, const int16_t* grad_y, int grad_stride,
                             const int16_t* dmv_scale_h, const int16_t* dmv_scale_v, uint8_t bidir);

typedef void (*BDOFSBFunction)(const int16_t* src0, int src0_stride,
                               const int16_t* src1, int src1_stride,
                               OVSample *dst, int dst_stride,
                               const int16_t *gradX0, const int16_t *gradX1,
                               const int16_t *gradY0, const int16_t *gradY1, int grad_stride,
                               int wgt_x, int wgt_y);

typedef void (*CIIPWeightedFuntion)(OVSample* dst, int dststride, const OVSample* src_intra,
                                    const OVSample* src_inter, int srcstride, int width, int height, int wt);

typedef void (*DFFilterFunction)(OVSample *src, const int stride, const int tc);


/**
 * The Context put together all functions used by strategies.
 */

struct MCFunctions{

    MCUniDirFunc unidir[4][8];

    MCBiDir0Func bidir0[4][8];
    MCBiDir1Func bidir1[4][8];

    MCUniDirWFunc unidir_w[4][8];
    MCBiDirWFunc bidir_w[4][8];

    MCBilinear    bilinear[4][8];
    MCRPRHor      rpr_h[2][8];
    MCRPRVerUni   rpr_v_uni[2][8];
    MCRPRVerBi    rpr_v_bi[2][8];
    MCRPRSum      rpr_sum;
    MCRPRWeighted rpr_w[8];

    void (*gpm_weighted)(OVSample* _dst, int _dststride, const int16_t* _src0,
                  int srcstride0, const int16_t* _src1, int srcstride1, int height,
                  int width, int step_x, int step_y, int16_t* weight);
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
  DCFunc pdpc[5][5];
};

struct PlanarFunctions
{
  PlanarFunc func;
  PlanarFunc pdpc[5][5];
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

struct OVRCNCtx;
struct MIPFunctions
{
  MIPUpSample upsample_h[2][4];
  MIPUpSample upsample_v[2][4];
  MIPMatMult matmult;

  void (*rcn_intra_mip)(const struct OVRCNCtx *const rcn_ctx,
                        uint8_t x0, uint8_t y0,
                        uint8_t log2_pb_w, uint8_t log2_pb_h,
                        uint8_t mip_opaque);
};

struct RectEntryInfo;
struct OVCTUDec;
struct RCNALF;

struct ALFFunctions{
    ALFClassifBlkFunc classif;
    ALFFilterBlkFunc luma[2];
    ALFChromaFilterBlkFunc chroma[2];
    CCALFFilterBlkFunc ccalf[2];
    void (*rcn_alf_filter_line)(struct OVCTUDec *const ctudec, const struct RectEntryInfo *const einfo, uint16_t ctb_y);

    void (*rcn_alf_reconstruct_coeff_APS)(struct RCNALF* alf, struct OVCTUDec *const ctudec, uint8_t luma_flag, uint8_t chroma_flag);

};

struct OVCTUDec;
struct SAOFunctions{
    SAOBandFilterFunc band;
    SAOEdgeFilterFunc edge[2];

    void (*rcn_sao_filter_line)(struct OVCTUDec *const ctudec,
                                const struct RectEntryInfo *const einfo,
                                uint16_t ctb_y);

    void (*rcn_sao_first_pix_rows)(struct OVCTUDec *const ctudec,
                                   const struct RectEntryInfo *const einfo,
                                   uint16_t ctb_y);
};

struct DMVRFunctions{
    DMVRSADFunc sad[2];
    DMVRComputeSADsFunc computeSB[2];
};

struct PROFFunctions{
    PROFGradFunction grad;
    PROFFunction rcn;
    void (*tmp_prof_mrg)(OVSample* _dst, ptrdiff_t _dststride,
                         const int16_t* _src0, ptrdiff_t _srcstride,
                         const int16_t* _src1, int height, intptr_t mx,
                         intptr_t my, int width);

    void (*tmp_prof_mrg_w)(OVSample* _dst, ptrdiff_t _dststride,
                           const int16_t* _src0, ptrdiff_t _srcstride,
                           const int16_t* _src1, int height, intptr_t mx,
                           intptr_t my, int width, int wt0, int wt1);

    void (*extend_prof_buff)(const OVSample *const src, int16_t *dst_prof, int16_t ref_stride,
                             uint8_t ext_x, uint8_t ext_y);


};

struct BDOFFunctions{
    PROFGradFunction grad;
    BDOFSBFunction subblock;

    void (*rcn_bdof)(struct BDOFFunctions *const bdof, OVSample *dst, int dst_stride,
                     const int16_t *ref_bdof0, const int16_t *ref_bdof1, int ref_stride,
                     const int16_t *grad_x0, const int16_t *grad_y0,
                     const int16_t *grad_x1, const int16_t *grad_y1,
                     int grad_stride, uint8_t pb_w, uint8_t pb_h);

    void (*extend_bdof_buff)(const OVSample *const src, int16_t *dst_prof,
                             int16_t ref_stride, int16_t pb_w, int16_t pb_h,
                             uint8_t ext_x, uint8_t ext_y);

};

struct CIIPFunctions{
    CIIPWeightedFuntion weighted;
};

struct DBFInfo;
struct DFFunctions{
    DFFilterFunction filter_h[11];
    DFFilterFunction filter_v[11];

    void (*rcn_dbf_ctu)(const struct OVRCNCtx  *const rcn_ctx, const struct DBFInfo *const dbf_info,
                        uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y);

    void (*rcn_dbf_truncated_ctu)(const struct OVRCNCtx  *const rcn_ctx, const struct DBFInfo *const dbf_info,
                                  uint8_t log2_ctu_s, uint8_t last_x, uint8_t last_y,
                                  uint8_t ctu_w, uint8_t ctu_h);
};

#include "rcn_dequant.h"
struct TMPBDCompat
{
    void (*filter_ref_samples)(const OVSample* const src, OVSample* const dst,
                               const OVSample* src2, int length);

    void (*fill_ref_left_0)(const OVSample* const src, int src_stride,
                            OVSample* const ref_left, uint64_t intra_map_cols,
                            uint64_t intra_map_rows, int8_t x0, int8_t y0, int log2_pb_w,
                            int log2_pb_h, int offset_y);

    void (*fill_ref_left_0_chroma)(const OVSample* const src, int src_stride,
                                   OVSample* const ref_left, uint64_t intra_map_cols,
                                   uint64_t intra_map_rows, int8_t x0, int8_t y0,
                                   int log2_pb_w, int log2_pb_h);

    void (*fill_ref_left_0_mref)(const OVSample* const src, int src_stride,
                                 OVSample* const ref_left, uint64_t intra_map_cols,
                                 uint64_t intra_map_rows, int mref_idx, int8_t x0,
                                 int8_t y0, int log2_pb_w, int log2_pb_h);

    void (*fill_ref_above_0)(const OVSample* const src, int src_stride,
                             OVSample* const ref_above, uint64_t intra_map_rows,
                             uint64_t intra_map_cols, int8_t x0, int8_t y0, int log2_pb_w,
                             int log2_pb_h, int offset_x);

    void (*fill_ref_above_0_chroma)(const OVSample* const src, int src_stride,
                                    OVSample* const ref_above, uint64_t intra_map_rows,
                                    uint64_t intra_map_cols, int8_t x0, int8_t y0,
                                    int log2_pb_w, int log2_pb_h);

    void (*fill_ref_above_0_mref)(const OVSample* const src, int src_stride,
                                  OVSample* const ref_above, uint64_t intra_map_rows,
                                  uint64_t intra_map_cols, int mref_idx, int8_t x0,
                                  int8_t y0, int log2_pb_w, int log2_pb_h);

    struct IQScale (*derive_dequant_sdh)(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h);

    struct IQScale (*derive_dequant_dpq)(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h);

    struct IQScale (*derive_dequant_ts)(int qp, uint8_t log2_tb_w, uint8_t log2_tb_h);

    void (*dequant_tb_4x4)(int16_t *dst, const int16_t *src, int scale, int shift,
                           uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map);

    void (*dequant_tb_4x4_neg)(int16_t *dst, const int16_t *src, int scale, int shift,
                               uint8_t log2_tb_w, uint8_t log2_tb_h, uint64_t sig_sb_map);

    void (*rcn_transform_tree)(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
                               uint8_t log2_tb_w, uint8_t log2_tb_h, uint8_t log2_max_tb_s,
                               uint8_t tr_depth, uint16_t cu_flags, const struct TUInfo *const tu_info);


    void (*rcn_tu_c)(OVCTUDec *const ctu_dec, uint8_t x0, uint8_t y0,
                     uint8_t log2_tb_w, uint8_t log2_tb_h,
                     uint16_t cu_flags, uint8_t cbf_mask,
                     const struct TUInfo *const tu_info);

    void (*rcn_tu_st)(OVCTUDec *const ctu_dec,
                      uint8_t x0, uint8_t y0,
                      uint8_t log2_tb_w, uint8_t log2_tb_h,
                      uint16_t cu_flags, uint8_t cbf_mask,
                      const struct TUInfo *const tu_info);

    void (*recon_isp_subtree_h)(OVCTUDec *const ctudec,
                                unsigned int x0, unsigned int y0,
                                unsigned int log2_cb_w, unsigned int log2_cb_h,
                                uint8_t intra_mode,
                                const struct ISPTUInfo *const tu_info);

    void (*recon_isp_subtree_v)(OVCTUDec *const ctudec,
                                unsigned int x0, unsigned int y0,
                                unsigned int log2_cb_w, unsigned int log2_cb_h,
                                uint8_t intra_mode,
                                const struct ISPTUInfo *const tu_info);
};

struct OVBuffInfo;

struct InterDRVCtx;
struct PROFInfo;
struct VVCGPM;

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
    LMCSReshapeFunc lmcs_reshape_forward;
    LMCSReshapeFunc lmcs_reshape_backward;

    void (*rcn_lmcs_compute_chroma_scale)(struct LMCSInfo *const lmcs_info,
                                          const struct CTUBitField *const progress_field,
                                          const OVSample *ctu_data_y, uint8_t x0, uint8_t y0);

    void (*rcn_init_lmcs)(struct LMCSInfo *lmcs_info, const struct OVLMCSData *const lmcs_data);


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

    const struct IntraAngularFunctions *intra_angular_gauss_h;
    const struct IntraAngularFunctions *intra_angular_gauss_v;
    const struct IntraAngularFunctions *intra_angular_cubic_h;
    const struct IntraAngularFunctions *intra_angular_cubic_v;
    const struct IntraAngularFunctions *intra_angular_nofrac_h;
    const struct IntraAngularFunctions *intra_angular_nofrac_v;
    const struct IntraAngularFunctions *intra_angular_c_h;
    const struct IntraAngularFunctions *intra_angular_c_v;
    const struct IntraMRLFunctions *intra_mrl;

    struct TMPBDCompat tmp;

    void (*intra_pred)(const struct OVRCNCtx *const rcn_ctx,
                       const struct OVBuffInfo* ctu_buff,
                       uint8_t intra_mode, int x0, int y0,
                       int log2_pb_w, int log2_pb_h, uint16_t cu_flags);

    void (*intra_pred_c)(const struct OVRCNCtx *const rcn_ctx,
                         uint8_t intra_mode, int x0, int y0,
                         int log2_pb_w, int log2_pb_h, uint16_t cu_flags);

    void (*intra_pred_isp)(const OVCTUDec *const ctudec,
                           OVSample *const src,
                           ptrdiff_t dst_stride,
                           uint8_t intra_mode,
                           int x0, int y0,
                           int log2_pb_w, int log2_pb_h,
                           int log2_cb_w, int log2_cb_h,
                           int offset_x, int offset_y);

    void (*intra_pred_mrl)(const OVCTUDec *const ctudec,
                           OVSample *const src,
                           ptrdiff_t dst_stride,
                           uint8_t intra_mode, int x0, int y0,
                           int log2_pb_w, int log2_pb_h,
                           int mrl_idx);

    void (*rcn_update_ctu_border)(struct OVRCNCtx *rcn_ctx, uint8_t log2_ctb_s);

    void (*rcn_write_ctu_to_frame)(const struct OVRCNCtx *const rcn_ctx, uint8_t log2_ctb_s);

    void (*rcn_intra_line_to_ctu)(const struct OVRCNCtx *const rcn_ctx, int x_l, uint8_t log2_ctb_s);

    void (*rcn_ctu_to_intra_line)(const struct OVRCNCtx *const rcn_ctx, int x_l, uint8_t log2_ctu_s);

    void (*rcn_write_ctu_to_frame_border)(const struct OVRCNCtx *const rcn_ctx,
                                          int last_ctu_w, int last_ctu_h);

    void (*rcn_buff_uninit)(struct OVRCNCtx *const rcn_ctx);

    void (*rcn_alloc_intra_line_buff)(struct OVRCNCtx *const rcn_ctx, int nb_ctu_w, uint8_t log2_ctb_s);

    void (*rcn_save_last_cols)(struct OVRCNCtx *const rcn_ctx, int x_pic_l, int y_pic_l, uint8_t is_border_rect);

    void (*rcn_save_last_rows)(struct OVRCNCtx *const rcn_ctx, OVSample** saved_rows, int x_l, int x_pic_l, int y_pic_l, uint8_t is_border_rect);

    void (*rcn_extend_filter_region)(struct OVRCNCtx *const rcn_ctx, OVSample** saved_rows, int x_l,
                                     int x_pic_l, int y_pic_l, uint8_t bnd_msk);

    void (*rcn_alloc_filter_buffers)(struct OVRCNCtx *const rcn_ctx, int nb_ctu_w, int margin, uint8_t log2_ctb_s);

    void (*rcn_attach_ctu_buff)(struct OVRCNCtx *const rcn_ctx);

    void (*rcn_attach_frame_buff)(struct OVRCNCtx *const rcn_ctx, const OVFrame *const f,
                                  const struct RectEntryInfo *const einfo, uint8_t log2_ctb_s);

    void (*rcn_next_buff_line)(struct OVRCNCtx *const rcn_ctx,  uint8_t log2_ctb_s);

    uint8_t (*rcn_dmvr_mv_refine)(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                                  uint8_t x0, uint8_t y0,
                                  uint8_t log2_pu_w, uint8_t log2_pu_h,
                                  OVMV *mv0, OVMV *mv1, uint8_t ref_idx0, uint8_t ref_idx1, uint8_t
                                  apply_bdof);

    void (*rcn_bdof_mcp_l)(OVCTUDec *const ctudec, struct OVBuffInfo dst,
                           uint8_t x0, uint8_t y0, uint8_t log2_pu_w, uint8_t log2_pu_h,
                           OVMV mv0, OVMV mv1, uint8_t ref_idx0, uint8_t ref_idx1);

    void (*rcn_mcp)(OVCTUDec *const ctudec, struct OVBuffInfo dst, int x0, int y0, int log2_pu_w, int log2_pu_h,
                    OVMV mv, uint8_t type, uint8_t ref_idx);

    void (*rcn_mcp_b)(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                      const OVPartInfo *const part_ctx,
                      const OVMV mv0, const OVMV mv1,
                      unsigned int x0, unsigned int y0,
                      unsigned int log2_pb_w, unsigned int log2_pb_h,
                      uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

    void (*rcn_mcp_b_l)(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                        const OVPartInfo *const part_ctx,
                        const OVMV mv0, const OVMV mv1,
                        unsigned int x0, unsigned int y0,
                        unsigned int log2_pb_w, unsigned int log2_pb_h,
                        uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

    void (*rcn_prof_mcp_b_l)(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                             const OVPartInfo *const part_ctx,
                             const OVMV mv0, const OVMV mv1,
                             unsigned int x0, unsigned int y0,
                             unsigned int log2_pb_w, unsigned int log2_pb_h,
                             uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1,
                             uint8_t prof_dir, const struct PROFInfo *const prof_info);

    void (*rcn_mcp_b_c)(OVCTUDec*const lc_ctx, struct OVBuffInfo dst, struct InterDRVCtx *const inter_ctx,
                        const OVPartInfo *const part_ctx,
                        const OVMV mv0, const OVMV mv1,
                        unsigned int x0, unsigned int y0,
                        unsigned int log2_pb_w, unsigned int log2_pb_h,
                        uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

    void (*rcn_ciip_b)(OVCTUDec*const ctudec,
                       const OVMV mv0, const OVMV mv1,
                       unsigned int x0, unsigned int y0,
                       unsigned int log2_pb_w, unsigned int log2_pb_h,
                       uint8_t inter_dir, uint8_t ref_idx0, uint8_t ref_idx1);

    void (*rcn_ciip)(OVCTUDec *const ctudec,
                     int x0, int y0, int log2_pb_w, int log2_pb_h,
                     OVMV mv, uint8_t ref_idx);

    void (*rcn_gpm_b)(OVCTUDec *const ctudec, struct VVCGPM* gpm_ctx,
                      int x0, int y0, int log2_pb_w, int log2_pb_h);
};


#endif
