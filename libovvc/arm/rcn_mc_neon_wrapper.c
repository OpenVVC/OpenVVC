#include "arm/rcn_neon.h"

#define SIZE_BLOCK_4 1
#define SIZE_BLOCK_8 2
#define SIZE_BLOCK_16 3
#define SIZE_BLOCK_32 4
#define SIZE_BLOCK_64 5
#define SIZE_BLOCK_128 6


const int8_t ov_mc_filters_neon[16][8] =
{
    {   0, 1,  -3, 63,  4,  -2,  1,  0 },
    {  -1, 2,  -5, 62,  8,  -3,  1,  0 },
    {  -1, 3,  -8, 60, 13,  -4,  1,  0 },
    {  -1, 4, -10, 58, 17,  -5,  1,  0 },
    {  -1, 4, -11, 52, 26,  -8,  3, -1 },
    {  -1, 3,  -9, 47, 31, -10,  4, -1 },
    {  -1, 4, -11, 45, 34, -10,  4, -1 },
    {  -1, 4, -11, 40, 40, -11,  4, -1 },
    {  -1, 4, -10, 34, 45, -11,  4, -1 },
    {  -1, 4, -10, 31, 47,  -9,  3, -1 },
    {  -1, 3,  -8, 26, 52, -11,  4, -1 },
    {   0, 1,  -5, 17, 58, -10,  4, -1 },
    {   0, 1,  -4, 13, 60,  -8,  3, -1 },
    {   0, 1,  -3,  8, 62,  -5,  2, -1 },
    {   0, 1,  -2,  4, 63,  -3,  1,  0 },

    //Hpel for amvr
    {  0, 3, 9, 20, 20, 9, 3, 0 }
};
void sink(){};

void ov_put_vvc_bi0_pel_pixels_10_4_neon();
void ov_put_vvc_bi0_pel_pixels_10_8_neon();
void ov_put_vvc_bi0_pel_pixels_10_16_neon();
void ov_put_vvc_bi0_pel_pixels_10_32_neon();

void ov_put_vvc_bi1_pel_pixels_10_4_neon();//untested
void ov_put_vvc_bi1_pel_pixels_10_8_neon();
void ov_put_vvc_bi1_pel_pixels_10_16_neon();
void ov_put_vvc_bi1_pel_pixels_10_32_neon();

void ov_put_vvc_uni_qpel_h_10_8_neon();
void ov_put_vvc_uni_qpel_h_10_16_neon();
void ov_put_vvc_uni_qpel_h_10_32_neon();

void ov_put_vvc_bi0_qpel_h_10_8_neon();
void ov_put_vvc_bi0_qpel_h_10_16_neon();
void ov_put_vvc_bi0_qpel_h_10_32_neon();

void ov_put_vvc_bi1_qpel_h_10_8_neon();
void ov_put_vvc_bi1_qpel_h_10_16_neon();
void ov_put_vvc_bi1_qpel_h_10_32_neon();

void ov_put_vvc_uni_qpel_h_10_8_neon_wrapper(uint16_t* _dst, ptrdiff_t _dststride,
                   const uint16_t* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width){
                     ov_put_vvc_uni_qpel_h_10_8_neon(_dst, _dststride, _src, _srcstride, height, ov_mc_filters_neon[mx-1], width);
                   }
void ov_put_vvc_uni_qpel_h_10_16_neon_wrapper(uint16_t* _dst, ptrdiff_t _dststride,
                   const uint16_t* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width){
                     ov_put_vvc_uni_qpel_h_10_16_neon(_dst, _dststride, _src, _srcstride, height, ov_mc_filters_neon[mx-1], width);
                   }
void ov_put_vvc_uni_qpel_h_10_32_neon_wrapper(uint16_t* _dst, ptrdiff_t _dststride,
                   const uint16_t* _src, ptrdiff_t _srcstride,
                   int height, intptr_t mx, intptr_t my, int width){
                     ov_put_vvc_uni_qpel_h_10_32_neon(_dst, _dststride, _src, _srcstride, height, ov_mc_filters_neon[mx-1], width);
                   }

void ov_put_vvc_bi0_qpel_h_10_8_neon_wrapper(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                  int height, intptr_t mx, intptr_t my, int width){
                     ov_put_vvc_bi0_qpel_h_10_8_neon(_dst, _src, _srcstride, height, ov_mc_filters_neon[mx-1], width);
                   }
void ov_put_vvc_bi0_qpel_h_10_16_neon_wrapper(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                  int height, intptr_t mx, intptr_t my, int width){
                    ov_put_vvc_bi0_qpel_h_10_16_neon(_dst, _src, _srcstride, height, ov_mc_filters_neon[mx-1], width);
                  }
void ov_put_vvc_bi0_qpel_h_10_32_neon_wrapper(int16_t* _dst, const uint16_t* _src, ptrdiff_t _srcstride,
                  int height, intptr_t mx, intptr_t my, int width){
                    ov_put_vvc_bi0_qpel_h_10_32_neon(_dst, _src, _srcstride, height, ov_mc_filters_neon[mx-1], width);
                  }

 void ov_put_vvc_bi1_qpel_h_10_8_neon_wrapper(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width){
                      ov_put_vvc_bi1_qpel_h_10_8_neon(_dst, _dststride, _src0, _srcstride, _src1, height, ov_mc_filters_neon[mx-1], width);
                    }
 void ov_put_vvc_bi1_qpel_h_10_16_neon_wrapper(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width){
                      ov_put_vvc_bi1_qpel_h_10_16_neon(_dst, _dststride, _src0, _srcstride, _src1, height, ov_mc_filters_neon[mx-1], width);
                    }
 void ov_put_vvc_bi1_qpel_h_10_32_neon_wrapper(uint16_t* _dst, ptrdiff_t _dststride, const uint16_t* _src0,
                   ptrdiff_t _srcstride, const int16_t* _src1, int height,
                   intptr_t mx, intptr_t my, int width){
                      ov_put_vvc_bi1_qpel_h_10_32_neon(_dst, _dststride, _src0, _srcstride, _src1, height, ov_mc_filters_neon[mx-1], width);
                    }

void
rcn_init_mc_functions_neon(struct RCNFunctions* const rcn_funcs)
{
  struct MCFunctions* const mc_l = &rcn_funcs->mc_l;
  struct MCFunctions* const mc_c = &rcn_funcs->mc_c;

  /* Luma functions */

  mc_l->bidir0[0][SIZE_BLOCK_8] = &ov_put_vvc_bi0_pel_pixels_10_8_neon;
  mc_l->bidir1[0][SIZE_BLOCK_8] = &ov_put_vvc_bi1_pel_pixels_10_8_neon;
  //
  mc_l->unidir[1][SIZE_BLOCK_8] = &ov_put_vvc_uni_qpel_h_10_8_neon_wrapper;
  mc_l->bidir0[1][SIZE_BLOCK_8] = &ov_put_vvc_bi0_qpel_h_10_8_neon_wrapper;
  mc_l->bidir1[1][SIZE_BLOCK_8] = &ov_put_vvc_bi1_qpel_h_10_8_neon_wrapper;
  //
  // mc_l->unidir[2][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir0[2][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir1[2][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_l->unidir[3][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir0[3][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir1[3][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  mc_l->bidir0[0][SIZE_BLOCK_16] = &ov_put_vvc_bi0_pel_pixels_10_16_neon;
  mc_l->bidir1[0][SIZE_BLOCK_16] = &ov_put_vvc_bi1_pel_pixels_10_16_neon;
  //
  mc_l->unidir[1][SIZE_BLOCK_16] = &ov_put_vvc_uni_qpel_h_10_16_neon_wrapper;
  mc_l->bidir0[1][SIZE_BLOCK_16] = &ov_put_vvc_bi0_qpel_h_10_16_neon_wrapper;
  mc_l->bidir1[1][SIZE_BLOCK_16] = &ov_put_vvc_bi1_qpel_h_10_16_neon_wrapper;
  //
  // mc_l->unidir[2][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir0[2][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir1[2][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_l->unidir[3][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir0[3][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir1[3][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  mc_l->bidir0[0][SIZE_BLOCK_32] = &ov_put_vvc_bi0_pel_pixels_10_32_neon;
  mc_l->bidir1[0][SIZE_BLOCK_32] = &ov_put_vvc_bi1_pel_pixels_10_32_neon;
  //
  mc_l->unidir[1][SIZE_BLOCK_32] = &ov_put_vvc_uni_qpel_h_10_32_neon_wrapper;
  mc_l->bidir0[1][SIZE_BLOCK_32] = &ov_put_vvc_bi0_qpel_h_10_32_neon_wrapper;
  mc_l->bidir1[1][SIZE_BLOCK_32] = &ov_put_vvc_bi1_qpel_h_10_32_neon_wrapper;
  //
  // mc_l->unidir[2][SIZE_BLOCK_32] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir0[2][SIZE_BLOCK_32] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir1[2][SIZE_BLOCK_32] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_l->unidir[3][SIZE_BLOCK_32] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir0[3][SIZE_BLOCK_32] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_l->bidir1[3][SIZE_BLOCK_32] = &dav1d_ff_hevc_transform_4x4_neon_8;

  /* Chroma functions */
  mc_c->bidir0[0][SIZE_BLOCK_4] = &ov_put_vvc_bi0_pel_pixels_10_4_neon;
  mc_c->bidir1[0][SIZE_BLOCK_4] = &ov_put_vvc_bi1_pel_pixels_10_4_neon;
  //
  // mc_c->unidir[1][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[1][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[1][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_c->unidir[2][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[2][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[2][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_c->unidir[3][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[3][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[3][SIZE_BLOCK_4] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  mc_c->bidir0[0][SIZE_BLOCK_8] = &ov_put_vvc_bi0_pel_pixels_10_8_neon;
  mc_c->bidir1[0][SIZE_BLOCK_8] = &ov_put_vvc_bi1_pel_pixels_10_8_neon;
  //
  // mc_c->unidir[1][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[1][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[1][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_c->unidir[2][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[2][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[2][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_c->unidir[3][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[3][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[3][SIZE_BLOCK_8] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  mc_c->bidir0[0][SIZE_BLOCK_16] = &ov_put_vvc_bi0_pel_pixels_10_16_neon;
  mc_c->bidir1[0][SIZE_BLOCK_16] = &ov_put_vvc_bi1_pel_pixels_10_16_neon;
  //
  // mc_c->unidir[1][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[1][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[1][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_c->unidir[2][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[2][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[2][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  //
  // mc_c->unidir[3][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir0[3][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
  // mc_c->bidir1[3][SIZE_BLOCK_16] = &dav1d_ff_hevc_transform_4x4_neon_8;
}
