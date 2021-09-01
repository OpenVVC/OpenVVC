#include "arm/rcn_neon.h"

#define SIZE_BLOCK_4 1
#define SIZE_BLOCK_8 2
#define SIZE_BLOCK_16 3
#define SIZE_BLOCK_32 4
#define SIZE_BLOCK_64 5
#define SIZE_BLOCK_128 6

void ov_put_vvc_bi0_pel_pixels_10_4_neon();
void ov_put_vvc_bi0_pel_pixels_10_8_neon();
void ov_put_vvc_bi0_pel_pixels_10_16_neon();
void ov_put_vvc_bi0_pel_pixels_10_32_neon();
// void ov_put_vvc_bi1_pel_pixels_10_4_neon();//untested
void ov_put_vvc_bi1_pel_pixels_10_8_neon();
void ov_put_vvc_bi1_pel_pixels_10_16_neon();
void ov_put_vvc_bi1_pel_pixels_10_32_neon();
void
rcn_init_mc_functions_neon(struct RCNFunctions* const rcn_funcs)
{
  struct MCFunctions* const mc_l = &rcn_funcs->mc_l;
  struct MCFunctions* const mc_c = &rcn_funcs->mc_c;

  /* Luma functions */

  mc_l->bidir0[0][SIZE_BLOCK_8] = &ov_put_vvc_bi0_pel_pixels_10_8_neon;
  mc_l->bidir1[0][SIZE_BLOCK_8] = &ov_put_vvc_bi1_pel_pixels_10_8_neon;
  //
  mc_l->bidir0[0][SIZE_BLOCK_16] = &ov_put_vvc_bi0_pel_pixels_10_16_neon;
  mc_l->bidir1[0][SIZE_BLOCK_16] = &ov_put_vvc_bi1_pel_pixels_10_16_neon;
  //
  mc_l->bidir0[0][SIZE_BLOCK_32] = &ov_put_vvc_bi0_pel_pixels_10_32_neon;
  mc_l->bidir1[0][SIZE_BLOCK_32] = &ov_put_vvc_bi1_pel_pixels_10_32_neon;
  //
  /* Chroma functions */
  mc_c->bidir0[0][SIZE_BLOCK_4] = &ov_put_vvc_bi0_pel_pixels_10_4_neon;
  // mc_c->bidir1[0][SIZE_BLOCK_4] = &ov_put_vvc_bi1_pel_pixels_10_4_neon;
  //
  mc_c->bidir0[0][SIZE_BLOCK_8] = &ov_put_vvc_bi0_pel_pixels_10_8_neon;
  mc_c->bidir1[0][SIZE_BLOCK_8] = &ov_put_vvc_bi1_pel_pixels_10_8_neon;
  //
  mc_c->bidir0[0][SIZE_BLOCK_16] = &ov_put_vvc_bi0_pel_pixels_10_16_neon;
  mc_c->bidir1[0][SIZE_BLOCK_16] = &ov_put_vvc_bi1_pel_pixels_10_16_neon;
  //
}
