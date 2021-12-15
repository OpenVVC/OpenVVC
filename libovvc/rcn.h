#ifndef RCN_H
#define RCN_H

#include "ovdefs.h"

struct RCNFunctions;

void rcn_init_gpm_params();

void rcn_init_functions(struct RCNFunctions *rcn_func, uint8_t ict_type, uint8_t lm_chroma_enabled,
                        uint8_t sps_chroma_vertical_collocated_flag, uint8_t lmcs_flag, uint8_t bitdepth);

void rcn_init_tr_functions(struct RCNFunctions *const rcn_funcs);

void rcn_init_ctu_buffs_10(struct RCNFunctions *rcn_func);

void rcn_init_intra_functions_10(struct RCNFunctions *rcn_func);

void rcn_init_inter_functions_10(struct RCNFunctions *rcn_func);

void rcn_init_cclm_functions_collocated_10(struct RCNFunctions *rcn_func);

void rcn_init_dc_planar_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_cclm_functions_10(struct RCNFunctions *rcn_func);

void rcn_init_lfnst_functions(struct RCNFunctions *rcn_func);

void rcn_init_mip_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_sao_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_dmvr_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_mc_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_prof_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_bdof_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_ciip_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_df_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_lmcs_function_10(struct RCNFunctions *rcn_func, uint8_t lmcs_flag);

void rcn_init_alf_functions_10(struct RCNFunctions *rcn_func);

void rcn_init_tr_functions_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_ict_functions_10(struct RCNFunctions *const rcn_funcs, uint8_t ict_type,
                               uint8_t bitdepth);
void rcn_init_intra_angular_functions_10_sse(struct RCNFunctions *rcn_func);

void rcn_init_fill_ref_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_dequant_10(struct RCNFunctions *rcn_funcs);

void rcn_init_transform_trees_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_ctu_buffs_8(struct RCNFunctions *rcn_func);

void rcn_init_intra_functions_8(struct RCNFunctions *rcn_func);

void rcn_init_inter_functions_8(struct RCNFunctions *rcn_func);

void rcn_init_cclm_functions_collocated_8(struct RCNFunctions *rcn_func);

void rcn_init_dc_planar_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_cclm_functions_8(struct RCNFunctions *rcn_func);

void rcn_init_lfnst_functions(struct RCNFunctions *rcn_func);

void rcn_init_mip_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_sao_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_mc_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_dmvr_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_prof_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_bdof_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_ciip_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_df_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_lmcs_function_8(struct RCNFunctions *rcn_func, uint8_t lmcs_flag);

void rcn_init_alf_functions_8(struct RCNFunctions *rcn_func);

void rcn_init_tr_functions_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_ict_functions_8(struct RCNFunctions *const rcn_funcs, uint8_t ict_type,
                              uint8_t bitdepth);

void rcn_init_fill_ref_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_transform_trees_8(struct RCNFunctions *const rcn_funcs);

void rcn_init_dequant_8(struct RCNFunctions *rcn_funcs);


void
vvc_add_residual(const int16_t *src, uint16_t *dst,
                 int log2_tb_w, int log2_tb_h,
                 int scale);

void
vvc_sub_residual(const int16_t *src, uint16_t *dst,
                int log2_tb_w, int log2_tb_h,
                int scale);

void
vvc_add_half_residual(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale);

void
vvc_sub_half_residual(const int16_t *src, uint16_t *dst,
                      int log2_tb_w, int log2_tb_h,
                      int scale);

void
vvc_scale_add_residual(const int16_t *src, uint16_t *dst,
                     int log2_tb_w, int log2_tb_h,
                     int scale);

void
vvc_scale_sub_residual(const int16_t *src, uint16_t *dst,
                     int log2_tb_w, int log2_tb_h,
                     int scale);

void
vvc_scale_add_half_residual(const int16_t *src, uint16_t *dst,
                          int log2_tb_w, int log2_tb_h,
                          int scale);

void
vvc_scale_sub_half_residual(const int16_t *src, uint16_t *dst,
                           int log2_tb_w, int log2_tb_h,
                           int scale);

#endif //RCN_H
