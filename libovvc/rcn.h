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

#ifndef RCN_H
#define RCN_H

#include "ovdefs.h"

struct RCNFunctions;

void rcn_init_gpm_params();

void rcn_init_functions(struct RCNFunctions *rcn_func, uint8_t ict_type, uint8_t lm_chroma_enabled,
                        uint8_t sps_chroma_vertical_collocated_flag, uint8_t lmcs_flag, uint8_t bitdepth,
                        uint8_t sh_dep_quant_used_flag);

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
                               
void rcn_init_ibc_10(struct RCNFunctions *rcn_funcs);

void rcn_init_fill_ref_10(struct RCNFunctions *const rcn_funcs);

void rcn_init_dequant_dpq_10(struct RCNFunctions *rcn_funcs);
void rcn_init_dequant_sdh_10(struct RCNFunctions *rcn_funcs);

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

void rcn_init_dequant_dpq_8(struct RCNFunctions *rcn_funcs);
void rcn_init_dequant_sdh_8(struct RCNFunctions *rcn_funcs);

void rcn_init_ibc_8(struct RCNFunctions *rcn_funcs);

void
vvc_add_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                 int log2_tb_w, int log2_tb_h,
                 int scale);

void
vvc_sub_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                int log2_tb_w, int log2_tb_h,
                int scale);

void
vvc_add_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                      int log2_tb_w, int log2_tb_h,
                      int scale);

void
vvc_sub_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                      int log2_tb_w, int log2_tb_h,
                      int scale);

void
vvc_scale_add_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                     int log2_tb_w, int log2_tb_h,
                     int scale);

void
vvc_scale_sub_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                     int log2_tb_w, int log2_tb_h,
                     int scale);

void
vvc_scale_add_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                          int log2_tb_w, int log2_tb_h,
                          int scale);

void
vvc_scale_sub_half_residual(const int16_t *src, uint16_t *dst, int16_t dst_stride,
                           int log2_tb_w, int log2_tb_h,
                           int scale);

#endif //RCN_H
