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

#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_private.h"

typedef struct OVGCI
{
    uint8_t gci_present_flag;
    /* general */
    uint8_t gci_intra_only_constraint_flag;
    uint8_t gci_all_layers_independent_constraint_flag;
    uint8_t gci_one_au_only_constraint_flag;
    /* picture format */
    uint8_t gci_sixteen_minus_max_bitdepth_constraint_idc;
    uint8_t gci_three_minus_max_chroma_format_constraint_idc;
    /* NAL unit type related */
    uint8_t gci_no_mixed_nalu_types_in_pic_constraint_flag;
    uint8_t gci_no_trail_constraint_flag;
    uint8_t gci_no_stsa_constraint_flag;
    uint8_t gci_no_rasl_constraint_flag;
    uint8_t gci_no_radl_constraint_flag;
    uint8_t gci_no_idr_constraint_flag;
    uint8_t gci_no_cra_constraint_flag;
    uint8_t gci_no_gdr_constraint_flag;
    uint8_t gci_no_aps_constraint_flag;
    uint8_t gci_no_idr_rpl_constraint_flag;
    /* tile, slice, subpicture partitioning */
    uint8_t gci_one_tile_per_pic_constraint_flag;
    uint8_t gci_pic_header_in_slice_header_constraint_flag;
    uint8_t gci_one_slice_per_pic_constraint_flag;
    uint8_t gci_no_rectangular_slice_constraint_flag;
    uint8_t gci_one_slice_per_subpic_constraint_flag;
    uint8_t gci_no_subpic_info_constraint_flag;
    /* CTU and block partitioning */
    uint8_t gci_three_minus_max_log2_ctu_size_constraint_idc;
    uint8_t gci_no_partition_constraints_override_constraint_flag;
    uint8_t gci_no_mtt_constraint_flag;
    uint8_t gci_no_qtbtt_dual_tree_intra_constraint_flag;
    /* intra */
    uint8_t gci_no_palette_constraint_flag;
    uint8_t gci_no_ibc_constraint_flag;
    uint8_t gci_no_isp_constraint_flag;
    uint8_t gci_no_mrl_constraint_flag;
    uint8_t gci_no_mip_constraint_flag;
    uint8_t gci_no_cclm_constraint_flag;
    /* inter */
    uint8_t gci_no_ref_pic_resampling_constraint_flag;
    uint8_t gci_no_res_change_in_clvs_constraint_flag;
    uint8_t gci_no_weighted_prediction_constraint_flag;
    uint8_t gci_no_ref_wraparound_constraint_flag;
    uint8_t gci_no_temporal_mvp_constraint_flag;
    uint8_t gci_no_sbtmvp_constraint_flag;
    uint8_t gci_no_amvr_constraint_flag;
    uint8_t gci_no_bdof_constraint_flag;
    uint8_t gci_no_smvd_constraint_flag;
    uint8_t gci_no_dmvr_constraint_flag;
    uint8_t gci_no_mmvd_constraint_flag;
    uint8_t gci_no_affine_motion_constraint_flag;
    uint8_t gci_no_prof_constraint_flag;
    uint8_t gci_no_bcw_constraint_flag;
    uint8_t gci_no_ciip_constraint_flag;
    uint8_t gci_no_gpm_constraint_flag;
    /* transform, quantization, residual */
    uint8_t gci_no_luma_transform_size_64_constraint_flag;
    uint8_t gci_no_transform_skip_constraint_flag;
    uint8_t gci_no_bdpcm_constraint_flag;
    uint8_t gci_no_mts_constraint_flag;
    uint8_t gci_no_lfnst_constraint_flag;
    uint8_t gci_no_joint_cbcr_constraint_flag;
    uint8_t gci_no_sbt_constraint_flag;
    uint8_t gci_no_act_constraint_flag;
    uint8_t gci_no_explicit_scaling_list_constraint_flag;
    uint8_t gci_no_dep_quant_constraint_flag;
    uint8_t gci_no_sign_data_hiding_constraint_flag;
    uint8_t gci_no_cu_qp_delta_constraint_flag;
    uint8_t gci_no_chroma_qp_offset_constraint_flag;
    /* loop filter */
    uint8_t gci_no_sao_constraint_flag;
    uint8_t gci_no_alf_constraint_flag;
    uint8_t gci_no_ccalf_constraint_flag;
    uint8_t gci_no_lmcs_constraint_flag;
    uint8_t gci_no_ladf_constraint_flag;
    uint8_t gci_no_virtual_boundaries_constraint_flag;
    /* unused */
    uint8_t gci_num_reserved_bits;
    uint8_t gci_reserved_zero_bit[256];
    uint8_t gci_alignment_zero_bit;
} OVGCI;

typedef struct OVPTL
{
    uint8_t general_profile_idc;
    uint8_t general_tier_flag;
    uint8_t general_level_idc;
    uint8_t ptl_frame_only_constraint_flag;
    uint8_t ptl_multilayer_enabled_flag;
    uint8_t ptl_sublayer_level_present_flag[64];
    uint8_t ptl_reserved_zero_bit;
    uint8_t sublayer_level_idc[64];
    uint8_t ptl_num_sub_profiles;
    uint8_t general_sub_profile_idc[64];
} OVPTL;

#if 0
profile_tier_level(profileTierPresentFlag, MaxNumSubLayersMinus1)
{
    if(profileTierPresentFlag) {
        ptl->general_profile_idc = nvcl_read_bits(rdr, 7);
        ptl->general_tier_flag = nvcl_read_flag(rdr);
    }

    ptl->general_level_idc = nvcl_read_bits(rdr, 8);
    ptl->ptl_frame_only_constraint_flag = nvcl_read_flag(rdr);
    ptl->ptl_multilayer_enabled_flag = nvcl_read_flag(rdr);

    if(profileTierPresentFlag) {
        general_constraints_info();
    }

    for (i = MaxNumSubLayersMinus1 - 1; i >= 0; i--) {
        ptl->ptl_sublayer_level_present_flag[i] = nvcl_read_flag(rdr);
    }

    while(!byte_aligned()) {
        ptl->ptl_reserved_zero_bit= nvcl_read_bits(rdr, 1);
    }

    for (i = MaxNumSubLayersMinus1 - 1; i >= 0; i--) {
        if(ptl_sublayer_level_present_flag[i]) {
            ptl->sublayer_level_idc[i]= nvcl_read_bits(rdr, 8);
        }
    }

    if(profileTierPresentFlag) {
        ptl->ptl_num_sub_profiles = nvcl_read_bits(rdr, 8);
        for (i = 0; i < ptl_num_sub_profiles; i++) {
            ptl->general_sub_profile_idc[i] = nvcl_read_bits(rdr, 32);
        }
    }
}
#endif

int
profile_tier_level_vps(OVNVCLReader *const rdr, uint8_t vps_pt_present_flag, uint8_t vps_ptl_max_tid)
{
    int i;
    struct OVPTL *ptl;
    struct OVPTL ptl_stck;
    ptl = &ptl_stck;
    if (vps_pt_present_flag) {
        ptl->general_profile_idc = nvcl_read_bits(rdr, 7);
        ptl->general_tier_flag = nvcl_read_flag(rdr);
    }

    ptl->general_level_idc = nvcl_read_bits(rdr, 8);
    ptl->ptl_frame_only_constraint_flag = nvcl_read_flag(rdr);
    ptl->ptl_multilayer_enabled_flag = nvcl_read_flag(rdr);

    if(vps_pt_present_flag) {
        general_constraints_info(rdr);
    }

    for (i = vps_ptl_max_tid - 1; i >= 0; i--) {
        ptl->ptl_sublayer_level_present_flag[i] = nvcl_read_flag(rdr);
    }

    #if 0
    while(!byte_aligned()) {
        ptl->ptl_reserved_zero_bit= nvcl_read_bits(rdr, 1);
    }
    #else
    nvcl_align(rdr);
    #endif

    for (i = vps_ptl_max_tid - 1; i >= 0; i--) {
        if(ptl->ptl_sublayer_level_present_flag[i]) {
            ptl->sublayer_level_idc[i]= nvcl_read_bits(rdr, 8);
        }
    }

    if (vps_pt_present_flag) {
        ptl->ptl_num_sub_profiles = nvcl_read_bits(rdr, 8);
        for (i = 0; i < ptl->ptl_num_sub_profiles; i++) {
            ptl->general_sub_profile_idc[i] = nvcl_read_bits(rdr, 32);
        }
    }
    return 0;
}

int
profile_tier_level_sps(OVNVCLReader *const rdr,  uint8_t sps_max_sublayers_minus1)
{
    int i;
    struct OVPTL *ptl;
    struct OVPTL ptl_stck;
    ptl = &ptl_stck;
    ptl->general_profile_idc = nvcl_read_bits(rdr, 7);
    ptl->general_tier_flag = nvcl_read_flag(rdr);

    #if 1
    ptl->general_level_idc = nvcl_read_bits(rdr, 8);
    ptl->ptl_frame_only_constraint_flag = nvcl_read_flag(rdr);
    ptl->ptl_multilayer_enabled_flag = nvcl_read_flag(rdr);
    #endif

    general_constraints_info(rdr);

    for (i = sps_max_sublayers_minus1 - 1; i >= 0; i--) {
        ptl->ptl_sublayer_level_present_flag[i] = nvcl_read_flag(rdr);
    }

    #if 0
    ptl->general_level_idc = nvcl_read_bits(rdr, 8);
    #endif
    /* FIXME check align function */
    #if 0
    while(!byte_aligned()) {
        ptl->ptl_reserved_zero_bit= nvcl_read_bits(rdr, 1);
    }
    #else
    nvcl_align(rdr);
    #endif

    for (i = sps_max_sublayers_minus1 - 1; i >= 0; i--) {
        if(ptl->ptl_sublayer_level_present_flag[i]) {
            ptl->sublayer_level_idc[i]= nvcl_read_bits(rdr, 8);
        }
    }

    ptl->ptl_num_sub_profiles = nvcl_read_bits(rdr, 8);
    for (i = 0; i < ptl->ptl_num_sub_profiles; i++) {
        ptl->general_sub_profile_idc[i] = nvcl_read_bits(rdr, 32);
    }
    nvcl_align(rdr);

    return 0;
}

int
profile_tier_level_dci(OVNVCLReader *const rdr)
{
    int i;
    struct OVPTL *ptl;
    struct OVPTL ptl_stck;
    ptl = &ptl_stck;

    ptl->general_profile_idc = nvcl_read_bits(rdr, 7);
    ptl->general_tier_flag = nvcl_read_flag(rdr);

    ptl->general_level_idc = nvcl_read_bits(rdr, 8);
    ptl->ptl_frame_only_constraint_flag = nvcl_read_flag(rdr);
    ptl->ptl_multilayer_enabled_flag = nvcl_read_flag(rdr);

    general_constraints_info(rdr);

    /* FIXME check align function */
    #if 0
    while(!byte_aligned()) {
        ptl->ptl_reserved_zero_bit= nvcl_read_bits(rdr, 1);
    }
    #else
    nvcl_align(rdr);
    #endif

    ptl->ptl_num_sub_profiles = nvcl_read_bits(rdr, 8);
    for (i = 0; i < ptl->ptl_num_sub_profiles; i++) {
        ptl->general_sub_profile_idc[i] = nvcl_read_bits(rdr, 32);
    }
    return 0;
}

int
general_constraints_info(OVNVCLReader *const rdr)
{
    /*FIXME read before function call */
    struct OVGCI *gci;
    struct OVGCI gci_stck;
    gci = &gci_stck;
    gci->gci_present_flag = nvcl_read_flag(rdr);
    if (gci->gci_present_flag) {
        int i;
        /* general */
        gci->gci_intra_only_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_all_layers_independent_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_one_au_only_constraint_flag = nvcl_read_flag(rdr);
        /* picture format */
        gci->gci_sixteen_minus_max_bitdepth_constraint_idc = nvcl_read_bits(rdr, 4);
        gci->gci_three_minus_max_chroma_format_constraint_idc = nvcl_read_bits(rdr, 2);
        /* NAL unit type related */
        gci->gci_no_mixed_nalu_types_in_pic_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_trail_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_stsa_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_rasl_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_radl_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_idr_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_cra_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_gdr_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_aps_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_idr_rpl_constraint_flag = nvcl_read_flag(rdr);
        /* tile, slice, subpicture partitioning */
        gci->gci_one_tile_per_pic_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_pic_header_in_slice_header_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_one_slice_per_pic_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_rectangular_slice_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_one_slice_per_subpic_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_subpic_info_constraint_flag = nvcl_read_flag(rdr);
        /* CTU and block partitioning */
        gci->gci_three_minus_max_log2_ctu_size_constraint_idc = nvcl_read_bits(rdr, 2);
        gci->gci_no_partition_constraints_override_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_mtt_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_qtbtt_dual_tree_intra_constraint_flag = nvcl_read_flag(rdr);
        /* intra */
        gci->gci_no_palette_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_ibc_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_isp_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_mrl_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_mip_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_cclm_constraint_flag = nvcl_read_flag(rdr);
        /* inter */
        gci->gci_no_ref_pic_resampling_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_res_change_in_clvs_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_weighted_prediction_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_ref_wraparound_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_temporal_mvp_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_sbtmvp_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_amvr_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_bdof_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_smvd_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_dmvr_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_mmvd_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_affine_motion_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_prof_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_bcw_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_ciip_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_gpm_constraint_flag = nvcl_read_flag(rdr);
        /* transform, quantization, residual */
        gci->gci_no_luma_transform_size_64_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_transform_skip_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_bdpcm_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_mts_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_lfnst_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_joint_cbcr_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_sbt_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_act_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_explicit_scaling_list_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_dep_quant_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_sign_data_hiding_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_cu_qp_delta_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_chroma_qp_offset_constraint_flag = nvcl_read_flag(rdr);
        /* loop filter */
        gci->gci_no_sao_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_alf_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_ccalf_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_lmcs_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_ladf_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_no_virtual_boundaries_constraint_flag = nvcl_read_flag(rdr);
        gci->gci_num_reserved_bits = nvcl_read_bits(rdr, 8);
        for (i = 0; i < gci->gci_num_reserved_bits; i++) {
            gci->gci_reserved_zero_bit[i] = nvcl_read_bits(rdr, 1);
        }
    }

    #if 0
    while( !byte_aligned()) {
        gci->gci_alignment_zero_bit = nvcl_read_bits(rdr, 1);
    }
    #else
    nvcl_align(rdr);
    #endif
    return 0;
}
