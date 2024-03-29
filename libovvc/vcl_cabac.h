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

#ifndef VCL_CABAC_H
#define VCL_CABAC_H
#include "stdint.h"

#define OVCABAC_NB_CTX 393

typedef struct OVCABACCtx OVCABACCtx;

/**
 * Offsets to first syntax element context in CABAC context tables
 */
enum SyntaxElemOffset
{
    SPLIT_FLAG_CTX_OFFSET = 0,
    SPLIT_QT_FLAG_CTX_OFFSET = 9,
    SPLIT_HV_FLAG_CTX_OFFSET = 15,
    SPLIT12_FLAG_CTX_OFFSET = 20,
    MODE_CONS_FLAG_CTX_OFFSET = 24,
    SKIP_FLAG_CTX_OFFSET = 26,
    MERGE_FLAG_CTX_OFFSET = 29,
    REGULAR_MERGE_FLAG_CTX_OFFSET = 30,
    MERGE_IDX_CTX_OFFSET = 32,
    MMVD_FLAG_CTX_OFFSET = 33,
    MMVD_MERGE_IDX_CTX_OFFSET = 34,
    MMVD_STEP_MVP_IDX_CTX_OFFSET = 35,
    PRED_MODE_CTX_OFFSET = 36,
    MULTI_REF_LINE_IDX_CTX_OFFSET = 38,
    INTRA_LUMA_MPM_FLAG_CTX_OFFSET = 40,
    INTRA_LUMA_PLANAR_FLAG_CTX_OFFSET = 41,
    CCLM_MODE_FLAG_CTX_OFFSET = 43,
    CCLM_MODE_IDX_CTX_OFFSET = 44,
    INTRA_CHROMA_PRED_MODE_CTX_OFFSET = 45,
    MIP_FLAG_CTX_OFFSET = 46,
    DELTA_QP_CTX_OFFSET = 50,
    INTER_DIR_CTX_OFFSET = 52,
    REF_PIC_CTX_OFFSET = 58,
    SUBBLOCK_MERGE_FLAG_CTX_OFFSET = 60,
    AFFINE_FLAG_CTX_OFFSET = 63,
    AFFINE_TYPE_CTX_OFFSET = 66,
    AFF_MERGE_IDX_CTX_OFFSET = 67,
    BCW_IDX_CTX_OFFSET = 68,
    MVD_CTX_OFFSET = 69,
    BDPCM_MODE_CTX_OFFSET = 71,
    QT_ROOT_CBF_CTX_OFFSET = 75,
    ACT_FLAG_CTX_OFFSET = 76,
    QT_CBF_CTX_OFFSET = 77,
    QT_CBF_CB_CTX_OFFSET = 81,
    QT_CBF_CR_CTX_OFFSET = 83,
    SIG_COEFF_GROUP_CTX_OFFSET = 86,
    SIG_COEFF_GROUP_C_CTX_OFFSET = 88,
    SIG_FLAG_CTX_OFFSET = 90,
    SIG_FLAG_C_CTX_OFFSET = 126,
    PAR_FLAG_CTX_OFFSET = 150,
    PAR_FLAG_C_CTX_OFFSET = 171,
    GT1_FLAG_CTX_OFFSET = 182,
    GT1_FLAG_C_CTX_OFFSET = 203,
    GT0_FLAG_CTX_OFFSET = 214,
    GT0_FLAG_C_CTX_OFFSET = 235,
    LAST_X_CTX_OFFSET = 246,
    LAST_X_C_CTX_OFFSET = 266,
    LAST_Y_CTX_OFFSET = 269,
    LAST_Y_C_CTX_OFFSET = 289,
    MVP_IDX_CTX_OFFSET = 292,
    SMVD_FLAG_CTX_OFFSET = 293,
    SAO_MERGE_FLAG_CTX_OFFSET = 294,
    SAO_TYPE_IDX_CTX_OFFSET = 295,
    LFNST_IDX_CTX_OFFSET = 296,
    PLT_FLAG_CTX_OFFSET = 299,
    ROTATION_FLAG_CTX_OFFSET = 300,
    RUN_TYPE_FLAG_CTX_OFFSET = 301,
    IDX_RUN_MODEL_CTX_OFFSET = 302,
    COPY_RUN_MODEL_CTX_OFFSET = 307,
    RDPCM_FLAG_CTX_OFFSET = 310,
    RDPCM_DIR_CTX_OFFSET = 312,
    TRANSFORM_SKIP_FLAG_CTX_OFFSET = 314,
    MTS_IDX_CTX_OFFSET = 316,
    ISP_MODE_CTX_OFFSET = 320,
    SBT_FLAG_CTX_OFFSET = 322,
    SBT_QUAD_FLAG_CTX_OFFSET = 324,
    SBT_HOR_FLAG_CTX_OFFSET = 325,
    SBT_POS_FLAG_CTX_OFFSET = 328,
    CROSS_COMP_PRED_CTX_OFFSET = 329,
    CHROMA_QP_ADJ_FLAG_CTX_OFFSET = 339,
    CHROMA_QP_ADJ_IDC_CTX_OFFSET = 340,
    IMV_FLAG_CTX_OFFSET = 341,
    CTB_ALF_FLAG_CTX_OFFSET = 346,
    CTB_ALF_ALTERNATIVE_CTX_OFFSET = 355,
    ALF_USE_TEMPORAL_FILT_CTX_OFFSET = 357,
    CC_ALF_FILTER_CONTROL_FLAG_CTX_OFFSET = 358,
    CIIP_FLAG_CTX_OFFSET = 364,
    IBC_FLAG_CTX_OFFSET = 365,
    JOINT_CB_CR_FLAG_CTX_OFFSET = 368,
    TS_SIG_COEFF_GROUP_CTX_OFFSET = 371,
    TS_SIG_FLAG_CTX_OFFSET = 374,
    TS_PAR_FLAG_CTX_OFFSET = 377,
    TS_GTX_FLAG_CTX_OFFSET = 378,
    TS_LRG1_FLAG_CTX_OFFSET = 383,
    TS_RESIDUAL_SIGN_CTX_OFFSET = 387
};

struct OVCABACCtx{
    const uint8_t *bytestream;
    const uint8_t *bytestream_start;
    const uint8_t *bytestream_end;
    uint32_t range;
    uint32_t low_b;
    uint64_t ctx_table[OVCABAC_NB_CTX];
};


int ovcabac_attach_entry(OVCABACCtx *const cabac_ctx, const uint8_t *const entry_point,
                         const uint8_t *const entry_end);

void ovcabac_init_slice_context_table(uint64_t ctx_table[], uint8_t slice_type, 
                                      int slice_qp);

int ovcabac_end_of_slice(OVCABACCtx *cabac_ctx);

#endif
