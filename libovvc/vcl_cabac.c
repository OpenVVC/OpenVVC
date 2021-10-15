#include "vcl_cabac.h"
#include "cabac_internal.h"

/* Default Window Size for CABAC rates */
#define DWS 8

#if 0
#define SPLIT_FLAG_CTX_OFFSET 0
#define SPLIT_QT_FLAG_CTX_OFFSET 9
#define SPLIT_HV_FLAG_CTX_OFFSET 15
#define SPLIT12_FLAG_CTX_OFFSET 20
#define MODE_CONS_FLAG_CTX_OFFSET 24
#define SKIP_FLAG_CTX_OFFSET 26
#define MERGE_FLAG_CTX_OFFSET 29
#define REGULAR_MERGE_FLAG_CTX_OFFSET 30
#define MERGE_IDX_CTX_OFFSET 32
#define MMVD_FLAG_CTX_OFFSET 33
#define MMVD_MERGE_IDX_CTX_OFFSET 34
#define MMVD_STEP_MVP_IDX_CTX_OFFSET 35
#define PRED_MODE_CTX_OFFSET 36
#define MULTI_REF_LINE_IDX_CTX_OFFSET 38
#define INTRA_LUMA_MPM_FLAG_CTX_OFFSET 40
#define INTRA_LUMA_PLANAR_FLAG_CTX_OFFSET 41
#define CCLM_MODE_FLAG_CTX_OFFSET 43
#define CCLM_MODE_IDX_CTX_OFFSET 44
#define INTRA_CHROMA_PRED_MODE_CTX_OFFSET 45
#define MIP_FLAG_CTX_OFFSET 46
#define DELTA_QP_CTX_OFFSET 50
#define INTER_DIR_CTX_OFFSET 52
#define REF_PIC_CTX_OFFSET 58
#define SUBBLOCK_MERGE_FLAG_CTX_OFFSET 60
#define AFFINE_FLAG_CTX_OFFSET 63
#define AFFINE_TYPE_CTX_OFFSET 66
#define AFF_MERGE_IDX_CTX_OFFSET 67
#define BCW_IDX_CTX_OFFSET 68
#define MVD_CTX_OFFSET 69
#define BDPCM_MODE_CTX_OFFSET 71
#define QT_ROOT_CBF_CTX_OFFSET 75
#define ACT_FLAG_CTX_OFFSET 76
#define QT_CBF_CTX_OFFSET 77
#define QT_CBF_CB_CTX_OFFSET 81
#define QT_CBF_CR_CTX_OFFSET 83
#define SIG_COEFF_GROUP_CTX_OFFSET 86
#define SIG_COEFF_GROUP_C_CTX_OFFSET 88
#define SIG_FLAG_CTX_OFFSET 90
#define SIG_FLAG_C_CTX_OFFSET 126
#define PAR_FLAG_CTX_OFFSET 150
#define PAR_FLAG_C_CTX_OFFSET 171
#define GT1_FLAG_CTX_OFFSET 182
#define GT1_FLAG_C_CTX_OFFSET 203
#define GT0_FLAG_CTX_OFFSET 214
#define GT0_FLAG_C_CTX_OFFSET 235
#define LAST_X_CTX_OFFSET 246
#define LAST_X_C_CTX_OFFSET 266
#define LAST_Y_CTX_OFFSET 269
#define LAST_Y_C_CTX_OFFSET 289
#define MVP_IDX_CTX_OFFSET 292
#define SMVD_FLAG_CTX_OFFSET 293
#define SAO_MERGE_FLAG_CTX_OFFSET 294
#define SAO_TYPE_IDX_CTX_OFFSET 295
#define LFNST_IDX_CTX_OFFSET 296
#define PLT_FLAG_CTX_OFFSET 299
#define ROTATION_FLAG_CTX_OFFSET 300
#define RUN_TYPE_FLAG_CTX_OFFSET 301
#define IDX_RUN_MODEL_CTX_OFFSET 302
#define COPY_RUN_MODEL_CTX_OFFSET 307
#define RDPCM_FLAG_CTX_OFFSET 310
#define RDPCM_DIR_CTX_OFFSET 312
#define TRANSFORM_SKIP_FLAG_CTX_OFFSET 314
#define MTS_IDX_CTX_OFFSET 316
#define ISP_MODE_CTX_OFFSET 320
#define SBT_FLAG_CTX_OFFSET 322
#define SBT_QUAD_FLAG_CTX_OFFSET 324
#define SBT_HOR_FLAG_CTX_OFFSET 325
#define SBT_POS_FLAG_CTX_OFFSET 328
#define CROSS_COMP_PRED_CTX_OFFSET 329
#define CHROMA_QP_ADJ_FLAG_CTX_OFFSET 339
#define CHROMA_QP_ADJ_IDC_CTX_OFFSET 340
#define IMV_FLAG_CTX_OFFSET 341
#define CTB_ALF_FLAG_CTX_OFFSET 346
#define CTB_ALF_ALTERNATIVE_CTX_OFFSET 355
#define ALF_USE_TEMPORAL_FILT_CTX_OFFSET 357
#define CC_ALF_FILTER_CONTROL_FLAG_CTX_OFFSET 358
#define CIIP_FLAG_CTX_OFFSET 364
#define IBC_FLAG_CTX_OFFSET 365
#define JOINT_CB_CR_FLAG_CTX_OFFSET 368
#define TS_SIG_COEFF_GROUP_CTX_OFFSET 371
#define TS_SIG_FLAG_CTX_OFFSET 374
#define TS_PAR_FLAG_CTX_OFFSET 377
#define TS_GTX_FLAG_CTX_OFFSET 378
#define TS_LRG1_FLAG_CTX_OFFSET 383
#define TS_RESIDUAL_SIGN_CTX_OFFSET 387
#endif

#define CNU 35

#define OV_ERROR -1

enum OVSliceType {
     OV_SLICE_B = 0,
     OV_SLICE_P = 1,
     OV_SLICE_I = 2,
};

#if 0
/* Values used for CABAC context table initialisation
 * for INTER B slice table
 */
static const uint8_t ctx_init_table_b_slice[OVCABAC_NB_CTX] = 
{
    18,  27,  15,  18,  28,  45,  26,   7,  23,   /*split_flag*/
    26,  36,  38,  18,  34,  21,   /*split_qt_flag*/
    43,  42,  37,  42,  44,   /*split_hv_flag*/
    28,  29,  28,  29,   /*split12_flag*/
    25,  20,   /*mode_cons_flag*/
    57,  60,  46,   /*skip_flag*/
    6,   /*merge_flag*/
    46,  15,   /*regular_merge_flag*/
    18,   /*merge_idx*/
    25,   /*mmvd_flag*/
    43,   /*mmvd_merge_idx*/
    59,   /*mmvd_step_mvp_idx*/
    40,  35,   /*pred_mode*/
    25,  59,   /*multi_ref_line_idx*/
    44,   /*intra_luma_mpm_flag*/
    13,   6,   /*intra_luma_planar_flag*/
    26,   /*cclm_mode_flag*/
    27,   /*cclm_mode_idx*/
    25,   /*intra_chroma_pred_mode*/
    56,  57,  50,  26,   /*mip_flag*/
    CNU, CNU,   /*delta_qp*/
    14,  13,   5,   4,   3,  40,   /*inter_dir*/
    5,  35,   /*ref_pic*/
    25,  58,  45,   /*subblock_merge_flag*/
    19,  13,   6,   /*affine_flag*/
    35,   /*affine_type*/
    4,   /*aff_merge_idx*/
    5,   /*bcw_idx*/
    51,  36,   /*mvd*/
    19,  21,   0,  28,   /*bdpcm_mode*/
    12,   /*qt_root_cbf*/
    46,   /*act_flag*/
    15,   6,   5,  14,   /*qt_cbf[]*/
    25,  37,   /*qt_cbf[]*/
    9,  36,  45,   /*qt_cbf[]*/
    25,  45,   /*sig_coeff_group[]*/
    25,  14,   /*sig_coeff_group[]*/
    17,  41,  49,  36,   1,  49,  50,  37,  48,  51,  58,  45,   /*sig_flag[]*/
    26,  45,  53,  46,  49,  54,  61,  39,  35,  39,  39,  39,   /*sig_flag[]*/
    19,  54,  39,  39,  50,  39,  39,  39,   0,  39,  39,  39,   /*sig_flag[]*/
    9,  49,  50,  36,  48,  59,  59,  38,   /*sig_flag_c[]*/
    34,  45,  38,  31,  58,  39,  39,  39,   /*sig_flag_c[]*/
    34,  38,  54,  39,  41,  39,  39,  39,   /*sig_flag_c[]*/
    33,  40,  25,  41,  26,  42,  25,  33,  26,  34,  27,  25,  41,  42,  42,  35,  33,  27,  35,  42,  43,   /*par_flag[]*/
    33,  25,  26,  34,  19,  27,  33,  42,  43,  35,  43,   /*par_flag_c[]*/
    25,   0,   0,  17,  25,  26,   0,   9,  25,  33,  19,   0,  25,  33,  26,  20,  25,  33,  27,  35,  22,   /*gt1_flag[]*/
    25,   1,  25,  33,  26,  12,  25,  33,  27,  28,  37,   /*gt1_flag_c[]*/
    0,   0,  33,  34,  35,  21,  25,  34,  35,  28,  29,  40,  42,  43,  29,  30,  49,  36,  37,  45,  38,   /*gt0_flag[]*/
    0,  40,  34,  43,  36,  37,  57,  52,  45,  38,  46,   /*gt0_flag_c[]*/
    6,   6,  12,  14,   6,   4,  14,   7,   6,   4,  29,   7,   6,   6,  12,  28,   7,  13,  13,  35,   /*last_x*/
    19,   5,   4,   /*last_x_c*/
    5,   5,  20,  13,  13,  19,  21,   6,  12,  12,  14,  14,   5,   4,  12,  13,   7,  13,  12,  41,   /*last_y*/
    11,   5,  27,   /*last_y_c*/
    34,   /*mvp_idx*/
    28,   /*smvd_flag*/
    2,   /*sao_merge_flag*/
    2,   /*sao_type_idx*/
    52,  37,  27,   /*lfnst_idx*/
    17,   /*plt_flag*/
    35,   /*rotation_flag*/
    50,   /*run_type_flag*/
    58,  45,  45,  30,  38,   /*idx_run_model*/
    45,  38,  46,   /*copy_run_model*/
    CNU, CNU,   /*rdpcm_flag*/
    CNU, CNU,   /*rdpcm_dir*/
    25,  17,   /*transform_skip_flag*/
    45,  25,  27,   0,   /*mts_idx*/
    33,  43,   /*isp_mode*/
    41,  57,   /*sbt_flag*/
    42,   /*sbt_quad_flag*/
    35,  51,  27,   /*sbt_hor_flag*/
    28,   /*sbt_pos_flag*/
    CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,   /*cross_comp_pred*/
    CNU,   /*chroma_qp_adj_flag*/
    CNU,   /*chroma_qp_adj_idc*/
    59,  26,  50,  60,  38,   /*imv_flag*/
    33,  52,  46,  25,  61,  54,  25,  61,  54,   /*ctb_alf_flag*/
    11,  26,   /*ctb_alf_alternative*/
    46,   /*alf_use_temporal_filt*/
    25,  35,  38,  25,  28,  38,   /*cc_alf_filter_control_flag*/
    57,   /*ciip_flag*/
    0,  43,  45,   /*ibc_flag*/
    42,  43,  52,   /*joint_cb_cr_flag*/
    18,  35,  45,   /*ts_sig_coeff_group*/
    25,  50,  37,   /*ts_sig_flag*/
    11,   /*ts_par_flag*/
    CNU,   3,   4,   4,   5,   /*ts_gtx_flag*/
    19,  11,   4,   6,   /*ts_lrg1_flag*/
    35,  25,  46,  28,  33,  38,   /*ts_residual_sign*/
};


/* Values used for CABAC context table initialisation
 * for INTER P slice table
 */
static const uint8_t ctx_init_table_p_slice[OVCABAC_NB_CTX] = 
{
     11,  35,  53,  12,   6,  30,  13,  15,  31,   /*split_flag*/
     20,  14,  23,  18,  19,   6,   /*split_qt_flag*/
     43,  35,  37,  34,  52,   /*split_hv_flag*/
     43,  37,  21,  22,   /*split12_flag*/
     25,  12,   /*mode_cons_flag*/
     57,  59,  45,   /*skip_flag*/
     21,   /*merge_flag*/
     38,   7,   /*regular_merge_flag*/
     20,   /*merge_idx*/
     26,   /*mmvd_flag*/
     43,   /*mmvd_merge_idx*/
     60,   /*mmvd_step_mvp_idx*/
     40,  35,   /*pred_mode*/
     25,  58,   /*multi_ref_line_idx*/
     36,   /*intra_luma_mpm_flag*/
     12,  20,   /*intra_luma_planar_flag*/
     34,   /*cclm_mode_flag*/
     27,   /*cclm_mode_idx*/
     25,   /*intra_chroma_pred_mode*/
     41,  57,  58,  26,   /*mip_flag*/
    CNU, CNU,   /*delta_qp*/
      7,   6,   5,  12,   4,  40,   /*inter_dir*/
     20,  35,   /*ref_pic*/
     48,  57,  44,   /*subblock_merge_flag*/
     12,  13,  14,   /*affine_flag*/
     35,   /*affine_type*/
      5,   /*aff_merge_idx*/
      4,   /*bcw_idx*/
     44,  43,   /*mvd*/
     40,  36,   0,  13,   /*bdpcm_mode*/
      5,   /*qt_root_cbf*/
     46,   /*act_flag*/
     23,   5,  20,   7,   /*qt_cbf[]*/
     25,  28,   /*qt_cbf[]*/
     25,  29,  45,   /*qt_cbf[]*/
     25,  30,   /*sig_coeff_group[]*/
     25,  45,   /*sig_coeff_group[]*/
     17,  41,  42,  29,  25,  49,  43,  37,  33,  58,  51,  30,   /*sig_flag[]*/
     19,  38,  38,  46,  34,  54,  54,  39,   6,  39,  39,  39,   /*sig_flag[]*/
     19,  39,  54,  39,  19,  39,  39,  39,  56,  39,  39,  39,   /*sig_flag[]*/
     17,  34,  35,  21,  41,  59,  60,  38,   /*sig_flag_c[]*/
     35,  45,  53,  54,  44,  39,  39,  39,   /*sig_flag_c[]*/
     34,  38,  62,  39,  26,  39,  39,  39,   /*sig_flag_c[]*/
     18,  17,  33,  18,  26,  42,  25,  33,  26,  42,  27,  25,  34,  42,  42,  35,  26,  27,  42,  20,  20,   /*par_flag[]*/
     25,  25,  26,  11,  19,  27,  33,  42,  35,  35,  43,   /*par_flag_c[]*/
     17,   0,   1,  17,  25,  18,   0,   9,  25,  33,  34,   9,  25,  18,  26,  20,  25,  18,  19,  27,  29,   /*gt1_flag[]*/
     17,   9,  25,  10,  18,   4,  17,  33,  19,  20,  29,   /*gt1_flag_c[]*/
      0,  17,  26,  19,  35,  21,  25,  34,  20,  28,  29,  33,  27,  28,  29,  22,  34,  28,  44,  37,  38,   /*gt0_flag[]*/
      0,  25,  19,  20,  13,  14,  57,  44,  30,  30,  23,   /*gt0_flag_c[]*/
      6,  13,  12,   6,   6,  12,  14,  14,  13,  12,  29,   7,   6,  13,  36,  28,  14,  13,   5,  26,   /*last_x*/
     12,   4,  18,   /*last_x_c*/
      5,   5,  12,   6,   6,   4,   6,  14,   5,  12,  14,   7,  13,   5,  13,  21,  14,  20,  12,  34,   /*last_y*/
     11,   4,  18,   /*last_y_c*/
     34,   /*mvp_idx*/
     28,   /*smvd_flag*/
     60,   /*sao_merge_flag*/
      5,   /*sao_type_idx*/
     37,  45,  27,   /*lfnst_idx*/
      0,   /*plt_flag*/
     42,   /*rotation_flag*/
     59,   /*run_type_flag*/
     51,  30,  30,  38,  23,   /*idx_run_model*/
     38,  53,  46,   /*copy_run_model*/
    CNU, CNU,   /*rdpcm_flag*/
    CNU, CNU,   /*rdpcm_dir*/
     25,   9,   /*transform_skip_flag*/
     45,  40,  27,   0,   /*mts_idx*/
     33,  36,   /*isp_mode*/
     56,  57,   /*sbt_flag*/
     42,   /*sbt_quad_flag*/
     20,  43,  12,   /*sbt_hor_flag*/
     28,   /*sbt_pos_flag*/
    CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,   /*cross_comp_pred*/
    CNU,   /*chroma_qp_adj_flag*/
    CNU,   /*chroma_qp_adj_idc*/
     59,  48,  58,  60,  60,   /*imv_flag*/
     13,  23,  46,   4,  61,  54,  19,  46,  54,   /*ctb_alf_flag*/
     20,  12,   /*ctb_alf_alternative*/
     46,   /*alf_use_temporal_filt*/
     18,  21,  38,  18,  21,  38,   /*cc_alf_filter_control_flag*/
     57,   /*ciip_flag*/
      0,  57,  44,   /*ibc_flag*/
     27,  36,  45,   /*joint_cb_cr_flag*/
     18,  12,  29,   /*ts_sig_coeff_group*/
     40,  35,  44,   /*ts_sig_flag*/
      3,   /*ts_par_flag*/
    CNU,   2,  10,   3,   3,   /*ts_gtx_flag*/
     18,  11,   4,  28,   /*ts_lrg1_flag*/
      5,  10,  53,  43,  25,  46,   /*ts_residual_sign*/
};

/* Values used for CABAC context table initialisation
 * for INTRA slice table
 */
static const uint8_t ctx_init_table_i_slice[OVCABAC_NB_CTX] = 
{
    19,  28,  38,  27,  29,  38,  20,  30,  31,   /*split_flag*/
    27,   6,  15,  25,  19,  37,   /*split_qt_flag*/
    43,  42,  29,  27,  44,   /*split_hv_flag*/
    36,  45,  36,  45,   /*split12_flag*/
    CNU, CNU,   /*mode_cons_flag*/
    0,  26,  28,   /*skip_flag*/
    26,   /*merge_flag*/
    CNU, CNU,   /*regular_merge_flag*/
    34,   /*merge_idx*/
    CNU,   /*mmvd_flag*/
    CNU,   /*mmvd_merge_idx*/
    CNU,   /*mmvd_step_mvp_idx*/
    CNU, CNU,   /*pred_mode*/
    25,  60,   /*multi_ref_line_idx*/
    45,   /*intra_luma_mpm_flag*/
    13,  28,   /*intra_luma_planar_flag*/
    59,   /*cclm_mode_flag*/
    27,   /*cclm_mode_idx*/
    34,   /*intra_chroma_pred_mode*/
    33,  49,  50,  25,   /*mip_flag*/
    CNU, CNU,   /*delta_qp*/
    CNU, CNU, CNU, CNU, CNU, CNU,   /*inter_dir*/
    CNU, CNU,   /*ref_pic*/
    CNU, CNU, CNU,   /*subblock_merge_flag*/
    CNU, CNU, CNU,   /*affine_flag*/
    CNU,   /*affine_type*/
    CNU,   /*aff_merge_idx*/
    CNU,   /*bcw_idx*/
    14,  45,   /*mvd*/
    19,  35,   1,  27,   /*bdpcm_mode*/
    6,   /*qt_root_cbf*/
    52,   /*act_flag*/
    15,  12,   5,   7,   /*qt_cbf[]*/
    12,  21,   /*qt_cbf[]*/
    33,  28,  36,   /*qt_cbf[]*/
    18,  31,   /*sig_coeff_group[]*/
    25,  15,   /*sig_coeff_group[]*/
    25,  19,  28,  14,  25,  20,  29,  30,  19,  37,  30,  38,   /*sig_flag[]*/
    11,  38,  46,  54,  27,  39,  39,  39,  44,  39,  39,  39,   /*sig_flag[]*/
    18,  39,  39,  39,  27,  39,  39,  39,   0,  39,  39,  39,   /*sig_flag[]*/
    25,  27,  28,  37,  34,  53,  53,  46,   /*sig_flag_c[]*/
    19,  46,  38,  39,  52,  39,  39,  39,   /*sig_flag_c[]*/
    11,  39,  39,  39,  19,  39,  39,  39,   /*sig_flag_c[]*/
    33,  25,  18,  26,  34,  27,  25,  26,  19,  42,  35,  33,  19,  27,  35,  35,  34,  42,  20,  43,  20,   /*par_flag[]*/
    33,  25,  26,  42,  19,  27,  26,  50,  35,  20,  43,   /*par_flag_c[]*/
    25,   1,  40,  25,  33,  11,  17,  25,  25,  18,   4,  17,  33,  26,  19,  13,  33,  19,  20,  28,  22,   /*gt1_flag[]*/
    40,   9,  25,  18,  26,  35,  25,  26,  35,  28,  37,   /*gt1_flag_c[]*/
    25,  25,  11,  27,  20,  21,  33,  12,  28,  21,  22,  34,  28,  29,  29,  30,  36,  29,  45,  30,  23,   /*gt0_flag[]*/
    40,  33,  27,  28,  21,  37,  36,  37,  45,  38,  46,   /*gt0_flag_c[]*/
    13,   5,   4,  21,  14,   4,   6,  14,  21,  11,  14,   7,  14,   5,  11,  21,  30,  22,  13,  42,   /*last_x*/
    12,   4,   3,   /*last_x_c*/
    13,   5,   4,   6,  13,  11,  14,   6,   5,   3,  14,  22,   6,   4,   3,   6,  22,  29,  20,  34,   /*last_y*/
    12,   4,   3,   /*last_y_c*/
    42,   /*mvp_idx*/
    CNU,   /*smvd_flag*/
    60,   /*sao_merge_flag*/
    13,   /*sao_type_idx*/
    28,  52,  42,   /*lfnst_idx*/
    25,   /*plt_flag*/
    42,   /*rotation_flag*/
    42,   /*run_type_flag*/
    50,  37,  45,  30,  46,   /*idx_run_model*/
    45,  38,  46,   /*copy_run_model*/
    CNU, CNU,   /*rdpcm_flag*/
    CNU, CNU,   /*rdpcm_dir*/
    25,   9,   /*transform_skip_flag*/
    29,   0,  28,   0,   /*mts_idx*/
    33,  43,   /*isp_mode*/
    CNU, CNU,   /*sbt_flag*/
    CNU,   /*sbt_quad_flag*/
    CNU, CNU, CNU,   /*sbt_hor_flag*/
    CNU,   /*sbt_pos_flag*/
    CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,   /*cross_comp_pred*/
    CNU,   /*chroma_qp_adj_flag*/
    CNU,   /*chroma_qp_adj_idc*/
    CNU,  34, CNU, CNU, CNU,   /*imv_flag*/
    62,  39,  39,  54,  39,  39,  31,  39,  39,   /*ctb_alf_flag*/
    11,  11,   /*ctb_alf_alternative*/
    46,   /*alf_use_temporal_filt*/
    18,  30,  31,  18,  30,  31,   /*cc_alf_filter_control_flag*/
    CNU,   /*ciip_flag*/
    17,  42,  36,   /*ibc_flag*/
    12,  21,  35,   /*joint_cb_cr_flag*/
    18,  20,  38,   /*ts_sig_coeff_group*/
    25,  28,  38,   /*ts_sig_flag*/
    11,   /*ts_par_flag*/
    CNU,  10,   3,   3,   3,   /*ts_gtx_flag*/
    11,   5,   5,  14,   /*ts_lrg1_flag*/
    12,  17,  46,  28,  25,  46,   /*ts_residual_sign*/
};
#endif
/* Values used for CABAC context table initialisation
 * for INTRA slice table
 */
static const uint8_t rate_init_table[OVCABAC_NB_CTX] = 
{
     12,  13,   8,   8,  13,  12,   5,   9,   9,   /*split_flag*/
      0,   8,   8,  12,  12,   8,   /*split_qt_flag*/
      9,   8,   9,   8,   5,   /*split_hv_flag*/
     12,  13,  12,  13,   /*split12_flag*/
      1,   0,   /*mode_cons_flag*/
      5,   4,   8,   /*skip_flag*/
      4,   /*merge_flag*/
      5,   5,   /*regular_merge_flag*/
      4,   /*merge_idx*/
      4,   /*mmvd_flag*/
     10,   /*mmvd_merge_idx*/
      0,   /*mmvd_step_mvp_idx*/
      5,   1,   /*pred_mode*/
      5,   8,   /*multi_ref_line_idx*/
      6,   /*intra_luma_mpm_flag*/
      1,   5,   /*intra_luma_planar_flag*/
      4,   /*cclm_mode_flag*/
      9,   /*cclm_mode_idx*/
      5,   /*intra_chroma_pred_mode*/
      9,  10,   9,   6,   /*mip_flag*/
    DWS, DWS,   /*delta_qp*/
      0,   0,   1,   4,   4,   0,   /*inter_dir*/
      0,   4,   /*ref_pic*/
      4,   4,   4,   /*subblock_merge_flag*/
      4,   0,   0,   /*affine_flag*/
      4,   /*affine_type*/
      0,   /*aff_merge_idx*/
      1,   /*bcw_idx*/
      9,   5,   /*mvd*/
      1,   4,   1,   0,   /*bdpcm_mode*/
      4,   /*qt_root_cbf*/
      1,   /*act_flag*/
      5,   1,   8,   9,   /*qt_cbf[]*/
      5,   0,   /*qt_cbf[]*/
      2,   1,   0,   /*qt_cbf[]*/
      8,   5,   /*sig_coeff_group[]*/
      5,   8,   /*sig_coeff_group[]*/
     12,   9,   9,  10,   9,   9,   9,  10,   8,   8,   8,  10,   /*sig_flag[]*/
      9,  13,   8,   8,   8,   8,   8,   5,   8,   0,   0,   0,   /*sig_flag[]*/
      8,   8,   8,   8,   8,   0,   4,   4,   0,   0,   0,   0,   /*sig_flag[]*/
     12,  12,   9,  13,   4,   5,   8,   9,   /*sig_flag_c[]*/
      8,  12,  12,   8,   4,   0,   0,   0,   /*sig_flag_c[]*/
      8,   8,   8,   8,   4,   0,   0,   0,   /*sig_flag_c[]*/
      8,   9,  12,  13,  13,  13,  10,  13,  13,  13,  13,  13,  13,  13,  13,  13,  10,  13,  13,  13,  13,   /*par_flag[]*/
      8,  12,  12,  12,  13,  13,  13,  13,  13,  13,  13,   /*par_flag_c[]*/
      1,   5,   9,   9,   9,   6,   5,   9,  10,  10,   9,   9,   9,   9,   9,   9,   6,   8,   9,   9,  10,   /*gt1_flag[]*/
      1,   5,   8,   8,   9,   6,   6,   9,   8,   8,   9,   /*gt1_flag_c[]*/
      9,   5,  10,  13,  13,  10,   9,  10,  13,  13,  13,   9,  10,  10,  10,  13,   8,   9,  10,  10,  13,   /*gt0_flag[]*/
      8,   8,   9,  12,  12,  10,   5,   9,   9,   9,  13,   /*gt0_flag_c[]*/
      8,   5,   4,   5,   4,   4,   5,   4,   1,   0,   4,   1,   0,   0,   0,   0,   1,   0,   0,   0,   /*last_x*/
      5,   4,   4,   /*last_x_c*/
      8,   5,   8,   5,   5,   4,   5,   5,   4,   0,   5,   4,   1,   0,   0,   1,   4,   0,   0,   0,   /*last_y*/
      6,   5,   5,   /*last_y_c*/
     12,   /*mvp_idx*/
      5,   /*smvd_flag*/
      0,   /*sao_merge_flag*/
      4,   /*sao_type_idx*/
      9,   9,  10,   /*lfnst_idx*/
      1,   /*plt_flag*/
      5,   /*rotation_flag*/
      9,   /*run_type_flag*/
      9,   6,   9,  10,   5,   /*idx_run_model*/
      0,   9,   5,   /*copy_run_model*/
    DWS, DWS,   /*rdpcm_flag*/
    DWS, DWS,   /*rdpcm_dir*/
      1,   1,   /*transform_skip_flag*/
      8,   0,   9,   0,   /*mts_idx*/
      9,   2,   /*isp_mode*/
      1,   5,   /*sbt_flag*/
     10,   /*sbt_quad_flag*/
      8,   4,   1,   /*sbt_hor_flag*/
     13,   /*sbt_pos_flag*/
    DWS, DWS, DWS, DWS, DWS, DWS, DWS, DWS, DWS, DWS,   /*cross_comp_pred*/
    DWS,   /*chroma_qp_adj_flag*/
    DWS,   /*chroma_qp_adj_idc*/
      0,   5,   0,   0,   4,   /*imv_flag*/
      0,   0,   0,   4,   0,   0,   1,   0,   0,   /*ctb_alf_flag*/
      0,   0,   /*ctb_alf_alternative*/
      0,   /*alf_use_temporal_filt*/
      4,   1,   4,   4,   1,   4,   /*cc_alf_filter_control_flag*/
      1,   /*ciip_flag*/
      1,   5,   8,   /*ibc_flag*/
      1,   1,   0,   /*joint_cb_cr_flag*/
      5,   8,   8,   /*ts_sig_coeff_group*/
     13,  13,   8,   /*ts_sig_flag*/
      6,   /*ts_par_flag*/
    DWS,   1,   1,   1,   1,   /*ts_gtx_flag*/
      4,   2,   1,   6,   /*ts_lrg1_flag*/
      1,   4,   4,   5,   8,   8,   /*ts_residual_sign*/
};

static const uint8_t init_values[4][OVCABAC_NB_CTX] =
{
    {
     18,  27,  15,  18,  28,  45,  26,   7,  23,   /*split_flag*/
     26,  36,  38,  18,  34,  21,   /*split_qt_flag*/
     43,  42,  37,  42,  44,   /*split_hv_flag*/
     28,  29,  28,  29,   /*split12_flag*/
     25,  20,   /*mode_cons_flag*/
     57,  60,  46,   /*skip_flag*/
      6,   /*merge_flag*/
     46,  15,   /*regular_merge_flag*/
     18,   /*merge_idx*/
     25,   /*mmvd_flag*/
     43,   /*mmvd_merge_idx*/
     59,   /*mmvd_step_mvp_idx*/
     40,  35,   /*pred_mode*/
     25,  59,   /*multi_ref_line_idx*/
     44,   /*intra_luma_mpm_flag*/
     13,   6,   /*intra_luma_planar_flag*/
     26,   /*cclm_mode_flag*/
     27,   /*cclm_mode_idx*/
     25,   /*intra_chroma_pred_mode*/
     56,  57,  50,  26,   /*mip_flag*/
    CNU, CNU,   /*delta_qp*/
     14,  13,   5,   4,   3,  40,   /*inter_dir*/
      5,  35,   /*ref_pic*/
     25,  58,  45,   /*subblock_merge_flag*/
     19,  13,   6,   /*affine_flag*/
     35,   /*affine_type*/
      4,   /*aff_merge_idx*/
      5,   /*bcw_idx*/
     51,  36,   /*mvd*/
     19,  21,   0,  28,   /*bdpcm_mode*/
     12,   /*qt_root_cbf*/
     46,   /*act_flag*/
     15,   6,   5,  14,   /*qt_cbf[]*/
     25,  37,   /*qt_cbf[]*/
      9,  36,  45,   /*qt_cbf[]*/
     25,  45,   /*sig_coeff_group[]*/
     25,  14,   /*sig_coeff_group[]*/
     17,  41,  49,  36,   1,  49,  50,  37,  48,  51,  58,  45,   /*sig_flag[]*/
     26,  45,  53,  46,  49,  54,  61,  39,  35,  39,  39,  39,   /*sig_flag[]*/
     19,  54,  39,  39,  50,  39,  39,  39,   0,  39,  39,  39,   /*sig_flag[]*/
      9,  49,  50,  36,  48,  59,  59,  38,   /*sig_flag_c[]*/
     34,  45,  38,  31,  58,  39,  39,  39,   /*sig_flag_c[]*/
     34,  38,  54,  39,  41,  39,  39,  39,   /*sig_flag_c[]*/
     33,  40,  25,  41,  26,  42,  25,  33,  26,  34,  27,  25,  41,  42,  42,  35,  33,  27,  35,  42,  43,   /*par_flag[]*/
     33,  25,  26,  34,  19,  27,  33,  42,  43,  35,  43,   /*par_flag_c[]*/
     25,   0,   0,  17,  25,  26,   0,   9,  25,  33,  19,   0,  25,  33,  26,  20,  25,  33,  27,  35,  22,   /*gt1_flag[]*/
     25,   1,  25,  33,  26,  12,  25,  33,  27,  28,  37,   /*gt1_flag_c[]*/
      0,   0,  33,  34,  35,  21,  25,  34,  35,  28,  29,  40,  42,  43,  29,  30,  49,  36,  37,  45,  38,   /*gt0_flag[]*/
      0,  40,  34,  43,  36,  37,  57,  52,  45,  38,  46,   /*gt0_flag_c[]*/
      6,   6,  12,  14,   6,   4,  14,   7,   6,   4,  29,   7,   6,   6,  12,  28,   7,  13,  13,  35,   /*last_x*/
     19,   5,   4,   /*last_x_c*/
      5,   5,  20,  13,  13,  19,  21,   6,  12,  12,  14,  14,   5,   4,  12,  13,   7,  13,  12,  41,   /*last_y*/
     11,   5,  27,   /*last_y_c*/
     34,   /*mvp_idx*/
     28,   /*smvd_flag*/
      2,   /*sao_merge_flag*/
      2,   /*sao_type_idx*/
     52,  37,  27,   /*lfnst_idx*/
     17,   /*plt_flag*/
     35,   /*rotation_flag*/
     50,   /*run_type_flag*/
     58,  45,  45,  30,  38,   /*idx_run_model*/
     45,  38,  46,   /*copy_run_model*/
    CNU, CNU,   /*rdpcm_flag*/
    CNU, CNU,   /*rdpcm_dir*/
     25,  17,   /*transform_skip_flag*/
     45,  25,  27,   0,   /*mts_idx*/
     33,  43,   /*isp_mode*/
     41,  57,   /*sbt_flag*/
     42,   /*sbt_quad_flag*/
     35,  51,  27,   /*sbt_hor_flag*/
     28,   /*sbt_pos_flag*/
    CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,   /*cross_comp_pred*/
    CNU,   /*chroma_qp_adj_flag*/
    CNU,   /*chroma_qp_adj_idc*/
     59,  26,  50,  60,  38,   /*imv_flag*/
     33,  52,  46,  25,  61,  54,  25,  61,  54,   /*ctb_alf_flag*/
     11,  26,   /*ctb_alf_alternative*/
     46,   /*alf_use_temporal_filt*/
     25,  35,  38,  25,  28,  38,   /*cc_alf_filter_control_flag*/
     57,   /*ciip_flag*/
      0,  43,  45,   /*ibc_flag*/
     42,  43,  52,   /*joint_cb_cr_flag*/
     18,  35,  45,   /*ts_sig_coeff_group*/
     25,  50,  37,   /*ts_sig_flag*/
     11,   /*ts_par_flag*/
    CNU,   3,   4,   4,   5,   /*ts_gtx_flag*/
     19,  11,   4,   6,   /*ts_lrg1_flag*/
     35,  25,  46,  28,  33,  38,   /*ts_residual_sign*/
    },
    {
     11,  35,  53,  12,   6,  30,  13,  15,  31,   /*split_flag*/
     20,  14,  23,  18,  19,   6,   /*split_qt_flag*/
     43,  35,  37,  34,  52,   /*split_hv_flag*/
     43,  37,  21,  22,   /*split12_flag*/
     25,  12,   /*mode_cons_flag*/
     57,  59,  45,   /*skip_flag*/
     21,   /*merge_flag*/
     38,   7,   /*regular_merge_flag*/
     20,   /*merge_idx*/
     26,   /*mmvd_flag*/
     43,   /*mmvd_merge_idx*/
     60,   /*mmvd_step_mvp_idx*/
     40,  35,   /*pred_mode*/
     25,  58,   /*multi_ref_line_idx*/
     36,   /*intra_luma_mpm_flag*/
     12,  20,   /*intra_luma_planar_flag*/
     34,   /*cclm_mode_flag*/
     27,   /*cclm_mode_idx*/
     25,   /*intra_chroma_pred_mode*/
     41,  57,  58,  26,   /*mip_flag*/
    CNU, CNU,   /*delta_qp*/
      7,   6,   5,  12,   4,  40,   /*inter_dir*/
     20,  35,   /*ref_pic*/
     48,  57,  44,   /*subblock_merge_flag*/
     12,  13,  14,   /*affine_flag*/
     35,   /*affine_type*/
      5,   /*aff_merge_idx*/
      4,   /*bcw_idx*/
     44,  43,   /*mvd*/
     40,  36,   0,  13,   /*bdpcm_mode*/
      5,   /*qt_root_cbf*/
     46,   /*act_flag*/
     23,   5,  20,   7,   /*qt_cbf[]*/
     25,  28,   /*qt_cbf[]*/
     25,  29,  45,   /*qt_cbf[]*/
     25,  30,   /*sig_coeff_group[]*/
     25,  45,   /*sig_coeff_group[]*/
     17,  41,  42,  29,  25,  49,  43,  37,  33,  58,  51,  30,   /*sig_flag[]*/
     19,  38,  38,  46,  34,  54,  54,  39,   6,  39,  39,  39,   /*sig_flag[]*/
     19,  39,  54,  39,  19,  39,  39,  39,  56,  39,  39,  39,   /*sig_flag[]*/
     17,  34,  35,  21,  41,  59,  60,  38,   /*sig_flag_c[]*/
     35,  45,  53,  54,  44,  39,  39,  39,   /*sig_flag_c[]*/
     34,  38,  62,  39,  26,  39,  39,  39,   /*sig_flag_c[]*/
     18,  17,  33,  18,  26,  42,  25,  33,  26,  42,  27,  25,  34,  42,  42,  35,  26,  27,  42,  20,  20,   /*par_flag[]*/
     25,  25,  26,  11,  19,  27,  33,  42,  35,  35,  43,   /*par_flag_c[]*/
     17,   0,   1,  17,  25,  18,   0,   9,  25,  33,  34,   9,  25,  18,  26,  20,  25,  18,  19,  27,  29,   /*gt1_flag[]*/
     17,   9,  25,  10,  18,   4,  17,  33,  19,  20,  29,   /*gt1_flag_c[]*/
      0,  17,  26,  19,  35,  21,  25,  34,  20,  28,  29,  33,  27,  28,  29,  22,  34,  28,  44,  37,  38,   /*gt0_flag[]*/
      0,  25,  19,  20,  13,  14,  57,  44,  30,  30,  23,   /*gt0_flag_c[]*/
      6,  13,  12,   6,   6,  12,  14,  14,  13,  12,  29,   7,   6,  13,  36,  28,  14,  13,   5,  26,   /*last_x*/
     12,   4,  18,   /*last_x_c*/
      5,   5,  12,   6,   6,   4,   6,  14,   5,  12,  14,   7,  13,   5,  13,  21,  14,  20,  12,  34,   /*last_y*/
     11,   4,  18,   /*last_y_c*/
     34,   /*mvp_idx*/
     28,   /*smvd_flag*/
     60,   /*sao_merge_flag*/
      5,   /*sao_type_idx*/
     37,  45,  27,   /*lfnst_idx*/
      0,   /*plt_flag*/
     42,   /*rotation_flag*/
     59,   /*run_type_flag*/
     51,  30,  30,  38,  23,   /*idx_run_model*/
     38,  53,  46,   /*copy_run_model*/
    CNU, CNU,   /*rdpcm_flag*/
    CNU, CNU,   /*rdpcm_dir*/
     25,   9,   /*transform_skip_flag*/
     45,  40,  27,   0,   /*mts_idx*/
     33,  36,   /*isp_mode*/
     56,  57,   /*sbt_flag*/
     42,   /*sbt_quad_flag*/
     20,  43,  12,   /*sbt_hor_flag*/
     28,   /*sbt_pos_flag*/
    CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,   /*cross_comp_pred*/
    CNU,   /*chroma_qp_adj_flag*/
    CNU,   /*chroma_qp_adj_idc*/
     59,  48,  58,  60,  60,   /*imv_flag*/
     13,  23,  46,   4,  61,  54,  19,  46,  54,   /*ctb_alf_flag*/
     20,  12,   /*ctb_alf_alternative*/
     46,   /*alf_use_temporal_filt*/
     18,  21,  38,  18,  21,  38,   /*cc_alf_filter_control_flag*/
     57,   /*ciip_flag*/
      0,  57,  44,   /*ibc_flag*/
     27,  36,  45,   /*joint_cb_cr_flag*/
     18,  12,  29,   /*ts_sig_coeff_group*/
     40,  35,  44,   /*ts_sig_flag*/
      3,   /*ts_par_flag*/
    CNU,   2,  10,   3,   3,   /*ts_gtx_flag*/
     18,  11,   4,  28,   /*ts_lrg1_flag*/
      5,  10,  53,  43,  25,  46,   /*ts_residual_sign*/
    },
    {
     19,  28,  38,  27,  29,  38,  20,  30,  31,   /*split_flag*/
     27,   6,  15,  25,  19,  37,   /*split_qt_flag*/
     43,  42,  29,  27,  44,   /*split_hv_flag*/
     36,  45,  36,  45,   /*split12_flag*/
    CNU, CNU,   /*mode_cons_flag*/
      0,  26,  28,   /*skip_flag*/
     26,   /*merge_flag*/
    CNU, CNU,   /*regular_merge_flag*/
     34,   /*merge_idx*/
    CNU,   /*mmvd_flag*/
    CNU,   /*mmvd_merge_idx*/
    CNU,   /*mmvd_step_mvp_idx*/
    CNU, CNU,   /*pred_mode*/
     25,  60,   /*multi_ref_line_idx*/
     45,   /*intra_luma_mpm_flag*/
     13,  28,   /*intra_luma_planar_flag*/
     59,   /*cclm_mode_flag*/
     27,   /*cclm_mode_idx*/
     34,   /*intra_chroma_pred_mode*/
     33,  49,  50,  25,   /*mip_flag*/
    CNU, CNU,   /*delta_qp*/
    CNU, CNU, CNU, CNU, CNU, CNU,   /*inter_dir*/
    CNU, CNU,   /*ref_pic*/
    CNU, CNU, CNU,   /*subblock_merge_flag*/
    CNU, CNU, CNU,   /*affine_flag*/
    CNU,   /*affine_type*/
    CNU,   /*aff_merge_idx*/
    CNU,   /*bcw_idx*/
     14,  45,   /*mvd*/
     19,  35,   1,  27,   /*bdpcm_mode*/
      6,   /*qt_root_cbf*/
     52,   /*act_flag*/
     15,  12,   5,   7,   /*qt_cbf[]*/
     12,  21,   /*qt_cbf[]*/
     33,  28,  36,   /*qt_cbf[]*/
     18,  31,   /*sig_coeff_group[]*/
     25,  15,   /*sig_coeff_group[]*/
     25,  19,  28,  14,  25,  20,  29,  30,  19,  37,  30,  38,   /*sig_flag[]*/
     11,  38,  46,  54,  27,  39,  39,  39,  44,  39,  39,  39,   /*sig_flag[]*/
     18,  39,  39,  39,  27,  39,  39,  39,   0,  39,  39,  39,   /*sig_flag[]*/
     25,  27,  28,  37,  34,  53,  53,  46,   /*sig_flag_c[]*/
     19,  46,  38,  39,  52,  39,  39,  39,   /*sig_flag_c[]*/
     11,  39,  39,  39,  19,  39,  39,  39,   /*sig_flag_c[]*/
     33,  25,  18,  26,  34,  27,  25,  26,  19,  42,  35,  33,  19,  27,  35,  35,  34,  42,  20,  43,  20,   /*par_flag[]*/
     33,  25,  26,  42,  19,  27,  26,  50,  35,  20,  43,   /*par_flag_c[]*/
     25,   1,  40,  25,  33,  11,  17,  25,  25,  18,   4,  17,  33,  26,  19,  13,  33,  19,  20,  28,  22,   /*gt1_flag[]*/
     40,   9,  25,  18,  26,  35,  25,  26,  35,  28,  37,   /*gt1_flag_c[]*/
     25,  25,  11,  27,  20,  21,  33,  12,  28,  21,  22,  34,  28,  29,  29,  30,  36,  29,  45,  30,  23,   /*gt0_flag[]*/
     40,  33,  27,  28,  21,  37,  36,  37,  45,  38,  46,   /*gt0_flag_c[]*/
     13,   5,   4,  21,  14,   4,   6,  14,  21,  11,  14,   7,  14,   5,  11,  21,  30,  22,  13,  42,   /*last_x*/
     12,   4,   3,   /*last_x_c*/
     13,   5,   4,   6,  13,  11,  14,   6,   5,   3,  14,  22,   6,   4,   3,   6,  22,  29,  20,  34,   /*last_y*/
     12,   4,   3,   /*last_y_c*/
     42,   /*mvp_idx*/
    CNU,   /*smvd_flag*/
     60,   /*sao_merge_flag*/
     13,   /*sao_type_idx*/
     28,  52,  42,   /*lfnst_idx*/
     25,   /*plt_flag*/
     42,   /*rotation_flag*/
     42,   /*run_type_flag*/
     50,  37,  45,  30,  46,   /*idx_run_model*/
     45,  38,  46,   /*copy_run_model*/
    CNU, CNU,   /*rdpcm_flag*/
    CNU, CNU,   /*rdpcm_dir*/
     25,   9,   /*transform_skip_flag*/
     29,   0,  28,   0,   /*mts_idx*/
     33,  43,   /*isp_mode*/
    CNU, CNU,   /*sbt_flag*/
    CNU,   /*sbt_quad_flag*/
    CNU, CNU, CNU,   /*sbt_hor_flag*/
    CNU,   /*sbt_pos_flag*/
    CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU, CNU,   /*cross_comp_pred*/
    CNU,   /*chroma_qp_adj_flag*/
    CNU,   /*chroma_qp_adj_idc*/
    CNU,  34, CNU, CNU, CNU,   /*imv_flag*/
     62,  39,  39,  54,  39,  39,  31,  39,  39,   /*ctb_alf_flag*/
     11,  11,   /*ctb_alf_alternative*/
     46,   /*alf_use_temporal_filt*/
     18,  30,  31,  18,  30,  31,   /*cc_alf_filter_control_flag*/
    CNU,   /*ciip_flag*/
     17,  42,  36,   /*ibc_flag*/
     12,  21,  35,   /*joint_cb_cr_flag*/
     18,  20,  38,   /*ts_sig_coeff_group*/
     25,  28,  38,   /*ts_sig_flag*/
     11,   /*ts_par_flag*/
    CNU,  10,   3,   3,   3,   /*ts_gtx_flag*/
     11,   5,   5,  14,   /*ts_lrg1_flag*/
     12,  17,  46,  28,  25,  46,   /*ts_residual_sign*/
    },
    {
     12,  13,   8,   8,  13,  12,   5,   9,   9,   /*split_flag*/
      0,   8,   8,  12,  12,   8,   /*split_qt_flag*/
      9,   8,   9,   8,   5,   /*split_hv_flag*/
     12,  13,  12,  13,   /*split12_flag*/
      1,   0,   /*mode_cons_flag*/
      5,   4,   8,   /*skip_flag*/
      4,   /*merge_flag*/
      5,   5,   /*regular_merge_flag*/
      4,   /*merge_idx*/
      4,   /*mmvd_flag*/
     10,   /*mmvd_merge_idx*/
      0,   /*mmvd_step_mvp_idx*/
      5,   1,   /*pred_mode*/
      5,   8,   /*multi_ref_line_idx*/
      6,   /*intra_luma_mpm_flag*/
      1,   5,   /*intra_luma_planar_flag*/
      4,   /*cclm_mode_flag*/
      9,   /*cclm_mode_idx*/
      5,   /*intra_chroma_pred_mode*/
      9,  10,   9,   6,   /*mip_flag*/
    DWS, DWS,   /*delta_qp*/
      0,   0,   1,   4,   4,   0,   /*inter_dir*/
      0,   4,   /*ref_pic*/
      4,   4,   4,   /*subblock_merge_flag*/
      4,   0,   0,   /*affine_flag*/
      4,   /*affine_type*/
      0,   /*aff_merge_idx*/
      1,   /*bcw_idx*/
      9,   5,   /*mvd*/
      1,   4,   1,   0,   /*bdpcm_mode*/
      4,   /*qt_root_cbf*/
      1,   /*act_flag*/
      5,   1,   8,   9,   /*qt_cbf[]*/
      5,   0,   /*qt_cbf[]*/
      2,   1,   0,   /*qt_cbf[]*/
      8,   5,   /*sig_coeff_group[]*/
      5,   8,   /*sig_coeff_group[]*/
     12,   9,   9,  10,   9,   9,   9,  10,   8,   8,   8,  10,   /*sig_flag[]*/
      9,  13,   8,   8,   8,   8,   8,   5,   8,   0,   0,   0,   /*sig_flag[]*/
      8,   8,   8,   8,   8,   0,   4,   4,   0,   0,   0,   0,   /*sig_flag[]*/
     12,  12,   9,  13,   4,   5,   8,   9,   /*sig_flag_c[]*/
      8,  12,  12,   8,   4,   0,   0,   0,   /*sig_flag_c[]*/
      8,   8,   8,   8,   4,   0,   0,   0,   /*sig_flag_c[]*/
      8,   9,  12,  13,  13,  13,  10,  13,  13,  13,  13,  13,  13,  13,  13,  13,  10,  13,  13,  13,  13,   /*par_flag[]*/
      8,  12,  12,  12,  13,  13,  13,  13,  13,  13,  13,   /*par_flag_c[]*/
      1,   5,   9,   9,   9,   6,   5,   9,  10,  10,   9,   9,   9,   9,   9,   9,   6,   8,   9,   9,  10,   /*gt1_flag[]*/
      1,   5,   8,   8,   9,   6,   6,   9,   8,   8,   9,   /*gt1_flag_c[]*/
      9,   5,  10,  13,  13,  10,   9,  10,  13,  13,  13,   9,  10,  10,  10,  13,   8,   9,  10,  10,  13,   /*gt0_flag[]*/
      8,   8,   9,  12,  12,  10,   5,   9,   9,   9,  13,   /*gt0_flag_c[]*/
      8,   5,   4,   5,   4,   4,   5,   4,   1,   0,   4,   1,   0,   0,   0,   0,   1,   0,   0,   0,   /*last_x*/
      5,   4,   4,   /*last_x_c*/
      8,   5,   8,   5,   5,   4,   5,   5,   4,   0,   5,   4,   1,   0,   0,   1,   4,   0,   0,   0,   /*last_y*/
      6,   5,   5,   /*last_y_c*/
     12,   /*mvp_idx*/
      5,   /*smvd_flag*/
      0,   /*sao_merge_flag*/
      4,   /*sao_type_idx*/
      9,   9,  10,   /*lfnst_idx*/
      1,   /*plt_flag*/
      5,   /*rotation_flag*/
      9,   /*run_type_flag*/
      9,   6,   9,  10,   5,   /*idx_run_model*/
      0,   9,   5,   /*copy_run_model*/
    DWS, DWS,   /*rdpcm_flag*/
    DWS, DWS,   /*rdpcm_dir*/
      1,   1,   /*transform_skip_flag*/
      8,   0,   9,   0,   /*mts_idx*/
      9,   2,   /*isp_mode*/
      1,   5,   /*sbt_flag*/
     10,   /*sbt_quad_flag*/
      8,   4,   1,   /*sbt_hor_flag*/
     13,   /*sbt_pos_flag*/
    DWS, DWS, DWS, DWS, DWS, DWS, DWS, DWS, DWS, DWS,   /*cross_comp_pred*/
    DWS,   /*chroma_qp_adj_flag*/
    DWS,   /*chroma_qp_adj_idc*/
      0,   5,   0,   0,   4,   /*imv_flag*/
      0,   0,   0,   4,   0,   0,   1,   0,   0,   /*ctb_alf_flag*/
      0,   0,   /*ctb_alf_alternative*/
      0,   /*alf_use_temporal_filt*/
      4,   1,   4,   4,   1,   4,   /*cc_alf_filter_control_flag*/
      1,   /*ciip_flag*/
      1,   5,   8,   /*ibc_flag*/
      1,   1,   0,   /*joint_cb_cr_flag*/
      5,   8,   8,   /*ts_sig_coeff_group*/
     13,  13,   8,   /*ts_sig_flag*/
      6,   /*ts_par_flag*/
    DWS,   1,   1,   1,   1,   /*ts_gtx_flag*/
      4,   2,   1,   6,   /*ts_lrg1_flag*/
      1,   4,   4,   5,   8,   8,   /*ts_residual_sign*/
    },
};

/* FIXME find a place for these two LUT if we inline cabac read functions
 * in a header file */
/* Look-Up Table used for Least probable Symbol
 * renormalisation process
 */
const uint8_t lps_renorm_table[64] =
{
  6,  5,  4,  4,
  3,  3,  3,  3,
  2,  2,  2,  2,
  2,  2,  2,  2,
  1,  1,  1,  1,
  1,  1,  1,  1,
  1,  1,  1,  1,
  1,  1,  1,  1,
  0,  0,  0,  0,
  0,  0,  0,  0,
  0,  0,  0,  0,
  0,  0,  0,  0,
  0,  0,  0,  0,
  0,  0,  0,  0,
  0,  0,  0,  0,
  0,  0,  0,  0,
};

/* Look-Up Table used avoid multiplication
 * in least probable symbol range computation
 */
const uint8_t range_lps_lut[512] =
{
    0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04,
    0x04, 0x04, 0x05, 0x05, 0x06, 0x06, 0x07, 0x07, 0x08, 0x08, 0x09, 0x09, 0x0a, 0x0a, 0x0b, 0x0b, 0x0c, 0x0c, 0x0d, 0x0d, 0x0e, 0x0e, 0x0f, 0x0f, 0x10, 0x10, 0x11, 0x11, 0x12, 0x12, 0x13, 0x13,
    0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x22, 0x23,
    0x04, 0x05, 0x07, 0x08, 0x0a, 0x0b, 0x0d, 0x0e, 0x10, 0x11, 0x13, 0x14, 0x16, 0x17, 0x19, 0x1a, 0x1c, 0x1d, 0x1f, 0x20, 0x22, 0x23, 0x25, 0x26, 0x28, 0x29, 0x2b, 0x2c, 0x2e, 0x2f, 0x31, 0x32,
    0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e, 0x10, 0x12, 0x14, 0x16, 0x18, 0x1a, 0x1c, 0x1e, 0x20, 0x22, 0x24, 0x26, 0x28, 0x2a, 0x2c, 0x2e, 0x30, 0x32, 0x34, 0x36, 0x38, 0x3a, 0x3c, 0x3e, 0x40, 0x42,
    0x04, 0x06, 0x09, 0x0b, 0x0e, 0x10, 0x13, 0x15, 0x18, 0x1a, 0x1d, 0x1f, 0x22, 0x24, 0x27, 0x29, 0x2c, 0x2e, 0x31, 0x33, 0x36, 0x38, 0x3b, 0x3d, 0x40, 0x42, 0x45, 0x47, 0x4a, 0x4c, 0x4f, 0x51,
    0x04, 0x07, 0x0a, 0x0d, 0x10, 0x13, 0x16, 0x19, 0x1c, 0x1f, 0x22, 0x25, 0x28, 0x2b, 0x2e, 0x31, 0x34, 0x37, 0x3a, 0x3d, 0x40, 0x43, 0x46, 0x49, 0x4c, 0x4f, 0x52, 0x55, 0x58, 0x5b, 0x5e, 0x61,
    0x04, 0x07, 0x0b, 0x0e, 0x12, 0x15, 0x19, 0x1c, 0x20, 0x23, 0x27, 0x2a, 0x2e, 0x31, 0x35, 0x38, 0x3c, 0x3f, 0x43, 0x46, 0x4a, 0x4d, 0x51, 0x54, 0x58, 0x5b, 0x5f, 0x62, 0x66, 0x69, 0x6d, 0x70,
    0x04, 0x08, 0x0c, 0x10, 0x14, 0x18, 0x1c, 0x20, 0x24, 0x28, 0x2c, 0x30, 0x34, 0x38, 0x3c, 0x40, 0x44, 0x48, 0x4c, 0x50, 0x54, 0x58, 0x5c, 0x60, 0x64, 0x68, 0x6c, 0x70, 0x74, 0x78, 0x7c, 0x80,
    0x04, 0x08, 0x0d, 0x11, 0x16, 0x1a, 0x1f, 0x23, 0x28, 0x2c, 0x31, 0x35, 0x3a, 0x3e, 0x43, 0x47, 0x4c, 0x50, 0x55, 0x59, 0x5e, 0x62, 0x67, 0x6b, 0x70, 0x74, 0x79, 0x7d, 0x82, 0x86, 0x8b, 0x8f,
    0x04, 0x09, 0x0e, 0x13, 0x18, 0x1d, 0x22, 0x27, 0x2c, 0x31, 0x36, 0x3b, 0x40, 0x45, 0x4a, 0x4f, 0x54, 0x59, 0x5e, 0x63, 0x68, 0x6d, 0x72, 0x77, 0x7c, 0x81, 0x86, 0x8b, 0x90, 0x95, 0x9a, 0x9f,
    0x04, 0x09, 0x0f, 0x14, 0x1a, 0x1f, 0x25, 0x2a, 0x30, 0x35, 0x3b, 0x40, 0x46, 0x4b, 0x51, 0x56, 0x5c, 0x61, 0x67, 0x6c, 0x72, 0x77, 0x7d, 0x82, 0x88, 0x8d, 0x93, 0x98, 0x9e, 0xa3, 0xa9, 0xae,
    0x04, 0x0a, 0x10, 0x16, 0x1c, 0x22, 0x28, 0x2e, 0x34, 0x3a, 0x40, 0x46, 0x4c, 0x52, 0x58, 0x5e, 0x64, 0x6a, 0x70, 0x76, 0x7c, 0x82, 0x88, 0x8e, 0x94, 0x9a, 0xa0, 0xa6, 0xac, 0xb2, 0xb8, 0xbe,
    0x04, 0x0a, 0x11, 0x17, 0x1e, 0x24, 0x2b, 0x31, 0x38, 0x3e, 0x45, 0x4b, 0x52, 0x58, 0x5f, 0x65, 0x6c, 0x72, 0x79, 0x7f, 0x86, 0x8c, 0x93, 0x99, 0xa0, 0xa6, 0xad, 0xb3, 0xba, 0xc0, 0xc7, 0xcd,
    0x04, 0x0b, 0x12, 0x19, 0x20, 0x27, 0x2e, 0x35, 0x3c, 0x43, 0x4a, 0x51, 0x58, 0x5f, 0x66, 0x6d, 0x74, 0x7b, 0x82, 0x89, 0x90, 0x97, 0x9e, 0xa5, 0xac, 0xb3, 0xba, 0xc1, 0xc8, 0xcf, 0xd6, 0xdd,
    0x04, 0x0b, 0x13, 0x1a, 0x22, 0x29, 0x31, 0x38, 0x40, 0x47, 0x4f, 0x56, 0x5e, 0x65, 0x6d, 0x74, 0x7c, 0x83, 0x8b, 0x92, 0x9a, 0xa1, 0xa9, 0xb0, 0xb8, 0xbf, 0xc7, 0xce, 0xd6, 0xdd, 0xe5, 0xec,
};

int
ovcabac_attach_entry(OVCABACCtx *const cabac_ctx, const uint8_t *const entry_point,
                     const uint8_t *const entry_end)
{
    cabac_ctx->bytestream_start = cabac_ctx->bytestream = entry_point;
    cabac_ctx->bytestream_end   = entry_end;

    cabac_ctx->low_b  =  (int)cabac_ctx->bytestream[0] << 18;
    cabac_ctx->low_b +=  (int)cabac_ctx->bytestream[1] << 10;
    cabac_ctx->bytestream += 2;

    if(((uintptr_t)cabac_ctx->bytestream & 1) == 0) {
        cabac_ctx->low_b |= 1 << 9;
    } else {
        cabac_ctx->low_b +=  (*cabac_ctx->bytestream++) << 2;
        cabac_ctx->low_b += 2;
    }

    cabac_ctx->range = 0x1FE;

    if ((cabac_ctx->range << (NB_CABAC_BITS + 1)) < cabac_ctx->low_b)
        return OV_ERROR;

    return 0;
}


void
ovcabac_init_slice_context_table(uint64_t cabac_states[], uint8_t slice_type, 
                                 int slice_qp)
{
    int i;
    uint64_t cabac_state;
    uint32_t ctx_table;
    uint32_t rates;
    uint16_t rate_0, rate_1;
    uint16_t state_0, state_1;
    const uint8_t *ctx_init_tab = init_values[slice_type];

    slice_qp = ov_clip(slice_qp, 0, 63);

    for (i = 0; i < OVCABAC_NB_CTX; i++) {
        int init_value    = ctx_init_tab[i];
        int log2_window_s = rate_init_table[i];
        int slope  =  (init_value >> 3) - 4;
        int offset = ((init_value & 0x7) * 18) + 1;
        int init_state  = ((slope * (slice_qp - 16)) >> 1) + offset;

        init_state = init_state < 1 ? 1 : init_state > 127 ? 127 : init_state;

        init_state <<= 8;

        state_0 = init_state & 0x7FE0;
        state_1 = init_state & 0x7FFE;

        /*FIXME we could set those values in the init tab instead
          of deriving them*/
        rate_0 = 2 + ((log2_window_s >> 2)   & 0x3);
        rate_1 = 3 + rate_0 + (log2_window_s & 0x3);

        ctx_table = ((uint32_t)state_0 << 16) + state_1;
        rates  = ((uint32_t)rate_0  << 16) + rate_1 ;

        cabac_state = ((uint64_t)ctx_table << 32) + rates;

        cabac_states[i] = cabac_state;
    }
}

/*FIXME call only on slice end
,*/
int
ovcabac_end_of_slice(OVCABACCtx *cabac_ctx)
{
    cabac_ctx->range -= 2;
    if (cabac_ctx->low_b < cabac_ctx->range << (NB_CABAC_BITS + 1)){
        uint8_t log2_renorm = lps_renorm_table[cabac_ctx->range >> 3];

        cabac_ctx->low_b <<= log2_renorm;
        cabac_ctx->range <<= log2_renorm;

        if (!(cabac_ctx->low_b & CABAC_MASK)){
            int num_bits = 0;
            int tmp_fill = -CABAC_MASK;

            tmp_fill += cabac_ctx->bytestream[0] << 9;
            tmp_fill += cabac_ctx->bytestream[1] << 1;
            cabac_ctx->low_b += tmp_fill << num_bits;

            if (cabac_ctx->bytestream < cabac_ctx->bytestream_end){
                cabac_ctx->bytestream += NB_CABAC_BITS >> 3;
            }
        }
        return 0;
    } else {
        /*FIXME use this result to check if end of stream is reached*/
        return cabac_ctx->bytestream - cabac_ctx->bytestream_start;
    }
}
