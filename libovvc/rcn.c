#include <stdint.h>
#include <string.h>
#include "nvcl_utils.h"
#include "ovutils.h"
#include "rcn_fill_ref.h"
#include "rcn_intra_angular.h"
#include "rcn_intra_dc_planar.h"
#include "data_rcn_angular.h"
#include "ctudec.h"
#include "rcn_intra_mip.h"
#include "rcn_structures.h"
#include "rcn_mc.h"
#include "rcn.h"
#include "rcn_lmcs.h"
#include "ovmem.h"
#include "ovconfig.h"
#include "drv.h"

#if ARCH_X86
  #if SSE_ENABLED
    #include "x86/rcn_sse.h"
  #elif AVX_ENABLED
    //Link AVX optims
  #else
    //Failover x86
  #endif
#elif ARCH_ARM
  #if NEON_ENABLED
    #include "arm/rcn_neon.h"
  #else
    //Failover ARM
  #endif
#else
  //Failover other arch
#endif

#include "dbf_utils.h"

#define BITDEPTH 10

void
rcn_residual(OVCTUDec *const ctudec,
             int16_t *const dst, int16_t *src,
             uint8_t x0, uint8_t y0,
             unsigned int log2_tb_w, unsigned int log2_tb_h,
             unsigned int lim_cg_w,
             uint8_t cu_mts_flag, uint8_t cu_mts_idx,
             uint8_t is_dc, uint8_t lfnst_flag, uint8_t is_mip, uint8_t lfnst_idx, uint8_t sbt)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;
    fill_bs_map(&ctudec->dbf_info.bs1_map, x0, y0, log2_tb_w, log2_tb_h);
    int shift_v = 6 + 1;
    int shift_h = (6 + 15 - 1) - BITDEPTH;

    DECLARE_ALIGNED(32, int16_t, tmp)[64*64];

    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag) {
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
        process_lfnst_luma(ctudec, src, ctudec->lfnst_subblock, log2_tb_w, log2_tb_h, x0, y0,
                           lfnst_idx);
        lim_cg_w = 8;
        is_dc = 0;
    }

    if (!is_mip && !cu_mts_flag && ctudec->mts_implicit && (log2_tb_w <= 4 || log2_tb_h <= 4) && !lfnst_flag) {
        /*FIXME condition on size in the if could be removed ?*/
        enum DCTType tr_h_idx = log2_tb_w <= 4 ? DST_VII : DCT_II;
        enum DCTType tr_v_idx = log2_tb_h <= 4 ? DST_VII : DCT_II;

        /* FIXME use coefficient zeroing in MTS */
        TRFunc->func[tr_v_idx][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, tb_w, shift_h);

    } else if (!cu_mts_flag) {

        if (is_dc) {

            TRFunc->dc(dst, log2_tb_w, log2_tb_h, src[0]);

        } else {
            int nb_row = OVMIN(lim_cg_w, 1 << log2_tb_w);
            int nb_col = OVMIN(lim_cg_w, 1 << log2_tb_h);

            TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, nb_row, nb_col, shift_v);
            TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, shift_h);
        }
    } else {
        enum DCTType tr_h_idx = cu_mts_idx  & 1;
        enum DCTType tr_v_idx = cu_mts_idx >> 1;

        TRFunc->func[tr_v_idx][log2_tb_h](src, tmp, tb_w, tb_w, tb_h, shift_v);
        TRFunc->func[tr_h_idx][log2_tb_w](tmp, dst, tb_h, tb_h, tb_w, shift_h);
    }
}

void
rcn_residual_c(OVCTUDec *const ctudec,
               int16_t *const dst, int16_t *src,
               uint8_t x0, uint8_t y0,
               uint8_t log2_tb_w, uint8_t log2_tb_h,
               uint16_t last_pos,
               uint8_t lfnst_flag, uint8_t lfnst_idx)
{
    struct TRFunctions *TRFunc = &ctudec->rcn_ctx.rcn_funcs.tr;

    const int shift_v = 6 + 1;
    const int shift_h = (6 + 15 - 1) - BITDEPTH;

    DECLARE_ALIGNED(32, int16_t, tmp)[32*32];

    int tb_w = 1 << log2_tb_w;
    int tb_h = 1 << log2_tb_h;

    memset(tmp, 0, sizeof(int16_t) << (log2_tb_w + log2_tb_h));

    if (lfnst_flag && log2_tb_w > 1 && log2_tb_h > 1) {
        /* FIXME separate lfnst mode derivation from lfnst reconstruction */
        int16_t lfnst_sb[16];
        memcpy(lfnst_sb     , &src[0], sizeof(int16_t) * 4);
        memcpy(lfnst_sb +  4, &src[1 << log2_tb_w], sizeof(int16_t) * 4);
        memcpy(lfnst_sb +  8, &src[2 << log2_tb_w], sizeof(int16_t) * 4);
        memcpy(lfnst_sb + 12, &src[3 << log2_tb_w], sizeof(int16_t) * 4);

        process_lfnst(ctudec, src, lfnst_sb, log2_tb_w, log2_tb_h,
                      x0, y0, lfnst_idx);
    }

    if (!last_pos && !lfnst_flag) {

        TRFunc->dc(dst, log2_tb_w, log2_tb_h, src[0]);

    } else {
        int lim_sb_s = ((((last_pos >> 8)) >> 2) + (((last_pos & 0xFF))>> 2) + 1) << 2;
        if (lfnst_flag && log2_tb_w > 1 && log2_tb_h > 1) lim_sb_s = 8;
        int nb_row =  OVMIN(lim_sb_s, 1 << log2_tb_w);
        int nb_col =  OVMIN(lim_sb_s, 1 << log2_tb_h);
        /*FIXME might be transform SKIP */

        TRFunc->func[DCT_II][log2_tb_h](src, tmp, tb_w, nb_row, nb_col, shift_v);
        TRFunc->func[DCT_II][log2_tb_w](tmp, dst, tb_h, tb_h, nb_row, shift_h);
    }
}

static void
rcn_init_intra_angular_functions(struct RCNFunctions *rcn_func, uint8_t bitdepth)
{
    switch (bitdepth) {
        case 8:
        {
            extern const struct IntraAngularFunctions angular_gauss_h_8;
            extern const struct IntraAngularFunctions angular_gauss_v_8;
            extern const struct IntraAngularFunctions angular_cubic_h_8;
            extern const struct IntraAngularFunctions angular_cubic_v_8;
            extern const struct IntraAngularFunctions angular_c_h_8;
            extern const struct IntraAngularFunctions angular_c_v_8;
            extern const struct IntraMRLFunctions mrl_func_8;
            extern const struct IntraAngularFunctions angular_nofrac_v_8;
            extern const struct IntraAngularFunctions angular_nofrac_h_8;

            rcn_func->intra_angular_gauss_h = &angular_gauss_h_8;
            rcn_func->intra_angular_gauss_v = &angular_gauss_v_8;
            rcn_func->intra_angular_cubic_h = &angular_cubic_h_8;
            rcn_func->intra_angular_cubic_v = &angular_cubic_v_8;
            rcn_func->intra_angular_c_h     = &angular_c_h_8;
            rcn_func->intra_angular_c_v     = &angular_c_v_8;
            rcn_func->intra_angular_nofrac_v = &angular_nofrac_v_8;
            rcn_func->intra_angular_nofrac_h = &angular_nofrac_h_8;
            rcn_func->intra_mrl = &mrl_func_8;
        }
            break;
        case 10:
        {
            extern const struct IntraAngularFunctions angular_gauss_h_10;
            extern const struct IntraAngularFunctions angular_gauss_v_10;
            extern const struct IntraAngularFunctions angular_cubic_h_10;
            extern const struct IntraAngularFunctions angular_cubic_v_10;
            extern const struct IntraAngularFunctions angular_c_h_10;
            extern const struct IntraAngularFunctions angular_c_v_10;
            extern const struct IntraMRLFunctions mrl_func_10;
            extern const struct IntraAngularFunctions angular_nofrac_v_10;
            extern const struct IntraAngularFunctions angular_nofrac_h_10;

            rcn_func->intra_angular_gauss_h = &angular_gauss_h_10;
            rcn_func->intra_angular_gauss_v = &angular_gauss_v_10;
            rcn_func->intra_angular_cubic_h = &angular_cubic_h_10;
            rcn_func->intra_angular_cubic_v = &angular_cubic_v_10;
            rcn_func->intra_angular_c_h     = &angular_c_h_10;
            rcn_func->intra_angular_c_v     = &angular_c_v_10;
            rcn_func->intra_angular_nofrac_v = &angular_nofrac_v_10;
            rcn_func->intra_angular_nofrac_h = &angular_nofrac_h_10;
            rcn_func->intra_mrl = &mrl_func_10;
        }
            break;
        default:
            {
            #if 0
            extern const struct IntraAngularFunctions angular_gauss_h;
            extern const struct IntraAngularFunctions angular_gauss_v;
            extern const struct IntraAngularFunctions angular_cubic_h;
            extern const struct IntraAngularFunctions angular_cubic_v;
            extern const struct IntraAngularFunctions angular_c_h;
            extern const struct IntraAngularFunctions angular_c_v;
            extern const struct IntraMRLFunctions mrl_func;
            extern const struct IntraAngularFunctions angular_nofrac_v;
            extern const struct IntraAngularFunctions angular_nofrac_h;

            rcn_func->intra_angular_gauss_h = &angular_gauss_h;
            rcn_func->intra_angular_gauss_v = &angular_gauss_v;
            rcn_func->intra_angular_cubic_h = &angular_cubic_h;
            rcn_func->intra_angular_cubic_v = &angular_cubic_v;
            rcn_func->intra_angular_c_h     = &angular_c_h;
            rcn_func->intra_angular_c_v     = &angular_c_v;
            rcn_func->intra_angular_nofrac_v = &angular_nofrac_v;
            rcn_func->intra_angular_nofrac_h = &angular_nofrac_h;
            rcn_func->intra_mrl = &mrl_func;
            #endif
            }
            break;
    }
}

void
rcn_init_functions(struct RCNFunctions *rcn_func, uint8_t ict_type, uint8_t lm_chroma_enabled,
                    uint8_t sps_chroma_vertical_collocated_flag, uint8_t lmcs_flag, uint8_t bitdepth)
{
  rcn_init_mc_functions(rcn_func);
  rcn_init_tr_functions(rcn_func);
  rcn_init_dc_planar_functions(rcn_func);
  rcn_init_ict_functions(rcn_func, ict_type);
  rcn_init_lfnst_functions(rcn_func);
  rcn_init_mip_functions(rcn_func);
  rcn_init_alf_functions(rcn_func);
  rcn_init_sao_functions(rcn_func);
  rcn_init_lmcs_function(rcn_func, lmcs_flag);
  rcn_init_dmvr_functions(rcn_func);
  rcn_init_prof_functions(rcn_func);
  rcn_init_bdof_functions(rcn_func);
  rcn_init_ciip_functions(rcn_func);
  rcn_init_df_functions(rcn_func);
  rcn_init_intra_angular_functions(rcn_func, bitdepth);
  if (lm_chroma_enabled) {
      /* FIXME add support vertical */
      if (sps_chroma_vertical_collocated_flag /*sps->sps_chroma_horizontal_collocated_flag*/) {
          rcn_init_cclm_functions_collocated(rcn_func);
      } else {
          rcn_init_cclm_functions(rcn_func);
      }
  }

  #if ARCH_X86
    #if SSE_ENABLED
      rcn_init_mc_functions_sse(rcn_func);
      rcn_init_tr_functions_sse(rcn_func);
      rcn_init_dc_planar_functions_sse(rcn_func);
      rcn_init_ict_functions_sse(rcn_func, ict_type);
      rcn_init_lfnst_functions_sse(rcn_func);
      rcn_init_mip_functions_sse(rcn_func);
      rcn_init_alf_functions_sse(rcn_func);
      rcn_init_sao_functions_sse(rcn_func);
      rcn_init_dmvr_functions_sse(rcn_func);
      rcn_init_prof_functions_sse(rcn_func);
      rcn_init_bdof_functions_sse(rcn_func);
      rcn_init_ciip_functions_sse(rcn_func);
      rcn_init_df_functions_sse(rcn_func);
      if (lm_chroma_enabled) {
          if (!sps_chroma_vertical_collocated_flag /*sps->sps_chroma_horizontal_collocated_flag*/) {
              rcn_init_cclm_functions_sse(rcn_func);
          }
      }
    #elif AVX_ENABLED
      //Link AVX optims
    #else
      //Failover x86
    #endif
  #elif ARCH_ARM
    #if NEON_ENABLED
      // rcn_init_mc_functions_sse(rcn_func);
      // rcn_init_tr_functions_sse(rcn_func);
      // rcn_init_dc_planar_functions_sse(rcn_func);
      // rcn_init_ict_functions_sse(rcn_func, ict_type);
      // rcn_init_lfnst_functions_sse(rcn_func);
      // rcn_init_mip_functions_sse(rcn_func);
      // rcn_init_alf_functions_sse(rcn_func);
      // rcn_init_sao_functions_sse(rcn_func);
      // rcn_init_dmvr_functions_sse(rcn_func);
      // rcn_init_prof_functions_sse(rcn_func);
      // rcn_init_bdof_functions_sse(rcn_func);
      // rcn_init_ciip_functions_sse(rcn_func);
      // rcn_init_df_functions_sse(rcn_func);
      // if (lm_chroma_enabled) {
      //     if (!sps_chroma_vertical_collocated_flag ) {
      //         rcn_init_cclm_functions_sse(rcn_func);
      //     }
      // }

      // to enable Assembly optimisation
      #if ARCH_AARCH64_ASSEMBLY
        rcn_init_mc_functions_neon(rcn_func);
      #endif
    #else
      //Failover ARM
    #endif
  #else
    //Failover other arch
  #endif
}
