#include <stdint.h>
#include <string.h>
#include "nvcl_utils.h"
#include "ovutils.h"
#include "rcn_intra_angular.h"
#include "data_rcn_angular.h"
#include "ctudec.h"
#include "rcn_structures.h"
#include "rcn_mc.h"
#include "rcn.h"
#include "rcn_lmcs.h"
#include "ovmem.h"
#include "ovconfig.h"
#include "drv.h"

#ifndef NO_SIMD
  #if HAVE_X86_OPTIM
    #if HAVE_SSE4_1
      #include "x86/rcn_sse.h"
    #endif
    #if HAVE_AVX2
      #include "x86/rcn_avx2.h"
    #endif
  #elif __ARM_ARCH
    #if __ARM_NEON
      #include "arm/rcn_neon.h"
    #else
      //Failover ARM
    #endif
  #else
    //Failover other arch
  #endif
#endif

#include "dbf_utils.h"

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
  rcn_init_ctu_buffs_10(rcn_func);
  rcn_init_mc_functions_10(rcn_func);
  rcn_init_tr_functions_10(rcn_func);
  rcn_init_dc_planar_functions_10(rcn_func);
  rcn_init_ict_functions_10(rcn_func, ict_type, bitdepth);
  rcn_init_intra_angular_functions(rcn_func, bitdepth);
  rcn_init_lfnst_functions(rcn_func);
  rcn_init_mip_functions_10(rcn_func);
  rcn_init_alf_functions_10(rcn_func);
  rcn_init_sao_functions_10(rcn_func);
  rcn_init_lmcs_function_10(rcn_func, lmcs_flag);
  rcn_init_dmvr_functions_10(rcn_func);
  rcn_init_prof_functions_10(rcn_func);
  rcn_init_bdof_functions_10(rcn_func);
  rcn_init_ciip_functions_10(rcn_func);
  rcn_init_df_functions_10(rcn_func);
  rcn_init_dequant_10(rcn_func);
  rcn_init_fill_ref_10(rcn_func);
  rcn_init_transform_trees_10(rcn_func);
  rcn_init_intra_functions_10(rcn_func);
  rcn_init_inter_functions_10(rcn_func);
  if (lm_chroma_enabled) {
      /* FIXME add support vertical */
      if (sps_chroma_vertical_collocated_flag /*sps->sps_chroma_horizontal_collocated_flag*/) {
          rcn_init_cclm_functions_collocated_10(rcn_func);
      } else {
          rcn_init_cclm_functions_10(rcn_func);
      }
  }
  if (bitdepth == 8) {
      rcn_init_ctu_buffs_8(rcn_func);
      rcn_init_mc_functions_8(rcn_func);
      rcn_init_tr_functions_8(rcn_func);
      rcn_init_dc_planar_functions_8(rcn_func);
      rcn_init_ict_functions_8(rcn_func, ict_type, bitdepth);
      rcn_init_intra_angular_functions(rcn_func, bitdepth);
      rcn_init_lfnst_functions(rcn_func);
      rcn_init_mip_functions_8(rcn_func);
      rcn_init_alf_functions_8(rcn_func);
      rcn_init_sao_functions_8(rcn_func);
      rcn_init_lmcs_function_8(rcn_func, lmcs_flag);
      rcn_init_dmvr_functions_8(rcn_func);
      rcn_init_prof_functions_8(rcn_func);
      rcn_init_bdof_functions_8(rcn_func);
      rcn_init_ciip_functions_8(rcn_func);
      rcn_init_df_functions_8(rcn_func);
      rcn_init_dequant_8(rcn_func);
      rcn_init_fill_ref_8(rcn_func);
      rcn_init_transform_trees_8(rcn_func);
      rcn_init_intra_functions_8(rcn_func);
      rcn_init_inter_functions_8(rcn_func);
      if (lm_chroma_enabled) {
          /* FIXME add support vertical */
          if (sps_chroma_vertical_collocated_flag /*sps->sps_chroma_horizontal_collocated_flag*/) {
              rcn_init_cclm_functions_collocated_8(rcn_func);
          } else {
              rcn_init_cclm_functions_8(rcn_func);
          }
      }
  }

  #ifndef NO_SIMD
    #if HAVE_X86_OPTIM
      #if HAVE_SSE4_1
      if (__builtin_cpu_supports("sse4.1") && bitdepth == 10) {
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
          rcn_init_intra_angular_functions_10_sse(rcn_func);
          rcn_init_dequant_sse(rcn_func);

          if (lm_chroma_enabled) {
            if (!sps_chroma_vertical_collocated_flag /*sps->sps_chroma_horizontal_collocated_flag*/) {
              rcn_init_cclm_functions_sse(rcn_func);
            }
          }
      }
      #endif
      #if HAVE_AVX2
        if (__builtin_cpu_supports("avx2") && bitdepth == 10) {
          rcn_init_alf_functions_avx2(rcn_func);
          rcn_init_sao_functions_avx2(rcn_func);
          rcn_init_ict_functions_avx2(rcn_func, ict_type);
          rcn_init_mip_functions_avx2(rcn_func);
          rcn_init_ciip_functions_avx2(rcn_func);
          rcn_init_mc_functions_avx2(rcn_func);
          rcn_init_dmvr_functions_avx2(rcn_func);
          rcn_init_prof_functions_avx2(rcn_func);
          rcn_init_bdof_functions_avx2(rcn_func);
          rcn_init_intra_angular_functions_10_avx2(rcn_func);
        }
      #endif
    #elif __ARM_ARCH
      #if __ARM_NEON
        #if ARM_SIMDE
          #ifndef EXCLUDE_FOR_CLANG
            rcn_init_dc_planar_functions_sse(rcn_func);
          #endif
          rcn_init_mc_functions_sse(rcn_func);
          rcn_init_tr_functions_sse(rcn_func);
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
          rcn_init_intra_angular_functions_10_sse(rcn_func);
          if (lm_chroma_enabled) {
              if (!sps_chroma_vertical_collocated_flag ) {
                  rcn_init_cclm_functions_sse(rcn_func);
              }
          }
          #endif
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
  #endif
}
