#include "cabac_internal.h"

static uint8_t
ovcabac_read_ae_alf(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
               uint8_t loop_filter_flags,
               uint8_t ctu_neighbour_flags,
               uint8_t  up_ctb_alf_flag, uint8_t left_ctb_alf_flag,
               unsigned int tile_group_num_aps,
               unsigned int num_alf_alternative)
{
    uint8_t ret_luma = 0;
    uint8_t ret_cb = 0;
    uint8_t ret_cr = 0;
    uint8_t ret;

    uint8_t ctx;
    uint8_t alf_idx;

    #if 0
    if(loop_filter_flags & VVC_ALF_LUMA_SLICE_FLAG){
        ctx  = ctu_neighbour_flags & VVC_CTU_LEFT_FLAG ? ((left_ctb_alf_flag & 4) >> 2) : 0;
        ctx += ctu_neighbour_flags & VVC_CTU_UP_FLAG   ? ((up_ctb_alf_flag   & 4) >> 2) : 0;
        ret_luma = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + ctx]);
        if(ret_luma){
            alf_idx = ovcabac_read_ae_alf_idx(cabac_ctx, cabac_state, tile_group_num_aps);
        }
    }

    //if(!(loop_filter_flags & VVC_ALF_CHROMA_CTB_FLAG)){
        if(loop_filter_flags & VVC_ALF_CB_SLICE_FLAG){
            int decoded = 0;

            ctx  = ctu_neighbour_flags & VVC_CTU_LEFT_FLAG ? ((left_ctb_alf_flag & 2) >> 1) : 0;
            ctx += ctu_neighbour_flags & VVC_CTU_UP_FLAG   ? ((up_ctb_alf_flag   & 2) >> 1) : 0;
            ret_cb = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + 3 + ctx]);
            while (ret_cb && decoded < num_alf_alternative - 1 && ovcabac_ae_read(cabac_ctx,
                                                                             &cabac_state[CTB_ALF_ALTERNATIVE_CTX_OFFSET])){
                ++decoded;
            }
        }
        if(loop_filter_flags & VVC_ALF_CR_SLICE_FLAG){
            int decoded = 0;
            ctx  = ctu_neighbour_flags & VVC_CTU_LEFT_FLAG ? (left_ctb_alf_flag & 1) : 0;
            ctx += ctu_neighbour_flags & VVC_CTU_UP_FLAG   ? (up_ctb_alf_flag   & 1) : 0;
            ret_cr = ovcabac_ae_read(cabac_ctx,&cabac_state[CTB_ALF_FLAG_CTX_OFFSET + 6 + ctx]);
            while (ret_cr && decoded < num_alf_alternative - 1 && ovcabac_ae_read(cabac_ctx,
                                                                             &cabac_state[CTB_ALF_ALTERNATIVE_CTX_OFFSET + 1])){
                ++decoded;
            }
        }
    //}
    #endif
    ret = (ret_luma << 2) | (ret_cb << 1) | ret_cr;
    return ret;
}

static uint8_t
ovcabac_read_ae_alf_idx(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                   unsigned int tile_group_num_aps)
{
    uint8_t filter_idx = 0;
    return filter_idx;
}

