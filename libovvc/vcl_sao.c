#include "cabac_internal.h"

static uint8_t
ovcabac_read_ae_sao_merge_type(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                          uint8_t neighbour_flags)
{
    uint8_t sao_merge_type = 0;

    #if 0
    if(neighbour_flags & VVC_CTU_LEFT_FLAG){
        sao_merge_type = ovcabac_ae_read(cabac_ctx,&cabac_state[SAO_MERGE_FLAG_CTX_OFFSET]);
    }

    if(!sao_merge_type && neighbour_flags & VVC_CTU_UP_FLAG ){
        sao_merge_type = ovcabac_ae_read(cabac_ctx,&cabac_state[SAO_MERGE_FLAG_CTX_OFFSET]);
        sao_merge_type = sao_merge_type << 1;
    }
    #endif

    return sao_merge_type;
}

static uint8_t
ovcabac_read_ae_sao_type_idx(OVCABACCtx *const cabac_ctx, uint64_t *const cabac_state,
                        uint8_t sao_flags, uint8_t num_bits_sao, uint8_t num_bits_sao_c)
{
    int sao_type_idx = 0;

    int offset[4]={1};
    int k;
    int i;

    #if 0
    if(sao_flags & VVC_SAO_LUMA_SLICE_FLAG){
        if(ovcabac_ae_read(cabac_ctx,&cabac_state[SAO_TYPE_IDX_CTX_OFFSET])){
            sao_type_idx = ovcabac_bypass_read(cabac_ctx) ? 2 : 1; //EO : BO
            //offsets
            for (i = 0; i < 4; i++){
                 //not 5
                for(  k = 0; k < num_bits_sao; k++ ) {
                    if( !ovcabac_bypass_read(cabac_ctx) ){
                        break;
                    }
                }
                offset[i] = k;
            }

            if(sao_type_idx & 1){
                for( k = 0; k < 4; k++ ) {
                  if(offset[k] && ovcabac_bypass_read(cabac_ctx)){
                    offset[k] = -offset[k];
                  }
                }//band position
                for(i=0; i < 5; i++)
                    ovcabac_bypass_read(cabac_ctx);
            } else {//edge
                ovcabac_bypass_read(cabac_ctx);
                ovcabac_bypass_read(cabac_ctx);
            }
        }
    }

    sao_type_idx = sao_type_idx << 2;

    if(sao_flags & VVC_SAO_CHROMA_SLICE_FLAG){
        if(ovcabac_ae_read(cabac_ctx,&cabac_state[SAO_TYPE_IDX_CTX_OFFSET])){
            sao_type_idx |= ovcabac_bypass_read(cabac_ctx) ? 2 : 1; //EO : BO
            //offsets
            for (i = 0; i < 4; i++){
                //not 5
                for(  k = 0; k < num_bits_sao_c; k++ ) {
                    if( !ovcabac_bypass_read(cabac_ctx) ){
                        break;
                    }
                }
                offset[i] = k;
            }//cb
            if(sao_type_idx & 1){
                for( k = 0; k < 4; k++ ) {
                  if(offset[k] && ovcabac_bypass_read(cabac_ctx)){
                    offset[k] = -offset[k];
                  }
                }//band position
                for(i=0; i < 5; i++)
                    ovcabac_bypass_read(cabac_ctx);
            } else {//edge
                ovcabac_bypass_read(cabac_ctx);
                ovcabac_bypass_read(cabac_ctx);
            }
            //cr
            for (i = 0; i < 4; i++){
                //not 5
                for(  k = 0; k < num_bits_sao_c; k++ ) {
                    if( !ovcabac_bypass_read(cabac_ctx) ){
                        break;
                    }
                }
                offset[i] = k;
            }
            if(sao_type_idx & 1){
                for( k = 0; k < 4; k++ ) {
                  if(offset[k] && ovcabac_bypass_read(cabac_ctx)){
                    offset[k] = -offset[k];
                  }
                }//band position
                for(i=0; i < 5; i++)
                    ovcabac_bypass_read(cabac_ctx);
            }
        }
    }
    #endif

    return sao_type_idx;
}

