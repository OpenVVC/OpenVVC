#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVDBPParams
{
    uint8_t dpb_max_dec_pic_buffering_minus1;
    uint8_t dpb_max_num_reorder_pics;
    uint8_t dpb_max_latency_increase_plus1;
} OVDBPParams;

/* sps_max_sublayers_minus1, sps_sublayer_dpb_params_flag */
/* vps_dpb_max_tid, vps_sublayer_dpb_params_present_flag */
dpb_parameters(OVNVCLReader *const rdr, int MaxSubLayersMinus1, int subLayerInfoFlag)
{
    /*FIXME loop outside of function */
    for( i = (subLayerInfoFlag ? 0 : MaxSubLayersMinus1) i <= MaxSubLayersMinus1; i++ ) {
        OVDPBParams *const dpb = dpb_list[i];
        dpb->dpb_max_dec_pic_buffering_minus1 = nvcl_read_u_expgolomb(rdr);
        dpb->dpb_max_num_reorder_pics         = nvcl_read_u_expgolomb(rdr);
        dpb->dpb_max_latency_increase_plus1   = nvcl_read_u_expgolomb(rdr);
    }
}
