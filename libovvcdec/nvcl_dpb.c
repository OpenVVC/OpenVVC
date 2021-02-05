#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_private.h"

/* sps_max_sublayers_minus1, sps_sublayer_dpb_params_flag */
/* vps_dpb_max_tid, vps_sublayer_dpb_params_present_flag */
int
dpb_parameters(OVNVCLReader *const rdr, int max_sub_layer_minus1, int sub_layer_info_flag)
{
    /*FIXME loop outside of function */
    int i;
    OVDPBParams dpb_list[64];
    for (i = (sub_layer_info_flag ? 0 : max_sub_layer_minus1); i <= max_sub_layer_minus1; ++i) {
        OVDPBParams *const dpb = &dpb_list[i];
        dpb->dpb_max_dec_pic_buffering_minus1 = nvcl_read_u_expgolomb(rdr);
        dpb->dpb_max_num_reorder_pics         = nvcl_read_u_expgolomb(rdr);
        dpb->dpb_max_latency_increase_plus1   = nvcl_read_u_expgolomb(rdr);
    }
    return 0;
}
