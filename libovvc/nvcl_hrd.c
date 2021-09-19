#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVSubLayerHRD
{
    uint8_t bit_rate_value_minus1;
    uint8_t cpb_size_value_minus1;
    uint8_t cpb_size_du_value_minus1;
    uint8_t bit_rate_du_value_minus1;
    uint8_t cbr_flag[0];
} OVSubLayerHRD;

typedef struct OVOLSHRDParams
{
    uint8_t fixed_pic_rate_general_flag;
    uint8_t fixed_pic_rate_within_cvs_flag;
    uint8_t elemental_duration_in_tc_minus1;
    uint8_t low_delay_hrd_flag;
} OVOLSHRDParams;

typedef struct OVGHRD
{
    uint32_t num_units_in_tick;
    uint32_t time_scale;
    uint8_t general_nal_hrd_params_present_flag;
    uint8_t general_vcl_hrd_params_present_flag;
    uint8_t general_same_pic_timing_in_all_ols_flag;
    uint8_t general_du_hrd_params_present_flag;
    uint8_t tick_divisor_minus2;
    uint8_t bit_rate_scale;
    uint8_t cpb_size_scale;
    uint8_t cpb_size_du_scale;
    uint8_t hrd_cpb_cnt_minus1;
} OVGHRD;

int
sublayer_hrd_parameters(OVNVCLReader *const rdr, int subLayerId)
{
    int j;
    OVGHRD *ghrd;
    OVSubLayerHRD **sl_hrd_list = &sl_hrd_list[subLayerId];
    for (j = 0; j <= ghrd->hrd_cpb_cnt_minus1; j++) {
        OVSubLayerHRD *sl_hrd = sl_hrd_list[j];

        sl_hrd->bit_rate_value_minus1 = nvcl_read_u_expgolomb(rdr);
        sl_hrd->cpb_size_value_minus1 = nvcl_read_u_expgolomb(rdr);

        if (ghrd->general_du_hrd_params_present_flag) {
            sl_hrd->cpb_size_du_value_minus1 = nvcl_read_u_expgolomb(rdr);
            sl_hrd->bit_rate_du_value_minus1 = nvcl_read_u_expgolomb(rdr);
        }

        sl_hrd->cbr_flag[j] = nvcl_read_flag(rdr);
    }
}

int
ols_timing_hrd_parameters(OVNVCLReader *const rdr, int firstSubLayer, int MaxSubLayersVal)
{
    int i;
    OVGHRD *ghrd;
    OVOLSHRDParams *ols_hrd_list[64] = {0};
    for (int i = firstSubLayer; i <= MaxSubLayersVal; i++) {
        OVOLSHRDParams *ols_hrd = ols_hrd_list[i];
        ols_hrd->fixed_pic_rate_general_flag = nvcl_read_flag(rdr);
        if (!ols_hrd->fixed_pic_rate_general_flag) {
            ols_hrd->fixed_pic_rate_within_cvs_flag = nvcl_read_flag(rdr);
        }

        if (ols_hrd->fixed_pic_rate_within_cvs_flag) {
            ols_hrd->elemental_duration_in_tc_minus1 = nvcl_read_u_expgolomb(rdr);
        } else if ((ghrd->general_nal_hrd_params_present_flag || ghrd->general_vcl_hrd_params_present_flag) && ghrd->hrd_cpb_cnt_minus1 == 0) {
            ols_hrd->low_delay_hrd_flag = nvcl_read_flag(rdr);
        }

        if (ghrd->general_nal_hrd_params_present_flag) {
            sublayer_hrd_parameters(rdr, i);
        }

        if (ghrd->general_vcl_hrd_params_present_flag) {
            sublayer_hrd_parameters(rdr, i);
        }
    }
}

int
general_timing_hrd_parameters(OVNVCLReader *const rdr)
{
    uint8_t ghrd_present;
    OVGHRD *ghrd;
    ghrd->num_units_in_tick = nvcl_read_bits(rdr, 32);
    ghrd->time_scale = nvcl_read_bits(rdr, 32);
    ghrd->general_nal_hrd_params_present_flag = nvcl_read_flag(rdr);
    ghrd->general_vcl_hrd_params_present_flag = nvcl_read_flag(rdr);

    ghrd_present  = ghrd->general_nal_hrd_params_present_flag;
    ghrd_present |= ghrd->general_vcl_hrd_params_present_flag;

    if (ghrd_present) {
        ghrd->general_same_pic_timing_in_all_ols_flag = nvcl_read_flag(rdr);
        ghrd->general_du_hrd_params_present_flag = nvcl_read_flag(rdr);
        if (ghrd->general_du_hrd_params_present_flag) {
            ghrd->tick_divisor_minus2 = nvcl_read_bits(rdr, 8);
        }
        ghrd->bit_rate_scale = nvcl_read_bits(rdr, 4);
        ghrd->cpb_size_scale = nvcl_read_bits(rdr, 4);
        if (ghrd->general_du_hrd_params_present_flag) {
            ghrd->cpb_size_du_scale = nvcl_read_bits(rdr, 4);
        }
        ghrd->hrd_cpb_cnt_minus1 = nvcl_read_u_expgolomb(rdr);
    }
}
