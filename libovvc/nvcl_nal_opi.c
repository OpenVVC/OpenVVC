#include "nvcl.h"
#include "nvcl_utils.h"

typedef struct OVOPI
{
    uint8_t opi_ols_info_present_flag;
    uint8_t opi_htid_info_present_flag;
    uint8_t opi_ols_idx;
    uint8_t opi_htid_plus1;
    uint8_t opi_extension_flag;
    uint8_t opi_extension_data_flag;

} OVOPI;

int
nvcl_opi_read(OVNVCLReader *const rdr, OVOPI *const opi,
                  OVNVCLCtx *const nvcl_ctx)
{
    opi->opi_ols_info_present_flag  = nvcl_read_flag(rdr);
    opi->opi_htid_info_present_flag = nvcl_read_flag(rdr);
    if(opi->opi_ols_info_present_flag) {
        opi->opi_ols_idx = nvcl_read_u_expgolomb(rdr);
    }

    if(opi->opi_htid_info_present_flag) {
        opi->opi_htid_plus1 = nvcl_read_bits(rdr, 3);
    }

    opi->opi_extension_flag = nvcl_read_flag(rdr);
    if (opi->opi_extension_flag) {
        while (more_rbsp_data()) {
            opi->opi_extension_data_flag = nvcl_read_flag(rdr);
        }
    }

    rbsp_trailing_bits();
}
