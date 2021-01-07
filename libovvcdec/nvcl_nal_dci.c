#include "nvcl.h"
#include "nvcl_utils.h"


typedef struct OVDCI
{
    uint8_t dci_reserved_zero_4bits;
    uint8_t dci_num_ptls_minus1;
    uint8_t dci_extension_flag;
    uint8_t dci_extension_data_flag;
} OVDCI;

int
nvcl_dci_read(OVNVCLReader *const rdr, OVDCI *const dci,
              OVNVCLCtx *const nvcl_ctx)
{
    int i;
    dci->dci_reserved_zero_4bits = nvcl_read_bits(rdr, 4);
    dci->dci_num_ptls_minus1 = nvcl_read_bits(rdr, 4);
    for(i = 0; i <= dci->dci_num_ptls_minus1; i++) {
        profile_tier_level(1, 0);
    }

    dci->dci_extension_flag = nvcl_read_flag;
    if (dci->dci_extension_flag) {
        while (more_rbsp_data()){
            dci->dci_extension_data_flag = nvcl_read_flag;
        }
    }

    rbsp_trailing_bits()
}
