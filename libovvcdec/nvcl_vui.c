#include "nvcl.h"
#include "nvcl_utils.h"

/* FIXME find other spec */

vui_payload( payloadSize ) {
    vui_parameters(payloadSize);
    if (more_data_in_payload()) {
        if (payload_extension_present()) {
            uint8_t vui_reserved_payload_extension_data;
        }
        uint8_t vui_payload_bit_equal_to_one;
        while(!byte_aligned()) {
            uint8_t vui_payload_bit_equal_to_zero;
        }
    }
}
