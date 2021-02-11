#include "nvcl.h"
#include "nvcl_utils.h"

/* FIXME find other spec */
sei_message(OVNVCLReader *const rdr) {
    int payloadType = 0;
    int payloadSize = 0;

    do {
        uint8_t payload_type_byte;
        payloadType += payload_type_byte;
    } while(payload_type_byte == 0xFF);

    do {
        uint8_t payload_size_byte;
        payloadSize += payload_size_byte;
    } while(payload_size_byte == 0xFF);

    sei_payload(payloadType, payloadSize);
}
