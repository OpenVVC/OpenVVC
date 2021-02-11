#ifndef OVANNEXB_H
#define OVANNEXB_H

#include <stdint.h>

#include "ovdmx.h"

int ovannexb_check_stc_or_epb(const uint8_t *byte);

int dmx_process_elem(OVVCDmx *const dmx,
                     const uint8_t *const bytestream,
                     uint64_t byte_pos,
                     int stc_or_epb);


#endif
