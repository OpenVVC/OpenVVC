#include "ovannexb.h"
#include "ovunits.h"


/* Check for start_code or emulation prevention byte
 * returns:
 *    0 if no start code was found
      1 if a valid start code was found
      2 if a valid zero byte emulation prevention was found
      FIXME: should we return 4 if byte2 was zero ? (might be useful to
      abort detection on filling data or detect a new sequence without EOS
      or EOB )
      negative value if something invalid emulation prevention
      or start code with an invalid NAL header
 */
int
ovannexb_check_stc_or_epb(const uint8_t *byte)
{
    /* we consider first byte has already been detected ?*/
    uint8_t byte1 = byte[1];
    uint8_t byte2 = byte[2];
    uint8_t byte3 = byte[3];

    if ((byte1 == 0) && !(byte2 & (~0x3))) {
        if (byte2 == 0x01) {
            uint8_t fbdn_zbit = byte3 & 0x01;
            if (fbdn_zbit) {
                goto invalid_data;
            }
            return 1;
        }

        /* Note this is only useful when removing rbsp
           and for probing*/
        if (byte2 == 0x03) {
            uint8_t invalid_emu = byte3 & (~0x3);
            if (invalid_emu) {
                goto invalid_data;
            }
            return 2;
        }
    } else if (byte1 == 0 && byte2 == 0) {
       /*FIXME check if this should be used as RBSP end detection*/
    }

    return 0;

invalid_data:
    /*TODO detect log what went wrong and return correct value
      Note this could be handled at a higher level based on return
      in case of failure we do not need to be fast*/
    return -1;
}

int
dmx_process_elem(OVVCDmx *const dmx, const uint8_t *const bytestream,
        uint64_t byte_pos, int stc_or_epb)
{
    /* TODO decide what to do with this this is supposed to become
       an internal function this could be used as a callback to do
       various things according to the demuxer context*/
    if (stc_or_epb == 1) {
        uint8_t byte4 = bytestream[4];
        enum OVNALUType nalu_type = (byte4 >> 3) & 0x1F;
        /* Start code */
        printf("STC at pos : %ld, NAL :%d\n", byte_pos - 8, nalu_type);
    } else {
        /* Emulation zero byte prevention */
        printf("EPB at pos : %ld\n", byte_pos - 8);
    }

    return 1;
}

