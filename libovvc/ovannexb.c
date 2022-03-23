/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

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

        /* Note this is only useful when removing rbsp and for probing*/
        if (byte2 == 0x03) {
            uint8_t invalid_emu = byte3 & (~0x3);
            if (invalid_emu) {
                goto invalid_data;
            }
            return 2;
        }
    } else if (byte1 == 0 && byte2 == 0) {
       /*FIXME check if this should be used as RBSP end detection*/
       return 3;
    }

    return 0;

invalid_data:
    /*TODO detect log what went wrong and return correct value
      Note this could be handled at a higher level based on return
      in case of failure we do not need to be fast*/
    return -1;
}

