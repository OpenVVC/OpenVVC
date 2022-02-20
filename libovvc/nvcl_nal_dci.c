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
