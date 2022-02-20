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

#include "ovdefs.h"
#include "nvcl.h"
#include "nvcl_utils.h"
#include "nvcl_structures.h"
#include "nvcl_private.h"

/* sps_max_sublayers_minus1, sps_sublayer_dpb_params_flag */
/* vps_dpb_max_tid, vps_sublayer_dpb_params_present_flag */
int
dpb_parameters(OVNVCLReader *const rdr, OVDPBParams *const dpb_list, int max_sub_layer_minus1, int sub_layer_info_flag)
{
    /*FIXME loop outside of function */
    int i;
    for (i = (sub_layer_info_flag ? 0 : max_sub_layer_minus1); i <= max_sub_layer_minus1; ++i) {
        OVDPBParams *const dpb = &dpb_list[i];
        dpb->dpb_max_dec_pic_buffering_minus1 = nvcl_read_u_expgolomb(rdr);
        dpb->dpb_max_num_reorder_pics         = nvcl_read_u_expgolomb(rdr);
        dpb->dpb_max_latency_increase_plus1   = nvcl_read_u_expgolomb(rdr);
    }
    return 0;
}
