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

#ifndef NVCL_PRIVATE_H
#define NVCL_PRIVATE_H
#include "ovdefs.h"

#include "nvcl.h"


int dpb_parameters(OVNVCLReader *const rdr, OVDPBParams *const dpb_list,
                   int max_sub_layer_minus1, int sub_layer_info_flag);

int profile_tier_level_vps(OVNVCLReader *const rdr, uint8_t vps_pt_present_flag, uint8_t vps_ptl_max_tid);

int profile_tier_level_sps(OVNVCLReader *const rdr,  uint8_t sps_max_sublayers_minus1);

int profile_tier_level_dci(OVNVCLReader *const rdr);

int general_constraints_info(OVNVCLReader *const rdr);

/* This one is called by PH/SH reader */
int nvcl_read_header_ref_pic_lists(OVNVCLReader *const rdr, OVHRPL *const rpl_h,
                                   const OVSPS *const sps, const OVPPS *pps);

int nvcl_read_sps_ref_pic_list(OVNVCLReader *const rdr, const OVSPS *const sps,
                               OVRPL *const rpl);


#endif
