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

#ifndef RCN_NEON_H
#define RCN_NEON_H
#include "rcn_structures.h"
#include "ovconfig.h"

#define SIMDE_ENABLE_NATIVE_ALIASES

void rcn_init_mc_functions_neon(struct RCNFunctions *const rcn_funcs);
void rcn_init_mc_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_tr_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_dc_planar_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_ict_functions_sse(struct RCNFunctions *rcn_func, uint8_t type);
void rcn_init_alf_functions_sse(struct RCNFunctions *rcn_func);
void rcn_init_cclm_functions_sse(struct RCNFunctions *rcn_func);
void rcn_init_lfnst_functions_sse(struct RCNFunctions *rcn_func);
void rcn_init_mip_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_sao_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_dmvr_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_prof_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_bdof_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_ciip_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_df_functions_sse(struct RCNFunctions *const rcn_funcs);
void rcn_init_intra_angular_functions_10_sse(struct RCNFunctions *rcn_func);

#endif//RCN_NEON_H
