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

#ifndef OV_NVCL_H
#define OV_NVCL_H
#include <stdint.h>

#include "ovdefs.h"
#include "ovunits.h"

#define OV_MAX_NUM_VPS 16
#define OV_MAX_NUM_SPS 16
#define OV_MAX_NUM_PPS 64
#define OV_MAX_NUM_APS 32

struct OVNVCLCtx
{
    /* TODO use an other typedef to store more info in
     * lists
     */
    struct HLSDataRef *vps_list[OV_MAX_NUM_VPS];
    struct HLSDataRef *sps_list[OV_MAX_NUM_SPS];
    struct HLSDataRef *pps_list[OV_MAX_NUM_PPS];
    OVAPS *lmcs_aps_list[OV_MAX_NUM_APS];
    OVAPS *alf_aps_list[OV_MAX_NUM_APS];
    OVAPS *scaling_list_aps_list[OV_MAX_NUM_APS];
    struct HLSDataRef *ph;
    struct HLSDataRef *sh;
    OVSEI *sei;
};

typedef union HLSData OVHLSData;


void nvcl_free_ctx(OVNVCLCtx *const nvcl_ctx);

/* Reading functions */
int nvcl_opi_read(OVNVCLReader *const rdr, OVOPI *const opi,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_dci_read(OVNVCLReader *const rdr, OVDCI *const dci,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_vps_read(OVNVCLReader *const rdr, OVHLSData *const hls_data,
                  const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

int nvcl_sps_read(OVNVCLReader *const rdr, OVHLSData *const sps,
                  const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

int nvcl_pps_read(OVNVCLReader *const rdr, OVHLSData *const pps,
                  const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

int nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_ph_read(OVNVCLReader *const rdr, OVHLSData *const ph,
                 const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

int nvcl_sei_read(OVNVCLReader *const rdr, OVSH *const sh,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sh_read(OVNVCLReader *const rdr, OVHLSData *const ph,
                 const OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

/* Decoding functions */
int nvcl_decode_nalu_hls_data(OVNVCLCtx *const nvcl_ctx, OVNALUnit *nal_unit);

int nvcl_decode_nalu_sh(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

int nvcl_decode_nalu_aps(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr, uint8_t nalu_type);

int nvcl_decode_nalu_sei(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr, uint8_t nalu_type);

void copy_sei_params(OVSEI **dst_p, OVSEI *src);

void nvcl_free_sei_params(OVSEI *sei);

#endif
