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

#ifndef OV_HLS_STRUCTURES_H
#define OV_HLS_STRUCTURES_H
#include "ovdefs.h"
#include "nvcl_structures.h"

union HLSData
{
    OVSH  sh;
    OVPH  ph;
    OVSPS sps;
    OVPPS pps;
    OVAPS aps;
    OVSEI sei;
};

struct HLSDataRef {
    struct HLSDataRef *data_ref;
    const union HLSData *data;
    void (*free)(struct HLSDataRef **ref, void *opaque);
    void *opaque;
    atomic_uint ref_count;
};

struct HLSReader
{
    const char *name;
    const size_t data_size;

    uint8_t (*probe_id)(OVNVCLReader *const rdr);

    struct HLSDataRef **(*find_storage)(OVNVCLReader *const rdr,
                                             OVNVCLCtx *const nvcl_ctx);

    int (*read)(OVNVCLReader *const rdr, OVHLSData *const hls_data,
                const OVNVCLCtx *const nvcl_ctx);

    int (*validate)(OVNVCLReader *rdr, const union HLSData *const hls_data);

    int (*replace)(const struct HLSReader *const manager,
                   struct HLSDataRef **storage,
                   const OVHLSData *const hls_data);

    void (*free)(const union HLSData *const hls_data);
};

void hlsdata_unref(struct HLSDataRef **dataref_p);

void hlsdata_ref_default_free(struct HLSDataRef **ref_p, void *opaque);

struct HLSDataRef * hlsdataref_create(union HLSData *data, void (*free)(struct HLSDataRef **ref, void *opaque), void *opaque);

int
hlsdata_newref(struct HLSDataRef **dst_p, struct HLSDataRef *src);
#endif
