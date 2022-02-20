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

#include <string.h>

#include "mempool.h"
#include "mempool_internal.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"

#include "ovframepool.h"
#include "nvcl.h"
#include "nvcl_structures.h"
#include "ovframe.h"
#include "ovdpb.h"

void
dpbpriv_uninit_framepool(struct DPBInternal *dpb_priv)
{
    ovframepool_uninit(&dpb_priv->frame_pool);
}

int
dpbpriv_init_framepool(struct DPBInternal *dpb_priv, const OVSPS *const sps)
{
    int ret;

    ret = ovframepool_init(&dpb_priv->frame_pool, sps->sps_chroma_format_idc,
                           sps->sps_bitdepth_minus8,
                           sps->sps_pic_width_max_in_luma_samples,
                           sps->sps_pic_height_max_in_luma_samples);
    if (ret < 0) {
        goto fail_init;
    }

    return 0;
fail_init:
    return ret;
}

int
dpbpriv_request_frame(struct DPBInternal *dpb_priv, OVFrame **frame_p)
{
    *frame_p = ovframepool_request_frame(dpb_priv->frame_pool);
    if (!*frame_p) {
        return OVVC_ENOMEM;
    }

    return 0;
}

