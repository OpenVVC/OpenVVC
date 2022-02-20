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

#include "overror.h"
#include "ovutils.h"

#include "ovframe.h"
#include "ovframepool.h"

int
ovframe_new_ref(OVFrame **dst, OVFrame *src)
{
     if (!src) {
         return -1;
     }

    unsigned ref_count = atomic_fetch_add_explicit(&src->internal.ref_count, 1, memory_order_acq_rel);
    ov_log(NULL, OVLOG_DEBUG, "NewRef Frame %p ref_count: %d\n", src, ref_count);

    *dst = src;

    return 0;
}

void
ovframe_unref(OVFrame **frame_p)
{
    if (!frame_p)
        return;

    if (!*frame_p){
        ov_log(NULL, OVLOG_ERROR, "Trying to unref NULL frame\n");
        return;
    }

    unsigned ref_count = atomic_fetch_add_explicit(&(*frame_p)->internal.ref_count, -1, memory_order_acq_rel);
    ov_log(NULL, OVLOG_DEBUG, "Unref Frame %p ref_count: %d\n", *frame_p, ref_count);

    if (!ref_count) {
        ovframepool_release_frame(frame_p);
    }

    *frame_p = NULL;
}
