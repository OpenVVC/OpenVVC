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

#ifndef OVVCDMX_H
#define OVVCDMX_H

/* Experimental raw video demuxer Annex B */

#include <stdio.h>
#include "ovio.h"
#include "ovunits.h"

typedef struct OVDemux OVVCDmx;
typedef OVVCDmx OVDemux;
typedef struct NALUnitsList NALUnitsList;

/* Initialize demuxer
 */
int ovdmx_init(OVDemux **ovdmx_p);

/* Close demuxer
 */
int ovdmx_close(OVDemux *ovdmx);

/* Attach an input stream to the demuxer
 */
int ovdmx_attach_stream(OVDemux *const ovdmx, OVIO *io);

/* Reinit the demuxer.
 */
void ovdmx_detach_stream(OVDemux *const ovdmx);

/* Create a OVPictureUnit from a NALUnitsList
 *
 * Note
 *    - OVNALUnit from src list are cleared by this function.
 *    - The NALUnitsList is assumed non empty.
 */
int ovdmx_init_pu_from_list(OVPictureUnit **ovpu_p, struct NALUnitsList *const src);

/* Extract a Picture Unit
 *
 * Note :
 *     - at the current time output does not correspond to
 *     a complete Picture Unit. The OVPictureUnit only contains
 *     one OVNALUnit.
 */
int ovdmx_extract_picture_unit(OVDemux *const ovdmx, OVPictureUnit **ovpu_p);

#endif

