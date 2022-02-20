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

#ifndef OVDEF_H
#define OVDEF_H

/* This file is used to make type definition available 
 * around the decoder so that we do not need to include
 * multiple headers from other header when declaring
 * functions.
 * This avoids including multiples headers from each others
 * when we only need a type definition leading to some functions
 * being available everywhere.
 * If knowledge of a structure are required from a .c file
 * the structure has to be defined in a specific header
 * such as in nvcl_structures.h or dec_structures.h.
 */

/* Generic Parameters Sets  */
typedef struct OVOPI OVOPI;
typedef struct OVDCI OVDCI;
typedef struct OVVPS OVVPS;
typedef struct OVSPS OVSPS;
typedef struct OVPPS OVPPS;
typedef struct OVAPS OVAPS;
typedef struct OVSEI OVSEI;
typedef struct OVPH OVPH;
typedef struct OVSH OVSH;

/* Generic Parameters Sets parts*/
typedef struct OVPTL OVPTL;
typedef struct OVGHRDTiming OVGHRDTiming;
typedef struct OOVOLSHRDTiming OVOLSHRDTiming;
typedef struct OVVUI OVVUI;
typedef struct OVDPBParams OVDPBParams;
typedef struct OVRPL OVRPL;
typedef struct OVHRPL OVHRPL;

/* Types related to NVCL reader */
typedef struct OVNVCLReader OVNVCLReader;

typedef struct OVNVCLCtx OVNVCLCtx;
typedef struct OVNVCLUnit OVNVCLUnit;

typedef struct OVVCDec OVVCDec;
typedef struct SubDec OVSubDec;

typedef struct OVPS OVPS;

/* Generic Frame Picture */
typedef struct Frame OVFrame;
typedef struct OVPicture OVPicture;

/* */
typedef struct DPB OVDPB;


typedef struct OVPartInfo OVPartInfo;
typedef struct VVCCU VVCCU;

/* Decoders related types */
typedef struct OVCTUDec OVCTUDec;

typedef struct OVSliceDec OVSliceDec;

typedef struct OVMV OVMV;
#endif
