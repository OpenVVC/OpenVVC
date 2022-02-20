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

#ifndef OVVC_ERROR_H
#define OVVC_ERROR_H

#include <stdint.h>

/*List of error*/
#define OVVC_ENOMEM           OVVC_ERROR_TAG('N','M','E','M')
#define OVVC_EINDATA          OVVC_ERROR_TAG('I','N','D','A')
#define OVVC_EUNSUPPORTED     OVVC_ERROR_TAG('U','N','S','P')
#define OVVC_EAGAIN           OVVC_ERROR_TAG('E','A','G','N')

#define OVVC_ERROR_TAG(a,b,c,d) -(((a)<<24)+((b)<<16)+((c)<<8)+(d))

#define OVVC_MAX_ERR_STRLEN 30

const char * ovvc_error_stringify(uint32_t error_code);

#endif/*OVVC_ERROR_H*/
