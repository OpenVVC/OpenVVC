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

#ifndef OVIO_H
#define OVIO_H

#include <stddef.h>
#include <stdio.h>
#include <stdint.h>

/*
 * This file contains wrappers for IO functions
 * names  should mimic stdio functions since they
 * are basically aiming to offer the same interface
 */

typedef struct OVIOStream OVIOStream;

typedef struct OVIO {
    int    (*const close)(struct OVIO*);
    size_t (*const read )(void *, struct OVIO*);
    int    (*const eof  )(struct OVIO*);
    size_t size;
} OVIO;

typedef struct OVFileIO {
    struct OVIO super;
    FILE* file;
} OVFileIO;

OVFileIO* ovio_new_fileio(const char* path, const char* mode);

OVIOStream *ovio_stream_open(OVIO *io);

void ovio_stream_close(OVIOStream *io_str);

size_t ovio_stream_read(const uint8_t **dst_buff, OVIOStream *const io_str);

int ovio_stream_eof(OVIOStream *const io_str);

size_t ovio_stream_buff_size(OVIOStream* const io_str);

#endif
