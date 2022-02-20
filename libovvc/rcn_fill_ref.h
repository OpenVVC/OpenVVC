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

#ifndef RCN_FILL_REF_H
#define RCN_FILL_REF_H

#include <stdint.h>

#define LOG2_UNIT_S 2

static inline uint64_t
non_available_units_map(uint64_t ref_map, uint8_t pos, uint8_t log2_pb_s)
{
    int pos_unit = pos >> LOG2_UNIT_S;
    int nb_ref_units = ((1 << (log2_pb_s + 1)) >> LOG2_UNIT_S) + 1;

    uint64_t req_unit_msk = (1llu << (nb_ref_units + 1)) - 1;
    uint64_t avl_map  = (ref_map >> pos_unit) & req_unit_msk;
    uint64_t navl_map = avl_map ^ req_unit_msk;

    return navl_map;
}

static inline uint64_t
available_units_map(uint64_t ref_map, uint8_t pos, uint8_t log2_pb_s)
{
    int pos_unit = pos >> LOG2_UNIT_S;
    int nb_ref_units = ((1 << (log2_pb_s + 1)) >> LOG2_UNIT_S) + 1;
    uint64_t req_unit_msk = (1llu << (nb_ref_units + 1)) - 1;

    uint64_t avl_map = (ref_map >> pos_unit) & req_unit_msk;

    return avl_map;
}

#endif // RCN_FILL_REF_H
