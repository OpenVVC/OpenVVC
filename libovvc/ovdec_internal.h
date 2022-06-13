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

#ifndef OVDEC_INTERNAL_H
#define OVDEC_INTERNAL_H

/* Private decoder functions and structures declarations and generic
 */

#include "ovdefs.h"
#include "dec_structures.h"
#include "mempool.h"

struct MVPlane
{

    struct TMVPMV *mvs;
    uint64_t *dirs;

    /* Pool elems */
    void *dir_elem;
    void *mv_elem;
};

struct MVPool
{
    MemPool *dir_pool;
    MemPool *mv_pool;
    /* FIXME dimension info ?*/
};

struct PicPartInfo;

int mvpool_init(struct MVPool **mv_pool_p, const struct PicPartInfo *const pinfo);

void mvpool_uninit(struct MVPool **mv_pool_p);

int mvpool_request_mv_plane(struct MVPool *mv_pool, struct MVPlane *mv_plane);

void mvpool_release_mv_plane(struct MVPlane *mv_plane);

#endif
