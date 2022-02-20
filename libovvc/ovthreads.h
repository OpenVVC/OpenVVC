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

#ifndef ovthread_H
#define ovthread_H

#include <stdint.h>

#include "slicedec.h"

#define USE_THREADS 1

struct SliceSynchro;
struct OVVCDec;
struct OVFrame;
/*
Functions for the threads decoding rectangular entries
*/
struct EntryThread
{
    struct MainThread *main_thread;
    pthread_t thread;
    pthread_mutex_t entry_mtx;
    pthread_cond_t  entry_cnd;

    /* CTU decoder associated to entry 
     * thread
     */
    OVCTUDec *ctudec;

    uint8_t state;
    uint8_t kill;
};

struct EntryJob{
    struct SliceSynchro *slice_sync;
    uint8_t entry_idx;
};

int ovthread_init_entry_thread(struct EntryThread *entry_th);

void ovthread_uninit_entry_thread(struct EntryThread *entry_th);

int ovthread_slice_add_entry_jobs(struct SliceSynchro *slice_sync, DecodeFunc decode_entry, int nb_entries);

int ovthread_slice_sync_init(struct SliceSynchro *slice_sync);

void ovthread_slice_sync_uninit(struct SliceSynchro *slice_sync);

#endif
