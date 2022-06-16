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

#include <pthread.h>
/* FIXME tmp*/
#include <stdatomic.h>

#include "slicedec.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"
#include "ovthreads.h"
#include "ovdpb.h"

static int
ovthread_decode_entry(struct EntryJob *entry_job, struct EntryThread *entry_th)
{   
    struct SliceSynchro* slice_sync = entry_job->slice_sync;
    uint8_t entry_idx               = entry_job->entry_idx;

    uint16_t nb_entries  = slice_sync->nb_entries;
    ov_log(NULL, OVLOG_DEBUG, "Decoder with POC %d, start entry %d\n", slice_sync->owner->pic->poc, entry_idx);
    
    OVCTUDec *const ctudec  = entry_th->ctudec;
    OVSliceDec *const sldec = slice_sync->owner;
    const OVPS *const prms  = &sldec->active_params;

    slice_sync->decode_entry(sldec, ctudec, prms, entry_idx);

    uint16_t nb_entries_decoded = atomic_fetch_add_explicit(&slice_sync->nb_entries_decoded, 1, memory_order_acq_rel);

    /* Last thread to exit loop will have entry_idx set to nb_entry - 1*/
    return nb_entries_decoded == nb_entries - 1 ;
}

static struct EntryJob *
fifo_pop_entry(struct EntriesFIFO *fifo)
{
    struct EntryJob *entry_job = NULL;

    if (fifo->first != fifo->last) {
        ptrdiff_t position = fifo->first - fifo->entries;
        uint16_t size = fifo->size;

        entry_job = fifo->first;

        fifo->first = &fifo->entries[(position + 1) % size];

    }
    return entry_job;
}

static void *
entry_thread_main_function(void *opaque)
{
    struct EntryThread *entry_th = (struct EntryThread *)opaque;
    struct MainThread* main_thread = entry_th->main_thread;

    pthread_mutex_lock(&entry_th->entry_mtx);
    entry_th->state = ACTIVE;
    pthread_mutex_unlock(&entry_th->entry_mtx);

    while (!entry_th->kill){
        struct EntriesFIFO *fifo = &main_thread->entries_fifo;

        pthread_mutex_lock(&main_thread->io_mtx);

        struct EntryJob entry_jobtmp;

        struct EntryJob *entry_job = NULL;
        struct EntryJob *entry_job2 = fifo_pop_entry(fifo);
        pthread_cond_signal(&main_thread->entry_threads_cnd);
        if (entry_job2) {
           entry_jobtmp = *entry_job2;
           entry_job = &entry_jobtmp;
        }

        pthread_mutex_unlock(&main_thread->io_mtx);

        if (entry_job) {
            slicedec_update_entry_decoder(entry_job->slice_sync->owner, entry_th->ctudec);
            pthread_mutex_lock(&entry_th->entry_mtx);
            pthread_mutex_unlock(&entry_th->entry_mtx);

            uint8_t is_last = ovthread_decode_entry(entry_job, entry_th);

            /* Check if the entry was the last of the slice */
            if (is_last) {
                slicedec_finish_decoding(entry_job->slice_sync->owner);
            }
        } else {
            pthread_mutex_lock(&main_thread->entry_threads_mtx);
            pthread_mutex_lock(&entry_th->entry_mtx);
            entry_th->state = IDLE;

            pthread_cond_signal(&main_thread->entry_threads_cnd);
            pthread_mutex_unlock(&main_thread->entry_threads_mtx);

            pthread_cond_wait(&entry_th->entry_cnd, &entry_th->entry_mtx);
            entry_th->state = ACTIVE;
            pthread_mutex_unlock(&entry_th->entry_mtx);
        }
    }
    return NULL;
}

int
ovthread_init_entry_thread(struct EntryThread *entry_th)
{
    entry_th->state = IDLE;
    entry_th->kill  = 0;

    int ret = ctudec_init(&entry_th->ctudec);
    if (ret < 0) {
        ov_log(NULL, OVLOG_ERROR, "Failed line decoder initialisation\n");
        ctudec_uninit(entry_th->ctudec);
        return OVVC_ENOMEM;
    }

    pthread_mutex_init(&entry_th->entry_mtx, NULL);
    pthread_cond_init(&entry_th->entry_cnd, NULL);
    pthread_mutex_lock(&entry_th->entry_mtx);

#if USE_THREADS
    if (pthread_create(&entry_th->thread, NULL, entry_thread_main_function, entry_th)) {
        pthread_mutex_unlock(&entry_th->entry_mtx);
        ov_log(NULL, OVLOG_ERROR, "Thread creation failed at decoder init\n");
        return OVVC_ENOMEM;
    }
#endif
    pthread_mutex_unlock(&entry_th->entry_mtx);
    return 1;
}

void
ovthread_uninit_entry_thread(struct EntryThread *entry_th)
{       
    pthread_mutex_destroy(&entry_th->entry_mtx);
    pthread_cond_destroy(&entry_th->entry_cnd);

    ctudec_uninit(entry_th->ctudec);
}

static int
fifo_push_entry(struct EntriesFIFO *fifo,
                struct SliceSynchro *slice_sync, int entry_idx)
{
    ptrdiff_t position = fifo->last - fifo->entries;
    ptrdiff_t next_pos = (position + 1) % fifo->size;
    do {
        if (&fifo->entries[next_pos] != fifo->first) {
            struct EntryJob *entry_job = &fifo->entries[position];

            entry_job->entry_idx  = entry_idx;
            entry_job->slice_sync = slice_sync;

            fifo->last = &fifo->entries[next_pos];
            return 0;
        } else {
            struct EntryThread *entry_threads_list = slice_sync->main_thread->entry_threads_list;
            for (int i = 0; i < slice_sync->main_thread->nb_entry_th; ++i){
                struct EntryThread *th_entry = &entry_threads_list[i];
                pthread_mutex_lock(&th_entry->entry_mtx);
                pthread_cond_signal(&th_entry->entry_cnd);
                pthread_mutex_unlock(&th_entry->entry_mtx);
            }
            pthread_cond_wait(&slice_sync->main_thread->entry_threads_cnd, &slice_sync->main_thread->io_mtx);
        }
    }while (1);
}

/*
Functions needed for the synchro of threads decoding the slice
*/
int
ovthread_slice_add_entry_jobs(struct SliceSynchro *slice_sync, DecodeFunc decode_entry, int nb_entries)
{
    struct MainThread* main_thread = slice_sync->main_thread;
    struct EntriesFIFO *entry_fifo = &main_thread->entries_fifo;

    slice_sync->nb_entries = nb_entries;
    slice_sync->decode_entry = decode_entry;

    atomic_store_explicit(&slice_sync->nb_entries_decoded, 0, memory_order_relaxed);

    pthread_mutex_lock(&main_thread->io_mtx);
    for (int i = 0; i < nb_entries; ++i) {
        fifo_push_entry(entry_fifo, slice_sync, i);
        ov_log(NULL, OVLOG_DEBUG, "Main adds POC %d entry %d\n", slice_sync->owner->pic->poc, i);
    }
    pthread_mutex_unlock(&main_thread->io_mtx);

    struct EntryThread *entry_threads_list = main_thread->entry_threads_list;
    for (int i = 0; i < main_thread->nb_entry_th; ++i){
        struct EntryThread *th_entry = &entry_threads_list[i];
        pthread_mutex_lock(&th_entry->entry_mtx);
        // ov_log(NULL, OVLOG_DEBUG,"main sign entry\n");
        pthread_cond_signal(&th_entry->entry_cnd);
        pthread_mutex_unlock(&th_entry->entry_mtx);
    }

    return 0;
}

int
ovthread_slice_sync_init(struct SliceSynchro *slice_sync)
{   
    atomic_init(&slice_sync->nb_entries_decoded, 0);

    pthread_mutex_init(&slice_sync->gnrl_mtx, NULL);
    pthread_cond_init(&slice_sync->gnrl_cnd,  NULL);

    return 0;
}

void
ovthread_slice_sync_uninit(struct SliceSynchro *slice_sync)
{   
    OVSliceDec * slicedec = slice_sync->owner;

    pthread_mutex_destroy(&slice_sync->gnrl_mtx);
    pthread_cond_destroy(&slice_sync->gnrl_cnd);

}

