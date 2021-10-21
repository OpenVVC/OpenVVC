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


int ovthread_init_entry_threads(OVVCDec *vvcdec, int nb_threads);

void ovthread_uninit_entry_threads(OVVCDec *vvcdec);

void ovthread_init_entry_jobs(OVVCDec *vvcdec, int nb_entry_th);

int ovthread_init_main_thread(OVVCDec *vvcdec);

int ovthread_uninit_main_thread(OVVCDec *vvcdec);


int ovthread_decode_entries(struct SliceSynchro *slice_sync, DecodeFunc decode_entry, int nb_entries);

int ovthread_add_entry_jobs(struct SliceSynchro *slice_sync, DecodeFunc decode_entry, int nb_entries);

int ovthread_slice_sync_init(struct SliceSynchro *slice_sync);

void ovthread_slice_sync_uninit(struct SliceSynchro *slice_sync);

#endif
