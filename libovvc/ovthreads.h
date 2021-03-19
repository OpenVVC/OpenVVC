#ifndef OVTHREADS_H
#define OVTHREADS_H

struct SliceThreads;

int ovthread_decode_entries(struct SliceThreads *th_info, DecodeFunc decode_entry, int nb_entries);

int init_entry_threads(struct SliceThreads *th_info, int nb_threads);

void uninit_entry_threads(struct SliceThreads *th_info);
#endif
