#ifndef OVTHREADS_H
#define OVTHREADS_H

#include <stdint.h>

#include "slicedec.h"

#define USE_THREADS 1

struct SliceThread;
struct OVVCDec;
struct OVFrame;

int ovthread_decode_entries(struct SliceThread *th_slice, DecodeFunc decode_entry, int nb_entries);

int init_entry_threads(struct SliceThread *th_slice, int nb_threads);

void uninit_entry_threads(struct SliceThread *th_slice);


int ovthread_slice_thread_init(struct SliceThread *th_slice, int nb_threads);

void ovthread_slice_thread_uninit(struct SliceThread *th_slice);


uint32_t write_decoded_frame_to_file(OVFrame *const frame, FILE *fp);

int ovthread_output_init(OVVCDec *dec, FILE* fout);

void ovthread_output_uninit(struct OutputThread* t_out);

#endif
