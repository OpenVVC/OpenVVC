#ifndef OVMEM_H
#define OVMEM_H
#include <stddef.h>

void *ov_malloc(size_t alloc_size);

void *ov_mallocz(size_t alloc_size);

void ov_free(void *ptr);

void ov_freep(void *ptr_ref);

#endif

