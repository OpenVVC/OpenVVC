#ifndef OVMEM_H
#define OVMEM_H
#include <stddef.h>

#define ALIGN 32

#define DECLARE_ALIGNED(n, t, v) t __attribute__ ((aligned (n))) v
#define ov_malloc_attrib   __attribute__((__malloc__))

void *ov_malloc(size_t alloc_size) ov_malloc_attrib;

void *ov_mallocz(size_t alloc_size) ov_malloc_attrib;

void ov_free(void *ptr);

void ov_freep(void *ptr_ref);

#endif
