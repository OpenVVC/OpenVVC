#include <stdlib.h>
#include <string.h>
#include "ovmem.h"

/* Wrappers around memory related function this will
   be useful if we plan later to adapt our allocations to
   to use aligned memory etc.  without changing it in
   every place we called malloc etc.*/

void *
ov_malloc(size_t alloc_size)
{
  void *ptr;

  //FIXME Should be checked on mutliple OS
#if HAVE_POSIX_MEMALIGN
  if (alloc_size){ //OS X on SDK 10.6 has a broken posix_memalign implementation
    if (posix_memalign(&ptr, ALIGN, alloc_size)){
      ptr = NULL;
    }
  } else {
      ptr = NULL;
  }
#elif HAVE_ALIGNED_MALLOC
  ptr = _aligned_malloc(alloc_size, ALIGN); //For windows
#elif HAVE_MEMALIGN
  #ifndef __DJGPP__
    ptr = memalign(ALIGN, alloc_size);
  #else
    ptr = memalign(alloc_size, ALIGN);
  #endif
#else
  ptr = malloc(alloc_size);
#endif

    // void *ptr = malloc(alloc_size);

    return ptr;
}

void *
ov_mallocz(size_t alloc_size)
{
    void *ptr = ov_malloc(alloc_size);

    if (ptr != NULL) {
        memset(ptr, 0, alloc_size);
    }

    return ptr;
}

void
ov_free(void *ptr)
{
    free(ptr);
}

/* This is inspired from FFmpeg av_freep function
 * so we can free and set a pointer to NULL
 */
void
ov_freep(void *ptr_ref)
{
     void *ptr_cpy;
    /* ptr_cpy = *ptr_ref = ptr*/
    memcpy(&ptr_cpy, ptr_ref, sizeof(ptr_cpy));

     /* *ptr_ref = &{NULL}  --> ptr = NULL */

    memcpy(ptr_ref, &(void *){NULL}, sizeof(ptr_cpy));

    ov_free(ptr_cpy);
}
