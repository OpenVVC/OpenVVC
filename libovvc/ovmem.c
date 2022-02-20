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

#include <stdlib.h>
#include <string.h>
#include "ovmem.h"
#include "ovconfig.h"


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
    if (alloc_size) { //OS X on SDK 10.6 has a broken posx_memalign implementation
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
#if HAVE_ALIGNED_MALLOC
    _aligned_free(ptr); //For windows
#else
    free(ptr);
#endif
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
