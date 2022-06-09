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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "ovmem.h"

#include "ovio.h"

/* Use Fixed size of demux read cache buffer to 64K */

#define OVIO_FILEIO_BUFF_SIZE (1 << 16)
#define OVIO_CACHE_PADDING 8

typedef struct OVReadBuff{
    uint8_t *bytestream;
    void *opaque;
}OVReadBuff;

struct OVIOStream {
    OVIO *io;
    const uint8_t *bytestream;
    OVReadBuff opaque_cache;
};

static int OVFileIOClose(OVIO* io)
{
    OVFileIO* file_io = (OVFileIO*) io;

    int ret = fclose(file_io->file);

    free(file_io);

    return ret;
}

static size_t OVFileIORead(void *ptr, OVIO* io)
{
    OVFileIO* file_io = (OVFileIO*) io;

    size_t nb_bytes_read = fread(ptr, 1, io->size, file_io->file);

    if (!nb_bytes_read) {
        nb_bytes_read = ftell(file_io->file);
    }

    return nb_bytes_read;
}

static int OVFileIOEOF(OVIO* io)
{
    OVFileIO* file_io = (OVFileIO*) io;

    return feof(file_io->file);
}

const OVFileIO defaultFileIO = {
  .super = { .close = OVFileIOClose, .read = OVFileIORead, .eof = OVFileIOEOF, .size = OVIO_FILEIO_BUFF_SIZE },
  .file = NULL
};

OVFileIO*
ovio_new_fileio(const char* path, const char* mode)
{
  OVFileIO* io = ov_malloc(sizeof(OVFileIO));
  if (!io) return NULL;
  memcpy(io, &defaultFileIO, sizeof(OVFileIO));

  io->file = fopen(path, mode);
  if (!io->file) goto erropen;

  return io;

erropen:
  ov_free(io);
  return NULL;

}

static int ovread_buff_init(struct OVReadBuff *const cache_buff,
                             size_t buff_size);

static void ovread_buff_close(struct OVReadBuff *const cache_buff);

OVIOStream *
ovio_stream_open(OVIO *io)
{
    OVIOStream *io_str;
    int ret;
    if (io == NULL) {
        return NULL;
    }

    io_str = ov_mallocz(sizeof(OVIOStream));
    if (io_str == NULL) {
        return io_str;
    }

    ret = ovread_buff_init(&io_str->opaque_cache, io->size);
    if (ret < 0) {
        ov_freep(&io_str);
        return io_str;
    }

    io_str->io = io;
    io_str->bytestream = io_str->opaque_cache.bytestream;

    return io_str;
}

void
ovio_stream_close(OVIOStream *io_str)
{
    ovread_buff_close(&io_str->opaque_cache);
    ov_free(io_str);
}

static int
ovread_buff_init(struct OVReadBuff *const cache_buff, size_t buff_size)
{
    uint8_t *byte_stream;

    /* Note we keep a 8 bytes left and 8 bytes padding at the right of
     * our buffer in order to check,to store the last bytes of the
     * previous chunk of the bytestream 8 bytes is because we might
     * want to use 64bit types in order to quickly probe for successive
     * zero bytes
     */
    byte_stream = ov_mallocz(OVIO_CACHE_PADDING + buff_size + OVIO_CACHE_PADDING);

    if (byte_stream == NULL) {
        return -1;
    }

    cache_buff->opaque = (void *)byte_stream;

    /* Bytestream will be cached from this position in allocated memory
     * the padding byte will be used to keep a copy of the last 8 bytes
     * in order to check for start codes overlapping between two chunks
     */
    cache_buff->bytestream = byte_stream + OVIO_CACHE_PADDING;

    /* last 16 bytes are set to 0xFF so we do not detect any zero byte
     * past the actual available data when checking for a start or emulation
     * prevention code.
     * Using. 0xFF should prevent patterns such as 0x000003 or 0x000001
     * we used 16 bytes so first copy from read function will copy 0XFF
     * bytes at the averlapping area reader will then ignore them since
     * it cannot be taken as start code.
     */
    memset(byte_stream + OVIO_CACHE_PADDING + buff_size, 0XFF, sizeof(*byte_stream) * OVIO_CACHE_PADDING);

    return 1;
}

void
ovread_buff_close(struct OVReadBuff *const cache_buff)
{
    cache_buff->bytestream = NULL;

    ov_freep(&cache_buff->opaque);
}

/* FIXME find a way to give a pointer to cache directly to bytestream */
size_t
ovio_stream_read(const uint8_t **dst_buff, OVIOStream *const io_str)
{
    OVIO *io = io_str->io;
    const size_t i_buff_size = io->size;
    uint8_t *cache_start = io_str->opaque_cache.bytestream;
    uint8_t *cache_end = cache_start + i_buff_size;
    size_t read_in_buf;

    memcpy(cache_start - OVIO_CACHE_PADDING, cache_end - OVIO_CACHE_PADDING, sizeof(*cache_start) * OVIO_CACHE_PADDING);

    read_in_buf = io->read(cache_start, io);

    *dst_buff = cache_start - OVIO_CACHE_PADDING;
    if (io->eof(io)) {
        read_in_buf += OVIO_CACHE_PADDING;
    }

    return read_in_buf;
}

/* FIXME decide of an EOF value check what is usual*/
int
ovio_stream_eof(OVIOStream *const io_str)
{
    OVIO *io = io_str->io;

    return io->eof(io);
}

size_t
ovio_stream_buff_size(OVIOStream* const io_str)
{
    OVIO *io = io_str->io;
    return io->size;
}
