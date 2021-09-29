#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "ovmem.h"

#include "ovio.h"

/* Use Fixed size of demux read cache buffer to 64K */

#define OVIO_FILEIO_BUFF_SIZE (1 << 16)

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
    int ret = 0;
    OVFileIO* file_io = (OVFileIO*) io;
    ret = fclose(file_io->file);
    free(file_io);
    return ret;
}

static size_t OVFileIORead(void *ptr, OVIO* io)
{
    int read = 0;
    OVFileIO* file_io = (OVFileIO*) io;
    read = fread(ptr, 1, io->size, file_io->file);
    if(!read)
        read = ftell(file_io->file);
    return  read;
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
  memcpy(io, &defaultFileIO, sizeof(OVFileIO));
  io->file = fopen(path, mode);
  return io;
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
    byte_stream = ov_mallocz(8 + buff_size + 8);

    if (byte_stream == NULL) {
        return -1;
    }

    cache_buff->opaque = (void *)byte_stream;

    /* Bytestream will be cached from this position in allocated memory
     * the padding byte will be used to keep a copy of the last 8 bytes
     * in order to check for start codes overlapping between two chunks
     */
    cache_buff->bytestream = byte_stream + 8;

    /* last 16 bytes are set to 0xFF so we do not detect any zero byte
     * past the actual available data when checking for a start or emulation
     * prevention code.
     * Using. 0xFF shoul prevent patterns such as 0x000003 or 0x000001
     * we used 16 bytes so first copy from read function will copy 0XFF
     * bytes at the averlapping area reader will then ignore them since
     * it cannot be taken as start code.
     */
    memset(byte_stream + 8 + buff_size, 0XFF, sizeof(*byte_stream) * 8);

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

    /* FIXME this might depend on the demux maybe this should
       be done somewhere else this force cache buffer */
    memcpy(cache_start - 8, cache_end - 8, sizeof(*cache_start) * 8);

    read_in_buf = io->read(cache_start, io);

    *dst_buff = cache_start;

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
