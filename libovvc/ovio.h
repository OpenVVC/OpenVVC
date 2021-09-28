#ifndef OVIO_H
#define OVIO_H

#include <stddef.h>
#include <stdio.h>

/*
 * This file contains wrappers for IO functions
 * names  should mimic stdio functions since they
 * are basically aiming to offer the same interface
 */

typedef struct OVIOStream OVIOStream;

/* TODO open / close */

typedef struct OVIO {
    struct OVIO* (*open)();
    int (*close)(struct OVIO*);
    size_t (*read)(void *, size_t, size_t, struct OVIO*);
    long (*tell)(struct OVIO*);
    int (*eof)(struct OVIO*);
    int (*error)(struct OVIO*);
} OVIO;

typedef struct OVFileIO{
    const struct OVIO super;
    FILE* file;
} OVFileIO;

struct OVFileIO* ovio_new_fileio(const char* path, const char* mode);

OVIOStream *ovio_stream_open(OVIO *io);

void ovio_stream_close(OVIOStream *io_str);

size_t ovio_stream_read(const uint8_t **dst_buff, size_t size, OVIOStream *const io_str);

int ovio_stream_eof(OVIOStream *const io_str);

int ovio_stream_error(OVIOStream *const io_str);

long int ovio_stream_tell(OVIOStream *const io_str);

#endif
