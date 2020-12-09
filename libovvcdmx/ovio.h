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
OVIOStream *ovio_stream_open(FILE *fstream);

void ovio_stream_close(OVIOStream *io_str);

size_t ovio_stream_read(const uint8_t **dst_buff, size_t size, OVIOStream *const io_str);

int ovio_stream_eof(OVIOStream *const io_str);

int ovio_stream_error(OVIOStream *const io_str);

long int ovio_stream_tell(OVIOStream *const io_str);

#endif
