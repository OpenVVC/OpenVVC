#include "nvcl.h"
#include "nvcl_utils.h"

int
nvcl_reader_init(OVNVCLReader *rdr, const uint8_t *bytestream_start,
              int bit_size)
{
    /* FIXME decide if stream check here */
    int buffer_size = (bit_size + 7) >> 3;

    rdr->bytestream       = bytestream_start;
    rdr->bytestream_end   = bytestream_start + buffer_size;
    rdr->size_in_bits     = bit_size;
    rdr->nb_bytes_read    = 0;

    fill_cache64(rdr);

    return 0;
}


