#include "ovmem.h"

#include "nvcl.h"
#include "nvcl_utils.h"

#define NB_ARRAY_ELEMS(x) sizeof(x)/sizeof(*(x))


/* FIXME give size in bytes instead and find SODB end
 * (RBSP stop bit) here instead
 */
int
nvcl_reader_init(OVNVCLReader *rdr, const uint8_t *bytestream_start,
                 int bit_size)
{
    int buffer_size = (bit_size + 7) >> 3;

    rdr->bytestream_start = bytestream_start;
    rdr->bytestream_end   = bytestream_start + buffer_size;
    rdr->bytestream       = bytestream_start;

    fill_cache64(rdr);

    /* FIXME properly read NAL Unit header */
    nvcl_skip_bits(rdr, 16);


    return 0;
}

void
nvcl_free_ctx(OVNVCLCtx *const nvcl_ctx)
{
    int i;
    int nb_elems = NB_ARRAY_ELEMS(nvcl_ctx->sps_list);
    for (i = 0; i < nb_elems; ++i) {
        if (nvcl_ctx->sps_list[i]) {
            ov_freep(&nvcl_ctx->sps_list[i]);
        }
    }

    nb_elems = NB_ARRAY_ELEMS(nvcl_ctx->pps_list);
    for (i = 0; i < nb_elems; ++i) {
        if (nvcl_ctx->pps_list[i]) {
            ov_freep(&nvcl_ctx->pps_list[i]);
        }
    }

    nb_elems = NB_ARRAY_ELEMS(nvcl_ctx->alf_aps_list);
    for (i = 0; i < nb_elems; ++i) {
        if (nvcl_ctx->alf_aps_list[i]) {
            ov_freep(&nvcl_ctx->alf_aps_list[i]);
        }
    }

    nb_elems = NB_ARRAY_ELEMS(nvcl_ctx->lmcs_aps_list);
    for (i = 0; i < nb_elems; ++i) {
        if (nvcl_ctx->lmcs_aps_list[i]) {
            ov_freep(&nvcl_ctx->lmcs_aps_list[i]);
        }
    }

    if (nvcl_ctx->ph) {
        ov_freep(&nvcl_ctx->ph);
    }

    if (nvcl_ctx->sh) {
        ov_freep(&nvcl_ctx->sh);
    }

    if (nvcl_ctx->sei) {
        nvcl_free_sei_params(nvcl_ctx->sei);
    }

}

uint32_t
nvcl_nb_bytes_read(const OVNVCLReader *const rdr)
{
    ptrdiff_t nb_bytes_read = rdr->bytestream - rdr->bytestream_start;
    return nb_bytes_read - ((rdr->nb_cached_bits + 7) >> 3);

}
