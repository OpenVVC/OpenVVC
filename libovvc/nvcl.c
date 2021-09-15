#include "ovmem.h"
#include "ovutils.h"

#include "nvcl.h"
#include "nvcl_utils.h"

#define NB_ARRAY_ELEMS(x) sizeof(x)/sizeof(*(x))


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

{
}
