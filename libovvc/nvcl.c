#include <string.h>
#include "ovmem.h"
#include "ovutils.h"

#include "nvcl.h"
#include "nvcl_utils.h"
#include "hls_structures.h"

#define NB_ARRAY_ELEMS(x) sizeof(x)/sizeof(*(x))


typedef int (*NALUnitAction)(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr,
                             uint8_t nalu_type);


static const char *nalu_name[32] =
{
    "TRAIL",
    "STSA",
    "RADL",
    "RASL",
    "RSVD_VCL",
    "RSVD_VCL",
    "RSVD_VCL",
    "IDR_W_RADL",
    "IDR_N_LP",
    "CRA",
    "GDR",
    "RSVD_IRAP_VCL",
    "OPI",
    "DCI",
    "VPS",
    "SPS",
    "PPS",
    "PREFIX_APS",
    "SUFFIX_APS",
    "PH",
    "AUD",
    "EOS",
    "EOB",
    "PREFIX_SEI",
    "SUFFIX_SEI",
    "FD",
    "RSVD_NVCL",
    "RSVD_NVCL",
    "UNSPEC",
    "UNSPEC",
    "UNSPEC",
    "UNSPEC"
};

static const struct HLSReader todo;
extern const struct HLSReader sps_manager;

static const struct HLSReader *nalu_reader[32] =
{
    &todo                , /* TRAIL */
    &todo                , /* STSA */
    &todo                , /* RADL */
    &todo                , /* RASL */
    &todo                , /* RSVD_VCL */
    &todo                , /* RSVD_VCL */
    &todo                , /* RSVD_VCL */
    &todo                , /* IDR_W_RADL */
    &todo                , /* IDR_N_LP */
    &todo                , /* CRA */
    &todo                , /* GDR */
    &todo                , /* RSVD_IRAP_VCL */
    &todo                , /* OPI */
    &todo                , /* DCI */
    &todo                , /* VPS */
    &sps_manager         , /* SPS */
    &todo                , /* PPS */
    &todo                , /* PREFIX_APS */
    &todo                , /* SUFFIX_APS */
    &todo                , /* PH */
    &todo                , /* AUD */
    &todo                , /* EOS */
    &todo                , /* EOB */
    &todo                , /* PREFIX_SEI */
    &todo                , /* SUFFIX_SEI */
    &todo                , /* FD */
    &todo                , /* RSVD_NVCL */
    &todo                , /* RSVD_NVCL */
    &todo                , /* UNSPEC */
    &todo                , /* UNSPEC */
    &todo                , /* UNSPEC */
    &todo                  /* UNSPEC */
};


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

static int
decode_nalu_hls_data(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr,
                     const struct HLSReader *const hls_hdl)
{
    const union HLSData **storage = hls_hdl->find_storage(rdr, nvcl_ctx);
    union HLSData data;
    int ret;

    if (*storage) {
        /* TODO compare RBSP data to avoid new read */
        uint8_t identical_rbsp = 0;
        if (identical_rbsp) goto duplicated;
    }

    memset(&data, 0, hls_hdl->data_size);

    ov_log(NULL, OVLOG_TRACE, "Reading new %s\n", hls_hdl->name);

    ret = hls_hdl->read(rdr, &data, nvcl_ctx);
    if (ret < 0)  goto failread;

    ret = hls_hdl->validate(rdr, &data);
    if (ret < 0)  goto invalid;

    ret = hls_hdl->replace(hls_hdl, storage, &data);

    return ret;

invalid:
    ov_log(NULL, OVLOG_ERROR, "Invalid %s\n", hls_hdl->name);
    return ret;

failread:
    ov_log(NULL, OVLOG_ERROR, "Error while reading %s\n", hls_hdl->name);
    return ret;

duplicated:
    ov_log(NULL, OVLOG_TRACE, "Ignored duplicated %s\n", hls_hdl->name);
    return 0;

}

static int warn_unsupported(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr,
                            uint8_t nalu_type)
{
    ov_log(NULL, OVLOG_WARNING, "Unsupported %s NAL unit.\n", nalu_name[nalu_type]);
    return 0;
}

static int warn_unspec(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr,
                       uint8_t nalu_type)
{
    ov_log(NULL, OVLOG_WARNING, "Unspec %s NAL unit.\n", nalu_name[nalu_type]);
    return 0;
}

static int log_ignored(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr,
                       uint8_t nalu_type)
{
    ov_log(NULL, OVLOG_TRACE, "Ignored %s NAL unit.\n", nalu_name[nalu_type]);
    return 0;
}

static int decode_nvcl_hls(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr,
                           uint8_t nalu_type)
{
    const struct HLSReader *const hls_reader = nalu_reader[nalu_type];
    if (hls_reader != &todo)
        return decode_nalu_hls_data(nvcl_ctx, rdr, hls_reader);
    return 0;
}

static const NALUnitAction nalu_action[32] =
{
    &log_ignored                , /* TRAIL */
    &log_ignored                , /* STSA */
    &log_ignored                , /* RADL */
    &log_ignored                , /* RASL */
    &warn_unspec                , /* RSVD_VCL */
    &warn_unspec                , /* RSVD_VCL */
    &warn_unspec                , /* RSVD_VCL */
    &log_ignored                , /* IDR_W_RADL */
    &log_ignored                , /* IDR_N_LP */
    &log_ignored                , /* CRA */
    &log_ignored                , /* GDR */
    &warn_unspec                , /* RSVD_IRAP_VCL */
    &warn_unsupported           , /* OPI */
    &warn_unsupported           , /* DCI */
    &warn_unsupported           , /* VPS */
    &decode_nvcl_hls            , /* SPS */
    &log_ignored                , /* PPS */
    &log_ignored                , /* PREFIX_APS */
    &log_ignored                , /* SUFFIX_APS */
    &log_ignored                , /* PH */
    &log_ignored                , /* AUD */
    &log_ignored                , /* EOS */
    &log_ignored                , /* EOB */
    &log_ignored                , /* PREFIX_SEI */
    &log_ignored                , /* SUFFIX_SEI */
    &log_ignored                , /* FD */
    &warn_unspec                , /* RSVD_NVCL */
    &warn_unspec                , /* RSVD_NVCL */
    &warn_unspec                , /* UNSPEC */
    &warn_unspec                , /* UNSPEC */
    &warn_unspec                , /* UNSPEC */
    &warn_unspec                  /* UNSPEC */
};

int
nvcl_decode_nalu_hls_data(OVNVCLCtx *const nvcl_ctx, OVNALUnit *const nalu)
{
    OVNVCLReader rdr;

    nvcl_reader_init(&rdr, nalu->rbsp_data, nalu->rbsp_size);

    /* FIXME properly read NAL Unit header */
    nvcl_skip_bits(&rdr, 16);

    uint8_t nalu_type = nalu->type & 0x1F;

    return nalu_action[nalu_type](nvcl_ctx, &rdr, nalu_type);
}
