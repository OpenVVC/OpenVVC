#include <stdlib.h>
#include <stdint.h>

#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"

#include "libovvcdmx/ovunits.h"

#include "nvcl.h"
#include "ovvcdec.h"

static const char *const decname = "Open VVC Decoder";

struct OVVCSubDec;

/* Main decoder structure */
struct OVVCDec
{
    const char *name;

    /* NAL Units to be decoded
     * Corresponding to a Picture Unit
     */
    OVPictureUnit *nalu_list;

    /* Paramters sets context */
    OVNVCLCtx *nvcl_ctx;

    /* List of Sub Decoders
     * Contains context for Tile / Slice / Picture / SubPicture
     * decoding
     */
    struct OVSubDec *subdec_list;

    /* Informations on decoder behaviour transmitted by user
     */
    struct {
        int opt1;
    }options;
};

static int nvcl_unsupported();

int (*nvcl_read_table[32])() =
{
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported,
    nvcl_unsupported
};

static int
nvcl_unsupported()
{
    enum OVNALUType nal_type = 31;

    ov_log(NULL, 2, "Ignored unsupported NAL type : %d \n", nal_type);
    return 0;
}

    #if 0
static int
decode_nal_unit(VVCContext *const s, const VVCNAL *const nal)
{
    VVCLocalContext *const lc = s->lc_ctx;
    GetBitContext   *const gb = &s->gb;
    int ctb_addr_ts, ret;
    VVCSliceHeaderData *const sh = &s->sh;
    int nalu_type = nal->type;

    /*FIXME do we really need to set a lc->gb */
    *gb            = nal->gb;

    switch (nalu_type) {
    case VVC_NAL_VPS:
        ret = ff_vvc_decode_nal_vps(gb, s->avctx, &s->ps);
        if (ret < 0) {
            goto fail;
        }
        break;
    case VVC_NAL_SPS:
        ret = ff_vvc_decode_nal_sps(gb, s->avctx, &s->ps);
        if (ret < 0) {
            goto fail;
        }
        break;
    case VVC_NAL_PPS:
        ret = ff_vvc_decode_nal_pps(gb, s->avctx, &s->ps);
        if (ret < 0) {
            goto fail;
        }
        break;
    case VVC_NAL_PH:
        ret = ff_vvc_decode_nal_ph(gb, s->avctx, &s->ps);
        if (ret < 0) {
            goto fail;
        }
        break;
    case VVC_NAL_SUFFIX_APS:
    case VVC_NAL_PREFIX_APS:
        ret = ff_vvc_decode_nal_aps(gb, s->avctx, &s->ps);
        if (ret < 0) {
            goto fail;
        }
        break;
    case VVC_NAL_PREFIX_SEI:
    case VVC_NAL_SUFFIX_SEI:
        ret = ff_vvc_decode_nal_sei(gb, s, &s->sei, &s->ps, nalu_type);
            /*FIXME would be safer to handle SEI process inside SEI
             instead of guessing the suffix is a hash*/
            if ((s->chk_sei_md5 || s->display_md5) && nalu_type == 24 && s->is_decoded)
                ret = vvc_check_picture_hash(s, s->active_pic->frame);
        if (ret < 0)
            goto fail;
        break;
    case VVC_NAL_TRAIL:
    case VVC_NAL_STSA:
    case VVC_NAL_IDR_W_RADL:
    case VVC_NAL_IDR_N_LP:
    case VVC_NAL_CRA:
    case VVC_NAL_RASL:
    case VVC_NAL_RADL:

        memset(sh, 0, sizeof(VVCSliceHeaderData));

        ret = parse_slice_header(s, gb, sh, &s->ps, nalu_type);
        if (ret < 0)
            return ret;

        /*FIXME At the current time we consider first slice in picture is the
         * slice with address 0 this is OK since we fail before if we find
         * multiple slice in same acces unit
         * we might want to start decoding if the first slice is missing*/
        if (sh->slice_address == 0) {
            ret = vvc_dpb_start_frame(s, nalu_type);
            if (ret < 0) {
                goto fail;
            }

            av_log(s->avctx, AV_LOG_DEBUG, "Successfull start of frame decoder with POC: %d\n",
                    s->active_pic->poc);

        } else if (!s->active_pic) {
            /* FIXME we could still init a frame decode it and just warn*/
            av_log(s->avctx, AV_LOG_ERROR, "First slice in a frame missing.\n");
            goto fail;
        }

        /*FIXME can be moved to dpb_frame_start ?*/
        ff_thread_finish_setup(s->avctx);

        ctb_addr_ts = hls_slice_data(s, sh, nal);

        if (ctb_addr_ts < 0) {
            ret = ctb_addr_ts;
            goto fail;
        }

        break;
    case VVC_NAL_EOS:
    case VVC_NAL_EOB:
        s->active_seq_id = (s->active_seq_id + 1) & 0xff;
        s->max_ra     = INT_MAX;
        break;
    case VVC_NAL_AUD:
        break;
    default:
        av_log(s->avctx, AV_LOG_TRACE,
               "Skipping unknown NAL unit type %d\n", nalu_type);
    }

    return 0;

fail:
    if (s->active_pic) {
        if (s->threads_type & FF_THREAD_FRAME) {
            ff_thread_report_progress(&s->active_pic->tf, INT_MAX, 0);
        }
        ff_vvc_unref_frame(s, s->active_pic, 0);
        s->active_pic = NULL;
    }
    return ret;
    return 0;
}
#endif


#if 0
static int
vvc_validate_access_unit(VVCContext *const s, const VVCPacket *const pkt)
{
    int eos_at_start = 1;
    int nb_vcl_nal = 0;
    int first_vcl_type = 0; /* Zero is non vcl Init value */
    int i;

    s->last_eos = s->eos;
    s->eos = 0;

    /*FIXME we could add a more in depth conformance check
     * e.g. for VCL and non VCL NALU order etc.
     * we coud also filter some NALUs out here if we 'd like
     * to ignore some NALUs type given decoder state
    */
    for (i = 0; i < pkt->nb_nals; i++) {
        int nal_type = pkt->nals[i].type;
        /*FIXME Check EOS is handled correctely afterwards
         * last_eos should trigger an increase in sequence_id
         * whereas eos is passed to next thread context
         */
        if (nal_type == VVC_NAL_EOB || nal_type == VVC_NAL_EOS) {
            if (eos_at_start) {
                s->last_eos = 1;
            } else {
                s->eos = 1;
            }
        } else {
            eos_at_start = 0;
        }

        if (nal_type <= 12) {
            if (!nb_vcl_nal)
                first_vcl_type = nal_type;
            else if (nal_type != first_vcl_type) {
                av_log(s, AV_LOG_ERROR, "Received NALUs of different types in same access unit\n");
                return AVERROR_INVALIDDATA;
            }
            ++nb_vcl_nal;
        }
    }

    if (nb_vcl_nal > 1) {
        /*FIXME we could probably decode first slice */
        av_log(s, AV_LOG_ERROR, "Multiple slice in same picture, this is unsupported\n");
        return AVERROR_INVALIDDATA;
    }

    return 0;
}
#endif

#if 0
static int
vvc_decode_access_unit(VVCContext *s, const uint8_t *buf, int length)
{
    int i, ret;
    VVCPacket *const pkt = &s->pkt;

    s->active_pic = NULL;

    /* split the input packet into NAL units, so we know the upper bound on the
     * number of slices in the frame */
    ret = ff_vvc_packet_split(pkt, buf, length, s->avctx, s->is_nalff,
                                s->nal_length_size, s->avctx->codec_id, 0);
    if (ret < 0) {
        av_log(s->avctx, AV_LOG_ERROR,
               "Error splitting the input into NAL units.\n");
        return ret;
    }

    ret = vvc_validate_access_unit(s, pkt);
    if (ret < 0) {
        return ret;
    }

    /* decode the NAL units */
    for (i = 0; i < pkt->nb_nals; i++) {
        const VVCNAL *const nal = &pkt->nals[i];

        ret = decode_nal_unit(s, nal);
        if (ret < 0) {
            av_log(s->avctx, AV_LOG_WARNING,
                   "Error parsing NAL unit #%d.\n", i);
            goto fail;
        }
    }

fail:
    if (s->active_pic && s->threads_type == FF_THREAD_FRAME)
        ff_thread_report_progress(&s->active_pic->tf, INT_MAX, 0);

    return ret;
}
#endif

int
ovdec_submit_picture_unit(OVVCDec *vvcdec, uint8_t *buff, size_t buff_size)
{
    int ret;

    #if 0
    if (!vvcdec->dmx) {
        OVVCDmx dmx;
        dmx = ovvcdmx_create();
        if (!dmx) {
            /* Failed to allocate AnnexB demuxer*/
            return -1;
        }
    }

    ret = ovcdmx_extract_pu_nal_units(OVVCDmx *dmx, OVVCPUPacket *pkt);
    if (ret < 0) {
        return ret;
    }

    ret = vvc_decode_picture_unit(vvcdec, pkt);
    #endif

    return ret;
}

int
ovdec_init(OVVCDec **vvcdec)
{
    *vvcdec = ov_mallocz(sizeof(OVVCDec));

    if (*vvcdec == NULL) goto fail;

    (*vvcdec)->name = decname;

    return 0;

fail:
    /* TODO proper error management (ENOMEM)*/
    return -1;
}

int
ovdec_close(OVVCDec *vvcdec)
{
    int not_dec;
    if (vvcdec != NULL) {

        not_dec = vvcdec->name != decname;

        if (not_dec) goto fail;

        ov_free(vvcdec);

        return 0;
    }

fail:
    ov_log(vvcdec, 3, "Trying to close a something not a decoder.\n");
    return -1;
}
