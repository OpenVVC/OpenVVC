#include <stdlib.h>
#include <stdint.h>

#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"
#include "libovvcutils/ovvcerror.h"

#include "libovvcdmx/ovunits.h"

#include "nvcl.h"
#include "ovvcdec.h"
#include "ctudec.h"
#include "decinit.h"
#include "slicedec.h"
#include "dec_structures.h"

static const char *const decname = "Open VVC Decoder";

struct OVVCSubDec;


/* Actions on unsupported NAL Unit types */
static int nalu_type_unsupported(enum OVNALUType nalu_type);

#if 0
/* Init / update vcl sub decoders based on slice header info */
static int init_vcl_decoder(const OVVCDec *dec, OVSubDec *sub_dec);
#endif

static int
nalu_type_unsupported(enum OVNALUType nalu_type)
{
    ov_log(NULL, 2, "Ignored unsupported NAL type : %d \n", nalu_type);

    return 0;
}

#if 0
/* TODO check Slice / Tiles / Sub Picture context */
static int
init_slice_decoder(OVVCDec *const dec, const OVNVCLCtx *const nvcl_ctx)
{
    const OVSH *const sh = nvcl_ctx->sh;
    int ret;

    #if 0
    init_slice_tools();
    #endif

    #if 0
    if (sh->nb_entry_points > 1){
        /* TILES or WPP */
        /* TODO loop other entry_points offsets and map entries to
         * tile / slice / sub picture decoder */
    } else {
        /* Only one slice in NAL Unit */
    }
    #endif
    return 0;
}
#endif

static int
init_subdec_list(OVVCDec *dec)
{
    int ret;

    /* TODO */
    if (!dec->subdec_list) {
        /* FIXME
         *    -at the current time we suppose only one slice decoder
         *    -select sub dec types according to active params sets
         *    -if parameters changed we might want to change sb_dec
         */
        ret = slicedec_init(&dec->subdec_list);
        if (ret < 0)
        {
            return OVVC_ENOMEM;
        }
    }

    return 0;
}


static OVSliceDec *
select_subdec(const OVVCDec *const dec)
{
    /* FIXME only one SLice dec at the current time */
    return dec->subdec_list;
}

static int
init_vcl_decoder(OVVCDec *const dec, const OVNVCLCtx *const nvcl_ctx)
{
    
    int ret;
    OVSliceDec *sldec;

    ret = decinit_update_params(dec, nvcl_ctx);
    if (ret < 0) {
        ov_log(dec, 3, "Failed to activate parameters\n");
        return ret;
    }

    #if 0
    if (dec->active_params.sh->sh_slice_address) {
    /* FIXME  dpb_request_picture */
    }
    #endif

    ret = init_subdec_list(dec);
    if (ret < 0) {
        return ret;
    }

    sldec = select_subdec(dec); 

    ret = slicedec_init_slice_tools(sldec, &dec->active_params);

    #if 0
    ret = decinti_update_subdec();
    if (ret < 0) {
        ov_log(dec, 3, "Failed to update sub decoder\n");
        return ret;
    }
    #endif

        /*TODO update DPB status */

        /* TODO init VCL decoder */

    #if 0
    ret = init_slice_decoder(dec, nvcl_ctx);
    #endif

    return 0;
}

static int
decode_nal_unit(OVVCDec *const vvcdec, const OVNALUnit *const nalu)
{
    OVNVCLReader rdr;
    OVNVCLCtx *const nvcl_ctx = &vvcdec->nvcl_ctx;
    enum OVNALUType nalu_type = nalu->type;
    /* FIXME add proper SH allocation */
    int ret;

    /* TODO init NVCLReader */
    nvcl_reader_init(&rdr, nalu->rbsp_data, (nalu->rbsp_size - 2) << 3);

    switch (nalu_type) {
    case OVNALU_TRAIL:
    case OVNALU_STSA:
    case OVNALU_RADL:
    case OVNALU_RASL:
    case OVNALU_IDR_W_RADL:
    case OVNALU_IDR_N_LP:
    case OVNALU_CRA:
    case OVNALU_GDR:

        ret = nvcl_decode_nalu_sh(&rdr, nvcl_ctx, nalu_type);

        if (ret < 0) {
            return ret;
        } else {
            #if 0
            int nb_sh_bytes = ret;
            #endif

           ret = init_vcl_decoder(vvcdec, nvcl_ctx);
           if (ret < 0) {
               goto fail;
           }

           #if 0
           ret = vcl_dec_decode_slice(vvcdec);
           #endif


        /* TODO start VCL decoder */
        }

        break;
    case OVNALU_VPS:
        ret = 0;
        if (ret < 0) {
            goto fail;
        }
        break;
    case OVNALU_SPS:
        ret = nvcl_decode_nalu_sps(&rdr, nvcl_ctx);
        if (ret < 0) {
            goto fail;
        }
        break;
    case OVNALU_PPS:
        ret = nvcl_decode_nalu_pps(&rdr, nvcl_ctx);
        if (ret < 0) {
            goto fail;
        }
        break;
    case OVNALU_PH:
        ret = nvcl_decode_nalu_ph(&rdr, nvcl_ctx);
        if (ret < 0) {
            goto fail;
        }
        break;
    case OVNALU_SUFFIX_APS:
    case OVNALU_PREFIX_APS:
        ret = 0;
        if (ret < 0) {
            goto fail;
        }
        break;
    case OVNALU_PREFIX_SEI:
    case OVNALU_SUFFIX_SEI:
        ret = 0;
        if (ret < 0)
            goto fail;
        break;
    case OVNALU_EOS:
    case OVNALU_EOB:
        /* TODO update DPB status (new cvs); 
         * call dpb_uninit
         */
        break;
    case OVNALU_AUD:
        /* AUD is ignored it should be the first NALU if
         * present
         */
        break;
    default:
        nalu_type_unsupported(nalu_type);
    }

    return 0;

fail:
    /* TODO error hanling */
    return ret;
}

static int
vvc_decode_picture_unit(OVVCDec *dec, const OVPictureUnit *pu)
{
    int i;
    int ret;
    for (i = 0; i < pu->nb_nalus; ++i) {
        ret = decode_nal_unit(dec, &pu->nalus[i]);
        if (ret < 0) {
            goto fail;
        }
    }
    return 0;

fail:
    /* Error processing if needed */
    return ret;
}

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
        if (nal_type == OVNALU_EOB || nal_type == OVNALU_EOS) {
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
ovdec_submit_picture_unit(OVVCDec *vvcdec, const OVPictureUnit *const pu)
{
    int ret = 0;

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
    #endif

    ret = vvc_decode_picture_unit(vvcdec, pu);

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

        nvcl_free_ctx(&vvcdec->nvcl_ctx);

        slicedec_uninit(&vvcdec->subdec_list);

        ov_free(vvcdec);

        return 0;
    }

fail:
    ov_log(vvcdec, 3, "Trying to close a something not a decoder.\n");
    return -1;
}
