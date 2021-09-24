#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "sys/resource.h"

#include "ovutils.h"
#include "ovmem.h"
#include "overror.h"

#include "ovunits.h"

#include "nvcl.h"
#include "ovdec.h"
#include "ctudec.h"
#include "decinit.h"
#include "slicedec.h"
#include "dec_structures.h"
#include "ovdec_internal.h"
#include "post_proc.h"
/* FIXME
 * To be removed includes
 */
#include "ovdpb.h"
#include "ovthreads.h"

static const char *const decname = "Open VVC Decoder";

static const char *option_names[OVDEC_NB_OPTIONS] =
{
    "frame threads",
    "entry threads",
    "display_output"
};

static int
ovdec_init_subdec_list(OVVCDec *dec)
{
    int ret;
    ov_log(NULL, OVLOG_TRACE, "Creating %d Slice decoders\n", dec->nb_frame_th);
    if (!dec->subdec_list) 
        dec->subdec_list = ov_mallocz(sizeof(OVSliceDec*) * dec->nb_frame_th);

    for (int i = 0; i < dec->nb_frame_th; ++i){
        dec->subdec_list[i] = ov_mallocz(sizeof(OVSliceDec));
        ret = slicedec_init(dec->subdec_list[i], dec->nb_entry_th);
        if (ret < 0) {
            return OVVC_ENOMEM;
        }
        dec->subdec_list[i]->th_slice.main_thread = &dec->main_thread;
        dec->subdec_list[i]->th_slice.output_thread = &dec->output_thread;
    }

    return 0;
}

static int
init_vcl_decoder(OVVCDec *const dec, OVSliceDec *sldec, const OVNVCLCtx *const nvcl_ctx,
                OVNALUnit * nalu, uint32_t nb_sh_bytes)
{

    int ret;

    ret = decinit_update_params(&dec->active_params, nvcl_ctx);
    if (ret < 0) {
        ov_log(dec, OVLOG_ERROR, "Failed to activate parameters\n");
        return ret;
    }

    if (!dec->dpb) {
         ret = ovdpb_init(&dec->dpb, &dec->active_params);
        dec->dpb->display_output = dec->display_output;
         if (ret < 0) {
             return ret;
         }
    }

    /* FIXME only if TMVP */
    //TODOpar: protect mv pool when more than one thread ?
    if (!dec->mv_pool) {
        ret = mvpool_init(&dec->mv_pool, &dec->active_params.pic_info);
    }

    //Temporary: copy active parameters
    slicedec_copy_params(sldec, &dec->active_params);

    if (!dec->active_params.sh->sh_slice_address) {
        ret = ovdpb_init_picture(dec->dpb, &sldec->pic, sldec->active_params, nalu->type, sldec, dec);
        if (ret < 0) {
            return ret;
        }
    } else {
        /* FIXME clean way on new slice with address 0 */
        ov_log(dec, OVLOG_ERROR, "Multiple slices in pic is not supported yet.\n");
    }

    ov_nalu_new_ref(&sldec->th_slice.slice_nalu, nalu);

    /*FIXME return checks */
    ret = slicedec_init_lines(sldec, sldec->active_params);

    ret = decinit_set_entry_points(sldec->active_params, nalu, nb_sh_bytes);

    ret = slicedec_update_entry_decoders(sldec, sldec->active_params);

    return 0;
}

OVSliceDec *
ovdec_select_subdec(OVVCDec *const dec)
{
    #if USE_THREADS
    OVSliceDec **sldec_list = dec->subdec_list;
    int nb_threads = dec->nb_frame_th;
    struct MainThread* th_main = &dec->main_thread;

    OVSliceDec * slicedec;
    struct SliceThread* th_slice;
    do{
        int min_idx_available = nb_threads;
        pthread_mutex_lock(&th_main->main_mtx);
        
        for(int i = nb_threads-1; i >= 0 ; i--){
            slicedec = sldec_list[i];
            th_slice = &slicedec->th_slice;
            
            //Unmark ref pict lists of decoded pics
            pthread_mutex_lock(&th_slice->gnrl_mtx);
            if(th_slice->active_state == DECODING_FINISHED){
                pthread_mutex_unlock(&th_slice->gnrl_mtx);
                min_idx_available = i;
                OVPicture *slice_pic = slicedec->pic;
                if(slice_pic && (slice_pic->flags & OV_IN_DECODING_PIC_FLAG)){
                    ov_log(NULL, OVLOG_TRACE, "Subdec %d Remove DECODING_PIC_FLAG POC: %d\n", min_idx_available, slice_pic->poc);
                    ovdpb_unref_pic(slice_pic, OV_IN_DECODING_PIC_FLAG);
                    ovdpb_unmark_ref_pic_lists(slicedec->slice_type, slice_pic);

                    pthread_mutex_lock(&th_slice->gnrl_mtx);
                    th_slice->active_state = IDLE;
                    pthread_mutex_unlock(&th_slice->gnrl_mtx);
                }
            } else if (th_slice->active_state == IDLE) {
                pthread_mutex_unlock(&th_slice->gnrl_mtx);
                min_idx_available = i;     
            } else {
                pthread_mutex_unlock(&th_slice->gnrl_mtx);
            }
        }

        if(min_idx_available < nb_threads){
            slicedec = sldec_list[min_idx_available];
            th_slice = &slicedec->th_slice;
        
            pthread_mutex_lock(&th_slice->gnrl_mtx);
            th_slice->active_state = ACTIVE;
            ov_log(NULL, OVLOG_TRACE, "Subdec %d selected\n", min_idx_available);
            pthread_mutex_unlock(&th_slice->gnrl_mtx);

            pthread_mutex_unlock(&th_main->main_mtx);
            return slicedec;
        }

        pthread_cond_wait(&th_main->main_cnd, &th_main->main_mtx);
        pthread_mutex_unlock(&th_main->main_mtx);

    } while(!th_main->kill);

    return NULL;

    #else
    OVSliceDec *sldec = vvcdec->subdec_list[0];
    #endif
}

static int
decode_nal_unit(OVVCDec *const vvcdec, OVNALUnit * nalu)
{
    OVNVCLReader rdr;
    OVNVCLCtx *const nvcl_ctx = &vvcdec->nvcl_ctx;
    enum OVNALUType nalu_type = nalu->type;

    int ret;

    nvcl_reader_init(&rdr, nalu->rbsp_data, nalu->rbsp_size);

    /* FIXME properly read NAL Unit header */
    nvcl_skip_bits(&rdr, 16);

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
            /* Select the first available subdecoder, or wait until one is available */
            OVSliceDec *sldec = ovdec_select_subdec(vvcdec);

            uint32_t nb_sh_bytes = nvcl_nb_bytes_read(&rdr);
                
            /* Beyond this point unref current picture on failure */
            ret = init_vcl_decoder(vvcdec, sldec, nvcl_ctx, nalu, nb_sh_bytes);

            if (ret < 0) {
                slicedec_finish_decoding(sldec);
                goto failvcl;
            }

            /* FIXME handle non rect entries later */
            ret = slicedec_decode_rect_entries(sldec, sldec->active_params);
        }

        break;
    case OVNALU_PREFIX_SEI:
    case OVNALU_SUFFIX_SEI:
    case OVNALU_SUFFIX_APS:
    case OVNALU_PREFIX_APS:
    case OVNALU_VPS:
    case OVNALU_SPS:
    case OVNALU_PPS:
    case OVNALU_PH:
    case OVNALU_EOS:
    case OVNALU_EOB:
    case OVNALU_AUD:
    default:
        ret = nvcl_decode_nalu_hls_data(nvcl_ctx, nalu);
        if (ret < 0) {
            goto fail;
        }
        break;
    }

    return 0;

fail:
    return ret;

failvcl:
    //TODOpar: change 0
    if (vvcdec->subdec_list[0]->pic) {
        ovdpb_unref_pic(vvcdec->subdec_list[0]->pic, ~0);
    }
    return ret;
}

static int
vvc_decode_picture_unit(OVVCDec *dec, const OVPictureUnit *pu)
{
    int i;
    int ret;
    for (i = 0; i < pu->nb_nalus; ++i) {
        ret = decode_nal_unit(dec, pu->nalus[i]);
        if (ret < 0) {
            goto fail;
        }
    }
    return 0;

fail:
    /* Error processing if needed */
    return ret;
}

int
ovdec_submit_picture_unit(OVVCDec *vvcdec, const OVPictureUnit *const pu)
{
    int ret = 0;

    ret = vvc_decode_picture_unit(vvcdec, pu);

    return ret;
}

int
ovdec_receive_picture(OVVCDec *dec, OVFrame **frame_p)
{
    /* FIXME this is temporary request output from DPB
     * instead
     */
    OVDPB *dpb = dec->dpb;
    uint16_t out_cvs_id;
    int ret = 0;

    if (!dpb) {
        ov_log(dec, OVLOG_ERROR, "No DPB on output request.\n");
        /* FIXME new return value */
        return 0;
    }

    OVPicture *pic = NULL;
    out_cvs_id = (dpb->cvs_id - 1) & 0xFF;
    ret = ovdpb_output_pic(dpb, &pic, out_cvs_id);

    if (pic) {
        *frame_p = pic->frame;
        pp_process_frame(pic->sei, dec->dpb, frame_p);

        //New ref if it is a frame already in a DPB pic
        if(*frame_p ==  pic->frame){
            ovframe_new_ref(frame_p, pic->frame);
        }
        /* we unref the picture even if ref failed the picture
         * will still be usable by the decoder if not bumped
         * */
        ovdpb_unref_pic(pic, OV_OUTPUT_PIC_FLAG | (pic->flags & OV_BUMPED_PIC_FLAG));
    } else {
      *frame_p = NULL;
    }

    return ret;
}

int
ovdec_drain_picture(OVVCDec *dec, OVFrame **frame_p)
{
    /* FIXME this is temporary request output from DPB
     * instead
     */
    OVDPB *dpb = dec->dpb;
    uint16_t out_cvs_id;
    int ret;

    ovdec_uninit_subdec_list(dec);

    if (!dpb) {
        ov_log(dec, OVLOG_ERROR, "No DPB on output request.\n");
        /* FIXME new return value */
        return OVVC_EINDATA;
    }

    OVPicture *pic = NULL;
    out_cvs_id = (dpb->cvs_id - 1) & 0xFF;
    ret = ovdpb_drain_frame(dpb, &pic, out_cvs_id);

    if (pic) {
        *frame_p = pic->frame;
        pp_process_frame(pic->sei, dec->dpb, frame_p);

        //New ref if it is a frame already in a DPB pic
        if(*frame_p ==  pic->frame){
            ovframe_new_ref(frame_p, pic->frame);
        }
        /* we unref the picture even if ref failed the picture
         * will still be usable by the decoder if not bumped
         * */
        ovdpb_unref_pic(pic, OV_OUTPUT_PIC_FLAG | (pic->flags & OV_BUMPED_PIC_FLAG));
    }

    return ret;
}

static int
set_display_output(OVVCDec *ovdec, int on_off)
{
    ovdec->display_output = !!on_off;
    return 0;
}

static int
set_nb_entry_threads(OVVCDec *ovdec, int nb_threads)
{
    ovdec->nb_entry_th = nb_threads;

    return 0;
}

static int
set_nb_frame_threads(OVVCDec *ovdec, int nb_threads)
{
    ovdec->nb_frame_th = nb_threads;

    return 0;
}

int
ovdec_set_option(OVVCDec *ovdec, enum OVOptions opt_id, int value)
{

    switch (opt_id){
        case OVDEC_NB_ENTRY_THREADS:
            set_nb_entry_threads(ovdec, value);
            break;
        case OVDEC_NB_FRAME_THREADS:
            set_nb_frame_threads(ovdec, value);
            break;
        case OVDEC_DISPLAY_OUTPUT:
            set_display_output(ovdec, value);
            break;
        default :
            if (opt_id < OVDEC_NB_OPTIONS) {
                ov_log(ovdec, OVLOG_ERROR, "Invalid option id %d.", opt_id);
                return OVVC_EINDATA;
            }
            break;
    }
    ov_log(ovdec, OVLOG_VERBOSE, "Option %s set to %d.\n", option_names[opt_id], value);

    return 0;
}

int
ovdec_init(OVVCDec **vvcdec, int display_output, int nb_frame_th, int nb_entry_th)
{

    if (nb_frame_th < 0) {
        nb_frame_th = 1;
    } else if (nb_frame_th == 0) {
        /* FIXME might not be available on every plateform */
        nb_frame_th = get_number_of_cores();
    }

    if (nb_entry_th < 1) {
        nb_entry_th = 1;
    }

    *vvcdec = ov_mallocz(sizeof(OVVCDec));

    if (*vvcdec == NULL) goto fail;

    (*vvcdec)->name = decname;

    ovdec_set_option(*vvcdec, OVDEC_NB_FRAME_THREADS, nb_frame_th);

    ovdec_set_option(*vvcdec, OVDEC_NB_ENTRY_THREADS, nb_entry_th);

    ovdec_set_option(*vvcdec, OVDEC_DISPLAY_OUTPUT, display_output);

    (*vvcdec)->display_output = !!display_output;

    ovdec_init_subdec_list(*vvcdec);

    ov_log(NULL, OVLOG_TRACE, "OpenVVC init at %p\n", *vvcdec);
    return 0;

fail:
    ov_log(NULL, OVLOG_ERROR, "Failed OpenVVC init\n");
    /* TODO proper error management (ENOMEM)*/
    return -1;
}

void
ovdec_uninit_subdec_list(OVVCDec *vvcdec)
{
    OVSliceDec *sldec;

    if (vvcdec != NULL)
    {
        if (vvcdec->subdec_list) {
            for (int i = 0; i < vvcdec->nb_frame_th; ++i){
                sldec = vvcdec->subdec_list[i];
                slicedec_uninit(&sldec);
                ov_log(NULL, OVLOG_INFO, "Main joined thread: %d\n", i);
            }
            ov_freep(&vvcdec->subdec_list);
        }

        // ovthread_output_uninit(&vvcdec->output_thread);
    }
}

int
ovdec_close(OVVCDec *vvcdec)
{
    int not_dec;
    if (vvcdec != NULL) {

        not_dec = vvcdec->name != decname;

        if (not_dec) goto fail;

        nvcl_free_ctx(&vvcdec->nvcl_ctx);

        ovdec_uninit_subdec_list(vvcdec);

        ovdpb_uninit(&vvcdec->dpb);

        if (vvcdec->mv_pool) {
            mvpool_uninit(&vvcdec->mv_pool);
        }

        pthread_mutex_destroy(&vvcdec->main_thread.main_mtx);
        pthread_cond_destroy(&vvcdec->main_thread.main_cnd);
        ov_free(vvcdec);

        return 0;
    }

fail:
    ov_log(vvcdec, 3, "Trying to close a something not a decoder.\n");
    return -1;
}

void ovdec_set_log_callback(void (*log_function)(void* ctx, int log_level, const char* log_content, va_list vl))
{
    set_log_callback(log_function);
}
