/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the 
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from 
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 * 
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *   
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 * 
 **/

#include <string.h>
#include <stddef.h>
#include <limits.h>
#include "ovutils.h"
#include "ovmem.h"
#include "nvcl_private.h"
#include "nvcl_structures.h"
#include "ovframe.h"
#include "ovunits.h"
#include "ovdpb.h"
#include "dec_structures.h"
#include "ovdpb_internal.h"
#include "overror.h"
#include "slicedec.h"
#include "ovdec_internal.h"

/* FIXME More global scope for this enum
 */
enum SliceType
{
   SLICE_B = 0,
   SLICE_P = 1,
   SLICE_I = 2,
};

const uint8_t SIZE_INT64 = 6;

static void tmvp_release_mv_planes(OVPicture *const pic);

static int dpb_init_params(OVDPB *dpb, OVDPBParams const *prm);

static void ovdpb_reset_decoded_ctus(OVPicture *const pic);

static void ovdpb_init_decoded_ctus(OVPicture *const pic, const OVPS *const ps);

static void ovdpb_uninit_decoded_ctus(OVPicture *const pic);

int
ovdpb_init(OVDPB **dpb_p, const OVPS *ps)
{
    #if 0
    OVDPB *dpb = *dpb_p;
    #endif
    int ret;

    *dpb_p = ov_mallocz(sizeof(**dpb_p));
    if (!*dpb_p) {
         ov_log(NULL, OVLOG_ERROR, "Failed DBP allocation attempt\n");
         return OVVC_ENOMEM;
    }

    ret = dpbpriv_init_framepool(&(*dpb_p)->internal, ps->sps);
    if (ret < 0) {
        goto failframepool;
    }

    /* FIXME handle temporal and sub layers*/
    dpb_init_params(*dpb_p, &ps->sps->dpb_parameters[ps->sps->sps_max_sublayers_minus1]);

    OVPicture *pic;
    int nb_dpb_pic = sizeof((*dpb_p)->pictures) / sizeof(*pic);
    for (int j = 0; j < nb_dpb_pic; j++) {
        pic = &(*dpb_p)->pictures[j];
        ovdpb_init_decoded_ctus(pic, ps);
    }


    return 0;

failframepool:
    ov_freep(dpb_p);
    return ret;
}

void
ovdpb_uninit(OVDPB **dpb_p)
{
    if (*dpb_p) {
        /* TODO
         * release all pics
         */
        ovdpb_flush_dpb(*dpb_p);
        dpbpriv_uninit_framepool(&(*dpb_p)->internal);

        ov_freep(dpb_p);
    }
}

/* Clean OVPicture from DPB and unref associated Frame */
static void
dpbpriv_release_pic(OVPicture *pic)
{
    if (pic->frame) {
        ov_log(NULL, OVLOG_TRACE, "Release Picture with POC %d\n", pic->poc);
        ovframe_unref(&pic->frame);

        /* FIXME better existence check */
        if (pic->mv_plane0.mvs){
            tmvp_release_mv_planes(pic);
        }

        pic->frame = NULL;

        if (pic->sei) {
            if (pic->sei->sei_fg) {
                ov_freep(&pic->sei->sei_fg);
            }
            if (pic->sei->sei_slhdr) {
                ov_freep(&pic->sei->sei_slhdr);
            }
            ov_freep(&pic->sei);
        }
    }
}

struct Rational {
    uint32_t num;
    uint32_t den;
};

static int
dpb_init_params(OVDPB *dpb, OVDPBParams const *prm)
{
    struct Rational framerate = {
        .num = 60,
        .den = 1
    };

    uint32_t time_scale         = 27000000;
    uint32_t nb_units_in_ticks  = framerate.den * time_scale / framerate.num;

    dpb->max_nb_dpb_pic       = prm->dpb_max_dec_pic_buffering_minus1 + 1;
    dpb->max_nb_reorder_pic   = prm->dpb_max_num_reorder_pics;
    dpb->max_latency_increase = prm->dpb_max_latency_increase_plus1 - 1;

    dpb->nb_units_in_ticks = nb_units_in_ticks;
    dpb->pts = 0;
    return 0;
}

static int
derive_poc(int poc_lsb, int log2_max_poc_lsb, int prev_poc)
{
    int max_poc_lsb  = 1 << log2_max_poc_lsb;
    int prev_poc_lsb = prev_poc & (max_poc_lsb - 1);
    int poc_msb = prev_poc - prev_poc_lsb;
    if((poc_lsb < prev_poc_lsb) &&
      ((prev_poc_lsb - poc_lsb) >= (max_poc_lsb >> 1))){
        poc_msb += max_poc_lsb;
    } else if((poc_lsb > prev_poc_lsb) &&
             ((poc_lsb - prev_poc_lsb) > (max_poc_lsb >> 1))){
        poc_msb -= max_poc_lsb;
    }
    return poc_msb + poc_lsb;
}

void
ovdpb_unref_pic(OVPicture *pic, int flags)
{
    /* pic->frame can be NULL if context init failed */
    if (!pic || (!pic->frame || !pic->frame->data[0]))
        return;

    pthread_mutex_lock(&pic->pic_mtx);
    pic->flags &= ~flags;

    pthread_mutex_unlock(&pic->pic_mtx);

    atomic_fetch_add_explicit(&pic->ref_count, -1, memory_order_acq_rel);
}

void
ovdpb_release_pic(OVDPB *dpb, OVPicture *pic)
{
    /* pic->frame can be NULL if context init failed */
    if (!pic->frame || !pic->frame->data[0])
        return;

    uint16_t ref_count = atomic_load(&pic->ref_count);

    /* If there is no more flags the picture can be
     * returned to the DPB;
     */
    pthread_mutex_lock(&pic->pic_mtx);
    if (!pic->flags && !ref_count) {
        /* Release TMVP  MV maps */
        dpbpriv_release_pic(pic);
    }
    pthread_mutex_unlock(&pic->pic_mtx);
}

int
ovdpb_new_ref_pic(OVPicture *pic, int flags)
{
    atomic_fetch_add_explicit(&pic->ref_count, 1, memory_order_acq_rel);

    pthread_mutex_lock(&pic->pic_mtx);
    pic->flags |= flags;

    pthread_mutex_unlock(&pic->pic_mtx);

    return 0;
}

static int
find_min_cvs_id(const OVDPB *const dpb)
{
    /* Count pictures in current output target Coded Video Sequence
     */
    int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    int min_cvs_id = INT_MAX;
    int i;

    for (i = 0; i < nb_dpb_pic; i++) {
        const OVPicture *pic = &dpb->pictures[i];
        if (pic->frame && pic->frame->data[0]) {
            if (pic->cvs_id < min_cvs_id) {
                min_cvs_id = pic->cvs_id;
            }
        }
    }

    if (min_cvs_id == 0) {
        int cvs_id_min1 = 0xFF;
        int got_pic = 0;
        do {
            for (i = 0; i < nb_dpb_pic; i++) {
                const OVPicture *pic = &dpb->pictures[i];
                if (pic->frame && pic->frame->data[0] && pic->cvs_id == cvs_id_min1) {
                    min_cvs_id = pic->cvs_id;
                    got_pic = 1;
                    cvs_id_min1--;
                    break;
                }
                got_pic = 0;
            }
        } while (got_pic);
    }

    return min_cvs_id;
}

/* Remove reference flags on all picture */
static void
ovdpb_clear_refs(OVDPB *dpb)
{
    //TODO: loop untill min_idx == nb_dpb_pic or nb_used_pic > dpb->max_nb_dpb_pic
    ov_log(NULL, OVLOG_DEBUG, "Release reference pictures\n");
    int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    int nb_used_pic = 0;
    int min_cvs_id = find_min_cvs_id(dpb);
    int i;

    /* Count pictures in current output target Coded Video Sequence
     */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];
        if (pic->frame && pic->frame->data[0]) {
            nb_used_pic++;
        }
    }

    while (nb_used_pic >= dpb->max_nb_dpb_pic) {
        int min_poc = INT_MAX;
        int min_idx = nb_dpb_pic;
        /* Determine the min POC among those pic
         */
        for (i = 0; i < nb_dpb_pic; i++) {
            OVPicture *pic = &dpb->pictures[i];
            uint16_t ref_count = atomic_load(&pic->ref_count);
            // is_output_cvs = pic->cvs_id == output_cvs_id;
            // if (is_output_cvs && pic->frame && pic->frame->data[0] && !ref_count) {
            if (pic->frame && pic->frame->data[0] && !ref_count && !pic->flags) {
                if (pic->poc < min_poc && pic->cvs_id == min_cvs_id) {
                    min_poc = pic->poc;
                    min_idx = i;
                }
            }
        }

        //TODOdpb: loop untill min_idx == nb_dpb_pic or nb_used_pic > dpb->max_nb_dpb_pic
        /* Try to release picture with POC == to min_poc
         */
        if (min_idx < nb_dpb_pic) {
            OVPicture *pic = &dpb->pictures[min_idx];
            ovdpb_release_pic(dpb, pic);
        } else {
            break;
        }

        nb_used_pic--;
    }
}


/* All pictures are removed from the DPB */
void
ovdpb_flush_dpb(OVDPB *dpb)
{
    int i;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    // const uint8_t flags = ~0;

    for (i = 0; i < nb_dpb_pic; i++) {
        dpb->pictures[i].flags = 0; 
        atomic_init( &dpb->pictures[i].ref_count, 0);
        ovdpb_release_pic(dpb, &dpb->pictures[i]);
    }
}

/*FIXME rename to request new picture */
static OVPicture *
alloc_frame(OVDPB *dpb, int poc)
{
    int i, ret;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];

        /* FIXME we could avoid checking for data in picture
         * if we are sure every unreference frame is set to NULL
         */
        if (pic->frame && pic->frame->data[0]) {
            continue;
        }

        ret = dpbpriv_request_frame(&dpb->internal, &pic->frame);
        
        if (ret < 0) {
            ov_log(NULL, OVLOG_ERROR, "Error while requesting picture from DPB\n");
            return NULL;
        }

        pic->flags = 0;
        atomic_init(&pic->ref_count, 0);

        ov_log(NULL, OVLOG_DEBUG, "Attached frame %p to Picture with POC %d\n", pic->frame, poc);

        return pic;
    }

    ov_log(NULL, OVLOG_ERROR, "DPB full\n");

    return NULL;
}

/* Allocate the current picture buffer */
static int
ovdpb_init_current_pic(OVDPB *dpb, OVPicture **pic_p, int poc, uint8_t ph_pic_output_flag)
{
    OVPicture *pic;
    int i;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);

    *pic_p = NULL;

    /* check that this POC doesn't already exist */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];

        if (pic->frame && pic->frame->data[0] && pic->cvs_id == dpb->cvs_id &&
            pic->poc == poc) {
            ov_log(NULL, OVLOG_ERROR, "Duplicate POC in a sequence: %d for cvs_id: %d.\n",
                   poc, pic->cvs_id);
            //*pic_p = pic;
            *pic_p = NULL;
            //ovdpb_report_decoded_frame(ref_pic);
            return OVVC_EINDATA;
        }
    }

    pic = alloc_frame(dpb, poc);

    if (!pic) {
        return OVVC_ENOMEM;
    }

    *pic_p = pic;

    ovdpb_reset_decoded_ctus(pic);

    if (ph_pic_output_flag) {
        ovdpb_new_ref_pic(pic, OV_OUTPUT_PIC_FLAG);
        ovdpb_new_ref_pic(pic, OV_IN_DECODING_PIC_FLAG);
    } else {
        ovdpb_new_ref_pic(pic, OV_IN_DECODING_PIC_FLAG);
    }

    pic->poc    = poc;
    pic->cvs_id = dpb->cvs_id;
    pic->frame->poc = poc;

    /* Copy display or conformance window properties */

    return 0;
}

static int
compute_ref_poc(const OVRPL *const rpl, struct RPLInfo *const rpl_info, int32_t poc, uint8_t weighted_pred)
{
    const int nb_refs = rpl->num_ref_entries;
    int i;
    rpl_info->nb_refs = nb_refs;
    rpl_info->nb_active_refs = rpl->num_ref_active_entries;
    int last_poc = poc;

    for (i = 0; i < nb_refs; ++i) {
        const struct RefPic *const rp = &rpl->rp_list[i];
        struct RefInfo *const rinfo = &rpl_info->ref_info[i];
        enum RefType ref_type = rp->st_ref_pic_flag ? ST_REF
                                                    : (rp->inter_layer_ref_pic_flag ? ILRP_REF
                                                                                    : LT_REF);
        int ref_poc = 0;
        rinfo->type = ref_type;


        switch (ref_type) {
        case ST_REF:
           ref_poc = !rp->strp_entry_sign_flag ? last_poc +  rp->abs_delta_poc_st + (!weighted_pred | !i)
                                               : last_poc - (rp->abs_delta_poc_st + (!weighted_pred | !i));
           rinfo->poc = ref_poc;

           last_poc = ref_poc;
        break;
        case LT_REF:
           /* FIXME
            *    - Handle when read in header (ltrp_in_header_flag)
            *    - Compute msb part of POC when needed
            */
           ref_poc = rp->rpls_poc_lsb_lt;

           rinfo->poc = ref_poc;
           ov_log(NULL, OVLOG_WARNING, "Partially supported Long Term Ref \n");

        break;
        case ILRP_REF:
           ov_log(NULL, OVLOG_ERROR, "Unsupported Inter Layer Ref \n");
           rinfo->poc = ref_poc;
        break;
        }

        if (poc == ref_poc) {
            goto self_ref;
        }
    }

    return 0;

self_ref:
    ov_log(NULL, OVLOG_ERROR, "Invalid self reference for picture with POC %d.\n", poc);
    return OVVC_EINDATA;
}


static int
vvc_mark_refs(OVDPB *dpb, const OVRPL *rpl, int32_t poc, OVPicture **dst_rpl, uint8_t weighted_pred)
{
    int i, j;
    uint8_t found, flag;
    int16_t ref_poc, ref_type;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);

    struct RPLInfo rpl_info;
    int ret = compute_ref_poc(rpl, &rpl_info, poc, weighted_pred);
    if (ret < 0) {
        return ret;
    }

    for (i = 0;  i < rpl->num_ref_active_entries; ++i){
        ref_poc  = rpl_info.ref_info[i].poc;
        ref_type = rpl_info.ref_info[i].type;
        flag = ref_type == ST_REF ? OV_ST_REF_PIC_FLAG : OV_LT_REF_PIC_FLAG;
        OVPicture *ref_pic;
        found = 0;
        for (j = 0; j < nb_dpb_pic; j++) {
            ref_pic = &dpb->pictures[j];
            if (ref_pic->poc == ref_poc && ref_pic->cvs_id == dpb->cvs_id){
                if(ref_pic->frame && ref_pic->frame->data[0]){
                    found = 1;
                    ov_log(NULL, OVLOG_TRACE, "Mark active reference %d for picture %d\n", ref_poc, dpb->poc);
                    ref_pic->flags &= ~(OV_LT_REF_PIC_FLAG | OV_ST_REF_PIC_FLAG);
                    ovdpb_new_ref_pic(ref_pic, flag);
                    dst_rpl[i] = ref_pic; 
                    break;
                }
            }
        }

        if (!found){
            /* If reference picture is not in the DPB we try create a new
             * Picture with requested POC ID in the DPB
             */
            ov_log(NULL, OVLOG_ERROR, "Generating missing reference %d for picture %d\n", ref_poc, dpb->poc);
            ref_pic = alloc_frame(dpb, ref_poc);

            if (ref_pic == NULL){
                return OVVC_ENOMEM;
            }

            ovdpb_report_decoded_frame(ref_pic);

            ref_pic->poc    = ref_poc;
            ref_pic->cvs_id = dpb->cvs_id;

            ref_pic->flags  = 0;

            ovdpb_new_ref_pic(ref_pic, OV_ST_REF_PIC_FLAG);

            /*FIXME  Set output / corrupt flag ? */
            dst_rpl[i] = ref_pic; 
        }
    }

    /* Mark non active refrences pictures as used for reference */
    for (; i < rpl->num_ref_entries; ++i) {
        ref_poc  = rpl_info.ref_info[i].poc;
        ref_type = rpl_info.ref_info[i].type;
        flag = ref_type == ST_REF ? OV_ST_REF_PIC_FLAG : OV_LT_REF_PIC_FLAG;
        OVPicture *ref_pic;
        found = 0;

        for (j = 0; j < nb_dpb_pic; j++) {
            ref_pic = &dpb->pictures[j];
            if (ref_pic->poc == ref_poc){
                if(ref_pic->frame && ref_pic->frame->data[0] && ref_pic->cvs_id == dpb->cvs_id){
                    found = 1;
                    ov_log(NULL, OVLOG_TRACE, "Mark non active reference %d for picture %d\n", ref_poc, dpb->poc);
                    ref_pic->flags &= ~(OV_LT_REF_PIC_FLAG | OV_ST_REF_PIC_FLAG);
                    ovdpb_new_ref_pic(ref_pic, flag);
                    dst_rpl[i] = ref_pic;
                    break;
                }
            }
        }

        if (!found){
            ov_log(NULL, OVLOG_TRACE, "Not found non active reference %d for picture %d\n", ref_poc, dpb->poc);
        }
    }

    return 0;
}

static void
vvc_unmark_refs(OVPicture * current_pic, OVPicture **dst_rpl, uint8_t nb_active_refs, uint8_t nb_refs)
{
    int i;
    for (i = 0;  i < nb_active_refs; ++i) {
        OVPicture *ref_pic = dst_rpl[i];
        if (ref_pic) {
            int16_t ref_poc  = ref_pic->poc;
            int16_t ref_type = ST_REF;
            uint8_t flag = ref_type == ST_REF ? OV_ST_REF_PIC_FLAG : OV_LT_REF_PIC_FLAG;
            ov_log(NULL, OVLOG_TRACE, "Unmark active reference %d from picture %d RPL\n", ref_poc, current_pic->poc);
            ovdpb_unref_pic(ref_pic, flag);
        }
    }

    for (;  i < nb_refs; ++i) {
        OVPicture *ref_pic = dst_rpl[i];
        if(ref_pic){
            int16_t ref_poc  = ref_pic->poc;
            int16_t ref_type = ST_REF;
            uint8_t flag = ref_type == ST_REF ? OV_ST_REF_PIC_FLAG : OV_LT_REF_PIC_FLAG;
            ov_log(NULL, OVLOG_TRACE, "Unmark non active reference %d from picture %d RPL\n", ref_poc, current_pic->poc);
            ovdpb_unref_pic(ref_pic, flag);
            dst_rpl[i] = NULL;
        }
     }
}

static void
dpb_pic_to_frame_ref(OVPicture *pic, OVFrame **dst, struct OVSEI **sei_p)
{
    if (pic->flags & OV_OUTPUT_PIC_FLAG) {
        ovframe_new_ref(dst, pic->frame);
    } else {
       *dst = NULL;
    }

    /* Move sei from pic to output
     * FIXME use ref/unref instead
     */
    *sei_p = pic->sei;
    pic->sei = NULL;

    ovdpb_unref_pic(pic, OV_OUTPUT_PIC_FLAG | (pic->flags & OV_BUMPED_PIC_FLAG));
}

int
ovdpb_drain_frame(OVDPB *dpb, OVFrame **out, OVSEI **sei_p)
{
    int output_cvs_id = find_min_cvs_id(dpb);

    do {
        const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
        int nb_output = 0;
        int min_poc   = INT_MAX;
        int min_idx   = INT_MAX;
        int i;

        /* Count pic marked for output in output cvs and find the min poc_id */
        for (i = 0; i < nb_dpb_pic; i++) {
            OVPicture *pic = &dpb->pictures[i];
            uint8_t output_flag = (pic->flags & OV_OUTPUT_PIC_FLAG);
            uint8_t is_output_cvs = pic->cvs_id == output_cvs_id;
            if (output_flag && is_output_cvs) {
                nb_output++;
                if (pic->poc < min_poc || nb_output == 1) {
                    min_poc = pic->poc;
                    min_idx = i;
                }
            }
        }

        /* If the number of pic to output is less than max_num_reorder_pics
         * in current cvs we wait for more pic before outputting any
         */
        if (nb_output) {
            OVPicture *pic = &dpb->pictures[min_idx];
            dpb->pts += dpb->nb_units_in_ticks;
            pic->frame->pts = dpb->pts;
            dpb_pic_to_frame_ref(pic, out, sei_p);
            return nb_output;
        }

        /* If no output pic found increase cvs_id and retry */
        if (output_cvs_id != dpb->cvs_id) {
            output_cvs_id = (output_cvs_id + 1) & 0xff;
        } else {
            break;
        }

    } while (1);

    ov_log(NULL, OVLOG_TRACE, "No picture to output\n");

    *out = NULL;

    return 0;
}

int
ovdpb_output_pic(OVDPB *dpb, OVFrame **out, OVSEI **sei_p)
{
    int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    int i;
    int output_cvs_id = find_min_cvs_id(dpb);

    int nb_output = 0;
    int min_poc   = INT_MAX;
    int min_idx   = nb_dpb_pic;
    uint8_t in_decoding;

    /* Count pic marked for output in output cvs and find the min poc_id */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];
        uint8_t output_flag = (pic->flags & OV_OUTPUT_PIC_FLAG);
        uint8_t is_output_cvs = pic->cvs_id == output_cvs_id;
        if (output_flag && is_output_cvs) {
            in_decoding = (pic->flags & OV_IN_DECODING_PIC_FLAG);
            if(!in_decoding){
                nb_output ++;
            }
            if (pic->poc < min_poc ){
                min_poc = pic->poc;
                if(!in_decoding){
                    min_idx = i;
                }
                else{
                    min_idx = nb_dpb_pic;
                }
            }
        }
    }

    /* If the number of pic to output is less than max_num_reorder_pics
     * in current cvs we wait for more pic before outputting any
     */
    if (output_cvs_id == dpb->cvs_id && nb_output <= dpb->max_nb_reorder_pic) {
        return 0;
    }

    if (min_idx < nb_dpb_pic) {
        OVPicture *pic = &dpb->pictures[min_idx];
        dpb->pts += dpb->nb_units_in_ticks;
        pic->frame->pts = dpb->pts;
        dpb_pic_to_frame_ref(pic, out, sei_p);
        return nb_output;
    }

    *out = NULL;

    ov_log(NULL, OVLOG_TRACE, "No picture to output\n");

    return 0;
}

void
ovdpb_unmark_ref_pic_lists(uint8_t slice_type, OVSliceDec *sldec)
{
    vvc_unmark_refs(sldec->pic, sldec->rpl0, sldec->nb_active_refs0, sldec->nb_refs0);

    if (slice_type == SLICE_B){
        vvc_unmark_refs(sldec->pic, sldec->rpl1, sldec->nb_active_refs1, sldec->nb_refs1);
    }
}

static int
mark_ref_pic_lists(OVDPB *const dpb, uint8_t slice_type, const struct OVRPL *const rpl0,
                   const struct OVRPL *const rpl1, OVSliceDec *const sldec)
{
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    uint32_t poc = dpb->poc;
    int i, ret;
    OVPicture *current_pic = NULL;

    /* This is the same as clear_refs except we do not remove
     * flags on current picture
     */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];

        if (pic->cvs_id == dpb->cvs_id && pic->poc == dpb->poc) {
            current_pic = pic;
            continue;
        }

        pic->flags &= ~(OV_LT_REF_PIC_FLAG | OV_ST_REF_PIC_FLAG);
    }
    uint8_t weighted_pred = sldec->active_params.sps->sps_weighted_pred_flag || sldec->active_params.sps->sps_weighted_bipred_flag;

    ret = vvc_mark_refs(dpb, rpl0, poc, sldec->rpl0, weighted_pred);

    sldec->nb_refs0 = rpl0->num_ref_entries;
    sldec->nb_active_refs0 = rpl0->num_ref_active_entries;

    if (slice_type == SLICE_B){
        ret |= vvc_mark_refs(dpb, rpl1, poc, sldec->rpl1, weighted_pred);
        sldec->nb_refs1 = rpl1->num_ref_entries;
        sldec->nb_active_refs1 = rpl1->num_ref_active_entries;
    } else {
        sldec->nb_active_refs1 = 0;
    }

    if ((slice_type != SLICE_I && !sldec->nb_active_refs0) || (!sldec->nb_active_refs1 && slice_type == SLICE_B)) {
         ret = OVVC_EINDATA;
    }

    if (ret < 0) {
        goto fail;
    }

    return 0;

fail:
    ovdpb_unmark_ref_pic_lists(slice_type, sldec);
    return ret;
}


static void
tmvp_release_mv_planes(OVPicture *const pic)
{
    /* FIXME check mv planes exists */
    if (pic->mv_plane0.dirs) {
        mvpool_release_mv_plane(&pic->mv_plane0);
    }

    if (pic->mv_plane1.dirs) {
        mvpool_release_mv_plane(&pic->mv_plane1);
    }
}

static int
tmvp_request_mv_plane(OVPicture *const pic, const OVVCDec *ovdec, uint8_t slice_type)
{
    struct MVPool *pool = ovdec->mv_pool;
    int ret;

    ret = mvpool_request_mv_plane(pool, &pic->mv_plane0);
    if (ret < 0) {
        return ret;
    }

    if (slice_type == SLICE_B) {
        ret = mvpool_request_mv_plane(pool, &pic->mv_plane1);
        if (ret < 0) {
            mvpool_release_mv_plane(&pic->mv_plane0);
            return ret;
        }
    }

    return 0;
}

static int
init_tmvp_info(OVPicture *const pic, const OVPS *const ps, const OVVCDec *ovdec)
{
    const OVSH *sh = ps->sh;
    const OVPH *ph = ps->ph;

    uint8_t slice_type = sh->sh_slice_type;

    /* The picture can contain inter slice thus Motions Vector */
    if (ph->ph_inter_slice_allowed_flag) {

        /* Request MV buffer to MV Pool */
        tmvp_request_mv_plane(pic, ovdec, slice_type);
    }

    return 0;
}

static int
update_rpl(const OVPPS *const pps,
           const OVSH  *const sh,
           const OVPH  *const ph,
           OVRPL *const rpl0,
           OVRPL *const rpl1,
           uint8_t slice_type)
{
    const OVHRPL *const hrpl = !pps->pps_rpl_info_in_ph_flag ? &sh->hrpl: &ph->hrpl;

    if (slice_type != 2) {
        /* copy rpl content */
        *rpl0 = hrpl->rpl_h0.rpl_data;
        *rpl1 = hrpl->rpl_h1.rpl_data;
        if (sh->sh_num_ref_idx_active_override_flag) {
            uint8_t nb_rpl_ref0 = hrpl->rpl_h0.rpl_data.num_ref_entries;

            if (nb_rpl_ref0 > 1) {
                uint8_t nb_active_ref0 = sh->sh_num_ref_idx_active_l0_minus1 + 1;
                rpl0->num_ref_active_entries = nb_active_ref0;
            } else {
                rpl0->num_ref_active_entries = nb_rpl_ref0;
            }

            if (slice_type == 0) {
                uint8_t nb_rpl_ref1 = hrpl->rpl_h1.rpl_data.num_ref_entries;

                if (nb_rpl_ref1 > 1) {
                    uint8_t nb_active_ref1 = sh->sh_num_ref_idx_active_l1_minus1 + 1;
                    rpl1->num_ref_active_entries = nb_active_ref1;
                } else {
                    rpl1->num_ref_active_entries = nb_rpl_ref1;
                }
            } else
                rpl1->num_ref_active_entries = 0;
        } else {
            uint8_t nb_active_ref0 = pps->pps_num_ref_idx_default_active_minus1[0] + 1;
            uint8_t nb_rpl_ref0    = hrpl->rpl_h0.rpl_data.num_ref_entries;

            rpl0->num_ref_active_entries = nb_active_ref0 < nb_rpl_ref0 ? nb_active_ref0 : nb_rpl_ref0;

            if (slice_type == 0) {
                uint8_t nb_active_ref1 = pps->pps_num_ref_idx_default_active_minus1[1] + 1;
                uint8_t nb_rpl_ref1    = hrpl->rpl_h1.rpl_data.num_ref_entries;

                rpl1->num_ref_active_entries = nb_active_ref1 < nb_rpl_ref1 ? nb_active_ref1 : nb_rpl_ref1;
            } else
                rpl1->num_ref_active_entries = 0;
        }
    } else {
        /* I slice might contain RPL info for non active refs
         * in this case all refs are non active so we still
         * copy rpl info.
         */
        *rpl0 = hrpl->rpl_h0.rpl_data;
        *rpl1 = hrpl->rpl_h1.rpl_data;
        rpl0->num_ref_active_entries = 0;
        rpl1->num_ref_active_entries = 0;
    }
    return 0;
}

int
ovdpb_init_picture(OVDPB *dpb, OVPicture **pic_p, const OVPS *const ps, uint8_t nalu_type,
                   OVSliceDec *const sldec, const OVVCDec *ovdec)
{

    const OVSH  *const sh  = ps->sh;
    const OVPH  *const ph  = ps->ph;
    int ret = 0;
    uint32_t poc = dpb->poc;
    uint8_t cra_flag = 0;
    uint8_t idr_flag = 0;

    idr_flag |= nalu_type == OVNALU_IDR_W_RADL;
    idr_flag |= nalu_type == OVNALU_IDR_N_LP;

    cra_flag |= nalu_type == OVNALU_CRA;
    cra_flag |= nalu_type == OVNALU_GDR;

    /* TODO move to dec init */
    if (idr_flag){
        /* New IDR involves a POC refresh and mark the start of
         * a new coded video sequence
         */
        dpb->cvs_id = (dpb->cvs_id + 1) & 0xFF;
        if (ps->ph->ph_poc_msb_cycle_present_flag) {
            uint8_t log2_max_poc_lsb = ps->sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4;
            poc = ps->ph->ph_poc_msb_cycle_val << log2_max_poc_lsb;
        } else {
            poc = 0;
        }
        poc += ps->ph->ph_pic_order_cnt_lsb;
    } else {
        uint32_t last_poc = dpb->poc;
        poc = derive_poc(ps->ph->ph_pic_order_cnt_lsb,
                         ps->sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4,
                         last_poc);
    }

    dpb->poc = poc;

    /* Find an available place in DPB and allocate/retrieve available memory
     * for the current picture data from the Frame Pool
     */
    ret = ovdpb_init_current_pic(dpb, pic_p, poc, ps->ph->ph_pic_output_flag);
    if (ret < 0) {
        goto fail;
    }
    //TODOrpr: put width and height as parameters of ovdpb_init_current_pic
    (*pic_p)->frame->width  = ps->pps->pps_pic_width_in_luma_samples;
    (*pic_p)->frame->height = ps->pps->pps_pic_height_in_luma_samples;

    if (ps->pps->pps_conformance_window_flag) {
        (*pic_p)->frame->output_window.offset_lft = ps->pps->pps_conf_win_left_offset;
        (*pic_p)->frame->output_window.offset_rgt = ps->pps->pps_conf_win_right_offset;
        (*pic_p)->frame->output_window.offset_abv = ps->pps->pps_conf_win_top_offset;
        (*pic_p)->frame->output_window.offset_blw = ps->pps->pps_conf_win_bottom_offset;
    } else {
        (*pic_p)->frame->output_window.offset_lft = ps->sps->sps_conf_win_left_offset;
        (*pic_p)->frame->output_window.offset_rgt = ps->sps->sps_conf_win_right_offset;
        (*pic_p)->frame->output_window.offset_abv = ps->sps->sps_conf_win_top_offset;
        (*pic_p)->frame->output_window.offset_blw = ps->sps->sps_conf_win_bottom_offset;
    }

    (*pic_p)->scale_info.scaling_win_left   = ps->pps->pps_scaling_win_left_offset;
    (*pic_p)->scale_info.scaling_win_right  = ps->pps->pps_scaling_win_right_offset;
    (*pic_p)->scale_info.scaling_win_top    = ps->pps->pps_scaling_win_top_offset;
    (*pic_p)->scale_info.scaling_win_bottom = ps->pps->pps_scaling_win_bottom_offset;
    (*pic_p)->scale_info.chroma_hor_col_flag = ps->sps->sps_chroma_horizontal_collocated_flag;
    (*pic_p)->scale_info.chroma_ver_col_flag = ps->sps->sps_chroma_vertical_collocated_flag;

    copy_sei_params(&(*pic_p)->sei, ovdec->active_params.sei);

    (*pic_p)->sei->upscale_flag = ovdec->upscale_flag;
    (*pic_p)->sei->scaling_info = (*pic_p)->scale_info;

    ov_log(NULL, OVLOG_TRACE, "DPB start new picture POC: %d\n", (*pic_p)->poc);

    /* If the picture is not an IDR Picture we set all flags to
     * FIXME in VVC we might still get some ref pic list in IDR
     * Slices it is not clear whether we should still mark them
     * or not
     */
    if (!idr_flag || ovdec->active_params.sps->sps_idr_rpl_present_flag) {
        const OVPPS *const pps = ps->pps;
        uint8_t slice_type = sh->sh_slice_type;
        uint8_t weighted_pred = ovdec->active_params.sps->sps_weighted_pred_flag || ovdec->active_params.sps->sps_weighted_bipred_flag;
        OVRPL rpl0, rpl1;
        update_rpl(pps, sh, ph, &rpl0, &rpl1, slice_type);
        if ((rpl0.num_ref_entries | rpl1.num_ref_entries) & ~0xF) {
            ov_log(NULL, OVLOG_ERROR, "Too many pictures in RPL for picture POC: %d\n", (*pic_p)->poc);
            ret = OVVC_EINDATA;
            goto failnoclear;
        }

        if ((rpl0.num_ref_active_entries > rpl0.num_ref_entries) ||
            (rpl1.num_ref_active_entries > rpl1.num_ref_entries)) {
            ret = OVVC_EINDATA;
            ov_log(NULL, OVLOG_ERROR, "Too many active pictures in RPL for picture POC: %d\n", (*pic_p)->poc);
            goto failnoclear;
        }

        ret = mark_ref_pic_lists(dpb, slice_type, &rpl0, &rpl1, sldec);
        if (ret < 0) {
            goto fail;
        }

        if (!idr_flag) {
            int last_dist = 0;
            for (int i = 0; i < rpl0.num_ref_entries; ++i) {
                const struct RefPic *const rp = &rpl0.rp_list[i];

                if (rp->st_ref_pic_flag) {
                    last_dist += rp->strp_entry_sign_flag ?  (rp->abs_delta_poc_st + (!weighted_pred | !i))
                                                          : -(rp->abs_delta_poc_st + (!weighted_pred | !i));

                    sldec->dist_ref_0[i] = ov_clip_intp2(last_dist, 8);
                } else {
                    sldec->dist_ref_0[i] = 0;
                }
            }

            last_dist = 0;
            for (int i = 0; i < rpl1.num_ref_entries; ++i) {
                const struct RefPic *const rp = &rpl1.rp_list[i];
                if (rp->st_ref_pic_flag) {
                    last_dist += rp->strp_entry_sign_flag ?  (rp->abs_delta_poc_st + (!weighted_pred | !i))
                                                          : -(rp->abs_delta_poc_st + (!weighted_pred | !i));

                    sldec->dist_ref_1[i] = ov_clip_intp2(last_dist, 8);
                } else {
                    sldec->dist_ref_1[i] = 0;
                }
            }
        }
    }

    ovdpb_clear_refs(dpb);

    /* Init picture TMVP info */
    if (ps->sps->sps_temporal_mvp_enabled_flag) {
        ret = init_tmvp_info(*pic_p, ps, ovdec);
    }

    return ret;

fail:
    ovdpb_clear_refs(dpb);

failnoclear:
    return ret;
}

static void
xctu_to_mask(uint64_t* mask, int mask_w, int xmin_ctu, int xmax_ctu)
{
    int sub_xmin_ctu, sub_xmax_ctu;
    for(int i = (xmin_ctu >> SIZE_INT64); i <= (xmax_ctu >> SIZE_INT64); i++) {
        mask[i] = 0;
        sub_xmin_ctu = xmin_ctu > (i<<SIZE_INT64)     ? xmin_ctu % (1<<SIZE_INT64) : 0;
        sub_xmax_ctu = xmax_ctu < ((i+1)<<SIZE_INT64) ? xmax_ctu % (1<<SIZE_INT64) : (1<<SIZE_INT64)-1;
        for(int ii = sub_xmin_ctu; ii <= sub_xmax_ctu; ii++) {
            mask[i] |= (uint64_t)1 << ii;
        }
    }
}

static void
ovdpb_no_synchro(const OVPicture *const ref_pic, int tl_ctu_x, int tl_ctu_y, int br_ctu_x, int br_ctu_y)
{
    return;
}

static void
ovdpb_synchro_ref_decoded_ctus(const OVPicture *const ref_pic, int tl_ctu_x, int tl_ctu_y, int br_ctu_x, int br_ctu_y)
{
    const struct PictureSynchro* sync = &ref_pic->sync;

    int mask_w = sync->map_w;
    uint64_t wanted_mask[5];
    xctu_to_mask(wanted_mask, mask_w, tl_ctu_x, br_ctu_x);

    uint8_t all_ctus_available;
    do {
        pthread_mutex_lock(sync->ref_mtx);

        all_ctus_available = 1;
        for (int ctu_y = tl_ctu_y; ctu_y <= br_ctu_y; ctu_y++ ) {
            const uint64_t *ctu_col = sync->decoded_ctus_map[ctu_y];
            for (int i = 0; i < mask_w; i++)
                all_ctus_available = all_ctus_available && ((ctu_col[i] & wanted_mask[i]) == wanted_mask[i]);
        }
        if (!all_ctus_available) {
            pthread_cond_wait(sync->ref_cnd, sync->ref_mtx);
        }
        pthread_mutex_unlock(sync->ref_mtx);

    } while (!all_ctus_available);
}

void
ovdpb_init_decoded_ctus(OVPicture *const pic, const OVPS *const ps)
{   
    int pic_w = ps->sps->sps_pic_width_max_in_luma_samples;
    int pic_h = ps->sps->sps_pic_height_max_in_luma_samples;
    uint8_t log2_ctb_s    = (ps->sps->sps_log2_ctu_size_minus5 + 5) & 0x7;
    uint16_t nb_ctb_pic_w = (pic_w + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;
    uint16_t nb_ctb_pic_h = (pic_h + ((1 << log2_ctb_s) - 1)) >> log2_ctb_s;

    struct PictureSynchro* sync = &pic->sync;
    sync->nb_ctu_h = nb_ctb_pic_h;
    sync->map_w = (nb_ctb_pic_w >> SIZE_INT64) + 1;

    sync->func = &sync->internal.sync_function;
    sync->ref_mtx = &sync->internal.ref_mtx;
    sync->ref_cnd = &sync->internal.ref_cnd;

    atomic_init(sync->func, (uintptr_t)&ovdpb_synchro_ref_decoded_ctus);
}

void
ovdpb_report_decoded_ctu_line(OVPicture *const pic, int y_ctu, int xmin_ctu, int xmax_ctu)
{
    struct PictureSynchro* sync = &pic->sync;
    int mask_w = sync->map_w;
    uint64_t mask[5];

    xctu_to_mask(mask, mask_w, xmin_ctu, xmax_ctu);

    pthread_mutex_lock(sync->ref_mtx);

    for (int i = 0; i < mask_w; i++) {
        sync->decoded_ctus_map[y_ctu][i] |= mask[i];
    }

    pthread_cond_broadcast(sync->ref_cnd);

    pthread_mutex_unlock(sync->ref_mtx);
}

void
ovdpb_report_decoded_frame(OVPicture *const pic)
{
    struct PictureSynchro* sync = &pic->sync;

    atomic_store(sync->func, (uintptr_t) &ovdpb_no_synchro);

    pthread_mutex_lock(sync->ref_mtx);
    for(int i = 0; i < sync->nb_ctu_h; i++){
        memset(sync->decoded_ctus_map[i], 0xFF, sync->map_w * sizeof(int64_t));
    }
    pthread_cond_broadcast(sync->ref_cnd);
    pthread_mutex_unlock(sync->ref_mtx);
}

static void
ovdpb_reset_decoded_ctus(OVPicture *const pic)
{
    struct PictureSynchro* sync = &pic->sync;
    pthread_mutex_lock(sync->ref_mtx);
    for(int i = 0; i < sync->nb_ctu_h; i++){
        memset(sync->decoded_ctus_map[i], 0, sync->map_w * sizeof(int64_t));
    }
    pthread_mutex_unlock(sync->ref_mtx);

    atomic_store(sync->func, (uintptr_t)&ovdpb_synchro_ref_decoded_ctus);
}

