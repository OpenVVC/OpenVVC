#include <string.h>
#include <stddef.h>
#include <limits.h>
#include "libovvcutils/ovvcutils.h"
#include "libovvcutils/ovmem.h"
#include "libovvcutils/ovvcerror.h"
#include "nvcl_private.h"
#include "nvcl_structures.h"
#include "ovframe.h"
#include "libovvcdmx/ovunits.h"
#if 0
#include "ovdpb.h"
#endif
#include "ovdpb_internal.h"


#if 1
#endif

#define OV_OUTPUT_PIC_FLAG (1 << 0)
#define OV_LT_REF_PIC_FLAG (1 << 1)
#define OV_ST_REF_PIC_FLAG (1 << 2)
#define OV_BUMPED_PIC_FLAG (1 << 3)

enum SliceType {
   SLICE_B = 0,
   SLICE_P = 1,
   SLICE_I = 2,
};

enum RefType {
    /* Short term reference Picture */
    ST_REF = 1,

    /* Long term reference Picture */
    LT_REF = 2,

    /* Inter Layer reference Picture */
    ILRP_REF = 3
};

struct RefInfo {
    enum RefType type;
    uint32_t poc;
};

struct RPLInfo {
   struct RefInfo ref_info[16];
   uint8_t nb_refs;
};

struct OVPicture{

   /* Associated frame */
   OVFrame *frame; 

   /* Flags used to mark Picture referenced by the 
    * active picture (current picture being decoded)
    * FIXME enum ?
    */
   uint8_t flags;

   /* Pointers to ref_pic_list */
   struct OVPicture *rpl1[16];

   struct TMVPInfo {
       OVPicture *collocated_ref;
       /* TODO tmvp scaling */
   }tmvp;

   uint32_t poc;

   /* Coded Video Sequence Id to which this Picture is
    * Associated : this avoid confusing ref with same POC
    * when the refresh period is shorter than DPB
    */
   uint16_t cvs_id;
};

struct DPB
{
   OVPicture pictures[64];
   uint8_t max_nb_dpb_pic;
   uint8_t max_nb_reorder_pic;
   uint8_t max_latency_increase;

   /* Coded Video Sequence Id 
    * It enables to reset the POC of pictures
    * when encountering IDR Picture Units
    */
   uint16_t cvs_id;

   /* DPB status info 
    * to be used to determine whether the decoder is waiting
    * for an IRAP picture or not
    */
   uint8_t state;

   struct DPBInternal internal;
};


/* Clean OVPicture from DPB and unref associated Frame */
static void
dpbpriv_release_pic(OVPicture *pic)
{
    if (pic->frame) {
        /* FIXME unref frame */
        ovframe_unref(&pic->frame);
        /* Do not delete frame the frame will delete itself 
         * when all its references are released
         */
        #if 0
        ov_freep(&pic->frame);
        #endif
    }
}

static int
dpb_init_params(OVDPB *dpb, OVDPBParams const *prm)
{
    dpb->max_nb_dpb_pic       = prm->dpb_max_dec_pic_buffering_minus1 + 1; 
    dpb->max_nb_reorder_pic   = prm->dpb_max_num_reorder_pics; 
    dpb->max_latency_increase = prm->dpb_max_latency_increase_plus1 - 1; 
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


void ovdpb_unref_frame(OVDPB *dpb, OVPicture *pic, int flags)
{
    /* pic->frame can be NULL if context init failed */
    if (!pic->frame || !pic->frame->data[0])
        return;

    pic->flags &= ~flags;

    /* If there is no more flags the picture can be
     * returned to the DPB;
     */
    if (!pic->flags) {
        dpbpriv_release_pic(pic);
    }
}

int ovdpb_ref_pic(OVDPB *dpb, OVPicture *dst, OVPicture *src)
{
    int ret;

    ret = ovframe_new_ref(&dst->frame, src->frame);
    if (ret < 0) {
        return ret;
    }

    /* TMVP */

    dst->poc     = src->poc;
    dst->flags   = src->flags;
    dst->cvs_id  = src->cvs_id;

    return 0;
}

static void 
vvc_clear_refs(OVDPB *dpb)
{
    int i;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    const uint8_t flags = OV_ST_REF_PIC_FLAG | OV_LT_REF_PIC_FLAG;

    for (i = 0; i < nb_dpb_pic; i++) {
        ovdpb_unref_frame(dpb, &dpb->pictures[i], flags);
    }
}

void ovdpb_flush_dpb(OVDPB *dpb)
{
    int i;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    const uint8_t flags = ~0;

    for (i = 0; i < nb_dpb_pic; i++) {
        ovdpb_unref_frame(dpb, &dpb->pictures[i], flags);
    }
}

static OVPicture *alloc_frame(OVDPB *dpb)
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
            return NULL;
        }

        return pic;
    }

    ov_log(NULL, 3, "DPB full\n");
    return NULL;
}

/* Allocate the current picture buffer */
int ovdpb_init_current_pic(OVDPB *dpb, OVFrame **frame, int poc)
{
    OVPicture *ref_pic;
    int i;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);

    /* check that this POC doesn't already exist */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];

        if (pic->frame->data[0] && pic->cvs_id == dpb->cvs_id &&
            pic->poc == poc) {
            ov_log(NULL, 3, "Duplicate POC in a sequence: %d.\n",
                   poc);
            return OVVC_EINDATA;
        }
    }

    ref_pic = alloc_frame(dpb);
    if (!ref_pic) {
        return OVVC_ENOMEM;
    }

    *frame = ref_pic->frame;

    #if 0
    dpb->active_pic = ref_pic;
    #endif

    #if 0
    if (dpb->ps.ph_data->ph_pic_output_flag) {
    #endif
        ref_pic->flags = OV_OUTPUT_PIC_FLAG | OV_ST_REF_PIC_FLAG;
    #if 0
    } else {
        ref_pic->flags = OV_ST_REF_PIC_FLAG;
    }
    #endif

    ref_pic->poc    = poc;
    ref_pic->cvs_id = dpb->cvs_id;

    /* Copy display or conformance window properties */

    return 0;
}

static int
derive_ref_poc(const OVRPL *const rpl, struct RPLInfo *const rpl_info, uint32_t poc)
{
    const int nb_refs = rpl->num_ref_entries;
    int i;
    for (i = 0; i < nb_refs; ++i) {
        const struct RefPic *const rp = &rpl->rp_list[i];
        struct RefInfo *const rinfo = &rpl_info->ref_info[i];
        enum RefType ref_type = rp->st_ref_pic_flag ? ST_REF
                                                    : (rp->inter_layer_ref_pic_flag ? ILRP_REF
                                                                                    : LT_REF);
        rinfo->type = ref_type;

        switch (ref_type) {
        int ref_poc = 0;
        case ST_REF:
           ref_poc = !rp->strp_entry_sign_flag ? poc +  rp->abs_delta_poc_st + 1
                                               : poc - (rp->abs_delta_poc_st + 1);
           rinfo->poc = ref_poc;

        break;
        case LT_REF:
           /* FIXME
            *    - Handle when read in header (ltrp_in_header_flag)
            *    - Compute msb part of POC when needed
            */
           ref_poc = rp->rpls_poc_lsb_lt;
                        
           rinfo->poc = ref_poc;
           ov_log(NULL, 2, "Partially supported Long Term Ref \n");

        break;
        case ILRP_REF:
           ov_log(NULL, 3, "Unsupported Inter Layer Ref \n");
           rinfo->poc = ref_poc;
        break;
        }
    }
    return 0;
}


static int
vvc_mark_refs(OVDPB *dpb, OVRPL *rpl, uint8_t poc)
{
    int i, j;
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    struct RPLInfo rpl_info;

    derive_ref_poc(rpl, &rpl_info, poc);

    for (i = 0;  i < rpl->num_ref_entries; ++i){
        int16_t ref_poc  = rpl_info.ref_info[i].poc;
        int16_t ref_type = rpl_info.ref_info[i].type;
        uint8_t flag = ref_type == ST_REF ? OV_ST_REF_PIC_FLAG : OV_LT_REF_PIC_FLAG;
        uint8_t found = 0;
        for (j = 0; j < nb_dpb_pic; j++) {
            OVPicture *ref_pic = &dpb->pictures[j];
            if (ref_pic->poc == ref_poc){
                if(ref_pic->frame->data[0]){
                    found = 1;
                    ref_pic->flags &= ~(OV_LT_REF_PIC_FLAG | OV_ST_REF_PIC_FLAG);
                    ref_pic->flags |= flag;
                }
            }
        }

        if (!found){
            /* If reference picture is not in the DPB we try create a new 
             * Picture with requested POC ID in the DPB
             */
            OVPicture *ref_pic = alloc_frame(dpb);

            if (ref_pic == NULL){
                return OVVC_ENOMEM;
            }

            ref_pic->poc    = ref_poc;
            ref_pic->cvs_id = dpb->cvs_id;

            ref_pic->flags  = 0;

            ref_pic->flags &= ~(OV_LT_REF_PIC_FLAG | OV_ST_REF_PIC_FLAG);
            ref_pic->flags |= flag;
            ov_log(NULL, 3, "Could not find ref %d for picture\n", ref_poc);
            return 0;
        }
    }

    return 0;
}

int ovdpb_output_frame(OVDPB *dpb, OVFrame *out, int flush)
{
    #if 0
    const VVCSPSData *const sps =  dpb->ps.sps_data;
    do {
        int nb_output = 0;
        int min_poc   = INT_MAX;
        int i, min_idx, ret;

        if (dpb->sh.no_output_of_prior_pics_flag == 1 && dpb->no_rasl_output_flag == 1) {
            /* Do not output previously decoded picture which are not already bumped 
             * Note that they can still be used by current pic
             */
            for (i = 0; i < nb_dpb_pic; i++) {
                OVPicture *pic = &dpb->DPB[i];
                uint8_t not_bumped = !(pic->flags & OV_BUMPED_PIC_FLAG);
                uint8_t not_current = pic->poc != dpb->poc;
                uint8_t is_output_cvs = pic->sequence == output_cvs_id;
                if (not_bumped && not_current && is_output_cvs) {
                    ovdpb_unref_frame(dpb, pic, OV_OUTPUT_PIC_FLAG);
                }
            }
        }

        /* Count pic marked for output in output cvs and find the min poc_id */
        for (i = 0; i < nb_dpb_pic; i++) {
            OVPicture *pic = &dpb->DPB[i];
            if ((pic->flags & OV_OUTPUT_PIC_FLAG) &&
                pic->sequence == output_cvs_id) {
                nb_output++;
                if (pic->poc < min_poc || nb_output == 1) {
                    min_poc = pic->poc;
                    min_idx = i;
                }
            }
        }

        /* If the number of pic to output is less than max_num_reorder_pics 
         * in current cvs we wait for more pic before outputing any
         * FIXME remove flush and provide a specialized function for flushing
         */
        if (!flush && output_cvs_id == dpb->cvs_id && dpb->ps.sps_data &&
            nb_output <= dpb->max_num_reorder_pics) {
            return 0;
        }

        if (nb_output) {
            OVPicture *pic = &dpb->pictures[min_idx];

            /* FIXME Add a reference on ouput_pic */
            ret = av_frame_ref(out, pic->frame);

            /* we unref the pic even if ref failed */
            if (pic->flags & VVC_FRAME_FLAG_BUMPING) {
                ovdpb_unref_frame(dpb, pic, OV_OUTPUT_PIC_FLAG | VVC_FRAME_FLAG_BUMPING);
            } else {
                ovdpb_unref_frame(dpb, pic, OV_OUTPUT_PIC_FLAG);
            }

            if (ret < 0) {
                return ret;
            }

            ov_log(NULL, 4, "Got ouput picture with POC %d.\n", pic->poc);

            return 1;
        }

        /* If no output pic found increase cvs_id and retry */
        if (output_cvs_id != dpb->active_seq_id) {
            output_cvs_id = (dpb->output_seq_id + 1) & 0xff;
        } else {
            break;
        }

    } while (1);
    ov_log(NULL, 4, "No picture to output\n");

    #endif
    return 0;
}

/*FIXME 
 *   There might be better ways instead of always looping over
 *   the whole DPB and check for POC and CVS.
 */
void ovdpb_bump_frame(OVDPB *dpb, uint32_t poc, uint16_t output_cvs_id)
{
    int nb_dpb_pic = 0;
    int min_poc = INT_MAX;
    int i;

    /* Count pictures in current output target Coded Video Sequence
     * which does not correspond to current picture
     */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];
        if ((pic->flags) && pic->cvs_id == output_cvs_id &&
            pic->poc != poc) {
            nb_dpb_pic++;
        }
    }

    if (nb_dpb_pic >= dpb->max_nb_dpb_pic) {
        /* Determine the min POC among those pic
         */
        for (i = 0; i < nb_dpb_pic; i++) {
            OVPicture *pic = &dpb->pictures[i];
            if ((pic->flags) &&
                pic->cvs_id == output_cvs_id &&
                pic->poc != poc) {
                if (pic->flags == OV_OUTPUT_PIC_FLAG && pic->poc < min_poc) {
                    min_poc = pic->poc;
                }
            }
        }

        /* Mark with bumping Flag picture with POC < to current
         */
        for (i = 0; i < nb_dpb_pic; i++) {
            OVPicture *pic = &dpb->pictures[i];
            if (pic->flags & OV_OUTPUT_PIC_FLAG &&
                pic->cvs_id == output_cvs_id &&
                pic->poc <= min_poc) {
                pic->flags |= OV_BUMPED_PIC_FLAG;
            }
        }
        nb_dpb_pic--;
    }
}


#if 0
static int
compute_tmvp_scale(int dist_current, int dist_colocated)
{
    int scale;
    if (dist_current == dist_colocated || !dist_colocated)
        return 256;

    dist_current = av_clip_int8(dist_current);
    dist_colocated = av_clip_int8(dist_colocated);

    scale = dist_current * ((0x4000 + abs(dist_colocated >> 1)) / dist_colocated);
    scale += 32;
    scale >>= 6;
    scale = av_clip(scale, -4096, 4095);
    return scale;
}

static void
compute_tmvp_scale_info()
{
    /*FIXME move those parts to function + test tmvp is enabled */
    tmvp_ctx.scale00 = compute_tmvp_scale(dpb->active_pic->dist_ref0,
                                          dpb->active_pic->collocated_ref->dist_ref0);

    tmvp_ctx.scale10 = compute_tmvp_scale(dpb->active_pic->dist_ref1,
                                          dpb->active_pic->collocated_ref->dist_ref0);

    tmvp_ctx.scale01 = compute_tmvp_scale(dpb->active_pic->dist_ref0,
                                          dpb->active_pic->collocated_ref->dist_ref1);

    tmvp_ctx.scale11 = compute_tmvp_scale(dpb->active_pic->dist_ref1,
                                          dpb->active_pic->collocated_ref->dist_ref1);
}
#endif

static int
mark_ref_pic_lists(OVDPB *const dpb, uint8_t slice_type, struct OVRPL *const rpl0,
                   struct OVRPL *const rpl1)
{
    const int nb_dpb_pic = sizeof(dpb->pictures) / sizeof(*dpb->pictures);
    uint32_t poc = 0;
    int i, ret;
    /* Reset all picture ref_flags except for current */
    /* FIXME difference with clear refs */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];

#if 0
        if (pic == dpb->active_pic)
            continue;
#endif

        pic->flags &= ~(OV_LT_REF_PIC_FLAG | OV_ST_REF_PIC_FLAG);
    }

    ret = vvc_mark_refs(dpb, rpl0, poc);

    if (ret < 0) {
        goto fail;
    }

    if (slice_type == SLICE_B){

        ret = vvc_mark_refs(dpb, rpl1, poc);
        if (ret < 0) {
            goto fail;
        }
    }

    /* FIXME TMVP here compute_tmvp_scale_info()
     */

    /* Unreference all non marked Picture */
    for (i = 0; i < nb_dpb_pic; i++) {
        OVPicture *pic = &dpb->pictures[i];
        ovdpb_unref_frame(dpb, pic, 0);
    }

    return 0;

fail:
    /* FIXME ref_marking failed but we allocated a new ref
     * to replace the missing one
     */
    return 1;
}

/* TODO rename to ovdpb_init_pic();*/
int
ovdpb_init_picture(OVDPB *dpb, const OVPH *const ph, const OVSPS *const sps, uint8_t nalu_type){

    int ret = 0;
    OVFrame *frame = NULL;
    uint32_t poc = 0;
    uint8_t cra_flag = 0;
    uint8_t idr_flag = 0;

    idr_flag |= nalu_type == OVNALU_IDR_W_RADL;
    idr_flag |= nalu_type == OVNALU_IDR_N_LP;

    /*FIXME Clarify how GDR NALO are supposed to be handled 
     * At the current time we consider it to have the same
     * effect as a CRA 
     */
    cra_flag |= nalu_type == OVNALU_CRA;
    cra_flag |= nalu_type == OVNALU_GDR;

    #if 0
    /* This means decoder was flushed or previous thread received EOS/EOB
     * thus sequence changed in previous thread thus we wait for a refresh
     * picture 
     * FIXME we need to handle this from DPB state 
     */
    if (dpb->max_ra == INT_MAX) {
        if (idr_cra_flag & (VVC_CRA_NAL_FLAG)) {
            dpb->max_ra = dpb->poc;
        } else if (idr_cra_flag & VVC_IDR_NAL_FLAG){
            dpb->max_ra = INT_MIN;
        }
    }

    dpb->no_rasl_output_flag = dpb->last_eos && (idr_cra_flag &
                              (VVC_IDR_NAL_FLAG | VVC_CRA_NAL_FLAG));
    #endif

    /* TODO move to dec init */
    if (idr_flag){
        /* New IDR involves a POC refresh and mark the start of
         * a new coded video sequence
         */
        dpb->cvs_id = (dpb->cvs_id + 1) & 0xFF;
        poc = 0;
    } else {
        /* FIXME arg should be last_poc */
        poc = derive_poc(ph->ph_pic_order_cnt_lsb,
                         sps->sps_log2_max_pic_order_cnt_lsb_minus4 + 4,
                         poc);
    }

    /* If the NALU is an Refresh Picture all previous pictures in DPB
     * can be unreferenced 
     */
    if (idr_flag | cra_flag) {
        vvc_clear_refs(dpb);
    }

    /* Find an available place in DPB and allocate/retrieve available memory
     * for the current picture data from the Frame Pool
     */
    ret = ovdpb_init_current_pic(dpb, &frame, poc);
    if (ret < 0) {
        goto fail;
    }

    /* If the picture is not an IDR Picture we set all flags to 
     * Note in VVC we might still get some ref pic list in IDR
     * Slices it is not clear whether we should still mark them
     * or not
     */
    if (!idr_flag) {
        OVRPL rpl0;
        OVRPL rpl1;
        uint8_t slice_type = SLICE_B;
        mark_ref_pic_lists(dpb, slice_type, &rpl0, &rpl1);
    }

    /* Mark previous pic for output */
    if (idr_flag | cra_flag) {
        /* FIXME */
        uint16_t out_cvs_id = (dpb->cvs_id - 1) & 0xFF;
        ovdpb_bump_frame(dpb, poc, out_cvs_id);
    }

    #if 0
    /* Place a frame on decoder output if we can */
    av_frame_unref(dpb->output_frame);
    ret = ovdpb_output_frame(dpb, dpb->output_frame, 0);
    #endif

    if (ret < 0) {
        goto fail;
    }

    return ret;

fail:
    #if 0
    if (dpb->active_pic){
        ovdpb_unref_frame(dpb, dpb->active_pic, ~0);
        dpb->active_pic = NULL;
    }
    #endif

    return ret;
}
