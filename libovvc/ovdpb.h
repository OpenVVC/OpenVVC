#ifndef OVDPB_H
#define OVDPB_H

#include <stdint.h>
#include "ovunits.h"
#include "ovdefs.h"
#include "ovdpb_internal.h"
#include "ovdec_internal.h"


#define OV_OUTPUT_PIC_FLAG (1 << 0)
#define OV_LT_REF_PIC_FLAG (1 << 1)
#define OV_ST_REF_PIC_FLAG (1 << 2)
#define OV_BUMPED_PIC_FLAG (1 << 3)
#define OV_IN_DECODING_PIC_FLAG (1 << 4)


struct DPBInternal;


/* OVDPB is intended to be in charge of Frame pool management
   Coded Video sequence switch and RPL list management */

enum RefType
{
    /* Short term reference Picture */
    ST_REF = 1,

    /* Long term reference Picture */
    LT_REF = 2,

    /* Inter Layer reference Picture */
    ILRP_REF = 3
};

struct RefInfo
{
    enum RefType type;
    int32_t poc;
};

struct RPLInfo
{
   struct RefInfo ref_info[16];
   uint8_t nb_refs;
   uint8_t nb_active_refs;
};

typedef void (*FrameSynchroFunction)(const OVPicture *const ref_pic, int tl_ctu_x, 
                                    int tl_ctu_y, int br_ctu_x, int br_ctu_y);

struct OVPicture
{
   /* Associated frame */
    OVFrame *frame;

    /* Flags used to mark Picture referenced by the
    * active picture (current picture being decoded)
    * FIXME enum ?
    */
    uint8_t flags;
    atomic_uint ref_count;
    pthread_mutex_t pic_mtx;

    //Map of decoded CTUs
    struct PicDecodedCtusInfo {
        uint64_t** mask;
        int mask_h;
        int mask_w;
        pthread_mutex_t ref_mtx;
        pthread_cond_t  ref_cnd;
    } decoded_ctus;

    atomic_uint idx_function;
    FrameSynchroFunction ovdpb_frame_synchro[2];

    /* Pointers to ref_pic_list */
    /* FIXME use frame directly ? */
    /* FIXME should be const */
    struct OVPicture *rpl0[16];
    struct OVPicture *rpl0_non_active[16];
    struct OVPicture *rpl1[16];
    struct OVPicture *rpl1_non_active[16];

    struct RPLInfo rpl_info0;
    struct RPLInfo rpl_info1;

    struct MVPlane mv_plane0;
    struct MVPlane mv_plane1;

    struct TMVPInfo {
       const struct OVPicture *collocated_ref;

       int16_t dist_ref_0[16];
       int16_t dist_ref_1[16];

       int16_t dist_col_0[16];
       int16_t dist_col_1[16];

        struct ColInfoPic {
            int8_t ref_idx_rpl0;
            int8_t ref_idx_rpl1;
        } col_info;

    } tmvp;

    OVSEI *sei;

    int32_t poc;

    /* Coded Video Sequence Id to which this Picture is
    * Associated : this avoid confusing ref with same POC
    * when the refresh period is shorter than DPB
    */
    uint16_t cvs_id;
};

/* Decoded Picture Buffer
 */
struct DPB
{
   OVPicture pictures[64];

   /* FIXME new struct cvs_info
    * could permit keeping track
    * of previous cvs
    */
   uint8_t max_nb_dpb_pic;
   uint8_t max_nb_reorder_pic;
   uint8_t max_latency_increase;

   /* Coded Video Sequence Id
    * It enables to reset the POC of pictures
    * when encountering IDR Picture Units
    */
   uint16_t cvs_id;
   uint32_t poc;

   /* DPB status info
    * to be used to determine whether the decoder is waiting
    * for an IRAP picture or not
    */
   uint8_t state;

    //Boolean that indicates if the video is displayed
    uint8_t display_output;
    
   struct DPBInternal internal;
   uint64_t nb_units_in_ticks;
   uint64_t pts;
};

int ovdpb_init(OVDPB **dpb_p, const OVPS *ps);

void ovdpb_uninit(OVDPB **dpb_p);

int ovdpb_init_current_pic(OVDPB *dpb, OVPicture **pic_p, int poc);

int ovdpb_init_picture(OVDPB *dpb, OVPicture **pic, const OVPS *const ps, uint8_t nalu_type, 
                   OVSliceDec *const sldec, const OVVCDec *ovdec);

void ovdpb_flush_dpb(OVDPB *dpb);

void ovdpb_unref_pic(OVPicture *pic, int flags);

void ovdpb_release_pic(OVDPB *dpb, OVPicture *pic);

int ovdpb_drain_frame(OVDPB *dpb, OVFrame **out);

int16_t tmvp_compute_scale(int32_t dist_current, int32_t dist_colocated);

int ovdpb_output_pic(OVDPB *dpb, OVFrame **out);

int ovdpb_unmark_ref_pic_lists(uint8_t slice_type, OVPicture * current_pic);

void ovdpb_report_decoded_ctu_line(OVPicture *const pic, int y_ctu, int xmin_ctu, int xmax_ctu);

void ovdpb_report_decoded_frame(OVPicture *const pic);

void ovdpb_get_lines_decoded_ctus(OVPicture *const pic, uint64_t* decoded, int y_start, int y_end );

#endif
