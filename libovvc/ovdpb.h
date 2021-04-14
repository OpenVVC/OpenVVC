#ifndef OVDPB_H
#define OVDPB_H

#include <stdint.h>
#include "ovunits.h"
#include "ovdefs.h"
#include "ovdpb_internal.h"
#include "ovdec_internal.h"

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
};


struct OVPicture
{
   /* Associated frame */
   OVFrame *frame;

   /* Flags used to mark Picture referenced by the
    * active picture (current picture being decoded)
    * FIXME enum ?
    */
   uint8_t flags;

   /* Pointers to ref_pic_list */
   /* FIXME use frame directly ? */
   const struct OVPicture *rpl0[16];
   const struct OVPicture *rpl1[16];

   /* FIXME Used only by TMPV? */
   #if 0
   uint32_t ref_poc0[16];
   uint32_t ref_poc1[16];
   #endif
   struct RPLInfo rpl_info0;
   struct RPLInfo rpl_info1;

   struct MVPlane mv_plane0;
   struct MVPlane mv_plane1;

   struct TMVPInfo {
       const struct OVPicture *collocated_ref;
       /* Per ref_idx Motion Scaling information */
       struct TMVPScale {
           int32_t scale;
       } scale_0[16];

       /* FIXME old compat  use sclae instead */
       int16_t scale00;
       int16_t scale01;
       int16_t scale10;
       int16_t scale11;
       /* Per CTU Bit fields for available Motion Vectors */
       void *mv_field0;
       void *mv_field1;
       /* TODO tmvp scaling */
   } tmvp;

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
   uint8_t poc;

   /* DPB status info
    * to be used to determine whether the decoder is waiting
    * for an IRAP picture or not
    */
   uint8_t state;

   struct DPBInternal internal;
};

int ovdpb_init(OVDPB **dpb_p, const OVPS *ps);

void ovdpb_uninit(OVDPB **dpb_p);

int ovdpb_init_current_pic(OVDPB *dpb, OVPicture **pic_p, int poc);

int ovdpb_init_picture(OVDPB *dpb, OVPicture **pic, const OVPS *const ps, uint8_t nalu_type, 
                   OVSliceDec *const sldec, const OVVCDec *ovdec);

void ovdpb_flush_dpb(OVDPB *dpb);

void ovdpb_unref_pic(OVDPB *dpb, OVPicture *pic, int flags);

int ovdpb_drain_frame(OVDPB *dpb, OVFrame **out, int output_cvs_id);

int ovdpb_output_frame(OVDPB *dpb, OVFrame **out, int output_cvs_id);
#endif
