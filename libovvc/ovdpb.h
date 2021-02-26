#ifndef OVDPB_H
#define OVDPB_H

#include <stdint.h>
#include "ovunits.h"
#include "ovdefs.h"
#include "ovdpb_internal.h"

/* OVDPB is intended to be in charge of Frame pool management
   Coded Video sequence switch and RPL list management */
#if 0
typedef struct OVDPB
{

  uint32_t current_poc_id;

  struct {
      uint16_t cvs_id;
      /* if last Picture Unit contained an EOS NAL Unit we need
         to increase cvs_id */
      uint8_t eos;
  }cvs_info;

} OVDPB;
#endif


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
                   OVSliceDec *const sldec);

void ovdpb_flush_dpb(OVDPB *dpb);

void ovdpb_unref_pic(OVDPB *dpb, OVPicture *pic, int flags);
#endif
