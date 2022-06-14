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

struct ScalingInfo
{
    uint16_t scaling_win_left;
    uint16_t scaling_win_right;
    uint16_t scaling_win_top;
    uint16_t scaling_win_bottom;
    uint8_t chroma_hor_col_flag;
    uint8_t chroma_ver_col_flag;
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
    atomic_uint ref_count;
    pthread_mutex_t pic_mtx;

    //Map of decoded CTUs
    struct PicDecodedCtusInfo {
        uint64_t mask[256][5];
        uint16_t mask_h;
        uint16_t mask_w;
        pthread_mutex_t *ref_mtx;
        pthread_cond_t  *ref_cnd;
    } decoded_ctus;

    atomic_uint *idx_function;
    FrameSynchroFunction ovdpb_frame_synchro[2];

    struct MVPlane mv_plane0;
    struct MVPlane mv_plane1;

    OVSEI *sei;

    struct ScalingInfo scale_info;
    struct PictureInternal
    {
        atomic_uint idx_function;
        pthread_mutex_t ref_mtx;
        pthread_cond_t  ref_cnd;
    } internal;

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

   struct DPBInternal internal;
   uint64_t nb_units_in_ticks;
   uint64_t pts;
};

int ovdpb_init(OVDPB **dpb_p, const OVPS *ps);

void ovdpb_uninit(OVDPB **dpb_p);

int ovdpb_init_picture(OVDPB *dpb, OVPicture **pic, const OVPS *const ps, uint8_t nalu_type, 
                   OVSliceDec *const sldec, const OVVCDec *ovdec);

void ovdpb_flush_dpb(OVDPB *dpb);

void ovdpb_unref_pic(OVPicture *pic, int flags);

void ovdpb_release_pic(OVDPB *dpb, OVPicture *pic);

int ovdpb_drain_frame(OVDPB *dpb, OVFrame **out, OVSEI **sei_p);

int ovdpb_output_pic(OVDPB *dpb, OVFrame **out, OVSEI **sei_p);

void ovdpb_unmark_ref_pic_lists(uint8_t slice_type, OVSliceDec *const sldec);

void ovdpb_report_decoded_ctu_line(OVPicture *const pic, int y_ctu, int xmin_ctu, int xmax_ctu);

void ovdpb_report_decoded_frame(OVPicture *const pic);

#endif
