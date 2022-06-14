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

#ifndef SLICEDEC_H
#define SLICEDEC_H
#include <stdint.h>
#include <pthread.h>
#include <stdatomic.h>

#include "ovdefs.h"
#include "ctudec.h"
#include "dec_structures.h"

struct EntryThread;

typedef int (*DecodeFunc)(OVSliceDec *sldec, OVCTUDec *const ctudec, const OVPS *const prms, uint16_t entry_idx);

struct SliceSynchro
{
    OVSliceDec *owner;
    // struct EntryThread *tdec;
    
    struct MainThread* main_thread;
    
    OVNALUnit* slice_nalu;

    //int nb_threads;
    uint8_t active_state;

    /* Information on current task */
    int nb_entries;
    atomic_uint nb_entries_decoded;

    DecodeFunc decode_entry;

    pthread_mutex_t gnrl_mtx;
    pthread_cond_t gnrl_cnd;
};

struct CCLines
{
    uint8_t *qt_depth_map_x;
    uint8_t *log2_cu_w_map_x;
    uint8_t *cu_mode_x;
    uint16_t nb_pb_w;
    /* TODO we could do the same for rows and allocate a
     * complete row instead of reset columns y buffers 
     * at each new line
     */
};

/* Structure used to retrieve above modes information for modes
 * derivation
 * FIXME realloc on picture width changes
 */
struct DRVLines
{
    /* Used for intra Most Probable Mode changes
     * Init value is set to PLANAR
     */
    uint8_t *intra_luma_x;

    /* Bit Field information on above reconstructed PU
     * Used for Intra reference construction
     * LSB correspond to first above PU
     */
    uint32_t *progress_map;

    /* Inter lines for Motion Vector Prediction */
    struct InterLines
    {
        /* Motion vectors of above line */
        OVMV *mv0;
        OVMV *mv1;

        /* Bit fields of above line */
        uint32_t *dir0;
        uint32_t *dir1;

        /* Bit fields for affine motion */
        uint32_t *affine;

        struct AffineInfo *aff_info;

    } inter_lines;

    struct DBFLines
    {
        /* QP Information for thresholds */
        int8_t *qp_x_map;
        int8_t *qp_x_map_cb;
        int8_t *qp_x_map_cr;

        /* Maps information */
        uint64_t *small_map;

        uint64_t *dbf_bs2_hor;
        uint64_t *dbf_bs2_hor_c;

        uint64_t *dbf_bs1_hor;
        uint64_t *dbf_bs1_hor_cb;
        uint64_t *dbf_bs1_hor_cr;

        uint64_t *dbf_affine;

        /* CU is large */
        uint64_t *large_map_c;
    } dbf_lines;

    struct IBCLines {
        IBCMV *mv;
        uint32_t *map;
    }ibc_lines;
    /*FIXME used */
    void *inter_data;
};


typedef struct OVSliceDec
{
   uint8_t slice_type;

   OVPS active_params;

   /* Lins for CABAC context derivation luma and chroma */
   struct CCLines cabac_lines[2];

   /* Lines used to retrieve local informations to be used 
    * by reconstructions such as MVs or intra modes
    */
   struct DRVLines drv_lines;

   /* Reference to current pic being decoded */
   OVPicture *pic;
   struct OVPicture *rpl0[16];
   struct OVPicture *rpl1[16];

   int16_t dist_ref_0[16];
   int16_t dist_ref_1[16];
   uint8_t nb_refs0;
   uint8_t nb_refs1;
   uint8_t nb_active_refs0;
   uint8_t nb_active_refs1;

   struct SliceSynchro slice_sync;

} OVSliceDec;

void slicedec_copy_params(OVSliceDec *sldec, struct OVPS* dec_params);

int slicedec_update_entry_decoder(OVSliceDec *sldec, OVCTUDec *ctudec);

int slicedec_decode_rect_entries(OVSliceDec *sldec, const OVPS *const prms, struct EntryThread* entry_th);

void slicedec_finish_decoding(OVSliceDec *sldec);

#if 0
int slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms);
#endif
int slicedec_init_lines(OVSliceDec *const sldec, const OVPS *const ps);

int slicedec_init(OVSliceDec *sldec);
void slicedec_uninit(OVSliceDec **sldec_p);
#endif
