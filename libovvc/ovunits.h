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

#ifndef OVUNITS_H
#define OVUNITS_H

#include <stdint.h>
#include <stddef.h>
#include <stdatomic.h>

enum OVNALUType
{
    /* VCL NALU Types*/
    OVNALU_TRAIL = 0,
    OVNALU_STSA  = 1,
    OVNALU_RADL  = 2,
    OVNALU_RASL  = 3,
    OVNALU_RSVD_VCL_4 = 4,
    OVNALU_RSVD_VCL_5 = 5,
    OVNALU_RSVD_VCL_6 = 6,

    OVNALU_IDR_W_RADL = 7,
    OVNALU_IDR_N_LP   = 8,
    OVNALU_CRA  =  9,
    OVNALU_GDR  = 10,
    OVNALU_RSVD_IRAP_VCL_11 = 11,

    /* Non VCL NALU Types*/
    OVNALU_OPI  = 12,
    OVNALU_DCI  = 13,
    OVNALU_VPS  = 14,
    OVNALU_SPS  = 15,
    OVNALU_PPS  = 16,
    OVNALU_PREFIX_APS  = 17,
    OVNALU_SUFFIX_APS  = 18,
    OVNALU_PH   = 19,
    OVNALU_AUD  = 20,
    OVNALU_EOS  = 21,
    OVNALU_EOB  = 22,
    OVNALU_PREFIX_SEI  = 23,
    OVNALU_SUFFIX_SEI  = 24,
    OVNALU_FD   = 25,
    OVNALU_RSVD_NVCL_26 = 26,
    OVNALU_RSVD_NVCL_27 = 27,
    OVNALU_UNSPEC_28 = 28,
    OVNALU_UNSPEC_29 = 29,
    OVNALU_UNSPEC_30 = 30,
    OVNALU_UNSPEC_31 = 31
};

/* NAL Units */
typedef struct OVNALUnit
{

  /* The type of current NAL Unit
   *
   * At the current time this needs to be set for the decoder to
   * select appropriate reader for NAL Unit data.
   */
  enum OVNALUType type;

  /* Associated Raw Byte Sequence Payload (RBSP)
   *
   * A pointer to a buffer containing NAL Unit payload.
   * Two bytes are required for the NALU Header.
   * The buffer requires an 8 bytes padding at the end of the RBSP data
   * so that the reader can process data using 64 bits fetches.
   * Note at the current time rbsp_data does not correspond to RBSP
   * but to the position of the NALU header so two bytes ahead of
   * the RBSP start. This might change in later versions to adjust
   * naming.
   */
  const uint8_t *rbsp_data;

  /* Size in bytes of the RBSP buffer plus 2 bytes
   *
   * Two bytes are required for the NALU Header.
   * Note this does not include the 8 padding bytes
   * at the end of the RBSP data buffer.
   */
  size_t rbsp_size;

  /* Emulation Prevention 0x3 Bytes position information
   *
   * If the NAL Unit byte stream was extracted from a raw video format
   * epb_pos should point to a table containing position offsets of
   * the removed 0x3 emulation prevention bytes in the rbsp_data buffer.
   */
  const uint32_t *epb_pos;

  /* Number of Emulation Prevention 0x3 Bytes removed by demuxer
   */
  int nb_epb;

  /* Reference counter
   *
   * This is used internally to know when the NAL Unit should be freed.
   * Do not set manually. See ov_nalu_ref_instead.
   */
  atomic_uint ref_count;

  /* Release function
   *
   * Callback function which will be called once OpenVVC does not need
   * a reference to this OVNALUnit anymore.
   */
  void (*release)(struct OVNALUnit **nalu_p);

  void *release_opaque;
} OVNALUnit;

/* Picture Unit */
typedef struct OVPictureUnit
{
    /* A table of OVNALUnits pointers
     *
     * Note :
     *    - must be at least nb_nalus wide
     */
    OVNALUnit **nalus;

    /* The number of NAL Units in this Picture Unit
    */
    uint8_t nb_nalus;

    /* Decoding Time Stamp (DTS) used for Presentation Time Stamps computation
     *
     * At the current time OpenVVC uses an internal clock with a
     * time scale value of 27 000 000 and 450 000 time units
     * per tick.
     *
     * Note:
     *     - currently unused for output PTS computation.
     */
    uint64_t dts;

    /* Reference counter
     *
     * This is used internally to know when the Picture Unit should be freed.
     * Do not set manually. See ovpu_ref_instead.
     */
    atomic_uint ref_count;

    /* Release function
     *
     * Callback function which will be called once OpenVVC does not need
     * a reference to this OVPictureUnit anymore.
     */
    void (*release)(struct OVNALUnit **nalu_p);
} OVPictureUnit;

/* Reference an OVNALUnit pointer
 *
 * Add a new reference to an OVNALUnit pointed by src and increase its
 * internal reference counter.
 *
 * Note:
 *     - Do not call if the value pointed by dst_p already
 *     references an OVNALUnit
 *     - src must be a valid OVNALUnit (i.e. from the decoder output)
 */
int ov_nalu_new_ref(OVNALUnit **dst_p, OVNALUnit *src);

/* Dereference an OVNALUnit pointer
 *
 * Decrements reference counter of the OVNALUnit pointed by ovnalu_p
 * and set ovnalu_p to NULL. If the OVNALUnit reference counter drops down
 * to zero OpenVVC will call the function pointed by release.
 */
void ov_nalu_unref(OVNALUnit **nalu_p);

int ov_nalu_init(OVNALUnit *nalu);


int ovnalu_init2(OVNALUnit **nalu);

/* Reference an OVPictureUnit pointer
 *
 * Add a new reference to an OVPictureUnit pointed by src and increase its
 * internal reference counter.
 *
 * Note:
 *     - Do not call if the value pointed by dst_p already
 *     references an OVPictureUnit
 *     - src must be a valid OVPictureUnit (i.e. from the decoder output)
 */
int ovpu_new_ref(OVPictureUnit **dst_p, OVPictureUnit *src);

/* Dereference an OVPictureUnit pointer
 *
 * Decrements reference counter of the OVPictureUnit pointed by ovpu_p
 * and set ovnalu_p to NULL. If the OVPictureUnit reference counter drops down
 * to zero OpenVVC will call the function pointed by release.
 */
void ovpu_unref(OVPictureUnit **ovpu_p);

void ov_free_pu(OVPictureUnit **pu);

int ovnalu_init(OVNALUnit *nalu, const uint8_t *rbsp_data, const uint32_t *epb_offset, size_t rbsp_size,
                uint32_t nb_epb, uint8_t nalu_type, void (*release_callback)(struct OVNALUnit **));

int ovpu_init(OVPictureUnit **pu, uint8_t nb_nalus);

#endif
