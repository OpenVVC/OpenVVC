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

#ifndef OV_fRAME_H
#define OV_fRAME_H

#include <stddef.h>
#include <stdint.h>
#include <stdatomic.h>
#include "ovdefs.h"

enum ChromaFmt
{
    OV_YUV_420_P8  = 0,
    OV_YUV_420_P10 = 1,
};

struct ColorDescription
{
    uint8_t colour_primaries;
    uint8_t transfer_characteristics;
    uint8_t matrix_coeffs;
    uint8_t full_range;
};

/* Miscelaneous information on Picture */
struct FrameInfo
{
    enum ChromaFmt chroma_format;
    struct ColorDescription color_desc;
};

struct FramePool;

/* OVFrame private data */
struct FrameInternal
{
    /* reference counter */
    atomic_uint ref_count;

    struct FramePool *frame_pool;
    void *felem;
    void *pool_elem[4];
};

struct Window
{
    uint16_t offset_lft;
    uint16_t offset_rgt;
    uint16_t offset_abv;
    uint16_t offset_blw;
};

struct Frame
{
    /* Pointer to Picture data planes per component
     */
    void *data[3];

    /* Per component line size in bytes
     */
    size_t linesize[3];

    /* Per component size in bytes of allocated plane data */
    size_t size[3];

    /* Picture width and height in samples
     */
    size_t width;
    size_t height;

    /* Picture conformance window information
     */
    struct Window output_window;

    /* Picture Order Count
     */
    uint32_t poc;

    /* Presentation Time Stamp
     * (Experimental)
     */
    uint64_t pts;

    struct FrameInfo frame_info;

    /* OVFrame Private data
     * Do not modify.
     */
    struct FrameInternal internal;

    /* Opaque Data */
    void *opaque;
};

/* Reference an OVFrame pointer
 *
 * Add a new reference to an OVFrame pointed by src and increase its
 * internal reference counter.
 *
 * Note:
 *     - Do not call if the value pointed by dst_p already
 *     references an OVFrame
 *     - src must be a valid OVFrame (i.e. from the decoder output)
 */
int ovframe_new_ref(OVFrame **dst_p, OVFrame *src);

/* Dereference an OVFrame pointer
 *
 * Decrements reference counter of the OVFrame pointed by ovframe_p
 * and set ovframe_p to NULL. If the reference counter drops down to zero
 * the OVFrame pointed by ovframe_p will be freed.
 */
void ovframe_unref(OVFrame **frame_p);

#endif
