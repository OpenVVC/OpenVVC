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

#ifndef OPENVVC_H
#define OPENVVC_H

#include <stdarg.h>
#include <stdint.h>
#include <stddef.h>

#include "ovunits.h"
#include "ovdefs.h"

enum OVOptions {
   /* Set the number of pictures being simultaneously decoded
    */
   OVDEC_NB_FRAME_THREADS = 0,

   /* Set the number of entries being simultaneously decoded
    *
    * Note:
    *    - Only tiles entries is currently supported.
    *    - The number of entries depends on the content of the
    *    coded video, so there will not be any actual speed up
    *    if the Picture Unit ony consists of one Slice and one
    *    entry.
    */
   OVDEC_NB_ENTRY_THREADS = 1,

   OVDEC_RPR_UPSCALE = 2,

   OVDEC_NB_OPTIONS,
};

typedef OVVCDec OVDec;

/**
 * Submit a new Picture Unit to the decoder.
 *
 * The decoder will process the Picture unit NAL Unit per NAL Unit.
 * It will try to read new Parameters Sets if any, ignoring unsupported
 * NAL Units types.
 * At each encountered VCL NAL Unit it will read the Slice Header, activate
 * the required Parameters sets and update the decoder status accordingly.
 * Once the decoder status updated it will wait until at least one sub decoder
 * thread is available before resuming to the processing of the remaining
 * NAL Units.
 *
 * return 0 on success,
 *        a negative number on failure.
 *
 * Notes:
 *    - Multiple Slices in same PU are currently unsupported.
 */
int ovdec_submit_picture_unit(OVDec *ovdec, const OVPictureUnit *pu);

/**
 * Request a reference to a picture from the decoder
 *
 * return the number of pictures available in the decoder output
 *         a negative number on failure
 *
 * Notes:
 *     - Once the user has finished working with current frame
 *     the frame must be unreferenced by calling ovframe_unref().
 */
int ovdec_receive_picture(OVDec *ovdec, OVFrame **frame_p);

/**
 * Drain the last pictures from the decoder output.
 *
 * Behaves the same as ovdec_receive_picture() except that it waits
 * for at least one picture to be decoded if any picture is still decoding
 * and forces the available picture out of the DPB without any regards to
 * the DPB status. See ovdec_receive_picture() for details.
 *
 * This is intended to retrieve the last pictures at the end of a sequence
 * if no EOS or EOB are explicitly signalled.
 *
 * Notes:
 *     - At the current time the decoder needs to be closed and restarted
 *     after a call to this function before attempting to resubmit a new
 *     Picture Unit. This is expected to change in future versions.
 */
int ovdec_drain_picture(OVDec *ovdec, OVFrame **frame_p);

/**
 * Initialise the OpenVVC Decoder
 *
 * return 0 on success,
 *        a negative number on failure.
 *
 * Notes:
 *    - Multiple Slices in same PU are currently unsupported.
 */
int ovdec_init(OVDec **ovdec_p);

int ovdec_config_threads(OVDec *ovdec, int nb_entry_th, int max_nb_frame_th);

int ovdec_start(OVDec *ovdec);

/**
 * Close OpenVVC decoder
 *
 * return 0 on success,
 *        a negative number on failure.
 *
 * Notes:
 *    - Multiple Slices in same PU are currently unsupported.
 */
int ovdec_close(OVDec *ovdec);

/**
 * Attempt to set a decoder options
 *
 * return 0 on success,
 *        a negative number on failure.
 *
 * Notes:
 *    - Setting some options such as the number of threads once the decoder
 *    received its first Picture Unit will not have any effect.
 *    See OVOptions for more details.
 */
int ovdec_set_option(OVDec *ovdec, enum OVOptions opt_id, int value);

void ovdec_set_log_callback(void (*log_function)(void* ctx, int log_level, const char* log_content, va_list vl));

const char* ovdec_version(void);

const char* ovdec_get_version(void);

#endif
