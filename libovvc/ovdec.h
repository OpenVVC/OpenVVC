#ifndef OPENVVC_H
#define OPENVVC_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#include "ovunits.h"

#include "ovdefs.h"


/**
 * Submit raw Annex B data corresponding to a Picture Unit
 * the decoder will then extract RBSP data of each NAL Unit
 * contained in the PU, sequentially read Non VCL data and
 * update the decoder status and attach the VCL data to a decoding
 * thread before returning
 * returns a negative number of failure, 0 otherwise
 */
int ovdec_submit_picture_unit(OVVCDec *vvcdec, const OVPictureUnit *pu);

/*
 * Request a reference to a picture from the decoder
 * and sets a new reference to the output frame
 * FIXME determine a more precise behaviour
 * returns the number of pictures available in the decoder
 * output;
 *         a negative number on failure
 * Once the user has finished working with current frame
 * the frame must be unreferenced by a calling
 * ovframe_unref()
 */
int ovdec_receive_picture(OVVCDec *dec, OVFrame **frame_p);

int ovdec_drain_picture(OVVCDec *vvcdec, OVFrame **frame);

int ovdec_init(OVVCDec **ovvcdec, int display_output, int nb_frame_th, int nb_entry_th);

int ovdec_close(OVVCDec *ovvcdec);

void ovdec_uninit_subdec_list(OVVCDec *vvcdec);

#endif
