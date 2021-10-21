#ifndef OVVCDMX_H
#define OVVCDMX_H

/* Experimental raw video demuxer Annex B */

#include <stdio.h>
#include "ovio.h"
#include "ovunits.h"

typedef struct OVVCDmx OVVCDmx;

/* Initialize demuxer
 */
int ovdmx_init(OVVCDmx **ovdmx_p);

/* Close demuxer
 */
int ovdmx_close(OVVCDmx *ovdmx);

/* Attach an input stream to the demuxer
 */
int ovdmx_attach_stream(OVVCDmx *const ovdmx, OVIO *io);

/* Reinit the demuxer.
 */
void ovdmx_detach_stream(OVVCDmx *const ovdmx);

/* Extract a Picture Unit
 *
 * Note :
 *     - at the current time output does not correspond to
 *     a complete Picture Unit. The OVPictureUnit only contains
 *     one OVNALUnit.
 */
int ovdmx_extract_picture_unit(OVVCDmx *const ovdmx, OVPictureUnit **ovpu_p);

#endif

