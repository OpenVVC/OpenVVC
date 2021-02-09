#ifndef SLICEDEC_H
#define SLICEDEC_H

#include "ovdefs.h"
#include "ctudec.h"

struct CCLines{
    uint8_t *qt_depth_map_x;
    uint8_t *log2_cu_w_map_x;
    uint8_t *cu_mode_x;
    /* TODO we could do the same for rows and allocate a
     * complete row instead of reset columns y buffers 
     * at each new line
     */
};

typedef struct OVSliceDec
{
   /* Lins for CABAC context derivation luma and chroma */
   struct CCLines cabac_lines[2];

   OVCTUDec *ctudec_list; 
   int nb_ctudec;
} OVSliceDec;

int slicedec_init_slice_tools(OVSliceDec *const sldec, const OVPS *const prms);

int slicedec_decode_rect_entries(OVSliceDec *sldec, const OVPS *const prms);

#if 0
int slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms);
#endif

int slicedec_init(OVSliceDec **dec);
void slicedec_uninit(OVSliceDec **sldec_p);
#endif
