#ifndef OVDPB_H
#define OVDPB_H

#include <stdint.h>
#include "libovvcdmx/ovunits.h"

typedef struct OVPicture OVPicture;
typedef struct DPB OVDPB;

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


#endif
