#ifndef NVCL_PRIVATE_H
#define NVCL_PRIVATE_H
#include "ovdefs.h"

#include "nvcl.h"


int dpb_parameters(OVNVCLReader *const rdr, OVDPBParams *const dpb_list,
                   int max_sub_layer_minus1, int sub_layer_info_flag);

int profile_tier_level_vps(OVNVCLReader *const rdr, uint8_t vps_pt_present_flag, uint8_t vps_ptl_max_tid);

int profile_tier_level_sps(OVNVCLReader *const rdr,  uint8_t sps_max_sublayers_minus1);

int profile_tier_level_dci(OVNVCLReader *const rdr);

int general_constraints_info(OVNVCLReader *const rdr);

/* This one is called by PH/SH reader */
int nvcl_read_header_ref_pic_lists(OVNVCLReader *const rdr, OVHRPL *const rpl_h,
                                   const OVSPS *const sps, const OVPPS *pps);

int nvcl_read_sps_ref_pic_list(OVNVCLReader *const rdr, const OVSPS *const sps,
                               OVRPL *const rpl);


#endif
