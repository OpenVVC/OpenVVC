#ifndef OVUNITS_H
#define OVUNITS_H

#include <stdint.h>
#include <stddef.h>

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

typedef struct OVVRBSPData{
    uint8_t *data;
    size_t size;
}OVRBSPData;

/* NAL Units */
typedef struct OVNALUnit {
  /* Associated Raw byte sequence stream payload
   */
  const uint8_t  *rbsp_data;
  const uint32_t *epb_pos;
  size_t rbsp_size;
  int nb_epb;
  enum OVNALUType type;
} OVNALUnit;

/* Picture Unit */
typedef struct OVPictureUnit
{
    /* A vector of NAL Units */
    OVNALUnit *nalus;

    /*
     * The number of NAL Units in this Picture Unit
     * int8_t should be sufficient it would be
     * surprising to get more than 127 NALUs in a PU
     */
    int8_t nb_nalus;

    /* Time stamps if needed later */
    uint64_t dts;
    uint64_t pts;
} OVPictureUnit;

/* Acces Unit Unit
   A list of Picitures Units, in mutli layer context,
   an Acces Unit may contain more than one Picture Unit
   It shoud not be used at the current time however
   it is defined here for later use
 */
typedef struct OVAccessUnit
{
    OVPictureUnit *picture_units;
} OVAccessUnit;

int ov_init_nalu(void);

void ov_free_pu(OVPictureUnit **pu);

#endif

