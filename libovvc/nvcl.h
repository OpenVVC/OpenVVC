#ifndef OV_NVCL_H
#define OV_NVCL_H
#include <stdint.h>

#include "ovdefs.h"
#include "ovunits.h"

#define OV_MAX_NUM_VPS 16
#define OV_MAX_NUM_SPS 16
#define OV_MAX_NUM_PPS 16
#define OV_MAX_NUM_APS 16

struct OVNVCLCtx
{
    /* TODO use an other typedef to store more info in
     * lists
     */
    OVSPS *sps_list[OV_MAX_NUM_SPS];
    OVPPS *pps_list[OV_MAX_NUM_PPS];
    OVAPS *lmcs_aps_list[OV_MAX_NUM_APS];
    OVAPS *alf_aps_list[OV_MAX_NUM_APS];
    OVPH *ph;
    OVSH *sh;
    OVSEI *sei;
};

typedef union HLSData OVHLSData;


void nvcl_free_ctx(OVNVCLCtx *const nvcl_ctx);

/* Reading functions */
int nvcl_opi_read(OVNVCLReader *const rdr, OVOPI *const opi,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_dci_read(OVNVCLReader *const rdr, OVDCI *const dci,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_vps_read(OVNVCLReader *const rdr, OVVPS *const vps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sps_read(OVNVCLReader *const rdr, OVHLSData *const sps,
                  const OVNVCLCtx *const nvcl_ctx);

int nvcl_pps_read(OVNVCLReader *const rdr, OVHLSData *const pps,
                  const OVNVCLCtx *const nvcl_ctx);

int nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_ph_read(OVNVCLReader *const rdr, OVHLSData *const ph,
                 const OVNVCLCtx *const nvcl_ctx);

int nvcl_sei_read(OVNVCLReader *const rdr, OVSH *const sh,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sh_read(OVNVCLReader *const rdr, OVSH *const sh,
                 OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

/* Decoding functions */
int nvcl_decode_nalu_hls_data(OVNVCLCtx *const nvcl_ctx, OVNALUnit *nal_unit);

int nvcl_decode_nalu_sh(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

int nvcl_decode_nalu_aps(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr, uint8_t nalu_type);

int nvcl_decode_nalu_sei(OVNVCLCtx *const nvcl_ctx, OVNVCLReader *const rdr, uint8_t nalu_type);

void copy_sei_params(OVSEI **dst_p, OVSEI *src);

void nvcl_free_sei_params(OVSEI *sei);

#endif
