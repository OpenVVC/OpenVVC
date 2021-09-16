#ifndef OV_HLS_STRUCTURES_H
#define OV_HLS_STRUCTURES_H
#include "ovdefs.h"
#include "nvcl_structures.h"

union HLSData
{
     OVSH  sh;
     OVPH  ph;
     OVSPS sps;
     OVPPS pps;
     OVAPS aps;
     OVSEI sei;
};

struct HLSReader
{
    const char *name;
    const size_t data_size;

    uint8_t (*probe_id)(OVNVCLReader *const rdr);

    const union HLSData **(*find_storage)(OVNVCLReader *const rdr,
                                          OVNVCLCtx *const nvcl_ctx);

    int (*read)(OVNVCLReader *const rdr, OVHLSData *const hls_data,
                const OVNVCLCtx *const nvcl_ctx);

    int (*validate)(OVNVCLReader *rdr, const union HLSData *const hls_data);

    int (*replace)(const struct HLSReader *const manager,
                   const union HLSData **storage,
                   const OVHLSData *const hls_data);

    void (*free)(const union HLSData *const hls_data);
};

#endif
