#ifndef OV_NVCL_H
#define OV_NVCL_H
#include <stdint.h>

#define OV_MAX_NUM_VPS 16
#define OV_MAX_NUM_SPS 16
#define OV_MAX_NUM_PPS 16

typedef struct OVNVCLReader OVNVCLReader;

typedef struct OVNVCLCtx OVNVCLCtx;
typedef struct OVNVCLUnit OVNVCLUnit;

typedef struct OVOPI OVOPI;
typedef struct OVDCI OVDCI;
typedef struct OVVPS OVVPS;
typedef struct OVSPS OVSPS;
typedef struct OVPPS OVPPS;
typedef struct OVAPS OVAPS;
typedef struct OVSEI OVSEI;
typedef struct OVPH OVPH;
typedef struct OVSH OVSH;

struct OVNVCLCtx
{
    /* TODO use an other typedef to store more info in
     * lists
     */
    OVSPS *sps_list[OV_MAX_NUM_SPS];
    OVPPS *pps_list[OV_MAX_NUM_PPS];
    OVPH *ph;
    OVSH *sh;
};

struct OVNVCLReader
{
    const uint8_t *bytestream;
    const uint8_t *bytestream_end;
    uint64_t cache;
    int nb_cached_bits;
    int nb_bytes_read;
    int size_in_bits;
};

/* Attach a RBSP to the NVCL Reader
 */
int nvcl_reader_init(OVNVCLReader *rdr, const uint8_t *bytestream_start,
                     int bit_size);

void nvcl_free_ctx(OVNVCLCtx *const nvcl_ctx);


/* Reading functions */
int nvcl_opi_read(OVNVCLReader *const rdr, OVOPI *const opi,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_dci_read(OVNVCLReader *const rdr, OVDCI *const dci,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_vps_read(OVNVCLReader *const rdr, OVVPS *const vps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sps_read(OVNVCLReader *const rdr, OVSPS *const sps,
                  const OVNVCLCtx *const nvcl_ctx);

int nvcl_pps_read(OVNVCLReader *const rdr, OVPPS *const pps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_ph_read(OVNVCLReader *const rdr, OVPH *const ph,
                 OVNVCLCtx *const nvcl_ctx);

int nvcl_sei_read(OVNVCLReader *const rdr, OVSH *const sh,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sh_read(OVNVCLReader *const rdr, OVSH *const sh,
                 OVNVCLCtx *const nvcl_ctx, uint8_t nalu_type);

/* Decoding functions */
int nvcl_decode_nalu_sps(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx);

int nvcl_decode_nalu_pps(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx);

int nvcl_decode_nalu_ph(OVNVCLReader *const rdr, OVNVCLCtx *const nvcl_ctx);
#endif
