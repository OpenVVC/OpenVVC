#ifndef OV_NVCL_H
#define OV_NVCL_H

typedef struct OVNVCLReader OVNVCLReader;
typedef struct OVNVCLCtx OVNVCLCtx;

typedef struct OVOPI OVOPI;
typedef struct OVDCI OVDCI;
typedef struct OVVPS OVVPS;
typedef struct OVSPS OVSPS;
typedef struct OVPPS OVPPS;
typedef struct OVAPS OVAPS;
typedef struct OVSEI OVSEI;
typedef struct OVPH OVPH;
typedef struct OVSH OVSH;

struct OVNVCLReader
{
    const uint8_t *bytestream;
    const uint8_t *bytestream_end;
    uint64_t cache;
    int nb_cached_bits;
    int nb_bytes_read;
    int size_in_bits;
};

/* Attach a bytestream to the NVCL Reader
 * return a negative number on error
 */
int nvcl_reader_init(OVNVCLReader *rdr, const uint8_t *bytestream_start,
                     int bit_size);


/* Reading functions */
int nvcl_opi_read(OVNVCLReader *const rdr, OVOPI *const opi,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_dci_read(OVNVCLReader *const rdr, OVDCI *const dci,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_vps_read(OVNVCLReader *const rdr, OVVPS *const vps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sps_read(OVNVCLReader *const rdr, OVSPS *const sps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_pps_read(OVNVCLReader *const rdr, OVPPS *const pps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_aps_read(OVNVCLReader *const rdr, OVAPS *const aps,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_ph_read(OVNVCLReader *const rdr, OVPH *const ph,
                 OVNVCLCtx *const nvcl_ctx);

int nvcl_sei_read(OVNVCLReader *const rdr, OVSH *const sh,
                  OVNVCLCtx *const nvcl_ctx);

int nvcl_sh_read(OVNVCLReader *const rdr, OVSH *const sh,
                 OVNVCLCtx *const nvcl_ctx);

#endif
