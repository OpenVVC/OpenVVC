#ifndef DEC_STRUCTURES_H
#define DEC_STRUCTURES_H

#include "ovdefs.h"
#include "nvcl.h"

struct OVPartInfo
{
    /**
     * Global limit on cb size
     */
    uint8_t log2_ctu_s;
    uint8_t log2_min_cb_s;

    /**
     *  Quad tree limits in depth an log2_min size
     */
    uint8_t log2_min_qt_s;

    /**
     *  Multi type tree limits in depth and max size
     */
    uint8_t max_mtt_depth;
    uint8_t log2_max_bt_s;
    uint8_t log2_max_tt_s;

    /**
     * Transform tree limits
     */
    uint8_t log2_max_tb_s;
};

struct QPInfo
{
    /* Slice QP used to CABAC tables initialisation */
    uint8_t slice_qp;

    /* Chroma QP tables used to derive chroma QP from slice QP
     * Note this is used when delta qp is on.
     */
    struct OVChromaQPTable{
        uint8_t qp[64];
    } chroma_qp_tables[3];

    /* TODO add chroma qp offsets lists */
};

/* Main decoder structure */
struct SPSInfo
{
    struct {
      uint8_t lfnst;
      uint8_t mts;
      uint8_t mrl;
      uint8_t mip;
      uint8_t isp;
      uint8_t cclm;
    } tool_flags;

    struct {
        uint8_t chroma_format;
        uint16_t pic_w;
        uint16_t pic_h;
        /*TODO conformance window */

    } pic_info;

    struct OVPartInfo part_info[3];

    #if 0
    struct {
    } DPBInfo;
    #endif

    /* Chroma QP tables */
    struct OVChromaQPTable qp_tables_c[3];

    uint8_t bitdepth;
};

struct PPSInfo
{

    /* Sub Picture / Tiles overrides */
    struct {
        uint16_t pic_w;
        uint16_t pic_h;
        uint8_t log2_ctb_s;
        /*TODO conformance and scaling window */

    } pic_info;

    struct {
      uint8_t pic_init_qp;
      /* TODO use offset tables */
      uint8_t qp_offset_c[3];
      struct {
        int8_t qp_offset[16];
      } offset_list[3];
    } pic_qp_info;

};

/* FIXME make PH/SH joint info for overrides ? */
struct PHInfo
{
    struct {
      uint8_t lfnst;
      uint8_t mts;
      uint8_t mrl;
      uint8_t mip;
      uint8_t isp;
      uint8_t cclm;
    } tool_flags;

    struct {
        uint8_t chroma_format;
        uint16_t pic_w;
        uint16_t pic_h;
        /*TODO conformance window */

    } pic_info;

    struct {
      uint8_t log2_ctb_s;
      uint8_t log2_min_cb_s;

    } part_info[2];

    #if 0
    struct {
    } DPBInfo;

    /* Chroma QP tables */
    struct {
    } qp_tables[3];
    #endif

    uint8_t bitdepth;
};

struct SHInfo
{

    /* Sub Picture / Tiles overrides */
    struct {
        uint16_t pic_w;
        uint16_t pic_h;
        uint8_t log2_ctb_s;
        /*TODO conformance and scaling window */

    } pic_info;

    struct {
      uint8_t pic_init_qp;
      /* TODO use offset tables */
      uint8_t qp_offset_c[3];
      struct {
        int8_t qp_offset[16];
      } offset_list[3];
    } pic_qp_info;

    /* Offset to entries from RBSP start */
    uint32_t entries_offsets[256];
    uint16_t nb_entries;
};

struct OVVCDec
{
    const char *name;

    /* NAL Units to be decoded
     * Corresponding to a Picture Unit
     */
    #if 0
    OVPictureUnit *nalu_list;
    #endif

    /* Paramters sets context */
    OVNVCLCtx nvcl_ctx;

    struct OVPS{
        /* Pointers to active parameter sets */
        const OVSPS *sps;
        const OVPPS *pps;
        const OVPH *ph;
        const OVSH *sh;

        /* Human readable information from active parameter sets */
        struct SPSInfo sps_info;
        struct PPSInfo pps_info;
        struct PHInfo ph_info;
        struct SHInfo sh_info;
    } active_params;

    struct OVDPB *dpb;

    /* List of Sub Decoders
     * Contains context for Tile / Slice / Picture / SubPicture
     * decoding
     */
    OVSliceDec *subdec_list;

    /* Informations on decoder behaviour transmitted by user
     */
    struct {
        int opt1;
    }options;
};

#endif
