#ifndef DEC_STRUCTURES_H
#define DEC_STRUCTURES_H

#include <stdio.h>
#include <pthread.h>

#include "ovdefs.h"
#include "nvcl.h"
#include "post_proc.h"

#define OV_BOUNDARY_LEFT_RECT      (1 << 1)
#define OV_BOUNDARY_UPPER_RECT     (1 << 3)
#define OV_BOUNDARY_RIGHT_RECT     (1 << 5)
#define OV_BOUNDARY_BOTTOM_RECT    (1 << 7)

struct MVPool;

struct RectEntryInfo {
    int tile_x;
    int tile_y;
    int ctb_x;
    int ctb_y;
    int nb_ctu_w;
    int nb_ctu_h;
    int nb_ctu_rect;
    const uint8_t *entry_start;
    const uint8_t *entry_end;
    uint8_t ngh_flag;
    uint8_t implicit_h;
    uint8_t implicit_w;
    int last_ctu_w;
    int last_ctu_h;
    int nb_ctb_pic_w;
};

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
        int8_t qp[64];
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

    struct OVPartInfo part_info[4];

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

    struct TileInfo {
        int16_t nb_ctu_w[16];
        int16_t nb_ctu_h[16];
        int16_t ctu_x[16];
        int16_t ctu_y[16];
        /*FIXME determine properly max_num_tiles */
        uint8_t nb_tile_rows;
        uint8_t nb_tile_cols;
    } tile_info;

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

    /* Entries points in  RBSP */
    const uint8_t *rbsp_entry[256];
    uint16_t nb_entries;
};

enum SAOType {
    SAO_NOT_APPLIED = 0,
    SAO_BAND,
    SAO_EDGE,
    SAO_APPLIED
};

enum SAOModeMergeTypes
{
    NUM_SAO_MERGE_TYPES=0,
    SAO_MERGE_LEFT,
    SAO_MERGE_ABOVE
};

typedef struct SAOParamsCtu {
    int offset_abs[3][4];   ///< sao_offset_abs
    int offset_sign[3][4];  ///< sao_offset_sign

    int16_t offset_val[3][5];   ///<SaoOffsetVal

    uint8_t band_position[3];   ///< sao_band_position
    uint8_t eo_class[3];        ///< sao_eo_class

    uint8_t type_idx[3];    ///< sao_type_idx
} SAOParamsCtu;


typedef struct ALFParamsCtu {
    uint8_t ctb_alf_flag;
    uint8_t ctb_alf_idx;
    uint8_t cb_alternative;
    uint8_t cr_alternative;
} ALFParamsCtu;

struct MainThread
{
    int kill;
    pthread_mutex_t main_mtx;
    pthread_cond_t main_cnd;
};

struct OVVCDec
{
    const char *name;

    /* Paramters sets context */
    OVNVCLCtx nvcl_ctx;

    struct OVPS{
        /* Pointers to active parameter sets */
        OVSPS *sps;
        OVPPS *pps;
        OVAPS *aps_alf[8];
        OVAPS *aps_alf_c;
        OVAPS *aps_cc_alf_cb;
        OVAPS *aps_cc_alf_cr;
        OVAPS *aps_lmcs;
        OVPH *ph;
        OVSH *sh;
        OVSEI *sei;

        /* Human readable information from active parameter sets */
        struct SPSInfo sps_info;
        struct PPSInfo pps_info;
        struct PHInfo ph_info;
        struct SHInfo sh_info;
        /* FIXME define this somewhere meaningful */

        struct PicPartInfo
        {
            uint16_t pic_w;
            uint16_t pic_h;

            uint16_t nb_ctb_w;
            uint16_t nb_ctb_h;

            uint16_t nb_pb_w;
            uint16_t nb_pb_h;

            uint8_t log2_min_cb_s;
            uint8_t log2_ctu_s;
        } pic_info;

    } active_params;

    //Boolean that indicates if the video is displayed
    uint8_t display_output;
    
    OVDPB *dpb;

    struct MVPool *mv_pool;

    /* List of Sub Decoders
     * Contains context for Tile / Slice / Picture / SubPicture
     * decoding
     */
    OVSliceDec **subdec_list;

    /* Number of available threads */
    int nb_frame_th;
    int nb_entry_th;

    struct MainThread main_thread;

    /* Informations on decoder behaviour transmitted by user
     */
    struct {
        int opt1;
    }options;

};

#endif
