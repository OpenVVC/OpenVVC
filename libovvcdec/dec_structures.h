#ifndef DEC_STRUCTURES_H
#define DEC_STRUCTURES_H

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
        uint8_t qp[64]
    } chroma_qp_tables[3];

    /* TODO add chroma qp offsets lists */
};

#endif
