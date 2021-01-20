#ifndef CTU_DEC_H
#define CTU_DEC_H

#include "ovvcdec.h"

/* Definitions of CTUDecoder structures */
struct InterDRVCtx
{
   struct RPLInfo *rpl_0;
   struct RPLInfo *rpl_1;
   /* Information to be used by motion vector predicition (MVP)
    * Contains Motion Vector and TMVP if needed
    */
   struct MVPCtx *mvp_context;
};

struct IntraDRVCtx;

struct FiltersDRVCtx
{
    struct OVDBFInfo *dbf_ctx;
    struct OVALFInfo *alf_ctx;
    struct OVSAOInfo *sao_ctx;
    struct OVLMCSInfo *lmcs_ctx;
};

typedef struct OVCTUDec {
    /**
     * Pointer to the parent frame decoder
     * Note this should be only be used only to retrieve
     * parent decoder error handling
     */
    const struct OVDec *parent_dec;
    /* TODO decide whether or not we should add thread info
     * here
     */

    /* Associated cabac_context
     * Contains context tables, cabac status and position in the
     * current entry
     */
    struct OVCABACCtx *cabac_context;

    /* List of flag/modes to match activated tools */
    struct {
      uint8_t flag;
    } tools_status;

    /* Derivations context according to activated tools
     */
    struct OVDrvCtx {
        /* Pointers to Information used by specific tools */
        struct InterDRVCtx *inter_ctx;
        struct IntraDRVCtx *intra_ctx;
        struct DeltaQPDRVCtx *delta_qp_ctx;
        struct FiltersDRVCtx *loop_filters_ctx;
    } drv_ctx;

    /* Reconstruction context */
    struct OVRCNCtx {
        /* Pointers to the first sample data of CTU in the current
         * picture
         */
        struct {
            uint16_t *data_y;
            uint16_t *data_cb;
            uint16_t *data_cr;
        } ctu_pos;

        /* Pointers to CTU reconstruction buffers to be used
         * to reconstruct current CTU.
         * These buffers will be written to the destination picture
         * before filtering operation
         */
        struct {
            uint16_t *data_y;
            uint16_t *data_cb;
            uint16_t *data_cr;
        } rcn_ctu_buff;

        /* Side Buffer to be used by reconstruction functions
         * when needed
         */
        struct {
            uint16_t *data;
        } tmp_buff;

        /* Bit fields corresponding to the decoding progress in
         * current CTU, and its borders those are used for example
         * in order to derive references samples for intra prediction
         */
        struct CTUBitField{
           uint64_t hfield[33];
           uint64_t vfield[33];
        } progress_field;
    } rcn_ctx;

    /* CTU neighbours availability flags
     * An aggregation of flag used to tell the decoder if
     * the CTU neighbours are supposed to be known from
     * the current CTU decoder.
     */
    uint8_t ctu_ngh_flags;

} OVCTUDec;

int ovdec_decode_ctu(OVVCDec *dec, OVCTUDec *ctu_dec);

/* TODO add hook functors for PU / TU init and update */

#endif
