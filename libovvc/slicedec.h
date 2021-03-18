#ifndef SLICEDEC_H
#define SLICEDEC_H
#include <pthread.h>
#include <stdatomic.h>

#include "ovdefs.h"
#include "ctudec.h"

struct EntryThread;

typedef int (*DecodeFunc)(OVSliceDec *sldec, OVCTUDec *const ctudec, const OVPS *const prms, uint16_t entry_idx);

struct SliceThreads
{
    OVSliceDec *owner;
    struct EntryThread *tdec;

    int nb_threads;
    int gnrl_state;

    /* Information on current task */
    int nb_task_threads;
    int nb_entries;

    DecodeFunc decode_entry;
    /* Pointers to functions arguments and returns */
    void *args;
    void *rets;

    atomic_uint first_job;
    atomic_uint last_entry_idx;

    /* Slice decoder thread to be used later if
     * multiple slices
     */
    pthread_t thread;
    pthread_mutex_t gnrl_mtx;
    pthread_cond_t gnrl_cnd;
};

struct CCLines
{
    uint8_t *qt_depth_map_x;
    uint8_t *log2_cu_w_map_x;
    uint8_t *cu_mode_x;
    /* TODO we could do the same for rows and allocate a
     * complete row instead of reset columns y buffers 
     * at each new line
     */
};

/* Structure used to retrieve above modes information for modes
 * derivation
 * FIXME realloc on picture width changes
 */
struct DRVLines
{
    /* Used for intra Most Probable Mode changes
     * Init value is set to PLANAR
     */
    uint8_t *intra_luma_x;

    /* Bit Field information on above reconstructed PU
     * Used for Intra reference construction
     * LSB correspond to first above PU
     */
    uint32_t *progress_map;

    /* Inter lines for Motion Vector Prediction */
    struct InterLines
    {
        /* Motion vectors of above line */
        OVMV *mv0;
        OVMV *mv1;

        /* Bit fields of above line */
        uint32_t *dir0;
        uint32_t *dir1;

    } inter_lines;

    struct DBFLines
    {
        /* QP Information for thresholds */
        int8_t *qp_x_map;
        int8_t *qp_x_map_cb;
        int8_t *qp_x_map_cr;
        int8_t *dbf_qp_ver;
        int8_t *dbf_qp_ver_cb;
        int8_t *dbf_qp_ver_cr;

        /* Maps information */
        uint64_t *dbf_edge_ver;
        uint64_t *dbf_edge_hor;

        uint64_t *dbf_bs2_ver;
        uint64_t *dbf_bs2_hor;

        uint64_t *dbf_bs1_ver;
        uint64_t *dbf_bs1_hor;
        uint64_t *dbf_bs1_ver_cb;
        uint64_t *dbf_bs1_hor_cb;
        uint64_t *dbf_bs1_ver_cr;
        uint64_t *dbf_bs1_hor_cr;

        /* CU is large */
        uint64_t *large_map_c;
    } dbf_lines;

    /*FIXME used */
    void *inter_data;
};

typedef struct OVSliceDec
{
   uint8_t slice_type;

   const OVPS *active_params;
   /* Lins for CABAC context derivation luma and chroma */
   struct CCLines cabac_lines[2];

   /* Lines used to retrieve local informations to be used 
    * by reconstructions such as MVs or intra modes
    */
   struct DRVLines drv_lines;

   /* Reference to current pic being decoded */
   OVPicture *pic;

   OVCTUDec **ctudec_list; 
   int nb_sbdec;

   struct SliceThreads th_info;

} OVSliceDec;

int slicedec_update_entry_decoders(OVSliceDec *sldec, const OVPS *const prms);

int slicedec_decode_rect_entries(OVSliceDec *sldec, const OVPS *const prms);

#if 0
int slicedec_decode_rect_entry(OVSliceDec *sldec, const OVPS *const prms);
#endif
int slicedec_init_lines(OVSliceDec *const sldec, const OVPS *const ps);

int slicedec_init(OVSliceDec **dec, int nb_ctudec);
void slicedec_uninit(OVSliceDec **sldec_p);
#endif
