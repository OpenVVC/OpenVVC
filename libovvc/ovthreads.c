#include <pthread.h>
/* FIXME tmp*/
#include <stdatomic.h>

#include "slicedec.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"
#include "ovthreads.h"
#include "ovdpb.h"

/*
Functions for the threads decoding rectangular entries
*/
struct EntryThread
{
    struct SliceThread *parent;
    pthread_t thread;
    pthread_mutex_t entry_mtx;
    pthread_cond_t  entry_cnd;


    /* CTU decoder associated to entry 
     * thread
     */
    OVCTUDec *ctudec;

    int state;
    int kill;
};

void uninit_entry_threads(struct SliceThread *th_slice);

static int
thread_decode_entries(struct SliceThread *th_slice, struct EntryThread *tdec)
{
    uint16_t nb_entries      = th_slice->nb_entries;
    uint16_t nb_entry_th = th_slice->nb_entry_th;

    unsigned first_job = atomic_fetch_add_explicit(&th_slice->first_job, 1, memory_order_acq_rel);
    unsigned entry_idx = first_job;
    ov_log(NULL, OVLOG_DEBUG, "Decoder with POC %d, start entry, nb_entries %d\n", th_slice->owner->pic->poc, nb_entries);
    do {
        OVCTUDec *const ctudec  = tdec->ctudec;
        OVSliceDec *const sldec = th_slice->owner;
        const OVPS *const prms  = sldec->active_params;

        th_slice->decode_entry(sldec, ctudec, prms, entry_idx);

        entry_idx = atomic_fetch_add_explicit(&th_slice->last_entry_idx, 1, memory_order_acq_rel);

    } while (entry_idx < nb_entries);

    /* Last thread to exit loop will have entry_idx set to nb_entry + nb_entry_th - 1*/
    return entry_idx == nb_entries + nb_entry_th - 1;
}

int
ovthread_decode_entries(struct SliceThread *th_slice, DecodeFunc decode_entry, int nb_entries)
{
    int i;
    int nb_entry_th = OVMIN(nb_entries, th_slice->nb_entry_th);

    th_slice->nb_entries = nb_entries;
    th_slice->decode_entry = decode_entry;

    atomic_store_explicit(&th_slice->first_job, 0, memory_order_relaxed);
    atomic_store_explicit(&th_slice->last_entry_idx, th_slice->nb_entry_th, memory_order_relaxed);

    /* Wake entry decoder threads by setting their state to 0 
     * and signaling on entry condition
     */
    for (i = 0; i < nb_entry_th; ++i) {
        struct EntryThread *tdec = &th_slice->tdec[i];
        pthread_mutex_lock(&tdec->entry_mtx);
        tdec->ctudec = th_slice->owner->ctudec_list[i];
        tdec->state = 0;
        pthread_cond_signal(&tdec->entry_cnd);
        pthread_mutex_unlock(&tdec->entry_mtx);
        ov_log(NULL, OVLOG_DEBUG, "Main launches POC %d entry %d\n", th_slice->owner->pic->poc, i);
    }

    //TODOpar: re-use when the entry entry will be launched by slice threads.  
    /* Main thread wait until all active_state has been set to 1 
     * by the last decoder thread
     */
    // pthread_mutex_lock(&th_slice->gnrl_mtx);
    // while (th_slice->active_state) {
    //     pthread_cond_wait(&th_slice->gnrl_cnd, &th_slice->gnrl_mtx);
    // }
    // // th_slice->active_state = 0;
    // pthread_mutex_unlock(&th_slice->gnrl_mtx);

    return 0;
}

static void *
entry_thread_main_function(void *opaque)
{
    struct EntryThread *tdec = (struct EntryThread *)opaque;

    pthread_mutex_lock(&tdec->entry_mtx);
    tdec->state = 1;
    pthread_cond_signal(&tdec->entry_cnd);

    while (!tdec->kill){
        do { 
            pthread_cond_wait(&tdec->entry_cnd, &tdec->entry_mtx);
            if (tdec->kill && tdec->state != 0) {
                pthread_mutex_unlock(&tdec->entry_mtx);
                return NULL;
            }
        /*FIXME determine state value to exit loop*/
        } while (tdec->state != 0 );
        tdec->state = 1;

        uint8_t is_last = thread_decode_entries(tdec->parent, tdec);

        /* Main thread is not last so we wake it
         * if its entry has already ended
         */
        if (is_last) {
            slicedec_finish_decoding(tdec->parent->owner);
        }
    }
    pthread_mutex_unlock(&tdec->entry_mtx);
    return NULL;
}

int
init_entry_threads(struct SliceThread *th_slice, int nb_threads)
{
    int i;
    ov_log(NULL, OVLOG_TRACE, "Creating %d entry threads\n", nb_threads);
    for (i = 0; i < nb_threads; ++i){
        struct EntryThread *tdec = &th_slice->tdec[i];

        tdec->state = 0;
        tdec->kill  = 0;

        tdec->parent = th_slice;

        pthread_mutex_init(&tdec->entry_mtx, NULL);
        pthread_cond_init(&tdec->entry_cnd, NULL);
        pthread_mutex_lock(&tdec->entry_mtx);

        if (pthread_create(&tdec->thread, NULL, entry_thread_main_function, tdec)) {
            pthread_mutex_unlock(&tdec->entry_mtx);
            ov_log(NULL, OVLOG_ERROR, "Thread creation failed at decoder init\n");
            goto failthread;
        }

        /* Wait until subdec is set */
        while (!tdec->state) {
            pthread_cond_wait(&tdec->entry_cnd, &tdec->entry_mtx);
        }
        pthread_mutex_unlock(&tdec->entry_mtx);
    }

    return 0;

failthread:
    ov_log(NULL, OVLOG_ERROR,  "Entry threads creation failed\n");
    uninit_entry_threads(th_slice);

    ov_freep(&th_slice->tdec);

    return OVVC_ENOMEM;
}

void
uninit_entry_threads(struct SliceThread *th_slice)
{
    int i;
    void *ret;
    ov_log(NULL, OVLOG_TRACE, "Deleting %d entry threads\n", th_slice->nb_entry_th);
    for (i = 0; i < th_slice->nb_entry_th; ++i){
        struct EntryThread *th_entry = &th_slice->tdec[i];
        pthread_mutex_lock(&th_entry->entry_mtx);
        th_entry->kill = 1;
        pthread_cond_signal(&th_entry->entry_cnd);
        pthread_mutex_unlock(&th_entry->entry_mtx);

        pthread_join(th_entry->thread, &ret);
        pthread_mutex_destroy(&th_entry->entry_mtx);
        pthread_cond_destroy(&th_entry->entry_cnd);
    }

    pthread_mutex_destroy(&th_slice->gnrl_mtx);
    pthread_cond_destroy(&th_slice->gnrl_cnd);
    ov_freep(&th_slice->tdec);

}


/*
Functions needed by the threads decoding an entire slice
*/
int
ovthread_slice_thread_init(struct SliceThread *th_slice, int nb_entry_th)
{   
    th_slice->nb_entry_th = nb_entry_th;
    th_slice->tdec = ov_mallocz(sizeof(struct EntryThread) * nb_entry_th);
    if (!th_slice->tdec) {
        goto failalloc;
    }

    atomic_init(&th_slice->first_job,      0);
    atomic_init(&th_slice->last_entry_idx, 0);

    pthread_mutex_init(&th_slice->gnrl_mtx, NULL);
    pthread_cond_init(&th_slice->gnrl_cnd,  NULL);

    // if (pthread_create(&tdec->thread, NULL, entry_thread_main_function, tdec)) {
    //     pthread_mutex_unlock(&tdec->entry_mtx);
    //     ov_log(NULL, OVLOG_ERROR, "Thread creation failed at decoder init\n");
    //     goto failthread;
    // }

    init_entry_threads(th_slice, nb_entry_th);

    return 0;

// failthread:
//     ov_freep(&th_slice->tdec);
//     return OVVC_ENOMEM;

failalloc:
    return OVVC_ENOMEM;
}


//TODO: function unecessary if the pthread in SliceThread is not launched.
void
ovthread_slice_thread_uninit(struct SliceThread *th_slice)
{   
    OVSliceDec * slicedec = th_slice->owner;
    //Unmark ref pict lists of decoded pics
    pthread_mutex_lock(&th_slice->gnrl_mtx);
    if (th_slice->active_state == DECODING_FINISHED) {
        pthread_mutex_unlock(&th_slice->gnrl_mtx);
        OVPicture *slice_pic = slicedec->pic;
        if(slice_pic && (slice_pic->flags & OV_IN_DECODING_PIC_FLAG)){
            ov_log(NULL, OVLOG_TRACE, "Remove DECODING_PIC_FLAG POC: %d\n", slice_pic->poc);
            ovdpb_unref_pic(slice_pic, OV_IN_DECODING_PIC_FLAG);
            ovdpb_unmark_ref_pic_lists(slicedec->slice_type, slice_pic);

            pthread_mutex_lock(&th_slice->gnrl_mtx);
            th_slice->active_state = IDLE;
            pthread_mutex_unlock(&th_slice->gnrl_mtx);
        }
    } else {
        pthread_mutex_unlock(&th_slice->gnrl_mtx);
    }

    uninit_entry_threads(th_slice);
}


/*
Functions needed by the thread writing the output frames
*/

uint32_t write_decoded_frame_to_file(OVFrame *const frame, FILE *fp){
  uint8_t component = 0;
  uint32_t ret = 0;
  for(component=0; component<3; component++){
    uint32_t frame_size = frame->height[component] * frame->linesize[component];
    ret +=fwrite(frame->data[component], frame_size, sizeof(uint8_t), fp);
  }
  return ret;
}

static void *
ovthread_out_frame_write(void *opaque)
{
    OVVCDec *dec = (struct OVVCDec *)opaque;
    OVFrame *frame;
    struct OutputThread* t_out = &dec->output_thread;
    FILE *fout = t_out->fout;
    int nb_pic = 0;
    t_out->write = 0;
    do {
        while(!t_out->write && !t_out->kill){
            pthread_cond_wait(&t_out->gnrl_cnd, &t_out->gnrl_mtx);
        }
        t_out->write = 0;

        do {
            frame = NULL;
            ovdec_receive_picture(dec, &frame);

            /* FIXME use ret instead of frame */
            if (frame) {
                write_decoded_frame_to_file(frame, fout);
                ++nb_pic;

                ov_log(NULL, OVLOG_DEBUG, "Got output picture with POC %d.\n", frame->poc);
            }
        } while (frame);
    } while (!t_out->kill);

    int ret = 1;
    while (ret > 0) {
        frame = NULL;
        ret = ovdec_receive_picture(dec, &frame);
        if (frame) {
            write_decoded_frame_to_file(frame, fout);
            ++nb_pic;
            ov_log(NULL, OVLOG_DEBUG, "Got output picture with POC %d.\n", frame->poc);

            ovframe_unref(&frame);
        }
    }

    ov_log(NULL, OVLOG_INFO, "Decoded %d pictures\n", nb_pic);
    return NULL;
}

int
ovthread_output_init(OVVCDec *dec, FILE* fout)
{
    struct OutputThread* th_out = &dec->output_thread;
    th_out->fout = fout;
    th_out->kill = 0;
    // th_out->state = 0;

    pthread_mutex_init(&th_out->gnrl_mtx, NULL);
    pthread_cond_init(&th_out->gnrl_cnd,  NULL);

    if (pthread_create(&th_out->thread, NULL, ovthread_out_frame_write, dec)) {
        pthread_mutex_unlock(&th_out->gnrl_mtx);
        ov_log(NULL, OVLOG_ERROR, "Thread creation failed for output frame init\n");
        goto failthread;
    }
    return 0;

failthread:
    return OVVC_ENOMEM;

}


void ovthread_output_uninit(struct OutputThread* t_out)
{
    //Signal output thread that main thread is finished
    pthread_mutex_lock(&t_out->gnrl_mtx);
    t_out->kill = 1;
    pthread_cond_signal(&t_out->gnrl_cnd);
    pthread_mutex_unlock(&t_out->gnrl_mtx);

    void *ret_join;
    pthread_join(t_out->thread, &ret_join);
    pthread_mutex_destroy(&t_out->gnrl_mtx);
    pthread_cond_destroy(&t_out->gnrl_cnd);
}