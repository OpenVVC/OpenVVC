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
    pthread_mutex_t task_mtx;
    pthread_cond_t  task_cnd;


    /* CTU decoder associated to entry 
     * thread
     */
    OVCTUDec *ctudec;

    int state;
    int kill;
};

void uninit_entry_threads(struct SliceThread *th_info);

static int
thread_decode_entries(struct SliceThread *th_info, struct EntryThread *tdec)
{
    uint16_t nb_entries      = th_info->nb_entries;
    uint16_t nb_task_threads = th_info->nb_task_threads;

    unsigned first_job = atomic_fetch_add_explicit(&th_info->first_job, 1, memory_order_acq_rel);
    unsigned entry_idx = first_job;
    ov_log(NULL, OVLOG_TRACE, "Decoder with POC %d, start entry, nb_entries %d\n", th_info->owner->pic->poc, nb_entries);
    do {
        OVCTUDec *const ctudec  = tdec->ctudec;
        OVSliceDec *const sldec = th_info->owner;
        const OVPS *const prms  = sldec->active_params;

        th_info->decode_entry(sldec, ctudec, prms, entry_idx);

        entry_idx = atomic_fetch_add_explicit(&th_info->last_entry_idx, 1, memory_order_acq_rel);

    } while (entry_idx < nb_entries);

    /* Last thread to exit loop will have entry_idx set to nb_entry + nb_task_threads - 1*/
    return entry_idx == nb_entries + nb_task_threads - 1;
}

int
ovthread_decode_entries(struct SliceThread *th_info, DecodeFunc decode_entry, int nb_entries)
{
    int i;
    int nb_task_threads = OVMIN(nb_entries, th_info->nb_threads);

    th_info->nb_task_threads = nb_task_threads;
    th_info->nb_entries = nb_entries;
    th_info->decode_entry = decode_entry;

    atomic_store_explicit(&th_info->first_job, 0, memory_order_relaxed);
    atomic_store_explicit(&th_info->last_entry_idx, nb_task_threads, memory_order_relaxed);

    /* Wake entry decoder threads by setting their state to 0 
     * and signaling on task condition
     */
    for (i = 0; i < nb_task_threads; ++i) {
        struct EntryThread *tdec = &th_info->tdec[i];
        pthread_mutex_lock(&tdec->task_mtx);
        tdec->ctudec = th_info->owner->ctudec_list[i];
        tdec->state = 0;
        pthread_cond_signal(&tdec->task_cnd);
        pthread_mutex_unlock(&tdec->task_mtx);
        ov_log(NULL, OVLOG_TRACE, "Main signals frame %d entrytask %d\n", th_info->owner->pic->poc, i);
    }


    //TODOpar: re-use when the entry task will be launched by slice threads.  
    /* Main thread wait until all gnrl_state has been set to 1 
     * by the last decoder thread
     */
    // pthread_mutex_lock(&th_info->gnrl_mtx);
    // while (th_info->gnrl_state) {
    //     pthread_cond_wait(&th_info->gnrl_cnd, &th_info->gnrl_mtx);
    // }
    // // th_info->gnrl_state = 0;
    // pthread_mutex_unlock(&th_info->gnrl_mtx);

    return 0;
}

static void *
thread_main_function(void *opaque)
{
    struct EntryThread *tdec = (struct EntryThread *)opaque;

    pthread_mutex_lock(&tdec->task_mtx);
    tdec->state = 1;
    pthread_cond_signal(&tdec->task_cnd);

    while (!tdec->kill){
        do { 
            pthread_cond_wait(&tdec->task_cnd, &tdec->task_mtx);
            if (tdec->kill && tdec->state != 0) {
                return NULL;
            }
        /*FIXME determine state value to exit loop*/
        } while (tdec->state != 0 );
        tdec->state = 1;

        uint8_t is_last = thread_decode_entries(tdec->parent, tdec);

        /* Main thread is not last so we wake it
         * if its task has already ended
         */
        if (is_last) {
            struct SliceThread *th_info = tdec->parent;
            //TODOpar: change location when using SliceThreads
            ov_nalu_unref(&th_info->slice_nalu);

            ov_log(NULL, OVLOG_TRACE, "Decoder with POC %d, finished frame \n", th_info->owner->pic->poc);
            ovdpb_unref_pic(th_info->owner->pic, OV_IN_DECODING_PIC_FLAG);
            ovdpb_unmark_ref_pic_lists(th_info->owner->slice_type, th_info->owner->pic);

            pthread_mutex_lock(&th_info->gnrl_mtx);
            th_info->gnrl_state = 0;
            // pthread_cond_signal(&th_info->gnrl_cnd);
            pthread_mutex_unlock(&th_info->gnrl_mtx);

            //Signal main thread that a slice thread is available
            struct MainThread* t_main = th_info->main_thread;
            if(t_main){
                pthread_mutex_lock(&t_main->main_mtx);
                pthread_cond_signal(&t_main->main_cnd);
                pthread_mutex_unlock(&t_main->main_mtx);
            }
        }
    }
    pthread_mutex_unlock(&tdec->task_mtx);
    return NULL;
}

int
init_entry_threads(struct SliceThread *th_info, int nb_threads)
{
    if(!th_info->tdec){
        th_info->nb_threads = nb_threads;

        th_info->tdec = ov_mallocz(sizeof(struct EntryThread) * nb_threads);

        if (!th_info->tdec) {
            goto failalloc;
        }

        atomic_init(&th_info->first_job,      0);
        atomic_init(&th_info->last_entry_idx, 0);

        pthread_mutex_init(&th_info->gnrl_mtx, NULL);
        pthread_cond_init(&th_info->gnrl_cnd,  NULL);
    }
    int i;
    for (i = 0; i < nb_threads; ++i){
        struct EntryThread *tdec = &th_info->tdec[i];

        tdec->state = 0;
        tdec->kill  = 0;

        tdec->parent = th_info;

        pthread_mutex_init(&tdec->task_mtx, NULL);
        pthread_cond_init(&tdec->task_cnd, NULL);

        pthread_mutex_lock(&tdec->task_mtx);

        if (pthread_create(&tdec->thread, NULL, thread_main_function, tdec)) {
            pthread_mutex_unlock(&tdec->task_mtx);
            ov_log(NULL, OVLOG_ERROR, "Thread creation failed at decoder init\n");
            goto failthread;
        }

        /* Wait until subdec is set */
        while (!tdec->state) {
            pthread_cond_wait(&tdec->task_cnd, &tdec->task_mtx);
        }
        pthread_mutex_unlock(&tdec->task_mtx);
    }

    return 0;

failthread:
    uninit_entry_threads(th_info);

    ov_freep(&th_info->tdec);

    return OVVC_ENOMEM;

failalloc:
    return OVVC_ENOMEM;
}

void
uninit_entry_threads(struct SliceThread *th_info)
{
    int i;
    void *ret;
    for (i = 0; i < th_info->nb_threads; ++i){
        struct EntryThread *th_dec = &th_info->tdec[i];
        pthread_mutex_lock(&th_dec->task_mtx);
        th_dec->kill = 1;
        pthread_cond_signal(&th_dec->task_cnd);
        pthread_mutex_unlock(&th_dec->task_mtx);

        pthread_join(th_dec->thread, &ret);
        pthread_mutex_destroy(&th_dec->task_mtx);
        pthread_cond_destroy(&th_dec->task_cnd);
    }

    pthread_mutex_destroy(&th_info->gnrl_mtx);
    pthread_cond_destroy(&th_info->gnrl_cnd);
    ov_freep(&th_info->tdec);

}


/*
Functions needed by the threads decoding an entire slice
*/
int
ovthread_slice_thread_init(struct SliceThread *th_slice, int nb_threads)
{   
    th_slice->nb_threads = nb_threads;
    th_slice->tdec = ov_mallocz(sizeof(struct EntryThread) * nb_threads);
    if (!th_slice->tdec) {
        goto failalloc;
    }

    atomic_init(&th_slice->first_job,      0);
    atomic_init(&th_slice->last_entry_idx, 0);

    pthread_mutex_init(&th_slice->gnrl_mtx, NULL);
    pthread_cond_init(&th_slice->gnrl_cnd,  NULL);

    // if (pthread_create(&tdec->thread, NULL, thread_main_function, tdec)) {
    //     pthread_mutex_unlock(&tdec->task_mtx);
    //     ov_log(NULL, OVLOG_ERROR, "Thread creation failed at decoder init\n");
    //     goto failthread;
    // }

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
    pthread_mutex_lock(&th_slice->gnrl_mtx);
    //TODO: need for a kill in SliceThread ?
    // th_slice->kill = 1;
    pthread_cond_signal(&th_slice->gnrl_cnd);
    pthread_mutex_unlock(&th_slice->gnrl_mtx);

    uninit_entry_threads(th_slice);

    // void *ret;
    // // pthread_join(th_slice->thread, &ret);
    // pthread_mutex_destroy(&th_slice->gnrl_mtx);
    // pthread_cond_destroy(&th_slice->gnrl_cnd);
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

                ov_log(NULL, OVLOG_DEBUG, "Got ouput picture with POC %d.\n", frame->poc);
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
            ov_log(NULL, OVLOG_DEBUG, "Got ouput picture with POC %d.\n", frame->poc);

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