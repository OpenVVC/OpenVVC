#include <pthread.h>
/* FIXME tmp*/
#include <stdatomic.h>
#include <stdio.h>

#include "slicedec.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"

struct EntryThread
{
    struct SliceThreads *parent;
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

void uninit_entry_threads(struct SliceThreads *th_info);

static int
thread_decode_entries(struct SliceThreads *th_info, struct EntryThread *tdec)
{
    uint16_t nb_entries      = th_info->nb_entries;
    uint16_t nb_task_threads = th_info->nb_task_threads;

    unsigned first_job = atomic_fetch_add_explicit(&th_info->first_job, 1, memory_order_acq_rel);
    unsigned entry_idx = first_job;

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
ovthread_decode_entries(struct SliceThreads *th_info, DecodeFunc decode_entry, int nb_entries)
{
    int i, is_last = 0;

    int nb_task_threads = OVMIN(nb_entries, th_info->nb_threads);

    th_info->nb_task_threads = nb_task_threads;;
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
    }

    /* Main thread wait until all gnrl_state has been set to 1 
     * by the last decoder thread
     */
    if (!is_last) {
        pthread_mutex_lock(&th_info->gnrl_mtx);
        while (!th_info->gnrl_state) {
            pthread_cond_wait(&th_info->gnrl_cnd, &th_info->gnrl_mtx);
        }
        th_info->gnrl_state = 0;
        pthread_mutex_unlock(&th_info->gnrl_mtx);
    }

    return 0;
}

static void *
thread_main_function(void *opaque)
{
    struct EntryThread *tdec = (struct EntryThread *)opaque;

    pthread_mutex_lock(&tdec->task_mtx);
    pthread_cond_signal(&tdec->task_cnd);

    tdec->state = 1;

    while (!tdec->kill){
        do {
            pthread_cond_wait(&tdec->task_cnd, &tdec->task_mtx);
            if (tdec->kill) {
                return NULL;
            }
        /*FIXME determine state value to exit loop*/
        } while (tdec->state != 0);

        if (!tdec->kill) {
            uint8_t is_last = thread_decode_entries(tdec->parent, tdec);

            /* Main thread is not last so we wake it
             * if its task has already ended
             */
            if (is_last) {
                struct SliceThreads *th_info = tdec->parent;
                pthread_mutex_lock(&th_info->gnrl_mtx);
                th_info->gnrl_state = 1;
                pthread_cond_signal(&th_info->gnrl_cnd);
                pthread_mutex_unlock(&th_info->gnrl_mtx);
            }
        }
    }
    return NULL;
}

int
init_entry_threads(struct SliceThreads *th_info, int nb_threads)
{
    int i;
    th_info->nb_threads = nb_threads;

    th_info->tdec = ov_mallocz(sizeof(struct EntryThread) * nb_threads);

    if (!th_info->tdec) {
        goto failalloc;
    }

    atomic_init(&th_info->first_job,      0);
    atomic_init(&th_info->last_entry_idx, 0);

    pthread_mutex_init(&th_info->gnrl_mtx, NULL);
    pthread_cond_init(&th_info->gnrl_cnd,  NULL);

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
uninit_entry_threads(struct SliceThreads *th_info)
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
