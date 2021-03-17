#include <pthread.h>
/* FIXME tmp*/
#include "stdio.h"

#include "slicedec.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"

struct EntryThread
{
    pthread_t thread;
    pthread_mutex_t task_mtx;
    pthread_cond_t  task_cnd;

    /* CTU decoder associated to entry 
     * thread
     */
    OVCTUDec *ctudec;

    int state;
    int end;
};

void uninit_entry_threads(struct SliceThreads *th_info);

static void *
t_function(void *opaque)
{
    struct EntryThread *tdec = (struct EntryThread *)opaque;

    pthread_mutex_lock(&tdec->task_mtx);
    pthread_cond_signal(&tdec->task_cnd);
    tdec->state = 1;

    while (!tdec->end){
        do {
            pthread_cond_wait(&tdec->task_cnd, &tdec->task_mtx);
            if (tdec->end) {
                return NULL;
            }
        /*FIXME determine state value to exit loop*/
        } while (tdec->state != 0);

        if (!tdec->end) {
            #if 0
            tdec->decode_entry_tile(tdec);
            #endif
        }
    }
    return NULL;
}

#if 0
static void
print_id(struct EntryThread *tdec)
{
    ov_log(NULL, OVLOG_DEBUG, "IDX : %d\n",tdec->idx);
}
#endif

int
init_entry_threads(struct SliceThreads *th_info, int nb_threads)
{
    int i;
    th_info->nb_threads = nb_threads;

    th_info->tdec = ov_mallocz(sizeof(struct EntryThread) * th_info->nb_threads);

    if (!th_info->tdec) goto failalloc;

    pthread_mutex_init(&th_info->gnrl_mtx, NULL);
    pthread_cond_init(&th_info->gnrl_cnd, NULL);

    for (i = 0; i < nb_threads; ++i){
        struct EntryThread *tdec = &th_info->tdec[i];
        #if 0
        tdec->idx = i;
        tdec->decode_entry_tile = &print_id;
        #endif
        tdec->state = 0;
        tdec->end = 0;

        pthread_mutex_init(&tdec->task_mtx, NULL);
        pthread_cond_init(&tdec->task_cnd, NULL);
        pthread_mutex_lock(&tdec->task_mtx);

        if (pthread_create(&tdec->thread, NULL, t_function, tdec)) {
            pthread_mutex_unlock(&tdec->task_mtx);
            ov_log(NULL, OVLOG_ERROR, "Thread creation failed at decoder init\n");
            goto failthread;
        }

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
        th_dec->end = 1;
        pthread_cond_signal(&th_dec->task_cnd);
        pthread_mutex_unlock(&th_dec->task_mtx);

        pthread_join(th_dec->thread, &ret);
        pthread_mutex_destroy(&th_dec->task_mtx);
        pthread_cond_destroy(&th_dec->task_cnd);
    }

    pthread_mutex_destroy(&th_info->gnrl_mtx);
    pthread_cond_destroy(&th_info->gnrl_cnd);
}
