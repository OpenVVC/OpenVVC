#include <pthread.h>
/* FIXME tmp*/
#include "stdio.h"

#include "dec_structures.h"
#include "overror.h"
#include "ovutils.h"
#include "ovmem.h"


void uninit_tiles_threads(struct ThreadInfo *th_info);

static void *
t_function(void *opaque)
{
    struct TileThread *tdec = (struct TileThread *)opaque;

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
            tdec->decode_entry_tile(tdec);
        }
    }
    return NULL;
}

static void
print_id(struct TileThread *tdec)
{
    ov_log(NULL, OVLOG_DEBUG, "IDX : %d\n",tdec->idx);
}

int
init_tiles_threads(struct ThreadInfo *th_info, int nb_threads)
{
    int i;
    th_info->nb_threads = nb_threads;

    th_info->tdec = ov_mallocz(sizeof(struct TileThread) * th_info->nb_threads);

    if (!th_info->tdec) goto failalloc;

    pthread_mutex_init(&th_info->gnrl_mtx, NULL);
    pthread_cond_init(&th_info->gnrl_cnd, NULL);

    for (i = 0; i < nb_threads; ++i){
        struct TileThread *tdec = &th_info->tdec[i];
        tdec->idx = i;
        tdec->state = 0;
        tdec->end = 0;
        tdec->decode_entry_tile = &print_id;

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
    uninit_tiles_threads(th_info);

    ov_freep(&th_info->tdec);

    return OVVC_ENOMEM;

failalloc:
    return OVVC_ENOMEM;
}

void
uninit_tiles_threads(struct ThreadInfo *th_info)
{
    int i;
    void *ret;
    for (i = 0; i < th_info->nb_threads; ++i){
        struct TileThread *th_dec = &th_info->tdec[i];
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
