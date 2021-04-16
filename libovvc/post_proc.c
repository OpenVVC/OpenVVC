#include "ovdpb.h"
#include "dec_structures.h"
#include "nvcl_structures.h"

#include "post_proc.h"

void
pp_init_functions(OVVCDec *dec, struct PostProcFunctions *const pp_funcs)
{
    pp_funcs->pp_apply_flag = 0;
    const OVSEI* sei = dec->active_params.sei;
    if (sei){
        //TODO: maybe not best way to know
        // ex: FG is applied on 1st frame but not 2nd
        if (sei->sei_fg){
            pp_funcs->pp_film_grain = fg_grain_apply_pic;
            pp_funcs->pp_apply_flag = 1;
        }
        else{
            pp_funcs->pp_film_grain = NULL;
        }
    }
}

int
pp_process_frame(OVVCDec *dec, OVDPB *dpb, OVFrame **frame_p)
{
    int ret=0;
    struct PostProcFunctions pp_funcs;
    pp_init_functions(dec, &pp_funcs);

    if (pp_funcs.pp_apply_flag){
        struct Frame* frame = *frame_p;
        struct Frame* frame_post_proc;
        ret = dpbpriv_request_frame(&dpb->internal, &frame_post_proc);

        int16_t* srcComp[3] = {(int16_t*)frame->data[0], (int16_t*)frame->data[1], (int16_t*)frame->data[2]};
        int16_t* dstComp[3] = {(int16_t*)frame_post_proc->data[0], (int16_t*)frame_post_proc->data[1], 
                                (int16_t*)frame_post_proc->data[2]};

        uint8_t enable_deblock = 1;
        pp_funcs.pp_film_grain(dstComp, srcComp, dec->active_params.sei->sei_fg, 
            frame->width[0], frame->height[0], frame->poc, 0, enable_deblock);

        *frame_p = frame_post_proc;

        //TODO: handle new ref and unref on post proc picture
        // ovdpb_unref_pic(dpb, vvcdec->subdec_list->pic, ~0);
    }
    return ret;
}
