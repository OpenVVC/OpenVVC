#include "ovdpb.h"
#include "ovconfig.h"

// #include "dec_structures.h"
#include "nvcl_structures.h"

#include "slicedec.h"
#include "post_proc.h"

#if ENABLE_SLHDR
#include "pp_wrapper_slhdr.h"
#endif
void
pp_init_functions(const OVSEI* sei, struct PostProcFunctions *const pp_funcs)
{
    pp_funcs->pp_apply_flag = 0;
    if (sei){
        //TODO: maybe not best way to know
        // ex: FG is applied on 1st frame but not 2nd
        if (sei->sei_fg){
            pp_funcs->pp_film_grain = fg_grain_apply_pic;
            pp_funcs->pp_apply_flag = 1;
        }
        else{
            pp_funcs->pp_film_grain = fg_grain_no_filter;
        }

#if ENABLE_SLHDR
        if (sei->sei_slhdr){
            pp_funcs->pp_sdr_to_hdr = pp_sdr_to_hdr;
            pp_funcs->pp_apply_flag = 1;
        }
        else{
            pp_funcs->pp_sdr_to_hdr = pp_slhdr_no_filter;
        }
#endif
    }
}

int
pp_process_frame(const OVSEI* sei, OVDPB *dpb, OVFrame **frame_p)
{
    int ret=0;
    struct PostProcFunctions pp_funcs;
    pp_init_functions(sei, &pp_funcs);

    //TODOpp: switch buffers src and dst when 2 or more post process are applied
    if (pp_funcs.pp_apply_flag){
        struct Frame* frame = *frame_p;
        struct Frame* frame_post_proc;
        ret = dpbpriv_request_frame(&dpb->internal, &frame_post_proc);

        int16_t* srcComp[3] = {(int16_t*)frame->data[0], (int16_t*)frame->data[1], (int16_t*)frame->data[2]};
        int16_t* dstComp[3] = {(int16_t*)frame_post_proc->data[0], (int16_t*)frame_post_proc->data[1], 
                                (int16_t*)frame_post_proc->data[2]};

        uint8_t enable_deblock = 1;
        pp_funcs.pp_film_grain(dstComp, srcComp, sei->sei_fg, 
            frame->width[0], frame->height[0], frame->poc, 0, enable_deblock);

#if ENABLE_SLHDR
        //TODOpp: redundant check with pp_init_functions
        if(sei->sei_slhdr){
            // uint8_t payload_example[87] = {0x06, 0x55, 0x04, 0x53, 0xb5, 0x00, 0x3a, 0x00, 0x01, 0x00, 0x30, 0x09, 0x00, 0x64, 0x00, 0x00, 0x21, 
            //     0x34, 0x9b, 0xaa, 0x19, 0x96, 0x08, 0xfc, 0x8a, 0x48, 0x39, 0x08, 0x3d, 0x13, 0x40, 0x42, 0x03, 0xe8, 0x00, 0x00, 0x03, 0x79, 0x01,
            //      0xd6, 0x01, 0x6e, 0x03, 0xe2, 0x00, 0x00, 0x06, 0x66, 0x00, 0x00, 0x00, 0x1a, 0x00, 0x58, 0xff, 0x00, 0x96, 0x1a, 0x15, 0x33, 
            //      0x30, 0x4d, 0x4e, 0x66, 0x69, 0x80, 0x84, 0x99, 0x9e, 0xb3, 0xb7, 0xcc, 0xd0, 0xe6, 0xe8, 0x66, 0x84, 0x9e, 0x88, 0xa6, 0x80, 
            //      0xae, 0x87, 0xc2, 0x96, 0xda, 0x87};
            // pp_funcs.pp_sdr_to_hdr(sei->sei_slhdr->slhdr_context, srcComp, dstComp, payload_example, frame->width[0], frame->height[0]);
            pp_funcs.pp_sdr_to_hdr(sei->sei_slhdr->slhdr_context, srcComp, dstComp, 
                                    sei->sei_slhdr->payload_array, frame->width[0], frame->height[0]);
        }
#endif

        *frame_p = frame_post_proc;
    }
    return ret;
}
