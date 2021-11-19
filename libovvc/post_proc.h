#ifndef RCN_POST_PROC_H
#define RCN_POST_PROC_H

#include <stdint.h>

struct OVSEIFGrain;
struct OVVCDec;
struct ScalingInfo;

typedef void (*FGFunc)(int16_t** dstComp, int16_t** srcComp, struct OVSEIFGrain* fgrain, 
                          int pic_w, int pic_h, int poc, uint8_t isIdrPic, uint8_t enableDeblocking);

typedef void (*SLHDRFunc)(void* slhdr_context, int16_t** sdr_pic, int16_t** hdr_pic, uint8_t* SEIPayload, int pic_width, int pic_height);

struct PostProcFunctions
{
    uint8_t pp_apply_flag;
    FGFunc pp_film_grain;
    SLHDRFunc pp_sdr_to_hdr;
};

int pp_process_frame(const OVSEI* sei, OVFrame **frame_p);


//TODO: change function names.
// void fg_data_base_generation(int8_t****  dataBase, uint8_t enableDeblocking)
void fg_data_base_generation(uint8_t enableDeblocking);

void fg_grain_apply_pic(int16_t** dstComp, int16_t** srcComp, struct OVSEIFGrain* fgrain, 
                          int pic_w, int pic_h, int poc, uint8_t isIdrPic, uint8_t enableDeblocking);

void fg_grain_no_filter(int16_t** dstComp, int16_t** srcComp, struct OVSEIFGrain* fgrain, 
                          int pic_w, int pic_h, int poc, uint8_t isIdrPic, uint8_t enableDeblocking);

void pp_sample_rate_conv(uint16_t* scaled_dst, uint16_t scaled_stride, int scaledWidth, int scaledHeight, 
                        uint16_t* orgSrc, uint16_t org_stride, int orgWidth, int orgHeight, 
                        struct ScalingInfo scale_info, uint8_t luma_flag );
#endif

