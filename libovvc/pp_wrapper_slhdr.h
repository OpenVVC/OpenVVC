#ifndef SLHDR_H
#define SLHDR_H

#ifdef __cplusplus
extern "C"
{
#endif

void
pp_slhdr_no_filter(void* slhdr_context, int16_t** sdr_pic, int16_t** hdr_pic, 
                uint8_t* SEIPayload, int pic_width, int pic_height);

void 
pp_sdr_to_hdr(void* slhdr_context, int16_t** sdr_pic, int16_t** hdr_pic, 
                uint8_t* SEIPayload, int pic_width, int pic_height);

void 
pp_init_slhdr_lib(void** pslhdr_context);

void 
pp_uninit_slhdr_lib(void* slhdr_context);

#ifdef __cplusplus
}
#endif

#endif //SLHDR_H