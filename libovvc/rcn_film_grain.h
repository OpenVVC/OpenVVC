#ifndef RCN_FILM_GRAIN_H
#define RCN_FILM_GRAIN_H


#define MAX_NUM_INTENSITIES                             256 // Maximum nuber of intensity intervals supported in FGC SEI
#define MAX_NUM_MODEL_VALUES                              6 // Maximum nuber of model values supported in FGC SEI

#define MAX_ALLOWED_MODEL_VALUES        3
#define MIN_LOG2SCALE_VALUE             2
#define MAX_LOG2SCALE_VALUE             7
#define FILM_GRAIN_MODEL_ID_VALUE       0
#define BLENDING_MODE_VALUE             0
#define MAX_STANDARD_DEVIATION          255
#define MIN_CUT_OFF_FREQUENCY           2
#define MAX_CUT_OFF_FREQUENCY           14
#define DEFAULT_HORZ_CUT_OFF_FREQUENCY  8
#define MAX_ALLOWED_COMP_MODEL_PAIRS    10

#define GRAIN_SCALE                     6
#define COLOUR_OFFSET_LUMA           0
#define COLOUR_OFFSET_CR             85
#define COLOUR_OFFSET_CB             170

#define DATA_BASE_SIZE               64
#define NUM_CUT_OFF_FREQ             13

#define MSB16(x) ((x&0xFFFF0000)>>16)
#define LSB16(x) (x&0x0000FFFF)
#define BIT0(x) (x&0x1)


//TODO: change function names.
// void dataBaseGen(int8_t****  dataBase, uint8_t enableDeblocking)
void dataBaseGen(uint8_t enableDeblocking);

void grainSynthesizeAndBlend(int16_t** decComp, struct OVSEIFGrain* fgrain, int pic_w, int pic_h, int poc, uint8_t isIdrPic, uint8_t enableDeblocking);


#endif