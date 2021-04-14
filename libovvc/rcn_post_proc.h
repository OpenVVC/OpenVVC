#ifndef RCN_POST_PROC_H
#define RCN_POST_PROC_H


//TODO: change function names.
// void dataBaseGen(int8_t****  dataBase, uint8_t enableDeblocking)
void dataBaseGen(uint8_t enableDeblocking);

void grainSynthesizeAndBlend(int16_t** dstComp, int16_t** srcComp, struct OVSEIFGrain* fgrain, 
    int pic_w, int pic_h, int poc, uint8_t isIdrPic, uint8_t enableDeblocking);


#endif