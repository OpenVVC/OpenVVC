#ifndef RCN_STRUCTURES_H
#define RCN_STRUCTURES_H
#include <stddef.h>
#include <stdint.h>

typedef void (*MCUniDirFunc)(uint16_t *_dst, ptrdiff_t _dststride,
                             const uint16_t *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBiDir0Func)(int16_t *_dst,
                             const uint16_t *_src, ptrdiff_t _srcstride,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*MCBiDir1Func)(uint16_t *_dst, ptrdiff_t _dststride,
                             const uint16_t *_src0, ptrdiff_t _srcstride,
                             const int16_t *_src1,
                             int height, intptr_t mx, intptr_t my, int width);

typedef void (*CCLMFunc)( const uint16_t* const src_luma, uint16_t* const dst_cb,
                          uint16_t* const dst_cr, int log2_pb_w, int log2_pb_h, int y0,
                          int up_available, int left_available);

typedef void (*MDLMFunc)(const uint16_t* const src_luma, uint16_t* const dst_cb,
                         uint16_t* const dst_cr, uint64_t intra_map_rows,
                         int log2_pb_w, int log2_pb_h, int x0, int y0,
                         uint8_t left_available, uint8_t up_available);

typedef void (*ResidualAddScaleFunc)(const int16_t *src, uint16_t *dst,
                                     int log2_tb_w, int log2_tb_h,
                                     int scale);

/**
 * The Context put together all functions used by strategies.
 */

struct MCFunctions{

    MCUniDirFunc unidir[4][8];

    MCBiDir0Func bidir0[4][8];
    MCBiDir1Func bidir1[4][8];
};

struct CCLMFunctions
{
    CCLMFunc cclm;
    MDLMFunc mdlm_left;
    MDLMFunc mdlm_top;
};

struct RCNFunctions
{
    /* Motion Compensation Luma */
    struct MCFunctions mc_l;

    /* Motion Compensation Chroma */
    struct MCFunctions mc_c;

    /* Motion Compensation Chroma */
    struct CCLMFunctions cclm;

    ResidualAddScaleFunc ict[3];
};


#endif
