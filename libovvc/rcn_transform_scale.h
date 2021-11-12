#ifndef RCN_TRANSFORM_SCALE_H
#define RCN_TRANSFORM_SCALE_H

#include <stdint.h>

struct ResidualScaleFunctions
{
    void (*scale_add_residual)(const int16_t *src, uint16_t *dst,
                               int log2_tb_w, int log2_tb_h,
                               int scale);

    void (*scale_sub_residual)(const int16_t *src, uint16_t *dst,
                               int log2_tb_w, int log2_tb_h,
                               int scale);

    void (*scale_add_half_residual)(const int16_t *src, uint16_t *dst,
                                    int log2_tb_w, int log2_tb_h,
                                    int scale);

    void (*scale_sub_half_residual)(const int16_t *src, uint16_t *dst,
                                    int log2_tb_w, int log2_tb_h,
                                    int scale);

    void (*add_residual)(const int16_t *src, uint16_t *dst,
                         int log2_tb_w, int log2_tb_h,
                         int scale);

    void (*sub_residual)(const int16_t *src, uint16_t *dst,
                         int log2_tb_w, int log2_tb_h,
                         int scale);

    void (*add_half_residual)(const int16_t *src, uint16_t *dst,
                              int log2_tb_w, int log2_tb_h,
                              int scale);

    void (*sub_half_residual)(const int16_t *src, uint16_t *dst,
                              int log2_tb_w, int log2_tb_h,
                              int scale);
};

#endif
