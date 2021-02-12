#ifndef RCN_STRUCT_H
#define RCN_STRUCT_H
#include <stddef.h>
#include <stdint.h>

enum DCTType
{
    DST_VII = 0,
    DCT_VIII = 1,
    DCT_II = 2,
    NB_TR_TYPES = 4
};

#define NB_TR_SIZES 7

struct TrFunc
{

    void (*transform)(const int16_t *src, int16_t *dst,
                      ptrdiff_t src_stride,
                      int num_lines, int num_columns, int shift);
};

struct TrRCN
{
    struct TrFunc func[NB_TR_SIZES];
};

#endif
