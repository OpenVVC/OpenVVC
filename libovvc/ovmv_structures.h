#ifndef OVMV_STRUCTURES_H
#define OVMV_STRUCTURES_H

#include <stdint.h>
#include "ctudec.h"

struct MVPDataP
{
    OVMV mvd;
    uint8_t ref_idx;
    uint8_t mvp_idx;
};

struct AffineMVPDataP
{
    struct AffineControlInfo mvd;
    uint8_t ref_idx;
    uint8_t mvp_idx;
};

struct MVPInfoP
{
    uint8_t cu_type;
    union {
        struct AffineMVPDataP cp_mvd_info;
        struct MVPDataP mvd_info;
    } data;
    uint8_t prec_amvr;
    uint8_t affine_type;
};


struct IBCMVPData
{
    IBCMV mvd;
    uint8_t mvp_idx;
};

enum MergeTypeP
{
    SB_MERGE,
    CIIP_MERGE,
    GPM_MERGE,
    MMVD_MERGE,
    DEFAULT_MERGE
};

struct MergeData
{
    enum MergeTypeP merge_type;
    uint8_t merge_idx;
};

struct SMVDData
{
    struct OVMV mvd;
    uint8_t mvp_idx0;
    uint8_t mvp_idx1;
};

struct AffineMVPDataB
{
    struct AffineControlInfo mvd0;
    struct AffineControlInfo mvd1;
    uint8_t ref_idx0;
    uint8_t mvp_idx0;
    uint8_t ref_idx1;
    uint8_t mvp_idx1;
};

struct MVPDataB
{
    OVMV mvd0;
    OVMV mvd1;
    uint8_t ref_idx0;
    uint8_t mvp_idx0;
    uint8_t ref_idx1;
    uint8_t mvp_idx1;
};

struct MVPInfoB
{
    uint8_t type;
    union
    {
        struct AffineMVPDataB aff_mvp;
        struct MVPDataB mvp;
        struct SMVDData smvd;
    } data;
    uint8_t prec_amvr;
    uint8_t bcw_idx;
    uint8_t affine_type;
};


#endif

