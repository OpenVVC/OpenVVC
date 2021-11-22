#ifndef RCN_FILL_REF_H
#define RCN_FILL_REF_H

#include <stdint.h>

#define LOG2_UNIT_S 2

static inline uint64_t
non_available_units_map(uint64_t ref_map, uint8_t pos, uint8_t log2_pb_s)
{
    int pos_unit = pos >> LOG2_UNIT_S;
    int nb_ref_units = ((1 << (log2_pb_s + 1)) >> LOG2_UNIT_S) + 1;

    uint64_t req_unit_msk = (1llu << (nb_ref_units + 1)) - 1;
    uint64_t avl_map  = (ref_map >> pos_unit) & req_unit_msk;
    uint64_t navl_map = avl_map ^ req_unit_msk;

    return navl_map;
}

static inline uint64_t
available_units_map(uint64_t ref_map, uint8_t pos, uint8_t log2_pb_s)
{
    int pos_unit = pos >> LOG2_UNIT_S;
    int nb_ref_units = ((1 << (log2_pb_s + 1)) >> LOG2_UNIT_S) + 1;
    uint64_t req_unit_msk = (1llu << (nb_ref_units + 1)) - 1;

    uint64_t avl_map = (ref_map >> pos_unit) & req_unit_msk;

    return avl_map;
}

#endif // RCN_FILL_REF_H
