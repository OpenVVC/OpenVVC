/**
 *
 *   OpenVVC is open-source real time software decoder compliant with the
 *   ITU-T H.266- MPEG-I - Part 3 VVC standard. OpenVVC is developed from
 *   scratch in C as a library that provides consumers with real time and
 *   energy-aware decoding capabilities under different OS including MAC OS,
 *   Windows, Linux and Android targeting low energy real-time decoding of
 *   4K VVC videos on Intel x86 and ARM platforms.
 *
 *   Copyright (C) 2020-2022  IETR-INSA Rennes :
 *
 *   Pierre-Loup CABARAT
 *   Wassim HAMIDOUCHE
 *   Guillaume GAUTIER
 *   Thomas AMESTOY
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *   USA
 *
 **/

#ifndef CU_UTILS_H
#define CU_UTILS_H

typedef uint16_t CUFlags;

#define FLG_STORE(name, dst) \
   dst |= -(!!name) & flg_##name;

#define DECL_FLG(name, pos)\
    flg_##name = (1llu << (pos))

enum VVCCUFlag
{
     DECL_FLG(cu_skip_flag,0),
     DECL_FLG(pred_mode_flag,1),
     DECL_FLG(mip_flag,2),
     DECL_FLG(isp_flag,3),
     DECL_FLG(mrl_flag,4),
     DECL_FLG(mpm_flag,5),
     DECL_FLG(cclm_flag,6),
     DECL_FLG(mpm_flag_c,7),
     DECL_FLG(intra_bdpcm_luma_flag,8),
     DECL_FLG(intra_bdpcm_chroma_flag,9),
     DECL_FLG(intra_bdpcm_luma_dir,10),
     DECL_FLG(intra_bdpcm_chroma_dir,11),
     DECL_FLG(ibc_flag,12),
     DECL_FLG(merge_flag,2),
     DECL_FLG(inter_dir,3),
};

enum CBSize
{
    CB_1x1     = 0,
    CB_1x2     = 1,
    CB_1x4     = 2,
    CB_1x8     = 3,
    CB_1x16    = 4,
    CB_1x32    = 5,
    CB_1x64    = 6,
    CB_1x128   = 7,
    CB_2x1     = 8,
    CB_2x2     = 9,
    CB_2x4     = 10,
    CB_2x8     = 11,
    CB_2x16    = 12,
    CB_2x32    = 13,
    CB_2x64    = 14,
    CB_2x128   = 15,
    CB_4x1     = 16,
    CB_4x2     = 17,
    CB_4x4     = 18,
    CB_4x8     = 19,
    CB_4x16    = 20,
    CB_4x32    = 21,
    CB_4x64    = 22,
    CB_4x128   = 23,
    CB_8x1     = 24,
    CB_8x2     = 25,
    CB_8x4     = 26,
    CB_8x8     = 27,
    CB_8x16    = 28,
    CB_8x32    = 29,
    CB_8x64    = 30,
    CB_8x128   = 31,
    CB_16x1    = 32,
    CB_16x2    = 33,
    CB_16x4    = 34,
    CB_16x8    = 35,
    CB_16x16   = 36,
    CB_16x32   = 37,
    CB_16x64   = 38,
    CB_16x128  = 39,
    CB_32x1    = 40,
    CB_32x2    = 41,
    CB_32x4    = 42,
    CB_32x8    = 43,
    CB_32x16   = 44,
    CB_32x32   = 45,
    CB_32x64   = 46,
    CB_32x128  = 47,
    CB_64x1    = 48,
    CB_64x2    = 49,
    CB_64x4    = 50,
    CB_64x8    = 51,
    CB_64x16   = 52,
    CB_64x32   = 53,
    CB_64x64   = 54,
    CB_64x128  = 55,
    CB_128x1   = 56,
    CB_128x2   = 57,
    CB_128x4   = 58,
    CB_128x8   = 59,
    CB_128x16  = 60,
    CB_128x32  = 61,
    CB_128x64  = 62,
    CB_128x128 = 63,
};

enum CUMode {
    OV_NA = 0xFF,
    OV_INTER = 1,
    OV_INTRA = 2,
    OV_INTER_SKIP = 3,
    OV_MIP = 4,
    OV_AFFINE = 5,
    OV_INTER_SKIP_AFFINE = 6,
    OV_IBC = 7,
    OV_IBC_SKIP = 8,
};

static inline uint8_t
log2_cb_size_2_idx(uint8_t log2_cb_w, uint8_t log2_cb_h)
{
    log2_cb_w &= 0x7;
    log2_cb_h &= 0x7;
    return (log2_cb_w << 3) | log2_cb_h;
}

#endif
