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
                        const struct ScalingInfo *const scale_info, uint8_t luma_flag );
#endif

