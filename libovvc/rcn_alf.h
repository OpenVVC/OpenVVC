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

#ifndef RCN_ALF_H
#define RCN_ALF_H

#include "nvcl_structures.h"

struct OVCTUDec;
struct RCNFunctions;
struct RectEntryInfo;

#define NUM_BITS  8
#define CLASSIFICATION_BLK_SIZE  32  //non-normative, local buffer size
#define ALF_VB_POS_ABOVE_CTUROW_LUMA   4
#define ALF_VB_POS_ABOVE_CTUROW_CHMA   2

#define ALF_FIXED_FILTER_NUM                           64
#define ALF_CTB_MAX_NUM_APS                             8
#define ALF_CTB_MAX_NUM_TRANSPOSE                       4
#define NUM_FIXED_FILTER_SETS                          16


typedef struct Area
{
  int x,y;
  int width,height;
} Area;

enum Direction
{
  HOR,
  VER,
  DIAG0,
  DIAG1,
  NUM_DIRECTIONS
};

typedef enum
{
  CHROMA_400        = 0,
  CHROMA_420        = 1,
  CHROMA_422        = 2,
  CHROMA_444        = 3,
  NUM_CHROMA_FORMAT = 4
}ChromaFormat;

typedef enum
{
  CHANNEL_TYPE_LUMA    = 0,
  CHANNEL_TYPE_CHROMA  = 1,
  MAX_NUM_CHANNEL_TYPE = 2
}ChannelType;

typedef struct RCNALF
{
  int16_t           filter_coeff_dec[NUM_FIXED_FILTER_SETS+ALF_CTB_MAX_NUM_APS][ALF_CTB_MAX_NUM_TRANSPOSE*MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t           filter_clip_dec[NUM_FIXED_FILTER_SETS+ALF_CTB_MAX_NUM_APS][ALF_CTB_MAX_NUM_TRANSPOSE*MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
  int16_t           chroma_coeff_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
  int16_t           chroma_clip_final[MAX_NUM_ALF_ALTERNATIVES_CHROMA][MAX_NUM_ALF_CHROMA_COEFF];
}RCNALF;

#endif
