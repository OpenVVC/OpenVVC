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

#include <stdio.h>
#include <pthread.h>
#include "ovlog.h"

static OVLogLevel ov_log_level = OVLOG_INFO;

static const char* vvctype = "VVCDec";

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RST "\x1B[0m"

static const char *OVLOG_COLORIFY[6] = { RED, YEL, BLU, CYN, GRN, MAG};
static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

void
ovlog_set_log_level(OVLogLevel log_level)
{
    ov_log_level = log_level;
}

static void
ov_log_default(void* ctx, int log_level, const char* log_content, va_list vl)
{
    if (log_level <= ov_log_level) {
        const char* type = "NULL";
        if (ctx != NULL) {
            type = vvctype;
        }
        pthread_mutex_lock(&mtx);
        fprintf(stderr, "%s", OVLOG_COLORIFY[log_level]);
        fprintf(stderr, "[%s @ %p] : ", type, ctx);
        vfprintf(stderr, log_content, vl);
        fprintf(stderr, "%s", RST);
        pthread_mutex_unlock(&mtx);
    }
}

static void (*ov_log_callback)(void* ctx, int log_level, const char* log_content, va_list vl) = ov_log_default;

void
ovlog_set_callback(void (*log_function)(void* ctx, int log_level, const char* log_content, va_list vl))
{
    ov_log_callback = log_function;
}

void
ov_log(void* ctx, int log_level, const char* log_content, ...)
{
    if (log_level <= ov_log_level) {
        va_list args;

        va_start(args, log_content);

        ov_log_callback(ctx, log_level, log_content, args);

        va_end(args);
    }
}

