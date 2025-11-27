//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of HPBLAS, a BLAS Implementation for Half-Precision
   Copyright (C) 2025 Martin Koehler

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
   */

#include <stdio.h>
#include <stdbool.h>

#include "lpf_internal.h"
#include "fp16_helper.h"



void LPF_GLOBAL(hrotg,HROTG)(lpf_float16_t *sa, lpf_float16_t *sb, lpf_float16_t *sc, lpf_float16_t *ss);


void check_hrotg(bool *ok)
{
#define NTEST 8
    lpf_float16_t SA, SB, SC, SS;
    lpf_blas_int_t i;

    lpf_float16_t DA1[NTEST] = {0.3f, 0.4f, -0.3f, -0.4f, -0.3f, 0.0f, 0.0f, 1.0f};
    lpf_float16_t DB1[NTEST] = {0.4f, 0.3f, 0.4f, 0.3f, -0.4f, 0.0f, 1.0f, 0.0f};
    lpf_float16_t DC1[NTEST] = {0.6f, 0.8f, -0.6f, 0.8f, 0.6f, 1.0f, 0.0f, 1.0f};
    lpf_float16_t DS1[NTEST] = {0.8f, 0.6f, 0.8f, -0.6f, 0.8f, 0.0f, 1.0f, 0.0f};
    lpf_float16_t DATRUE[NTEST] = {0.5f, 0.5f, 0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 1.0f};
    lpf_float16_t DBTRUE[NTEST] = {1.0f/0.6f, 0.6f, -1.0f/0.6f, -0.6f, 1.0f/0.6f, 0.0f, 1.0f, 0.0f};

    for ( i = 0 ; i < NTEST; i++) {
        bool lok = true;
        SA = DA1[i];
        SB = DB1[i];
        LPF_GLOBAL(hrotg,HROTG)(&SA, &SB, &SC, &SS);

        if ( diff_fp16( SA, DATRUE[i]) > __FLT16_EPSILON__ * abs_fp16(DATRUE[i])* 10.0) {
            printf("HROTG   -- FAIL -- Testcase %2d: SA : Result = %10f, Expected = %10f\n", i+1, (float)SA, (float)DATRUE[i]);
            lok = false;
        }
        if ( diff_fp16( SB, DBTRUE[i]) > __FLT16_EPSILON__ * abs_fp16(DBTRUE[i])* 10.0) {
            printf("HROTG   -- FAIL -- Testcase %2d: SB : Result = %10f, Expected = %10f\n", i+1, (float)SB, (float)DBTRUE[i]);
            lok = false;
        }
        if ( diff_fp16( SC, DC1[i]) > __FLT16_EPSILON__ * abs_fp16(DC1[i])* 10.0) {
            printf("HROTG   -- FAIL -- Testcase %2d: SC : Result = %10f, Expected = %10f\n", i+1, (float)SC, (float)DC1[i]);
            lok = false;
        }
        if ( diff_fp16( SS, DS1[i]) > __FLT16_EPSILON__ * abs_fp16(DS1[i])* 10.0) {
            printf("HROTG   -- FAIL -- Testcase %2d: SS : Result = %10f, Expected = %10f\n", i+1, (float)SS, (float)DS1[i]);
            lok = false;
        }

        if ( lok ) {
            printf("HROTG   -- PASS -- Testcase %2d.\n", i+1);
        } else {
            *ok = false;
        }
    }
}


int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;

    bool ok = true;
    check_hrotg(&ok);
    return (ok)?0:1;
}
