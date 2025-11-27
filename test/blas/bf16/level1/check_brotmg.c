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
#include "bf16_helper.h"



void LPF_GLOBAL(brotmg,BROTMG)(lpf_bfloat16_t *d1, lpf_bfloat16_t *d2, lpf_bfloat16_t *d3, lpf_bfloat16_t *d4, lpf_bfloat16_t *d5);

void check_brotmg(bool *ok)
{
#define NTEST 4
    lpf_blas_int_t i,k;

    lpf_bfloat16_t DAB[NTEST][4] = {
        /* 1 */ {0.1f, 0.3f, 1.2f, 0.2f},
        /* 2 */ {0.7f, 0.2f, 0.6f, 4.2f},
        /* 3 */ {0.0f, 0.0f, 0.0f, 0.0f},
        /* 4 */ {4.0f, -1.0f, 2.0f,4.0f},
    };
    lpf_bfloat16_t DTEMP[9];
    lpf_bfloat16_t DTRUE[NTEST][9] = {
        /* 1 */ { 0.0e0f, 0.0e0f,   1.30e0f,     .20e0f,  0.0e0f, 0.0e0f,     0.0e0f,   .50e0f, 0.0e0f},
        /* 2 */ { 0.0e0f, 0.0e0f,   4.50e0f,    4.20e0f,  1.0e0f, .50e0f,     0.0e0f,   0.0e0f, 0.0e0f},
        /* 3 */ { 0.0e0f, 0.0e0f,   0.0e0f,      0.0e0f, -2.0e0f, 0.0e0f,     0.0e0f,   0.0e0f, 0.0e0f},
        /* 4 */ { 0.0e0f, 0.0e0f,   0.0e0f,      4.0e0f, -1.0e0f, 0.0e0f,     0.0e0f,   0.0e0f, 0.0e0f},
    };
    // Initialize computed values for DTRUE
    DTRUE[0][0] = 12.0f / 130.0f;
    DTRUE[0][1] = 36.0f / 130.0f;
    DTRUE[0][6] = -1.0f / 6.0f;
    DTRUE[1][0] = 14.0f / 75.0f;
    DTRUE[1][1] = 49.0f / 75.0f;
    DTRUE[1][8] = 1.0f / 7.0f;


    for ( i = 0 ; i < NTEST; i++) {
        bool lok = true;

        for (k = 0; k < 4; k++) {
            DTEMP[k] = DAB[i][k];
            DTEMP[k + 4] = 0.0f;
        }
        DTEMP[8] = 0.0f;
        LPF_GLOBAL(brotmg, BROTMG)(DTEMP, DTEMP + 1, DTEMP + 2, DTEMP + 3, DTEMP + 4);

        for(k = 0; k < 9; k++) {
            lpf_bfloat16_t comp = DTEMP[k];
            lpf_bfloat16_t true_val = DTRUE[i][k];
            lpf_bfloat16_t mag = abs_bf16(DTRUE[i][k]);
            if ( diff_bf16( comp, true_val) > __FLT16_EPSILON__ * mag * 10.0) {
                printf("BROTMG  -- FAIL -- Testcase %2d: Value %d : Result = %10f, Expected = %10f\n", i+1, k, (float)comp, (float)true_val);
                lok = false;
            }
        }
        if ( lok ) {
            printf("BROTMG  -- PASS -- Testcase %2d.\n", i+1);
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
    check_brotmg(&ok);
    return (ok)?0:1;
}
