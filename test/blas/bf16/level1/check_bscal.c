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



void LPF_GLOBAL(bscal,BSCAL)(lpf_blas_int_t *, lpf_bfloat16_t *, lpf_bfloat16_t *, lpf_blas_int_t *);

// Testfunktion check_bnrm2
void check_bscal( bool *ok) {
#define NTEST 6
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 1, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 1, 1};
    lpf_bfloat16_t data[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,0,0,0,0,0,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9}
    };
    lpf_bfloat16_t sa [ NTEST ] = { 0.1, 0.0625, 1, -4, 10, 0};
    lpf_bfloat16_t bscal_expected[NTEST][10] = {
        {0.1, -0.2, 0.3, -0.4, 0.5, 0,0,0,0,0},
        {0.09375, -0.15625, 0.21875,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0},
        {-4, 2, 4, 2, -4, 2, 4, 2, 1, 2 },
        {100,0,0,0,0,0,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9}

    };

    for ( int i = 0; i < NTEST; i++) {
        bool lok = true;
        LPF_GLOBAL(bscal, BSCAL)(&n[i], &sa[i], data[i], &incx[i]);

        for ( int k = 0;  k < 10; k++) {
            if ( diff_bf16(data[i][k], bscal_expected[i][k]) > __FLT16_EPSILON__ * 10) {
                printf("BSCAL   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, K = %2d, Result = %10f, Expected = %10f\n", i+1, n[i], incx[i], k+1, (float)data[i][k], (float)bscal_expected[i][k]);
                *ok = false;
                lok = false;
            }
        }
        if (lok) {
            printf("BSCAL   -- PASS -- Testcase %2d: N = %3d, INCX = %3d, SA = %10f\n", i+1, n[i], incx[i], (float)sa[i]);
        }
    }
}
int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_bscal(&ok);
    return (ok)?0:1;
}
