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



void LPF_GLOBAL(hswap,HSWAP)(lpf_blas_int_t *, lpf_float16_t *, lpf_blas_int_t *, lpf_float16_t *, lpf_blas_int_t *);

void check_hswap( bool *ok) {
#define NTEST 8
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 4, 3, 3, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 2, 1, 1, 1};
    lpf_blas_int_t incy[NTEST] = { 1, 1, 1, 2, 1, 2, 1, 1};

    lpf_float16_t x[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {3,-4, 0,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,4,-4,-5,0.4,12,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 }
    };
    lpf_float16_t y[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {0,3,3,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {0,0,0,0,0 , -9, -9, -9, -9, -9},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,4,-4,-5,0.4,12,0,0,0,0 }
    };

    lpf_float16_t ytemp[10];
    lpf_float16_t xtemp[10];
    lpf_blas_int_t k;

    for ( int i = 0; i < NTEST; i++) {
        bool lok = true;
        for (k = 0; k < 10; ++k) {
            ytemp[k] = y[i][k];
            xtemp[k] = x[i][k];
        }
        LPF_GLOBAL(hswap, HSWAP)(&n[i], xtemp, &incx[i], ytemp, &incy[i]);

        for ( k = 0;  k < n[i]; k++) {
            if ( xtemp[k*incx[i]] != y[i][k*incy[i]] || ytemp[k*incy[i]] != x[i][k*incx[i]] ) {
                printf("HCOPY   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d, K = %2d \n", i+1, n[i], incx[i], incy[i], k);
                *ok = false;
                lok = false;
            }
        }
        if (lok) {
            printf("HCOPY   -- PASS -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d\n", i+1, n[i], incx[i], incy[i]);
        }

    }
}

int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_hswap(&ok);
    return (ok)?0:1;
}
