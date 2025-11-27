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

void LPF_GLOBAL(baxpy,BAXPY)(lpf_blas_int_t *, lpf_bfloat16_t *, lpf_bfloat16_t *, lpf_blas_int_t *, lpf_bfloat16_t *, lpf_blas_int_t *);

void check_baxpy( bool *ok) {
#define NTEST 8
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 4, 3, 3, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 2, 1, 1, 1};
    lpf_blas_int_t incy[NTEST] = { 1, 1, 1, 2, 1, 2, 1, 1};
    lpf_bfloat16_t alpha[NTEST] = { 1, 2, sqrt_bf16(2.0f), 0, -1, -2, 0, 2};


    lpf_bfloat16_t x[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {3,-4, 0,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,4,-4,-5,0.4,12,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 }
    };
    lpf_bfloat16_t y[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {0,3,3,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,4,-4,-5,0.4,12,0,0,0,0 }
    };
    lpf_bfloat16_t baxpy_expected[NTEST][10] = {
        {2,-4,6,-8,10,0,0,0,0 },
        {4.5,-7.5,10.5,0,0,0,0,0,0,0 },
        {5.242188,-3.656250,-1,2,1,2,-1,2,1,2 },
        {0,3,3,0,0,0,0,0,0,0 },
        {-9,6,-1.400391,2,1,2,-1,2,1,2 },
        {-9, 9, -9, 9 , -9, 9, 9, 9, 9, 9 },
        {1,2,-1,2, 1, 2, -1, 2, 1, 2 },
        {10, 4, -4, -5, 0.399902,12,0,0,0,0 }
    };

    lpf_bfloat16_t ytemp[10];
    lpf_blas_int_t k;

    for ( int i = 0; i < NTEST; i++) {
        bool lok = true;
        for (k = 0; k < 10; ++k) {
            ytemp[k] = y[i][k];
        }
        LPF_GLOBAL(baxpy, BAXPY)(&n[i], &alpha[i], x[i], &incx[i], ytemp, &incy[i]);

        for ( k = 0;  k < 10; k++) {
            lpf_bfloat16_t mag = abs_bf16(baxpy_expected[i][k]);
            if ( diff_bf16(ytemp[k], baxpy_expected[i][k]) > __FLT16_EPSILON__ * mag * 10 ) {
                printf("BAXPY   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d, K = %2d, Result = %10f, Expected = %10f\n", i+1, n[i], incx[i], incy[i], k, (float)ytemp[k], (float)baxpy_expected[i][k]);
                *ok = false;
                lok = false;
            }
        }
        if (lok) {
            printf("BAXPY   -- PASS -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d\n", i+1, n[i], incx[i], incy[i]);
        }

    }
}

int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_baxpy(&ok);
    return (ok)?0:1;
}
