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



int16_t LPF_GLOBAL(hshdot,HSHDOT)(lpf_blas_int_t *, lpf_float16_t *, lpf_float16_t *, lpf_blas_int_t *, lpf_float16_t *, lpf_blas_int_t *);

static lpf_float16_t from_u16(uint16_t val)
{
    lpf_fp16_u16_t r;
    r.u16 = val;
    return r.f16;
}


void check_hshdot( bool *ok) {
#define NTEST 8
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 4, 3, 3, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 2, 1, 1, 1};
    lpf_blas_int_t incy[NTEST] = { 1, 1, 1, 2, 1, 2, 1, 1};
    lpf_float16_t b = sqrt_fp16(2.0f);

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
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,4,-4,-5,0.4,12,0,0,0,0 }
    };
    lpf_float16_t hshdot_expected[NTEST] = {55 + b , 20.75 + b , -5 + b, -3 + b, 1.599609f + b, 243 + b, 6 + b, b};

    for ( int i = 0; i < NTEST; i++) {
        lpf_float16_t result = from_u16(LPF_GLOBAL(hshdot, HSHDOT)(&n[i], &b, x[i], &incx[i], y[i], &incy[i]));
        if ( diff_fp16(result, hshdot_expected[i]) < __FLT16_EPSILON__ * 10) {
            printf("HSHDOT -- PASS -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d, Result = %10f, Expected = %10f\n", (int)i+1,(int) n[i],(int) incx[i],(int) incy[i], (float)result, (float)hshdot_expected[i]);
        } else{
            printf("HSHDOT -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d, Result = %10f, Expected = %10f\n", (int)i+1,(int) n[i],(int) incx[i],(int) incy[i], (float)result, (float)hshdot_expected[i]);
            *ok = false;
        }
    }
}

int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_hshdot(&ok);
    return (ok)?0:1;
}
