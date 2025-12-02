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



int16_t LPF_GLOBAL(hnrm2,HNRM2)(lpf_blas_int_t *, lpf_float16_t *, lpf_blas_int_t *);

static lpf_float16_t from_u16(uint16_t val)
{
    lpf_fp16_u16_t r;
    r.u16 = val;
    return r.f16;
}


// Testfunktion check_hnrm2
void check_hnrm2( bool *ok) {
#define NTEST 6
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 1, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 1, 1};
    lpf_float16_t data[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0},
        {1, 2, 1, 2, 1, 2, 1, 2, 1, 2 },
        {10,0,0,0,0,0,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9}
    };
    lpf_float16_t hnrm2_expected[NTEST] = {7.416198, 4.5552, 0, 2.0, 10, 0};

    for ( int i = 0; i < NTEST; i++) {
        lpf_float16_t result = from_u16(LPF_GLOBAL(hnrm2, HNRM2)(&n[i], data[i], &incx[i]));
        if ( diff_fp16(result, hnrm2_expected[i]) < __FLT16_EPSILON__ * 10) {
            printf("HNRM2   -- OK   -- Testcase %2d: N = %3d, INCX = %3d, Result = %10f, Expected = %10f\n",(int) i+1,(int) n[i], (int)incx[i], (float)result, (float)hnrm2_expected[i]);
        } else{
            printf("HNRM2   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, Result = %10f, Expected = %10f\n",(int) i+1,(int) n[i], (int)incx[i], (float)result, (float)hnrm2_expected[i]);
            *ok = false;
        }
    }
}
int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_hnrm2(&ok);
    return (ok)?0:1;
}
