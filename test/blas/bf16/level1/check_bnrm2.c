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



int16_t LPF_GLOBAL(bnrm2,BNRM2)(lpf_blas_int_t *, lpf_bfloat16_t *, lpf_blas_int_t *);

static lpf_bfloat16_t from_u16(uint16_t val)
{
    lpf_bf16_u16_t r;
    r.u16 = val;
    return r.f16;
}


// Testfunktion check_bnrm2
void check_bnrm2( bool *ok) {
#define NTEST 6
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 1, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 1, 1};
    lpf_bfloat16_t data[NTEST][10] =  {
        {1.0, -2.0, 3.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -2.5, 3.5,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0},
        {1, 2, 1, 2, 1, 2, 1, 2, 1, 2 },
        {10,0,0,0,0,0,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9}
    };
    lpf_bfloat16_t bnrm2_expected[NTEST] = {7.437500, 4.5552, 0, 2.0, 10, 0};

    for ( int i = 0; i < NTEST; i++) {
        lpf_bfloat16_t result = from_u16(LPF_GLOBAL(bnrm2, BNRM2)(&n[i], data[i], &incx[i]));
        if ( diff_bf16(result, bnrm2_expected[i]) < __FLT16_EPSILON__ * 10) {
            printf("BNRM2   -- OK   -- Testcase %2d: N = %3d, INCX = %3d, Result = %10f, Expected = %10f\n", i+1, n[i], incx[i], (float)result, (float)bnrm2_expected[i]);
        } else{
            printf("BNRM2   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, Result = %10f, Expected = %10f\n", i+1, n[i], incx[i], (float)result, (float)bnrm2_expected[i]);
            *ok = false;
        }
    }
}
int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_bnrm2(&ok);
    return (ok)?0:1;
}
