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



lpf_blas_int_t LPF_GLOBAL(ibamax,IBAMAX)(lpf_blas_int_t *, lpf_bfloat16_t *, lpf_blas_int_t *);

// Testfunktion check_bnrm2
void check_ibamax( bool *ok) {
#define NTEST 6
    lpf_blas_int_t n[NTEST]    = { 5, 3, 4, 4, 1, 0};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 1, 1};
    lpf_bfloat16_t data[NTEST][10] =  {
        {1.0, -2.0, 33.0, -4.0, 5.0, 0,0,0,0,0},
        {1.5, -4.5, 3.5,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0},
        {1, 2, -1, 2, 1, 2, -1, 2, 1, 2 },
        {10,0,0,0,0,0,0,0,0,0 },
        {9, 9, 9, 9, 9 , 9, 9, 9, 9, 9}
    };
    lpf_blas_int_t ibamax_expected[NTEST] = {3, 2, 1, 1, 1, 0};

    for ( int i = 0; i < NTEST; i++) {
        lpf_blas_int_t result = LPF_GLOBAL(ibamax, IBAMAX)(&n[i], data[i], &incx[i]);
        if ( result- ibamax_expected[i] == 0) {
            printf("IBAMAX  -- PASS -- Testcase %2d: N = %3d, INCX = %3d, Result = %2d, Expected = %2d\n", (int) i+1, (int) n[i], (int) incx[i], (int)result, (int)ibamax_expected[i]);
        } else{
            printf("IBAMAX  -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, Result = %2d, Expected = %2d\n", (int) i+1, (int) n[i],  (int)incx[i], (int)result, (int)ibamax_expected[i]);
            *ok = false;
        }
    }
}

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;

    bool ok = true;
    check_ibamax(&ok);
    return (ok)?0:1;
}
