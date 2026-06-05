//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of LPF-BLAS, a BLAS Implementation for Half-Precision
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
#include "fp16_helper.h"
#include "lpf_internal.h"
#include <stdint.h>

/**
 * @brief Half-precision Euclidean norm.
 *
 * Performs the operation: result := ||x||_2
 *
 * @param[in] n Number of elements in vector x.
 * @param[in] x Vector x.
 * @param[in] incx Increment for the elements of x.
 * @return The Euclidean norm of vector x.
 */
int16_t LPF_GLOBAL(hnrm2, HNRM2)(int64_t* n, lpf_float16_t* x, int64_t* incx)
{

    int64_t i__1, i__2;
    lpf_float16_t ret_val, r__1;

    int64_t ix;
    lpf_float16_t ssq, norm, scale, absxi;

    --x;

    if (*n < 1 || *incx < 1)
    {
        norm = 0.f;
    }
    else if (*n == 1)
    {
        norm = abs_fp16(x[1]);
    }
    else
    {
        scale = 0.f;
        ssq = 1.f;

        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2)
        {
            if (x[ix] != 0.f)
            {
                r__1 = x[ix];
                absxi = abs_fp16(x[ix]);
                if (scale < absxi)
                {

                    r__1 = scale / absxi;
                    ssq = ssq * (r__1 * r__1) + 1.f;
                    scale = absxi;
                }
                else
                {

                    r__1 = absxi / scale;
                    ssq += r__1 * r__1;
                }
            }
        }
        norm = scale * sqrt_fp16(ssq);
    }

    ret_val = norm;

    RETURN_FP16(ret_val);
}

#include <ISO_Fortran_binding.h>

lpf_ffloat16_t lpf_blas_hnrm2_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _x,
                                                  int64_t* incx)
{
    lpf_float16_t* x = _x->base_addr;
    lpf_float16_t r_val = LPF_GLOBAL(hnrm2, HNRM2)(n, x, incx);
    lpf_ffloat16_t r;
    r.value = r_val;
    return r;
}

lpf_ffloat16_t lpf_blas_hnrm2_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _x,
                                                  int32_t* incx)
{
    lpf_float16_t* x = _x->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    lpf_float16_t r_val = LPF_GLOBAL(hnrm2, HNRM2)(&_n, x, &_incx);
    lpf_ffloat16_t r;
    r.value = r_val;
    return r;
}
