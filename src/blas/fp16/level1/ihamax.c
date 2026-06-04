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

int64_t LPF_GLOBAL(ihamax, IHAMAX)(int64_t* n, lpf_float16_t* sx, int64_t* incx)
{

    int64_t ret_val, i__1;
    lpf_float16_t r__1;

    int64_t i__, ix;
    lpf_float16_t smax;

    --sx;

    ret_val = 0;
    if (*n < 1 || *incx <= 0)
    {
        return ret_val;
    }
    ret_val = 1;
    if (*n == 1)
    {
        return ret_val;
    }
    if (*incx == 1)
    {

        smax = abs_fp16(sx[1]);
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__)
        {
            r__1 = sx[i__];
            if (abs_fp16(r__1) > smax)
            {
                ret_val = i__;
                smax = abs_fp16(r__1);
            }
        }
    }
    else
    {

        ix = 1;
        smax = abs_fp16(sx[1]);
        ix += *incx;
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__)
        {
            r__1 = abs_fp16(sx[ix]);
            if (r__1 > smax)
            {
                ret_val = i__;
                smax = r__1;
            }
            ix += *incx;
        }
    }
    return ret_val;
}

#include <ISO_Fortran_binding.h>

int64_t lpf_blas_ihamax_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _sx,
                                                   int64_t* incx)
{
    lpf_float16_t* sx = _sx->base_addr;
    int64_t res = LPF_GLOBAL(ihamax, IHAMAX)(n, sx, incx);
    return (int64_t)res;
}

int32_t lpf_blas_ihamax_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _sx,
                                                   int32_t* incx)
{
    lpf_float16_t* sx = _sx->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t res = LPF_GLOBAL(ihamax, IHAMAX)(&_n, sx, &_incx);
    return (int32_t)res;
}
