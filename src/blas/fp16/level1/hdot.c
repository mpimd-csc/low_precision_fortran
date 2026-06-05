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
#include "lpf_internal.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

/**
 * @brief Half-precision Dot product.
 *
 * Performs the operation: result := x^T * y
 *
 * @param[in] n Number of elements in vectors x and y.
 * @param[in] sx Vector x.
 * @param[in] incx Increment for the elements of x.
 * @param[in] sy Vector y.
 * @param[in] incy Increment for the elements of y.
 * @return The result of the dot product.
 */
int16_t LPF_GLOBAL(hdot, HDOT)(int64_t* n, lpf_float16_t* sx, int64_t* incx,
                               lpf_float16_t* sy, int64_t* incy)
{

    int64_t i__1;
    lpf_float16_t ret_val;

    int64_t i__, m, ix, iy, mp1;
    lpf_float16_t stemp;

    --sy;
    --sx;

    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0)
    {
        RETURN_FP16(ret_val);
    }
    if (*incx == 1 && *incy == 1)
    {

        m = *n % 5;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                stemp += sx[i__] * sy[i__];
            }
            if (*n < 5)
            {
                ret_val = stemp;
                RETURN_FP16(ret_val);
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 5)
        {
            stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] +
                    sx[i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] +
                    sx[i__ + 4] * sy[i__ + 4];
        }
    }
    else
    {

        ix = 1;
        iy = 1;
        if (*incx < 0)
        {
            ix = (-(*n) + 1) * *incx + 1;
        }
        if (*incy < 0)
        {
            iy = (-(*n) + 1) * *incy + 1;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            stemp += sx[ix] * sy[iy];
            ix += *incx;
            iy += *incy;
        }
    }
    ret_val = stemp;

    RETURN_FP16(ret_val);
}

#include <ISO_Fortran_binding.h>

lpf_ffloat16_t lpf_blas_hdot_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _sx,
                                                 int64_t* incx,
                                                 CFI_cdesc_t* _sy,
                                                 int64_t* incy)
{
    lpf_float16_t* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    lpf_ffloat16_t r;
    r.value = LPF_GLOBAL(hdot, HDOT)(n, sx, incx, sy, incy);
    return r;
}

lpf_ffloat16_t lpf_blas_hdot_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _sx,
                                                 int32_t* incx,
                                                 CFI_cdesc_t* _sy,
                                                 int32_t* incy)
{
    lpf_float16_t* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    lpf_ffloat16_t r;
    r.value = LPF_GLOBAL(hdot, HDOT)(&_n, sx, &_incx, sy, &_incy);
    return r;
}
