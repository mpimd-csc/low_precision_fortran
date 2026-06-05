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
 * @brief Bfloat16 AXPY operation (BAXPY).
 *
 * Performs the operation: y := alpha * x + y
 *
 * @param[in] n Number of elements in vectors x and y.
 * @param[in] sa Scalar multiplier alpha.
 * @param[in] sx Vector x.
 * @param[in] incx Increment for the elements of x.
 * @param[in,out] sy Vector y.
 * @param[in] incy Increment for the elements of y.
 */
void LPF_GLOBAL(baxpy, BAXPY)(int64_t* n, lpf_bfloat16_t* sa,
                              lpf_bfloat16_t* sx, int64_t* incx,
                              lpf_bfloat16_t* sy, int64_t* incy)
{

    int64_t i__1;

    int64_t i__, m, ix, iy, mp1;

    --sy;
    --sx;

    if (*n <= 0)
    {
        return;
    }
    if (*sa == 0.f)
    {
        return;
    }
    if (*incx == 1 && *incy == 1)
    {

        m = *n % 4;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                sy[i__] += *sa * sx[i__];
            }
        }
        if (*n < 4)
        {
            return;
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 4)
        {
            sy[i__] += *sa * sx[i__];
            sy[i__ + 1] += *sa * sx[i__ + 1];
            sy[i__ + 2] += *sa * sx[i__ + 2];
            sy[i__ + 3] += *sa * sx[i__ + 3];
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
            sy[iy] += *sa * sx[ix];
            ix += *incx;
            iy += *incy;
        }
    }
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_baxpy_fortran_dyn_rank_64(int64_t* n, lpf_fbfloat16_t* sa,
                                        CFI_cdesc_t* _sx, int64_t* incx,
                                        CFI_cdesc_t* _sy, int64_t* incy)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    lpf_bfloat16_t* sy = _sy->base_addr;
    LPF_GLOBAL(baxpy, BAXPY)(n, (lpf_bfloat16_t*)sa, (lpf_bfloat16_t*)sx, incx,
                             (lpf_bfloat16_t*)sy, incy);
}

void lpf_blas_baxpy_fortran_dyn_rank_32(int32_t* n, lpf_fbfloat16_t* sa,
                                        CFI_cdesc_t* _sx, int32_t* incx,
                                        CFI_cdesc_t* _sy, int32_t* incy)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    lpf_bfloat16_t* sy = _sy->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(baxpy, BAXPY)(&_n, (lpf_bfloat16_t*)sa, (lpf_bfloat16_t*)sx,
                             &_incx, (lpf_bfloat16_t*)sy, &_incy);
}
