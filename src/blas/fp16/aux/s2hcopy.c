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

void LPF_GLOBAL(s2hcopy, S2HCOPY)(int64_t* n, float* sx, int64_t* incx,
                                  lpf_float16_t* sy, int64_t* incy)
{

    int64_t i__1;

    int64_t i__, m, ix, iy, mp1;

    --sy;
    --sx;

    if (*n <= 0)
    {
        return;
    }
    if (*incx == 1 && *incy == 1)
    {

        m = *n % 7;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                sy[i__] = (lpf_float16_t)sx[i__];
            }
            if (*n < 7)
            {
                return;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 7)
        {
            sy[i__] = (lpf_float16_t)sx[i__];
            sy[i__ + 1] = (lpf_float16_t)sx[i__ + 1];
            sy[i__ + 2] = (lpf_float16_t)sx[i__ + 2];
            sy[i__ + 3] = (lpf_float16_t)sx[i__ + 3];
            sy[i__ + 4] = (lpf_float16_t)sx[i__ + 4];
            sy[i__ + 5] = (lpf_float16_t)sx[i__ + 5];
            sy[i__ + 6] = (lpf_float16_t)sx[i__ + 6];
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
            sy[iy] = (lpf_float16_t)sx[ix];
            ix += *incx;
            iy += *incy;
        }
    }
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_s2hcopy_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _sx,
                                          int64_t* incx, CFI_cdesc_t* _sy,
                                          int64_t* incy)
{
    float* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    LPF_GLOBAL(s2hcopy, S2HCOPY)(n, sx, incx, (lpf_float16_t*)sy, incy);
}

void lpf_blas_s2hcopy_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _sx,
                                          int32_t* incx, CFI_cdesc_t* _sy,
                                          int32_t* incy)
{
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;

    float* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    LPF_GLOBAL(s2hcopy, S2HCOPY)(&_n, sx, &_incx, (lpf_float16_t*)sy, &_incy);
}
