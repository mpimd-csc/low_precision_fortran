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

void LPF_GLOBAL(hswap, HSWAP)(int64_t* n, lpf_float16_t* sx, int64_t* incx,
                              lpf_float16_t* sy, int64_t* incy)
{

    int64_t i__1;

    int64_t i__, m, ix, iy, mp1;
    lpf_float16_t stemp;

    --sy;
    --sx;

    if (*n <= 0)
    {
        return;
    }
    if (*incx == 1 && *incy == 1)
    {

        m = *n % 3;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                stemp = sx[i__];
                sx[i__] = sy[i__];
                sy[i__] = stemp;
            }
            if (*n < 3)
            {
                return;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 3)
        {
            stemp = sx[i__];
            sx[i__] = sy[i__];
            sy[i__] = stemp;
            stemp = sx[i__ + 1];
            sx[i__ + 1] = sy[i__ + 1];
            sy[i__ + 1] = stemp;
            stemp = sx[i__ + 2];
            sx[i__ + 2] = sy[i__ + 2];
            sy[i__ + 2] = stemp;
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
            stemp = sx[ix];
            sx[ix] = sy[iy];
            sy[iy] = stemp;
            ix += *incx;
            iy += *incy;
        }
    }
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hswap_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _sx,
                                        int64_t* incx, CFI_cdesc_t* _sy,
                                        int64_t* incy)
{
    lpf_float16_t* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    LPF_GLOBAL(hswap, HSWAP)(n, sx, incx, sy, incy);
}

void lpf_blas_hswap_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _sx,
                                        int32_t* incx, CFI_cdesc_t* _sy,
                                        int32_t* incy)
{
    lpf_float16_t* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(hswap, HSWAP)(&_n, sx, &_incx, sy, &_incy);
}
