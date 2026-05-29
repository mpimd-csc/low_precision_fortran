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

void LPF_GLOBAL(hrotm, HROTM)(int64_t* n, lpf_float16_t* sx, int64_t* incx,
                              lpf_float16_t* sy, int64_t* incy,
                              lpf_float16_t* sparam)
{

    const lpf_float16_t zero = 0.f;
    const lpf_float16_t two = 2.f;

    int64_t i__1, i__2;

    int64_t i__;
    lpf_float16_t w, z__;
    int64_t kx, ky;
    lpf_float16_t sh11, sh12, sh21, sh22, sflag;
    int64_t nsteps;

    --sparam;
    --sy;
    --sx;

    sflag = sparam[1];
    if (*n <= 0 || sflag + two == zero)
    {
        return;
    }
    if (*incx == *incy && *incx > 0)
    {

        nsteps = *n * *incx;
        if (sflag < zero)
        {
            sh11 = sparam[2];
            sh12 = sparam[4];
            sh21 = sparam[3];
            sh22 = sparam[5];
            i__1 = nsteps;
            i__2 = *incx;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
            {
                w = sx[i__];
                z__ = sy[i__];
                sx[i__] = w * sh11 + z__ * sh12;
                sy[i__] = w * sh21 + z__ * sh22;
            }
        }
        else if (sflag == zero)
        {
            sh12 = sparam[4];
            sh21 = sparam[3];
            i__2 = nsteps;
            i__1 = *incx;
            for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
            {
                w = sx[i__];
                z__ = sy[i__];
                sx[i__] = w + z__ * sh12;
                sy[i__] = w * sh21 + z__;
            }
        }
        else
        {
            sh11 = sparam[2];
            sh22 = sparam[5];
            i__1 = nsteps;
            i__2 = *incx;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
            {
                w = sx[i__];
                z__ = sy[i__];
                sx[i__] = w * sh11 + z__;
                sy[i__] = -w + sh22 * z__;
            }
        }
    }
    else
    {
        kx = 1;
        ky = 1;
        if (*incx < 0)
        {
            kx = (1 - *n) * *incx + 1;
        }
        if (*incy < 0)
        {
            ky = (1 - *n) * *incy + 1;
        }

        if (sflag < zero)
        {
            sh11 = sparam[2];
            sh12 = sparam[4];
            sh21 = sparam[3];
            sh22 = sparam[5];
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                w = sx[kx];
                z__ = sy[ky];
                sx[kx] = w * sh11 + z__ * sh12;
                sy[ky] = w * sh21 + z__ * sh22;
                kx += *incx;
                ky += *incy;
            }
        }
        else if (sflag == zero)
        {
            sh12 = sparam[4];
            sh21 = sparam[3];
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                w = sx[kx];
                z__ = sy[ky];
                sx[kx] = w + z__ * sh12;
                sy[ky] = w * sh21 + z__;
                kx += *incx;
                ky += *incy;
            }
        }
        else
        {
            sh11 = sparam[2];
            sh22 = sparam[5];
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                w = sx[kx];
                z__ = sy[ky];
                sx[kx] = w * sh11 + z__;
                sy[ky] = -w + sh22 * z__;
                kx += *incx;
                ky += *incy;
            }
        }
    }
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hrotm_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _sx,
                                        int64_t* incx, CFI_cdesc_t* _sy,
                                        int64_t* incy, CFI_cdesc_t* _sparam)
{
    lpf_float16_t* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    lpf_float16_t* sparam = _sparam->base_addr;
    LPF_GLOBAL(hrotm, HROTM)(n, sx, incx, sy, incy, sparam);
}

void lpf_blas_hrotm_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _sx,
                                        int32_t* incx, CFI_cdesc_t* _sy,
                                        int32_t* incy, CFI_cdesc_t* _sparam)
{
    lpf_float16_t* sx = _sx->base_addr;
    lpf_float16_t* sy = _sy->base_addr;
    lpf_float16_t* sparam = _sparam->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(hrotm, HROTM)(&_n, sx, &_incx, sy, &_incy, sparam);
}
