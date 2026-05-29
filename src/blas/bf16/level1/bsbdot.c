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

int16_t LPF_GLOBAL(bsbdot, BSBDOT)(int64_t* n, lpf_bfloat16_t* sb,
                                   lpf_bfloat16_t* sx, int64_t* incx,
                                   lpf_bfloat16_t* sy, int64_t* incy)
{

    int64_t i__1, i__2;
    lpf_bfloat16_t ret_val;

    int64_t i__, ns, kx, ky;
    float dbdot;

    --sy;
    --sx;

    dbdot = *sb;
    if (*n <= 0)
    {
        ret_val = dbdot;
        RETURN_BF16(ret_val);
    }
    if (*incx == *incy && *incx > 0)
    {

        ns = *n * *incx;
        i__1 = ns;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            dbdot += (float)sx[i__] * (float)sy[i__];
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
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__)
        {
            dbdot += (float)sx[kx] * (float)sy[ky];
            kx += *incx;
            ky += *incy;
        }
    }
    ret_val = dbdot;

    RETURN_BF16(ret_val);
}

#include <ISO_Fortran_binding.h>

lpf_fbfloat16_t
lpf_blas_bsbdot_fortran_dyn_rank_64(int64_t* n, lpf_fbfloat16_t* sb,
                                    CFI_cdesc_t* _sx, int64_t* incx,
                                    CFI_cdesc_t* _sy, int64_t* incy)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    lpf_bfloat16_t* sy = _sy->base_addr;
    lpf_fbfloat16_t r;
    r.value =
        LPF_GLOBAL(bsbdot, BSBDOT)(n, (lpf_bfloat16_t*)sb, sx, incx, sy, incy);
    return r;
}

lpf_fbfloat16_t
lpf_blas_bsbdot_fortran_dyn_rank_32(int32_t* n, lpf_fbfloat16_t* sb,
                                    CFI_cdesc_t* _sx, int32_t* incx,
                                    CFI_cdesc_t* _sy, int32_t* incy)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    lpf_bfloat16_t* sy = _sy->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    lpf_fbfloat16_t r;
    r.value = LPF_GLOBAL(bsbdot, BSBDOT)(&_n, (lpf_bfloat16_t*)sb, sx, &_incx,
                                         sy, &_incy);
    return r;
}
