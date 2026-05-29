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
#include <string.h>

void LPF_GLOBAL(hgemv, HGEMV)(char* trans, int64_t* m, int64_t* n,
                              lpf_float16_t* alpha, lpf_float16_t* a,
                              int64_t* lda, lpf_float16_t* x, int64_t* incx,
                              lpf_float16_t* beta, lpf_float16_t* y,
                              int64_t* incy, lpf_fortran_strlen_t trans_len)
{

    int64_t a_dim1, a_offset, i__1, i__2;

    int64_t i__, j, ix, iy, jx, jy, kx, ky, info;
    lpf_float16_t temp;
    int64_t lenx, leny;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    (void)trans_len;

    info = 0;
    if (!(strncasecmp(trans, "N", 1) == 0) &&
        !(strncasecmp(trans, "T", 1) == 0) &&
        !(strncasecmp(trans, "C", 1) == 0))
    {
        info = 1;
    }
    else if (*m < 0)
    {
        info = 2;
    }
    else if (*n < 0)
    {
        info = 3;
    }
    else if (*lda < LPF_MAX(1, *m))
    {
        info = 6;
    }
    else if (*incx == 0)
    {
        info = 8;
    }
    else if (*incy == 0)
    {
        info = 11;
    }
    if (info != 0)
    {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HGEMV ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*m == 0 || *n == 0 || (*alpha == 0.f && *beta == 1.f))
    {
        return;
    }

    if (strncasecmp(trans, "N", 1) == 0)
    {
        lenx = *n;
        leny = *m;
    }
    else
    {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (leny - 1) * *incy;
    }

    if (*beta != 1.f)
    {
        if (*incy == 1)
        {
            if (*beta == 0.f)
            {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    y[i__] = 0.f;
                }
            }
            else
            {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    y[i__] = *beta * y[i__];
                }
            }
        }
        else
        {
            iy = ky;
            if (*beta == 0.f)
            {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    y[iy] = 0.f;
                    iy += *incy;
                }
            }
            else
            {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    y[iy] = *beta * y[iy];
                    iy += *incy;
                }
            }
        }
    }
    if (*alpha == 0.f)
    {
        return;
    }
    if (strncasecmp(trans, "N", 1) == 0)
    {

        jx = kx;
        if (*incy == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp = *alpha * x[jx];
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    y[i__] += temp * a[i__ + j * a_dim1];
                }
                jx += *incx;
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp = *alpha * x[jx];
                iy = ky;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    y[iy] += temp * a[i__ + j * a_dim1];
                    iy += *incy;
                }
                jx += *incx;
            }
        }
    }
    else
    {

        jy = ky;
        if (*incx == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp = 0.f;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    temp += a[i__ + j * a_dim1] * x[i__];
                }
                y[jy] += *alpha * temp;
                jy += *incy;
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp = 0.f;
                ix = kx;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    temp += a[i__ + j * a_dim1] * x[ix];
                    ix += *incx;
                }
                y[jy] += *alpha * temp;
                jy += *incy;
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hgemv_fortran_dyn_rank_64(char* trans, int64_t* m, int64_t* n,
                                        lpf_ffloat16_t* alpha, CFI_cdesc_t* a,
                                        int64_t* lda, CFI_cdesc_t* x,
                                        int64_t* incx, lpf_ffloat16_t* beta,
                                        CFI_cdesc_t* y, int64_t* incy)
{
    lpf_float16_t* a_ptr = a->base_addr;
    lpf_float16_t* x_ptr = x->base_addr;
    lpf_float16_t* y_ptr = y->base_addr;

    LPF_GLOBAL(hgemv, HGEMV)(trans, m, n, (lpf_float16_t*)alpha,
                             (lpf_float16_t*)a_ptr, lda, (lpf_float16_t*)x_ptr,
                             incx, (lpf_float16_t*)beta, (lpf_float16_t*)y_ptr,
                             incy, 1);
}

void lpf_blas_hgemv_fortran_dyn_rank_32(char* trans, int32_t* m, int32_t* n,
                                        lpf_ffloat16_t* alpha, CFI_cdesc_t* a,
                                        int32_t* lda, CFI_cdesc_t* x,
                                        int32_t* incx, lpf_ffloat16_t* beta,
                                        CFI_cdesc_t* y, int32_t* incy)
{
    lpf_float16_t* a_ptr = a->base_addr;
    lpf_float16_t* x_ptr = x->base_addr;
    lpf_float16_t* y_ptr = y->base_addr;

    int64_t _m = *m;
    int64_t _n = *n;
    int64_t _lda = *lda;
    int64_t _incx = *incx;
    int64_t _incy = *incy;

    LPF_GLOBAL(hgemv, HGEMV)(
        trans, &_m, &_n, (lpf_float16_t*)alpha, (lpf_float16_t*)a_ptr, &_lda,
        (lpf_float16_t*)x_ptr, &_incx, (lpf_float16_t*)beta,
        (lpf_float16_t*)y_ptr, &_incy, 1);
}
