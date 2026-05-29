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

void LPF_GLOBAL(hsyr2, HSYR2)(char* uplo, int64_t* n, lpf_float16_t* alpha,
                              lpf_float16_t* x, int64_t* incx, lpf_float16_t* y,
                              int64_t* incy, lpf_float16_t* a, int64_t* lda,
                              lpf_fortran_strlen_t uplo_len)
{

    int64_t a_dim1, a_offset, i__1, i__2;

    int64_t i__, j, ix, iy, jx, jy, kx, ky, info;
    lpf_float16_t temp1, temp2;

    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    (void)uplo_len;

    info = 0;
    if (!(strncasecmp(uplo, "U", 1) == 0) && !(strncasecmp(uplo, "L", 1) == 0))
    {
        info = 1;
    }
    else if (*n < 0)
    {
        info = 2;
    }
    else if (*incx == 0)
    {
        info = 5;
    }
    else if (*incy == 0)
    {
        info = 7;
    }
    else if (*lda < LPF_MAX(1, *n))
    {
        info = 9;
    }
    if (info != 0)
    {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HSYR2 ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*n == 0 || *alpha == 0.f)
    {
        return;
    }

    if (*incx != 1 || *incy != 1)
    {
        if (*incx > 0)
        {
            kx = 1;
        }
        else
        {
            kx = 1 - (*n - 1) * *incx;
        }
        if (*incy > 0)
        {
            ky = 1;
        }
        else
        {
            ky = 1 - (*n - 1) * *incy;
        }
        jx = kx;
        jy = ky;
    }

    if (strncasecmp(uplo, "U", 1) == 0)
    {

        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (x[j] != 0.f || y[j] != 0.f)
                {
                    temp1 = *alpha * y[j];
                    temp2 = *alpha * x[j];
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        a[i__ + j * a_dim1] = a[i__ + j * a_dim1] +
                                              x[i__] * temp1 + y[i__] * temp2;
                    }
                }
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (x[jx] != 0.f || y[jy] != 0.f)
                {
                    temp1 = *alpha * y[jy];
                    temp2 = *alpha * x[jx];
                    ix = kx;
                    iy = ky;
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        a[i__ + j * a_dim1] =
                            a[i__ + j * a_dim1] + x[ix] * temp1 + y[iy] * temp2;
                        ix += *incx;
                        iy += *incy;
                    }
                }
                jx += *incx;
                jy += *incy;
            }
        }
    }
    else
    {

        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (x[j] != 0.f || y[j] != 0.f)
                {
                    temp1 = *alpha * y[j];
                    temp2 = *alpha * x[j];
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        a[i__ + j * a_dim1] = a[i__ + j * a_dim1] +
                                              x[i__] * temp1 + y[i__] * temp2;
                    }
                }
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (x[jx] != 0.f || y[jy] != 0.f)
                {
                    temp1 = *alpha * y[jy];
                    temp2 = *alpha * x[jx];
                    ix = jx;
                    iy = jy;
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        a[i__ + j * a_dim1] =
                            a[i__ + j * a_dim1] + x[ix] * temp1 + y[iy] * temp2;
                        ix += *incx;
                        iy += *incy;
                    }
                }
                jx += *incx;
                jy += *incy;
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hsyr2_fortran_dyn_rank_64(char* uplo, int64_t* n,
                                        lpf_ffloat16_t* alpha, CFI_cdesc_t* x,
                                        int64_t* incx, CFI_cdesc_t* y,
                                        int64_t* incy, CFI_cdesc_t* a,
                                        int64_t* lda)
{
    lpf_float16_t* x_ptr = x->base_addr;
    lpf_float16_t* y_ptr = y->base_addr;
    lpf_float16_t* a_ptr = a->base_addr;

    LPF_GLOBAL(hsyr2, HSYR2)(uplo, n, (lpf_float16_t*)alpha,
                             (lpf_float16_t*)x_ptr, incx, (lpf_float16_t*)y_ptr,
                             incy, (lpf_float16_t*)a_ptr, lda, 1);
}

void lpf_blas_hsyr2_fortran_dyn_rank_32(char* uplo, int32_t* n,
                                        lpf_ffloat16_t* alpha, CFI_cdesc_t* x,
                                        int32_t* incx, CFI_cdesc_t* y,
                                        int32_t* incy, CFI_cdesc_t* a,
                                        int32_t* lda)
{
    lpf_float16_t* x_ptr = x->base_addr;
    lpf_float16_t* y_ptr = y->base_addr;
    lpf_float16_t* a_ptr = a->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    int64_t _lda = *lda;

    LPF_GLOBAL(hsyr2, HSYR2)(
        uplo, &_n, (lpf_float16_t*)alpha, (lpf_float16_t*)x_ptr, &_incx,
        (lpf_float16_t*)y_ptr, &_incy, (lpf_float16_t*)a_ptr, &_lda, 1);
}
