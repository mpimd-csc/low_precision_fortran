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

#include <string.h>

void LPF_GLOBAL(bsbmv, BSBMV)(char* uplo, int64_t* n, int64_t* k,
                              lpf_bfloat16_t* alpha, lpf_bfloat16_t* a,
                              int64_t* lda, lpf_bfloat16_t* x, int64_t* incx,
                              lpf_bfloat16_t* beta, lpf_bfloat16_t* y,
                              int64_t* incy, lpf_fortran_strlen_t uplo_len)
{

    int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int64_t i__, j, l, ix, iy, jx, jy, kx, ky, info;
    lpf_bfloat16_t temp1, temp2;
    int64_t kplus1;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

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
    else if (*k < 0)
    {
        info = 3;
    }
    else if (*lda < *k + 1)
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
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSBMV ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*n == 0 || (*alpha == 0.f && *beta == 1.f))
    {
        return;
    }

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

    if (*beta != 1.f)
    {
        if (*incy == 1)
        {
            if (*beta == 0.f)
            {
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    y[i__] = 0.f;
                }
            }
            else
            {
                i__1 = *n;
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
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__)
                {
                    y[iy] = 0.f;
                    iy += *incy;
                }
            }
            else
            {
                i__1 = *n;
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
    if (strncasecmp(uplo, "U", 1) == 0)
    {

        kplus1 = *k + 1;
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                l = kplus1 - j;

                i__2 = 1, i__3 = j - *k;
                i__4 = j - 1;
                for (i__ = LPF_MAX(i__2, i__3); i__ <= i__4; ++i__)
                {
                    y[i__] += temp1 * a[l + i__ + j * a_dim1];
                    temp2 += a[l + i__ + j * a_dim1] * x[i__];
                }
                y[j] = y[j] + temp1 * a[kplus1 + j * a_dim1] + *alpha * temp2;
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp1 = *alpha * x[jx];
                temp2 = 0.f;
                ix = kx;
                iy = ky;
                l = kplus1 - j;

                i__4 = 1, i__2 = j - *k;
                i__3 = j - 1;
                for (i__ = LPF_MAX(i__4, i__2); i__ <= i__3; ++i__)
                {
                    y[iy] += temp1 * a[l + i__ + j * a_dim1];
                    temp2 += a[l + i__ + j * a_dim1] * x[ix];
                    ix += *incx;
                    iy += *incy;
                }
                y[jy] = y[jy] + temp1 * a[kplus1 + j * a_dim1] + *alpha * temp2;
                jx += *incx;
                jy += *incy;
                if (j > *k)
                {
                    kx += *incx;
                    ky += *incy;
                }
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
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                y[j] += temp1 * a[j * a_dim1 + 1];
                l = 1 - j;

                i__4 = *n, i__2 = j + *k;
                i__3 = LPF_MIN(i__4, i__2);
                for (i__ = j + 1; i__ <= i__3; ++i__)
                {
                    y[i__] += temp1 * a[l + i__ + j * a_dim1];
                    temp2 += a[l + i__ + j * a_dim1] * x[i__];
                }
                y[j] += *alpha * temp2;
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp1 = *alpha * x[jx];
                temp2 = 0.f;
                y[jy] += temp1 * a[j * a_dim1 + 1];
                l = 1 - j;
                ix = jx;
                iy = jy;

                i__4 = *n, i__2 = j + *k;
                i__3 = LPF_MIN(i__4, i__2);
                for (i__ = j + 1; i__ <= i__3; ++i__)
                {
                    ix += *incx;
                    iy += *incy;
                    y[iy] += temp1 * a[l + i__ + j * a_dim1];
                    temp2 += a[l + i__ + j * a_dim1] * x[ix];
                }
                y[jy] += *alpha * temp2;
                jx += *incx;
                jy += *incy;
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_bsbmv_fortran_dyn_rank_64(char* uplo, int64_t* n, int64_t* k,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int64_t* lda, CFI_cdesc_t* _x,
                                        int64_t* incx, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _y, int64_t* incy)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    LPF_GLOBAL(bsbmv, BSBMV)(uplo, n, k, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)a, lda, (lpf_bfloat16_t*)x, incx,
                             (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)y, incy,
                             1);
}

void lpf_blas_bsbmv_fortran_dyn_rank_32(char* uplo, int32_t* n, int32_t* k,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int32_t* lda, CFI_cdesc_t* _x,
                                        int32_t* incx, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _y, int32_t* incy)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    int64_t _n = *n;
    int64_t _k = *k;
    int64_t _lda = *lda;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(bsbmv, BSBMV)(uplo, &_n, &_k, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)a, &_lda, (lpf_bfloat16_t*)x,
                             &_incx, (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)y,
                             &_incy, 1);
}
