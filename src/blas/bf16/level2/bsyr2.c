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

/**
 * @brief Bfloat16 Symmetric Rank-2 Update (BSYR2).
 *
 * Performs the operation: A := A + alpha * (x * y^T + y * x^T)
 *
 * @param[in] uplo Character specifying the upper or lower triangle of A.
 * @param[in] n Order of matrix A.
 * @param[in] alpha Scalar multiplier.
 * @param[in] x Vector X.
 * @param[in] incx Increment for the elements of x.
 * @param[in] y Vector Y.
 * @param[in] incy Increment for the elements of y.
 * @param[in,out] a Symmetric matrix A.
 * @param[in] lda Leading dimension of matrix A.
 */
void LPF_GLOBAL(bsyr2, BSYR2)(char* uplo, int64_t* n, lpf_bfloat16_t* alpha,
                              lpf_bfloat16_t* x, int64_t* incx,
                              lpf_bfloat16_t* y, int64_t* incy,
                              lpf_bfloat16_t* a, int64_t* lda,
                              lpf_fortran_strlen_t uplo_len)
{

    int64_t a_dim1, a_offset, i__1, i__2;

    int64_t i__, j, ix, iy, jx, jy, kx, ky, info;
    lpf_bfloat16_t temp1, temp2;

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
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSYR2 ", &infox,
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

void lpf_blas_bsyr2_fortran_dyn_rank_64(char* uplo, int64_t* n,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _x,
                                        int64_t* incx, CFI_cdesc_t* _y,
                                        int64_t* incy, CFI_cdesc_t* _a,
                                        int64_t* lda)
{
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    lpf_bfloat16_t* a = _a->base_addr;
    LPF_GLOBAL(bsyr2, BSYR2)(uplo, n, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)x, incx, (lpf_bfloat16_t*)y, incy,
                             (lpf_bfloat16_t*)a, lda, 1);
}

void lpf_blas_bsyr2_fortran_dyn_rank_32(char* uplo, int32_t* n,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _x,
                                        int32_t* incx, CFI_cdesc_t* _y,
                                        int32_t* incy, CFI_cdesc_t* _a,
                                        int32_t* lda)
{
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    lpf_bfloat16_t* a = _a->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    int64_t _lda = *lda;
    LPF_GLOBAL(bsyr2, BSYR2)(uplo, &_n, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)x, &_incx, (lpf_bfloat16_t*)y,
                             &_incy, (lpf_bfloat16_t*)a, &_lda, 1);
}
