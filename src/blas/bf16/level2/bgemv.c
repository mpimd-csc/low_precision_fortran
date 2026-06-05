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
 * @brief Bfloat16 General Matrix-Vector Multiplication (BGEMV).
 *
 * Performs the operation: y := alpha * A * x + beta * y
 *
 * @param[in] trans Character specifying the transpose. 'N' for no transpose, 'T' for transpose.
 * @param[in] m Number of rows of matrix A.
 * @param[in] n Number of columns of matrix A.
 * @param[in] alpha Scalar multiplier for the matrix-vector product.
 * @param[in] a Matrix A.
 * @param[in] lda Leading dimension of matrix A.
 * @param[in] x Vector X.
 * @param[in] incx Increment for the elements of x.
 * @param[in] beta Scalar multiplier for vector y.
 * @param[in,out] y Vector Y.
 * @param[in] incy Increment for the elements of y.
 */
void LPF_GLOBAL(bgemv, BGEMV)(char* trans, int64_t* m, int64_t* n,
                              lpf_bfloat16_t* alpha, lpf_bfloat16_t* a,
                              int64_t* lda, lpf_bfloat16_t* x, int64_t* incx,
                              lpf_bfloat16_t* beta, lpf_bfloat16_t* y,
                              int64_t* incy, lpf_fortran_strlen_t trans_len)
{

    int64_t a_dim1, a_offset, i__1, i__2;

    int64_t i__, j, ix, iy, jx, jy, kx, ky, info;
    lpf_bfloat16_t temp;
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
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BGEMV ", &infox,
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

void lpf_blas_bgemv_fortran_dyn_rank_64(char* trans, int64_t* m, int64_t* n,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int64_t* lda, CFI_cdesc_t* _x,
                                        int64_t* incx, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _y, int64_t* incy)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    LPF_GLOBAL(bgemv, BGEMV)(trans, m, n, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)a, lda, (lpf_bfloat16_t*)x, incx,
                             (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)y, incy,
                             1);
}

void lpf_blas_bgemv_fortran_dyn_rank_32(char* trans, int32_t* m, int32_t* n,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int32_t* lda, CFI_cdesc_t* _x,
                                        int32_t* incx, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _y, int32_t* incy)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    int64_t _m = *m;
    int64_t _n = *n;
    int64_t _lda = *lda;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(bgemv, BGEMV)(trans, &_m, &_n, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)a, &_lda, (lpf_bfloat16_t*)x,
                             &_incx, (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)y,
                             &_incy, 1);
}
