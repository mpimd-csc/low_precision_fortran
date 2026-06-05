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
 * @brief Bfloat16 Symmetric Positive Definite Banded Matrix-Vector Multiplication (BSPMV).
 *
 * Performs the operation: y := alpha * A * x + beta * y
 * where A is a symmetric positive definite banded matrix.
 *
 * @param[in] uplo Character specifying the triangular part of the matrix A to be used: 'U' for upper, 'L' for lower.
 * @param[in] n Order of matrix A.
 * @param[in] alpha Scalar multiplier.
 * @param[in] ap The symmetric positive definite banded matrix A in packed storage.
 * @param[in] x Vector x.
 * @param[in] incx Increment for the elements of x.
 * @param[in] beta Scalar multiplier for vector y.
 * @param[in,out] y Vector y.
 * @param[in] incy Increment for the elements of y.
 */
void LPF_GLOBAL(bspmv, BSPMV)(char* uplo, int64_t* n, lpf_bfloat16_t* alpha,
                              lpf_bfloat16_t* ap, lpf_bfloat16_t* x,
                              int64_t* incx, lpf_bfloat16_t* beta,
                              lpf_bfloat16_t* y, int64_t* incy,
                              lpf_fortran_strlen_t uplo_len)
{

    int64_t i__1, i__2;

    int64_t i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    lpf_bfloat16_t temp1, temp2;

    --y;
    --x;
    --ap;
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
        info = 6;
    }
    else if (*incy == 0)
    {
        info = 9;
    }
    if (info != 0)
    {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSPMV ", &infox,
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
    kk = 1;
    if (strncasecmp(uplo, "U", 1) == 0)
    {

        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                k = kk;
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    y[i__] += temp1 * ap[k];
                    temp2 += ap[k] * x[i__];
                    ++k;
                }
                y[j] = y[j] + temp1 * ap[kk + j - 1] + *alpha * temp2;
                kk += j;
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
                i__2 = kk + j - 2;
                for (k = kk; k <= i__2; ++k)
                {
                    y[iy] += temp1 * ap[k];
                    temp2 += ap[k] * x[ix];
                    ix += *incx;
                    iy += *incy;
                }
                y[jy] = y[jy] + temp1 * ap[kk + j - 1] + *alpha * temp2;
                jx += *incx;
                jy += *incy;
                kk += j;
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
                y[j] += temp1 * ap[kk];
                k = kk + 1;
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__)
                {
                    y[i__] += temp1 * ap[k];
                    temp2 += ap[k] * x[i__];
                    ++k;
                }
                y[j] += *alpha * temp2;
                kk += *n - j + 1;
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
                y[jy] += temp1 * ap[kk];
                ix = jx;
                iy = jy;
                i__2 = kk + *n - j;
                for (k = kk + 1; k <= i__2; ++k)
                {
                    ix += *incx;
                    iy += *incy;
                    y[iy] += temp1 * ap[k];
                    temp2 += ap[k] * x[ix];
                }
                y[jy] += *alpha * temp2;
                jx += *incx;
                jy += *incy;
                kk += *n - j + 1;
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_bspmv_fortran_dyn_rank_64(char* uplo, int64_t* n,
                                        lpf_fbfloat16_t* alpha,
                                        CFI_cdesc_t* _ap, CFI_cdesc_t* _x,
                                        int64_t* incx, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _y, int64_t* incy)
{
    lpf_bfloat16_t* ap = _ap->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    LPF_GLOBAL(bspmv, BSPMV)(uplo, n, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)ap, (lpf_bfloat16_t*)x, incx,
                             (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)y, incy,
                             1);
}

void lpf_blas_bspmv_fortran_dyn_rank_32(char* uplo, int32_t* n,
                                        lpf_fbfloat16_t* alpha,
                                        CFI_cdesc_t* _ap, CFI_cdesc_t* _x,
                                        int32_t* incx, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _y, int32_t* incy)
{
    lpf_bfloat16_t* ap = _ap->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* y = _y->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(bspmv, BSPMV)(uplo, &_n, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)ap, (lpf_bfloat16_t*)x, &_incx,
                             (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)y, &_incy,
                             1);
}
