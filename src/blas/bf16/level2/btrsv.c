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
 * @brief Bfloat16 Triangular System Solve (BTRSV).
 *
 * Solves a system of linear equations:
 * A * x = b
 * where A is a triangular matrix.
 *
 * @param[in] uplo Character specifying the triangular part of the matrix A to be used: 'U' for upper, 'L' for lower.
 * @param[in] trans Character specifying the transpose. 'N' for no transpose, 'T' for transpose.
 * @param[in] diag Character specifying whether the diagonal of A is unit: 'U' for unit, 'N' for non-unit.
 * @param[in] n Order of matrix A.
 * @param[in] a The triangular matrix A.
 * @param[in] lda Leading dimension of matrix A.
 * @param[in,out] x Vector X (input `b`, output `x`).
 * @param[in] incx Increment for the elements of x.
 */
void LPF_GLOBAL(btrsv, BTRSV)(char* uplo, char* trans, char* diag, int64_t* n,
                              lpf_bfloat16_t* a, int64_t* lda,
                              lpf_bfloat16_t* x, int64_t* incx,
                              lpf_fortran_strlen_t uplo_len,
                              lpf_fortran_strlen_t trans_len,
                              lpf_fortran_strlen_t diag_len)
{

    int64_t a_dim1, a_offset, i__1, i__2;

    int64_t i__, j, ix, jx, kx, info;
    lpf_bfloat16_t temp;
    lpf_logical_t nounit;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    (void)diag_len;
    (void)trans_len;
    (void)uplo_len;

    info = 0;
    if (!(strncasecmp(uplo, "U", 1) == 0) && !(strncasecmp(uplo, "L", 1) == 0))
    {
        info = 1;
    }
    else if (!(strncasecmp(trans, "N", 1) == 0) &&
             !(strncasecmp(trans, "T", 1) == 0) &&
             !(strncasecmp(trans, "C", 1) == 0))
    {
        info = 2;
    }
    else if (!(strncasecmp(diag, "U", 1) == 0) &&
             !(strncasecmp(diag, "N", 1) == 0))
    {
        info = 3;
    }
    else if (*n < 0)
    {
        info = 4;
    }
    else if (*lda < LPF_MAX(1, *n))
    {
        info = 6;
    }
    else if (*incx == 0)
    {
        info = 8;
    }
    if (info != 0)
    {
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BTRSV ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*n == 0)
    {
        return;
    }

    nounit = strncasecmp(diag, "N", 1) == 0;

    if (*incx <= 0)
    {
        kx = 1 - (*n - 1) * *incx;
    }
    else if (*incx != 1)
    {
        kx = 1;
    }

    if (strncasecmp(trans, "N", 1) == 0)
    {

        if (strncasecmp(uplo, "U", 1) == 0)
        {
            if (*incx == 1)
            {
                for (j = *n; j >= 1; --j)
                {
                    if (x[j] != 0.f)
                    {
                        if (nounit)
                        {
                            x[j] /= a[j + j * a_dim1];
                        }
                        temp = x[j];
                        for (i__ = j - 1; i__ >= 1; --i__)
                        {
                            x[i__] -= temp * a[i__ + j * a_dim1];
                        }
                    }
                }
            }
            else
            {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j)
                {
                    if (x[jx] != 0.f)
                    {
                        if (nounit)
                        {
                            x[jx] /= a[j + j * a_dim1];
                        }
                        temp = x[jx];
                        ix = jx;
                        for (i__ = j - 1; i__ >= 1; --i__)
                        {
                            ix -= *incx;
                            x[ix] -= temp * a[i__ + j * a_dim1];
                        }
                    }
                    jx -= *incx;
                }
            }
        }
        else
        {
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (x[j] != 0.f)
                    {
                        if (nounit)
                        {
                            x[j] /= a[j + j * a_dim1];
                        }
                        temp = x[j];
                        i__2 = *n;
                        for (i__ = j + 1; i__ <= i__2; ++i__)
                        {
                            x[i__] -= temp * a[i__ + j * a_dim1];
                        }
                    }
                }
            }
            else
            {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (x[jx] != 0.f)
                    {
                        if (nounit)
                        {
                            x[jx] /= a[j + j * a_dim1];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__2 = *n;
                        for (i__ = j + 1; i__ <= i__2; ++i__)
                        {
                            ix += *incx;
                            x[ix] -= temp * a[i__ + j * a_dim1];
                        }
                    }
                    jx += *incx;
                }
            }
        }
    }
    else
    {

        if (strncasecmp(uplo, "U", 1) == 0)
        {
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    temp = x[j];
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        temp -= a[i__ + j * a_dim1] * x[i__];
                    }
                    if (nounit)
                    {
                        temp /= a[j + j * a_dim1];
                    }
                    x[j] = temp;
                }
            }
            else
            {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    temp = x[jx];
                    ix = kx;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        temp -= a[i__ + j * a_dim1] * x[ix];
                        ix += *incx;
                    }
                    if (nounit)
                    {
                        temp /= a[j + j * a_dim1];
                    }
                    x[jx] = temp;
                    jx += *incx;
                }
            }
        }
        else
        {
            if (*incx == 1)
            {
                for (j = *n; j >= 1; --j)
                {
                    temp = x[j];
                    i__1 = j + 1;
                    for (i__ = *n; i__ >= i__1; --i__)
                    {
                        temp -= a[i__ + j * a_dim1] * x[i__];
                    }
                    if (nounit)
                    {
                        temp /= a[j + j * a_dim1];
                    }
                    x[j] = temp;
                }
            }
            else
            {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j)
                {
                    temp = x[jx];
                    ix = kx;
                    i__1 = j + 1;
                    for (i__ = *n; i__ >= i__1; --i__)
                    {
                        temp -= a[i__ + j * a_dim1] * x[ix];
                        ix -= *incx;
                    }
                    if (nounit)
                    {
                        temp /= a[j + j * a_dim1];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                }
            }
        }
    }

    return;
}

void lpf_blas_btrsv_fortran_dyn_rank_64(char* uplo, char* trans, char* diag,
                                        int64_t* n, CFI_cdesc_t* _a,
                                        int64_t* lda, CFI_cdesc_t* _x,
                                        int64_t* incx)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    LPF_GLOBAL(btrsv, BTRSV)(uplo, trans, diag, n, (lpf_bfloat16_t*)a, lda,
                             (lpf_bfloat16_t*)x, incx, 1, 1, 1);
}

void lpf_blas_btrsv_fortran_dyn_rank_32(char* uplo, char* trans, char* diag,
                                        int32_t* n, CFI_cdesc_t* _a,
                                        int32_t* lda, CFI_cdesc_t* _x,
                                        int32_t* incx)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    int64_t _n = *n;
    int64_t _lda = *lda;
    int64_t _incx = *incx;
    LPF_GLOBAL(btrsv, BTRSV)(uplo, trans, diag, &_n, (lpf_bfloat16_t*)a, &_lda,
                             (lpf_bfloat16_t*)x, &_incx, 1, 1, 1);
}
