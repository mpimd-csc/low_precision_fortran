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
 * @brief Half-precision Symmetric Rank-k Update (HSYRK).
 *
 * Performs the operation: C := alpha * (A * B^T + B * A^T) + beta * C
 * where C is a symmetric matrix.
 *
 * @param[in] uplo Character specifying the triangular part of the matrix C to be used: 'U' for upper, 'L' for lower.
 * @param[in] trans Character specifying the transpose of matrices A and B: 'N' for no transpose, 'T' for transpose.
 * @param[in] n Order of matrix C.
 * @param[in] k Common dimension of matrices A and B.
 * @param[in] alpha Scalar multiplier.
 * @param[in] a Matrix A.
 * @param[in] lda Leading dimension of matrix A.
 * @param[in] b Matrix B.
 * @param[in] ldb Leading dimension of matrix B.
 * @param[in] beta Scalar multiplier for matrix C.
 * @param[in,out] c Symmetric matrix C.
 * @param[in] ldc Leading dimension of matrix C.
 */
void LPF_GLOBAL(hsyr2k, HSYR2K)(char* uplo, char* trans, int64_t* n, int64_t* k,
                                lpf_float16_t* alpha, lpf_float16_t* a,
                                int64_t* lda, lpf_float16_t* b, int64_t* ldb,
                                lpf_float16_t* beta, lpf_float16_t* c__,
                                int64_t* ldc, lpf_fortran_strlen_t uplo_len,
                                lpf_fortran_strlen_t trans_len)
{

    int64_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
        i__3;

    int64_t i__, j, l, info;
    lpf_float16_t temp1, temp2;
    int64_t nrowa;
    lpf_logical_t upper;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    (void)uplo_len;
    (void)trans_len;

    if (strncasecmp(trans, "N", 1) == 0)
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    upper = strncasecmp(uplo, "U", 1) == 0;

    info = 0;
    if (!upper && !(strncasecmp(uplo, "L", 1) == 0))
    {
        info = 1;
    }
    else if (!(strncasecmp(trans, "N", 1) == 0) &&
             !(strncasecmp(trans, "T", 1) == 0) &&
             !(strncasecmp(trans, "C", 1) == 0))
    {
        info = 2;
    }
    else if (*n < 0)
    {
        info = 3;
    }
    else if (*k < 0)
    {
        info = 4;
    }
    else if (*lda < LPF_MAX(1, nrowa))
    {
        info = 7;
    }
    else if (*ldb < LPF_MAX(1, nrowa))
    {
        info = 9;
    }
    else if (*ldc < LPF_MAX(1, *n))
    {
        info = 12;
    }
    if (info != 0)
    {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HSYR2K", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*n == 0 || ((*alpha == 0.f || *k == 0) && *beta == 1.f))
    {
        return;
    }

    if (*alpha == 0.f)
    {
        if (upper)
        {
            if (*beta == 0.f)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
        else
        {
            if (*beta == 0.f)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
        return;
    }

    if (strncasecmp(trans, "N", 1) == 0)
    {

        if (upper)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (*beta == 0.f)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    if (a[j + l * a_dim1] != 0.f || b[j + l * b_dim1] != 0.f)
                    {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = j;
                        for (i__ = 1; i__ <= i__3; ++i__)
                        {
                            c__[i__ + j * c_dim1] =
                                c__[i__ + j * c_dim1] +
                                a[i__ + l * a_dim1] * temp1 +
                                b[i__ + l * b_dim1] * temp2;
                        }
                    }
                }
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (*beta == 0.f)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    if (a[j + l * a_dim1] != 0.f || b[j + l * b_dim1] != 0.f)
                    {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = *n;
                        for (i__ = j; i__ <= i__3; ++i__)
                        {
                            c__[i__ + j * c_dim1] =
                                c__[i__ + j * c_dim1] +
                                a[i__ + l * a_dim1] * temp1 +
                                b[i__ + l * b_dim1] * temp2;
                        }
                    }
                }
            }
        }
    }
    else
    {

        if (upper)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    temp1 = 0.f;
                    temp2 = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l)
                    {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
                    }
                    if (*beta == 0.f)
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * temp2;
                    }
                    else
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] +
                                                *alpha * temp1 + *alpha * temp2;
                    }
                }
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = *n;
                for (i__ = j; i__ <= i__2; ++i__)
                {
                    temp1 = 0.f;
                    temp2 = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l)
                    {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
                    }
                    if (*beta == 0.f)
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * temp2;
                    }
                    else
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] +
                                                *alpha * temp1 + *alpha * temp2;
                    }
                }
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hsyr2k_fortran_dyn_rank_64(char* uplo, char* trans, int64_t* n,
                                         int64_t* k, lpf_ffloat16_t* alpha,
                                         CFI_cdesc_t* _a, int64_t* lda,
                                         CFI_cdesc_t* _b, int64_t* ldb,
                                         lpf_ffloat16_t* beta, CFI_cdesc_t* _c,
                                         int64_t* ldc)
{
    lpf_float16_t* a = _a->base_addr;
    lpf_float16_t* b = _b->base_addr;
    lpf_float16_t* c = _c->base_addr;

    LPF_GLOBAL(hsyr2k, HSYR2K)(uplo, trans, n, k, (lpf_float16_t*)alpha,
                               (lpf_float16_t*)a, lda, (lpf_float16_t*)b, ldb,
                               (lpf_float16_t*)beta, (lpf_float16_t*)c, ldc, 1,
                               1);
}

void lpf_blas_hsyr2k_fortran_dyn_rank_32(char* uplo, char* trans, int32_t* n,
                                         int32_t* k, lpf_ffloat16_t* alpha,
                                         CFI_cdesc_t* _a, int32_t* lda,
                                         CFI_cdesc_t* _b, int32_t* ldb,
                                         lpf_ffloat16_t* beta, CFI_cdesc_t* _c,
                                         int32_t* ldc)
{
    int64_t _n = *n;
    int64_t _k = *k;
    int64_t _lda = *lda;
    int64_t _ldb = *ldb;
    int64_t _ldc = *ldc;

    lpf_float16_t* a = _a->base_addr;
    lpf_float16_t* b = _b->base_addr;
    lpf_float16_t* c = _c->base_addr;

    LPF_GLOBAL(hsyr2k, HSYR2K)(uplo, trans, &_n, &_k, (lpf_float16_t*)alpha,
                               (lpf_float16_t*)a, &_lda, (lpf_float16_t*)b,
                               &_ldb, (lpf_float16_t*)beta, (lpf_float16_t*)c,
                               &_ldc, 1, 1);
}
