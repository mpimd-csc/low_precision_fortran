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
 * @brief Half-precision General Matrix-Matrix Multiplication with Triangular Result (HGEMMTR).
 *
 * Performs the operation: C := alpha * A * B + beta * C
 * where C is a triangular matrix.
 *
 * @param[in] uplo Character specifying the triangular part of the matrix C to be used: 'U' for upper, 'L' for lower.
 * @param[in] transa Character specifying the transpose of matrix A: 'N' for no transpose, 'T' for transpose.
 * @param[in] transb Character specifying the transpose of matrix B: 'N' for no transpose, 'T' for transpose.
 * @param[in] n Order of matrix C.
 * @param[in] k Common dimension of matrices A and B.
 * @param[in] alpha Scalar multiplier.
 * @param[in] a Matrix A.
 * @param[in] lda Leading dimension of matrix A.
 * @param[in] b Matrix B.
 * @param[in] ldb Leading dimension of matrix B.
 * @param[in] beta Scalar multiplier for matrix C.
 * @param[in,out] c Triangular matrix C.
 * @param[in] ldc Leading dimension of matrix C.
 */
void LPF_GLOBAL(hgemmtr, HGEMMTR)(char* uplo, char* transa, char* transb,
                                  int64_t* n, int64_t* k, lpf_float16_t* alpha,
                                  lpf_float16_t* a, int64_t* lda,
                                  lpf_float16_t* b, int64_t* ldb,
                                  lpf_float16_t* beta, lpf_float16_t* c__,
                                  int64_t* ldc, lpf_fortran_strlen_t uplo_len,
                                  lpf_fortran_strlen_t transa_len,
                                  lpf_fortran_strlen_t transb_len)
{

    int64_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
        i__3;

    int64_t i__, j, l, info;
    lpf_logical_t nota, notb;
    lpf_float16_t temp;
    int64_t nrowa, nrowb;
    lpf_logical_t upper;
    int64_t istop;
    int64_t istart;

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
    (void)transa_len;
    (void)transb_len;

    nota = strncasecmp(transa, "N", 1) == 0;
    notb = strncasecmp(transb, "N", 1) == 0;
    if (nota)
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    if (notb)
    {
        nrowb = *k;
    }
    else
    {
        nrowb = *n;
    }
    upper = strncasecmp(uplo, "U", 1) == 0;

    info = 0;
    if (!upper && !(strncasecmp(uplo, "L", 1) == 0))
    {
        info = 1;
    }
    else if (!nota && !(strncasecmp(transa, "C", 1) == 0) &&
             !(strncasecmp(transa, "T", 1) == 0))
    {
        info = 2;
    }
    else if (!notb && !(strncasecmp(transb, "C", 1) == 0) &&
             !(strncasecmp(transb, "T", 1) == 0))
    {
        info = 3;
    }
    else if (*n < 0)
    {
        info = 4;
    }
    else if (*k < 0)
    {
        info = 5;
    }
    else if (*lda < LPF_MAX(1, nrowa))
    {
        info = 8;
    }
    else if (*ldb < LPF_MAX(1, nrowb))
    {
        info = 10;
    }
    else if (*ldc < LPF_MAX(1, *n))
    {
        info = 13;
    }
    if (info != 0)
    {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HGEMMTR", &infox,
                                                     (lpf_fortran_strlen_t)7);
        return;
    }

    if (*n == 0)
    {
        return;
    }

    if (*alpha == 0.f)
    {
        if (*beta == 0.f)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (upper)
                {
                    istart = 1;
                    istop = j;
                }
                else
                {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__)
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
                if (upper)
                {
                    istart = 1;
                    istop = j;
                }
                else
                {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__)
                {
                    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                }
            }
        }
        return;
    }

    if (notb)
    {
        if (nota)
        {

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (upper)
                {
                    istart = 1;
                    istop = j;
                }
                else
                {
                    istart = j;
                    istop = *n;
                }
                if (*beta == 0.f)
                {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    temp = *alpha * b[l + j * b_dim1];
                    i__3 = istop;
                    for (i__ = istart; i__ <= i__3; ++i__)
                    {
                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                    }
                }
            }
        }
        else
        {

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (upper)
                {
                    istart = 1;
                    istop = j;
                }
                else
                {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__)
                {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l)
                    {
                        temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                    }
                    if (*beta == 0.f)
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    }
                    else
                    {
                        c__[i__ + j * c_dim1] =
                            *alpha * temp + *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
    }
    else
    {
        if (nota)
        {

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (upper)
                {
                    istart = 1;
                    istop = j;
                }
                else
                {
                    istart = j;
                    istop = *n;
                }
                if (*beta == 0.f)
                {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    temp = *alpha * b[j + l * b_dim1];
                    i__3 = istop;
                    for (i__ = istart; i__ <= i__3; ++i__)
                    {
                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                    }
                }
            }
        }
        else
        {

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (upper)
                {
                    istart = 1;
                    istop = j;
                }
                else
                {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__)
                {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l)
                    {
                        temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
                    }
                    if (*beta == 0.f)
                    {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    }
                    else
                    {
                        c__[i__ + j * c_dim1] =
                            *alpha * temp + *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hgemmtr_fortran_dyn_rank_64(
    char* uplo, char* transa, char* transb, int64_t* n, int64_t* k,
    lpf_ffloat16_t* alpha, CFI_cdesc_t* _a, int64_t* lda, CFI_cdesc_t* _b,
    int64_t* ldb, lpf_ffloat16_t* beta, CFI_cdesc_t* _c, int64_t* ldc)
{
    lpf_float16_t* a = _a->base_addr;
    lpf_float16_t* b = _b->base_addr;
    lpf_float16_t* c = _c->base_addr;

    LPF_GLOBAL(hgemmtr, HGEMMTR)(uplo, transa, transb, n, k,
                                 (lpf_float16_t*)alpha, (lpf_float16_t*)a, lda,
                                 (lpf_float16_t*)b, ldb, (lpf_float16_t*)beta,
                                 (lpf_float16_t*)c, ldc, 1, 1, 1);
}

void lpf_blas_hgemmtr_fortran_dyn_rank_32(
    char* uplo, char* transa, char* transb, int32_t* n, int32_t* k,
    lpf_ffloat16_t* alpha, CFI_cdesc_t* _a, int32_t* lda, CFI_cdesc_t* _b,
    int32_t* ldb, lpf_ffloat16_t* beta, CFI_cdesc_t* _c, int32_t* ldc)
{
    int64_t _n = *n;
    int64_t _k = *k;
    int64_t _lda = *lda;
    int64_t _ldb = *ldb;
    int64_t _ldc = *ldc;

    lpf_float16_t* a = _a->base_addr;
    lpf_float16_t* b = _b->base_addr;
    lpf_float16_t* c = _c->base_addr;

    LPF_GLOBAL(hgemmtr, HGEMMTR)(
        uplo, transa, transb, &_n, &_k, (lpf_float16_t*)alpha,
        (lpf_float16_t*)a, &_lda, (lpf_float16_t*)b, &_ldb,
        (lpf_float16_t*)beta, (lpf_float16_t*)c, &_ldc, 1, 1, 1);
}
