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
 * @brief Bfloat16 Triangular Packed Matrix-Vector Multiplication (BTPMV).
 *
 * Performs the operation: x := A * x
 * where A is a triangular packed matrix.
 *
 * @param[in] uplo Character specifying the triangular part of the matrix A to be used: 'U' for upper, 'L' for lower.
 * @param[in] trans Character specifying the transpose. 'N' for no transpose, 'T' for transpose.
 * @param[in] diag Character specifying whether the diagonal of A is unit: 'U' for unit, 'N' for non-unit.
 * @param[in] n Order of matrix A.
 * @param[in] ap The triangular packed matrix A.
 * @param[in,out] x Vector X.
 * @param[in] incx Increment for the elements of x.
 */
void LPF_GLOBAL(btpmv, BTPMV)(char* uplo, char* trans, char* diag, int64_t* n,
                              lpf_bfloat16_t* ap, lpf_bfloat16_t* x,
                              int64_t* incx, lpf_fortran_strlen_t uplo_len,
                              lpf_fortran_strlen_t trans_len,
                              lpf_fortran_strlen_t diag_len)
{

    int64_t i__1, i__2;

    int64_t i__, j, k, kk, ix, jx, kx, info;
    lpf_bfloat16_t temp;
    lpf_logical_t nounit;

    --x;
    --ap;

    (void)uplo_len;
    (void)trans_len;
    (void)diag_len;

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
    else if (*incx == 0)
    {
        info = 7;
    }
    if (info != 0)
    {
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BTPMV ", &infox,
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
            kk = 1;
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (x[j] != 0.f)
                    {
                        temp = x[j];
                        k = kk;
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            x[i__] += temp * ap[k];
                            ++k;
                        }
                        if (nounit)
                        {
                            x[j] *= ap[kk + j - 1];
                        }
                    }
                    kk += j;
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
                        temp = x[jx];
                        ix = kx;
                        i__2 = kk + j - 2;
                        for (k = kk; k <= i__2; ++k)
                        {
                            x[ix] += temp * ap[k];
                            ix += *incx;
                        }
                        if (nounit)
                        {
                            x[jx] *= ap[kk + j - 1];
                        }
                    }
                    jx += *incx;
                    kk += j;
                }
            }
        }
        else
        {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1)
            {
                for (j = *n; j >= 1; --j)
                {
                    if (x[j] != 0.f)
                    {
                        temp = x[j];
                        k = kk;
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__)
                        {
                            x[i__] += temp * ap[k];
                            --k;
                        }
                        if (nounit)
                        {
                            x[j] *= ap[kk - *n + j];
                        }
                    }
                    kk -= *n - j + 1;
                }
            }
            else
            {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j)
                {
                    if (x[jx] != 0.f)
                    {
                        temp = x[jx];
                        ix = kx;
                        i__1 = kk - (*n - (j + 1));
                        for (k = kk; k >= i__1; --k)
                        {
                            x[ix] += temp * ap[k];
                            ix -= *incx;
                        }
                        if (nounit)
                        {
                            x[jx] *= ap[kk - *n + j];
                        }
                    }
                    jx -= *incx;
                    kk -= *n - j + 1;
                }
            }
        }
    }
    else
    {

        if (strncasecmp(uplo, "U", 1) == 0)
        {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1)
            {
                for (j = *n; j >= 1; --j)
                {
                    temp = x[j];
                    if (nounit)
                    {
                        temp *= ap[kk];
                    }
                    k = kk - 1;
                    for (i__ = j - 1; i__ >= 1; --i__)
                    {
                        temp += ap[k] * x[i__];
                        --k;
                    }
                    x[j] = temp;
                    kk -= j;
                }
            }
            else
            {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j)
                {
                    temp = x[jx];
                    ix = jx;
                    if (nounit)
                    {
                        temp *= ap[kk];
                    }
                    i__1 = kk - j + 1;
                    for (k = kk - 1; k >= i__1; --k)
                    {
                        ix -= *incx;
                        temp += ap[k] * x[ix];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    kk -= j;
                }
            }
        }
        else
        {
            kk = 1;
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    temp = x[j];
                    if (nounit)
                    {
                        temp *= ap[kk];
                    }
                    k = kk + 1;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__)
                    {
                        temp += ap[k] * x[i__];
                        ++k;
                    }
                    x[j] = temp;
                    kk += *n - j + 1;
                }
            }
            else
            {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    temp = x[jx];
                    ix = jx;
                    if (nounit)
                    {
                        temp *= ap[kk];
                    }
                    i__2 = kk + *n - j;
                    for (k = kk + 1; k <= i__2; ++k)
                    {
                        ix += *incx;
                        temp += ap[k] * x[ix];
                    }
                    x[jx] = temp;
                    jx += *incx;
                    kk += *n - j + 1;
                }
            }
        }
    }

    return;
}

void lpf_blas_btpmv_fortran_dyn_rank_64(char* uplo, char* trans, char* diag,
                                        int64_t* n, CFI_cdesc_t* _ap,
                                        CFI_cdesc_t* _x, int64_t* incx)
{
    lpf_bfloat16_t* ap = _ap->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    LPF_GLOBAL(btpmv, BTPMV)(uplo, trans, diag, n, (lpf_bfloat16_t*)ap,
                             (lpf_bfloat16_t*)x, incx, 1, 1, 1);
}

void lpf_blas_btpmv_fortran_dyn_rank_32(char* uplo, char* trans, char* diag,
                                        int32_t* n, CFI_cdesc_t* _ap,
                                        CFI_cdesc_t* _x, int32_t* incx)
{
    lpf_bfloat16_t* ap = _ap->base_addr;
    lpf_bfloat16_t* x = _x->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    LPF_GLOBAL(btpmv, BTPMV)(uplo, trans, diag, &_n, (lpf_bfloat16_t*)ap,
                             (lpf_bfloat16_t*)x, &_incx, 1, 1, 1);
}
