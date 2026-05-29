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

void LPF_GLOBAL(btrsm, BTRSM)(
    char* side, char* uplo, char* transa, char* diag, int64_t* m, int64_t* n,
    lpf_bfloat16_t* alpha, lpf_bfloat16_t* a, int64_t* lda, lpf_bfloat16_t* b,
    int64_t* ldb, lpf_fortran_strlen_t side_len, lpf_fortran_strlen_t uplo_len,
    lpf_fortran_strlen_t transa_len, lpf_fortran_strlen_t diag_len)
{

    int64_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    int64_t i__, j, k, info;
    lpf_bfloat16_t temp;
    lpf_logical_t lside;
    int64_t nrowa;
    lpf_logical_t upper;
    lpf_logical_t nounit;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    (void)diag_len;
    (void)side_len;
    (void)uplo_len;
    (void)transa_len;

    lside = strncasecmp(side, "L", 1) == 0;
    if (lside)
    {
        nrowa = *m;
    }
    else
    {
        nrowa = *n;
    }
    nounit = strncasecmp(diag, "N", 1) == 0;
    upper = strncasecmp(uplo, "U", 1) == 0;

    info = 0;
    if (!lside && !(strncasecmp(side, "R", 1) == 0))
    {
        info = 1;
    }
    else if (!upper && !(strncasecmp(uplo, "L", 1) == 0))
    {
        info = 2;
    }
    else if (!(strncasecmp(transa, "N", 1) == 0) &&
             !(strncasecmp(transa, "T", 1) == 0) &&
             !(strncasecmp(transa, "C", 1) == 0))
    {
        info = 3;
    }
    else if (!(strncasecmp(diag, "U", 1) == 0) &&
             !(strncasecmp(diag, "N", 1) == 0))
    {
        info = 4;
    }
    else if (*m < 0)
    {
        info = 5;
    }
    else if (*n < 0)
    {
        info = 6;
    }
    else if (*lda < LPF_MAX(1, nrowa))
    {
        info = 9;
    }
    else if (*ldb < LPF_MAX(1, *m))
    {
        info = 11;
    }
    if (info != 0)
    {
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BTRSM ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*m == 0 || *n == 0)
    {
        return;
    }

    if (*alpha == 0.f)
    {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__)
            {
                b[i__ + j * b_dim1] = 0.f;
            }
        }
        return;
    }

    if (lside)
    {
        if (strncasecmp(transa, "N", 1) == 0)
        {

            if (upper)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    for (k = *m; k >= 1; --k)
                    {
                        if (b[k + j * b_dim1] != 0.f)
                        {
                            if (nounit)
                            {
                                b[k + j * b_dim1] /= a[k + k * a_dim1];
                            }
                            i__2 = k - 1;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                b[i__ + j * b_dim1] -=
                                    b[k + j * b_dim1] * a[i__ + k * a_dim1];
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
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    i__2 = *m;
                    for (k = 1; k <= i__2; ++k)
                    {
                        if (b[k + j * b_dim1] != 0.f)
                        {
                            if (nounit)
                            {
                                b[k + j * b_dim1] /= a[k + k * a_dim1];
                            }
                            i__3 = *m;
                            for (i__ = k + 1; i__ <= i__3; ++i__)
                            {
                                b[i__ + j * b_dim1] -=
                                    b[k + j * b_dim1] * a[i__ + k * a_dim1];
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
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        temp = *alpha * b[i__ + j * b_dim1];
                        i__3 = i__ - 1;
                        for (k = 1; k <= i__3; ++k)
                        {
                            temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
                        }
                        if (nounit)
                        {
                            temp /= a[i__ + i__ * a_dim1];
                        }
                        b[i__ + j * b_dim1] = temp;
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    for (i__ = *m; i__ >= 1; --i__)
                    {
                        temp = *alpha * b[i__ + j * b_dim1];
                        i__2 = *m;
                        for (k = i__ + 1; k <= i__2; ++k)
                        {
                            temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
                        }
                        if (nounit)
                        {
                            temp /= a[i__ + i__ * a_dim1];
                        }
                        b[i__ + j * b_dim1] = temp;
                    }
                }
            }
        }
    }
    else
    {
        if (strncasecmp(transa, "N", 1) == 0)
        {

            if (upper)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    i__2 = j - 1;
                    for (k = 1; k <= i__2; ++k)
                    {
                        if (a[k + j * a_dim1] != 0.f)
                        {
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__)
                            {
                                b[i__ + j * b_dim1] -=
                                    a[k + j * a_dim1] * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (nounit)
                    {
                        temp = 1.f / a[j + j * a_dim1];
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        }
                    }
                }
            }
            else
            {
                for (j = *n; j >= 1; --j)
                {
                    if (*alpha != 1.f)
                    {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1];
                        }
                    }
                    i__1 = *n;
                    for (k = j + 1; k <= i__1; ++k)
                    {
                        if (a[k + j * a_dim1] != 0.f)
                        {
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                b[i__ + j * b_dim1] -=
                                    a[k + j * a_dim1] * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (nounit)
                    {
                        temp = 1.f / a[j + j * a_dim1];
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        }
                    }
                }
            }
        }
        else
        {

            if (upper)
            {
                for (k = *n; k >= 1; --k)
                {
                    if (nounit)
                    {
                        temp = 1.f / a[k + k * a_dim1];
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                        }
                    }
                    i__1 = k - 1;
                    for (j = 1; j <= i__1; ++j)
                    {
                        if (a[j + k * a_dim1] != 0.f)
                        {
                            temp = a[j + k * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__)
                            {
                                b[i__ + j * b_dim1] -=
                                    temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (*alpha != 1.f)
                    {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__)
                        {
                            b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1];
                        }
                    }
                }
            }
            else
            {
                i__1 = *n;
                for (k = 1; k <= i__1; ++k)
                {
                    if (nounit)
                    {
                        temp = 1.f / a[k + k * a_dim1];
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                        }
                    }
                    i__2 = *n;
                    for (j = k + 1; j <= i__2; ++j)
                    {
                        if (a[j + k * a_dim1] != 0.f)
                        {
                            temp = a[j + k * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__)
                            {
                                b[i__ + j * b_dim1] -=
                                    temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    if (*alpha != 1.f)
                    {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__)
                        {
                            b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1];
                        }
                    }
                }
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_btrsm_fortran_dyn_rank_64(char* side, char* uplo, char* transa,
                                        char* diag, int64_t* m, int64_t* n,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int64_t* lda, CFI_cdesc_t* _b,
                                        int64_t* ldb)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* b = _b->base_addr;
    LPF_GLOBAL(btrsm, BTRSM)(side, uplo, transa, diag, m, n,
                             (lpf_bfloat16_t*)alpha, (lpf_bfloat16_t*)a, lda,
                             (lpf_bfloat16_t*)b, ldb, 1, 1, 1, 1);
}

void lpf_blas_btrsm_fortran_dyn_rank_32(char* side, char* uplo, char* transa,
                                        char* diag, int32_t* m, int32_t* n,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int32_t* lda, CFI_cdesc_t* _b,
                                        int32_t* ldb)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* b = _b->base_addr;
    int64_t _m = *m;
    int64_t _n = *n;
    int64_t _lda = *lda;
    int64_t _ldb = *ldb;
    LPF_GLOBAL(btrsm, BTRSM)(side, uplo, transa, diag, &_m, &_n,
                             (lpf_bfloat16_t*)alpha, (lpf_bfloat16_t*)a, &_lda,
                             (lpf_bfloat16_t*)b, &_ldb, 1, 1, 1, 1);
}
