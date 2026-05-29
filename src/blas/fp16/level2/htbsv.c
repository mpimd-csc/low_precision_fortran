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

void LPF_GLOBAL(htbsv, HTBSV)(char* uplo, char* trans, char* diag, int64_t* n,
                              int64_t* k, lpf_float16_t* a, int64_t* lda,
                              lpf_float16_t* x, int64_t* incx,
                              lpf_fortran_strlen_t uplo_len,
                              lpf_fortran_strlen_t trans_len,
                              lpf_fortran_strlen_t diag_len)
{

    int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int64_t i__, j, l, ix, jx, kx, info;
    lpf_float16_t temp;
    int64_t kplus1;
    lpf_logical_t nounit;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

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
    else if (*k < 0)
    {
        info = 5;
    }
    else if (*lda < *k + 1)
    {
        info = 7;
    }
    else if (*incx == 0)
    {
        info = 9;
    }
    if (info != 0)
    {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HTBSV ", &infox,
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
            kplus1 = *k + 1;
            if (*incx == 1)
            {
                for (j = *n; j >= 1; --j)
                {
                    if (x[j] != 0.f)
                    {
                        l = kplus1 - j;
                        if (nounit)
                        {
                            x[j] /= a[kplus1 + j * a_dim1];
                        }
                        temp = x[j];

                        i__2 = 1, i__3 = j - *k;
                        i__1 = LPF_MAX(i__2, i__3);
                        for (i__ = j - 1; i__ >= i__1; --i__)
                        {
                            x[i__] -= temp * a[l + i__ + j * a_dim1];
                        }
                    }
                }
            }
            else
            {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j)
                {
                    kx -= *incx;
                    if (x[jx] != 0.f)
                    {
                        ix = kx;
                        l = kplus1 - j;
                        if (nounit)
                        {
                            x[jx] /= a[kplus1 + j * a_dim1];
                        }
                        temp = x[jx];

                        i__2 = 1, i__3 = j - *k;
                        i__1 = LPF_MAX(i__2, i__3);
                        for (i__ = j - 1; i__ >= i__1; --i__)
                        {
                            x[ix] -= temp * a[l + i__ + j * a_dim1];
                            ix -= *incx;
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
                        l = 1 - j;
                        if (nounit)
                        {
                            x[j] /= a[j * a_dim1 + 1];
                        }
                        temp = x[j];

                        i__3 = *n, i__4 = j + *k;
                        i__2 = LPF_MIN(i__3, i__4);
                        for (i__ = j + 1; i__ <= i__2; ++i__)
                        {
                            x[i__] -= temp * a[l + i__ + j * a_dim1];
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
                    kx += *incx;
                    if (x[jx] != 0.f)
                    {
                        ix = kx;
                        l = 1 - j;
                        if (nounit)
                        {
                            x[jx] /= a[j * a_dim1 + 1];
                        }
                        temp = x[jx];

                        i__3 = *n, i__4 = j + *k;
                        i__2 = LPF_MIN(i__3, i__4);
                        for (i__ = j + 1; i__ <= i__2; ++i__)
                        {
                            x[ix] -= temp * a[l + i__ + j * a_dim1];
                            ix += *incx;
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
            kplus1 = *k + 1;
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    temp = x[j];
                    l = kplus1 - j;

                    i__2 = 1, i__3 = j - *k;
                    i__4 = j - 1;
                    for (i__ = LPF_MAX(i__2, i__3); i__ <= i__4; ++i__)
                    {
                        temp -= a[l + i__ + j * a_dim1] * x[i__];
                    }
                    if (nounit)
                    {
                        temp /= a[kplus1 + j * a_dim1];
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
                    l = kplus1 - j;

                    i__4 = 1, i__2 = j - *k;
                    i__3 = j - 1;
                    for (i__ = LPF_MAX(i__4, i__2); i__ <= i__3; ++i__)
                    {
                        temp -= a[l + i__ + j * a_dim1] * x[ix];
                        ix += *incx;
                    }
                    if (nounit)
                    {
                        temp /= a[kplus1 + j * a_dim1];
                    }
                    x[jx] = temp;
                    jx += *incx;
                    if (j > *k)
                    {
                        kx += *incx;
                    }
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
                    l = 1 - j;

                    i__1 = *n, i__3 = j + *k;
                    i__4 = j + 1;
                    for (i__ = LPF_MIN(i__1, i__3); i__ >= i__4; --i__)
                    {
                        temp -= a[l + i__ + j * a_dim1] * x[i__];
                    }
                    if (nounit)
                    {
                        temp /= a[j * a_dim1 + 1];
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
                    l = 1 - j;

                    i__4 = *n, i__1 = j + *k;
                    i__3 = j + 1;
                    for (i__ = LPF_MIN(i__4, i__1); i__ >= i__3; --i__)
                    {
                        temp -= a[l + i__ + j * a_dim1] * x[ix];
                        ix -= *incx;
                    }
                    if (nounit)
                    {
                        temp /= a[j * a_dim1 + 1];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    if (*n - j >= *k)
                    {
                        kx -= *incx;
                    }
                }
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_htbsv_fortran_dyn_rank_64(char* uplo, char* trans, char* diag,
                                        int64_t* n, int64_t* k, CFI_cdesc_t* a,
                                        int64_t* lda, CFI_cdesc_t* x,
                                        int64_t* incx)
{
    lpf_float16_t* a_ptr = a->base_addr;
    lpf_float16_t* x_ptr = x->base_addr;

    LPF_GLOBAL(htbsv, HTBSV)(uplo, trans, diag, n, k, (lpf_float16_t*)a_ptr,
                             lda, (lpf_float16_t*)x_ptr, incx, 1, 1, 1);
}

void lpf_blas_htbsv_fortran_dyn_rank_32(char* uplo, char* trans, char* diag,
                                        int32_t* n, int32_t* k, CFI_cdesc_t* a,
                                        int32_t* lda, CFI_cdesc_t* x,
                                        int32_t* incx)
{
    lpf_float16_t* a_ptr = a->base_addr;
    lpf_float16_t* x_ptr = x->base_addr;
    int64_t _n = *n;
    int64_t _k = *k;
    int64_t _lda = *lda;
    int64_t _incx = *incx;

    LPF_GLOBAL(htbsv, HTBSV)(uplo, trans, diag, &_n, &_k, (lpf_float16_t*)a_ptr,
                             &_lda, (lpf_float16_t*)x_ptr, &_incx, 1, 1, 1);
}
