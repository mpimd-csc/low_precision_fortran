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
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef BLAS_IS_MKL
#include <mkl.h>
#endif

void LPF_GLOBAL(bgemm, BGEMM)(char* transa, char* transb, int64_t* m,
                              int64_t* n, int64_t* k, lpf_bfloat16_t* alpha,
                              lpf_bfloat16_t* a, int64_t* lda,
                              lpf_bfloat16_t* b, int64_t* ldb,
                              lpf_bfloat16_t* beta, lpf_bfloat16_t* c__,
                              int64_t* ldc, lpf_fortran_strlen_t transa_len,
                              lpf_fortran_strlen_t transb_len)
{

    int64_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
        i__3;

    int64_t i__, j, l, info;
    lpf_logical_t nota, notb;
    lpf_bfloat16_t temp;
    int64_t nrowa, nrowb;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
#ifndef USE_BLAS_BGEMM
    a -= a_offset;
    b -= b_offset;
    c__ -= c_offset;
#endif

    (void)transa_len;
    (void)transb_len;

    nota = strncasecmp(transa, "N", 1) == 0;
    notb = strncasecmp(transb, "N", 1) == 0;
    if (nota)
    {
        nrowa = *m;
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

    info = 0;
    if (!nota && !(strncasecmp(transa, "C", 1) == 0) &&
        !(strncasecmp(transa, "T", 1) == 0))
    {
        info = 1;
    }
    else if (!notb && !(strncasecmp(transb, "C", 1) == 0) &&
             !(strncasecmp(transb, "T", 1) == 0))
    {
        info = 2;
    }
    else if (*m < 0)
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
    else if (*ldc < LPF_MAX(1, *m))
    {
        info = 13;
    }
    if (info != 0)
    {
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BGEMM ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*m == 0 || *n == 0 || ((*alpha == 0.f || *k == 0) && *beta == 1.f))
    {
        return;
    }

#ifdef USE_BLAS_BGEMM
    CBLAS_TRANSPOSE ta, tb;
    union
    {
        lpf_bfloat16_t f;
        MKL_F16 i;
    } malpha, mbeta;

    if (tolower(transa[0]) == 'n')
    {
        ta = CblasNoTrans;
    }
    else
    {
        ta = CblasTrans;
    }
    if (tolower(transb[0]) == 'n')
    {
        tb = CblasNoTrans;
    }
    else
    {
        tb = CblasTrans;
    }

    malpha.f = *alpha;
    mbeta.f = *beta;

    cblas_bgemm(CblasColMajor, ta, tb, *m, *n, *k, malpha.i, (MKL_F16*)a, *lda,
                (MKL_F16*)b, *ldb, mbeta.i, (MKL_F16*)c__, *ldc);

    return;

#else

    if (*alpha == 0.f)
    {
        if (*beta == 0.f)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = *m;
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
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
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
                if (*beta == 0.f)
                {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    temp = *alpha * b[l + j * b_dim1];
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__)
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
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
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
                if (*beta == 0.f)
                {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = 0.f;
                    }
                }
                else if (*beta != 1.f)
                {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l)
                {
                    temp = *alpha * b[j + l * b_dim1];
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__)
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
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__)
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

#endif
}

void lpf_blas_bgemm_fortran_dyn_rank_64(char* transa, char* transb, int64_t* m,
                                        int64_t* n, int64_t* k,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int64_t* lda, CFI_cdesc_t* _b,
                                        int64_t* ldb, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _c, int64_t* ldc)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* b = _b->base_addr;
    lpf_bfloat16_t* c = _c->base_addr;
    LPF_GLOBAL(bgemm, BGEMM)(transa, transb, m, n, k, (lpf_bfloat16_t*)alpha,
                             (lpf_bfloat16_t*)a, lda, (lpf_bfloat16_t*)b, ldb,
                             (lpf_bfloat16_t*)beta, (lpf_bfloat16_t*)c, ldc, 1,
                             1);
}

void lpf_blas_bgemm_fortran_dyn_rank_32(char* transa, char* transb, int32_t* m,
                                        int32_t* n, int32_t* k,
                                        lpf_fbfloat16_t* alpha, CFI_cdesc_t* _a,
                                        int32_t* lda, CFI_cdesc_t* _b,
                                        int32_t* ldb, lpf_fbfloat16_t* beta,
                                        CFI_cdesc_t* _c, int32_t* ldc)
{
    lpf_bfloat16_t* a = _a->base_addr;
    lpf_bfloat16_t* b = _b->base_addr;
    lpf_bfloat16_t* c = _c->base_addr;
    int64_t _m = *m;
    int64_t _n = *n;
    int64_t _k = *k;
    int64_t _lda = *lda;
    int64_t _ldb = *ldb;
    int64_t _ldc = *ldc;
    LPF_GLOBAL(bgemm, BGEMM)(transa, transb, &_m, &_n, &_k,
                             (lpf_bfloat16_t*)alpha, (lpf_bfloat16_t*)a, &_lda,
                             (lpf_bfloat16_t*)b, &_ldb, (lpf_bfloat16_t*)beta,
                             (lpf_bfloat16_t*)c, &_ldc, 1, 1);
}
