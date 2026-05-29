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

void LPF_GLOBAL(bspr, BSPR)(char* uplo, int64_t* n, lpf_bfloat16_t* alpha,
                            lpf_bfloat16_t* x, int64_t* incx,
                            lpf_bfloat16_t* ap, lpf_fortran_strlen_t uplo_len)
{

    int64_t i__1, i__2;

    int64_t i__, j, k, kk, ix, jx, kx, info;
    lpf_bfloat16_t temp;

    --ap;
    --x;
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
    if (info != 0)
    {
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSPR  ", &infox,
                                                     (lpf_fortran_strlen_t)6);
        return;
    }

    if (*n == 0 || *alpha == 0.f)
    {
        return;
    }

    if (*incx <= 0)
    {
        kx = 1 - (*n - 1) * *incx;
    }
    else if (*incx != 1)
    {
        kx = 1;
    }

    kk = 1;
    if (strncasecmp(uplo, "U", 1) == 0)
    {

        if (*incx == 1)
        {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                if (x[j] != 0.f)
                {
                    temp = *alpha * x[j];
                    k = kk;
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        ap[k] += x[i__] * temp;
                        ++k;
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
                    temp = *alpha * x[jx];
                    ix = kx;
                    i__2 = kk + j - 1;
                    for (k = kk; k <= i__2; ++k)
                    {
                        ap[k] += x[ix] * temp;
                        ix += *incx;
                    }
                }
                jx += *incx;
                kk += j;
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
                    temp = *alpha * x[j];
                    k = kk;
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__)
                    {
                        ap[k] += x[i__] * temp;
                        ++k;
                    }
                }
                kk = kk + *n - j + 1;
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
                    temp = *alpha * x[jx];
                    ix = jx;
                    i__2 = kk + *n - j;
                    for (k = kk; k <= i__2; ++k)
                    {
                        ap[k] += x[ix] * temp;
                        ix += *incx;
                    }
                }
                jx += *incx;
                kk = kk + *n - j + 1;
            }
        }
    }

    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_bspr_fortran_dyn_rank_64(char* uplo, int64_t* n,
                                       lpf_fbfloat16_t* alpha, CFI_cdesc_t* _x,
                                       int64_t* incx, CFI_cdesc_t* _ap)
{
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* ap = _ap->base_addr;
    LPF_GLOBAL(bspr, BSPR)(uplo, n, (lpf_bfloat16_t*)alpha, (lpf_bfloat16_t*)x,
                           incx, (lpf_bfloat16_t*)ap, 1);
}

void lpf_blas_bspr_fortran_dyn_rank_32(char* uplo, int32_t* n,
                                       lpf_fbfloat16_t* alpha, CFI_cdesc_t* _x,
                                       int32_t* incx, CFI_cdesc_t* _ap)
{
    lpf_bfloat16_t* x = _x->base_addr;
    lpf_bfloat16_t* ap = _ap->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    LPF_GLOBAL(bspr, BSPR)(uplo, &_n, (lpf_bfloat16_t*)alpha,
                           (lpf_bfloat16_t*)x, &_incx, (lpf_bfloat16_t*)ap, 1);
}
