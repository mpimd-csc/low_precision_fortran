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


#include <math.h>
#include <ctype.h>
#ifdef BLAS_IS_ILP64
#define MKL_ILP64
#endif
#include <mkl.h>

#include "lpf_fp8_e4m3_emu.h"

void __lpf_blas_gemm_e4m3_fp32 (
        char *_transa, char *_transb, int64_t m, int64_t n, int64_t k,
        float alpha, fp8_e4m3_t *a, int64_t lda, fp8_e4m3_t *b, int64_t ldb,
        float beta, float *c, int64_t ldc
        )
{
    MKL_INT _m = m;
    MKL_INT _n = n;
    MKL_INT _k = k;
    MKL_INT _lda = lda;
    MKL_INT _ldb = ldb;
    MKL_INT _ldc = ldc;

    char transa  = tolower(_transa[0]);
    char transb  = tolower(_transb[0]);

    CBLAS_TRANSPOSE ctransa, ctransb;

    if ( transa ==  'n' || transa == 'N')
    {
        ctransa = CblasNoTrans;
    }
    else if ( transa ==  't' || transa == 'T' || transa == 'c' || transa == 'C')
    {
        ctransa = CblasTrans;
    }

    if ( transb ==  'n' || transb == 'N')
    {
        ctransb = CblasNoTrans;
    }
    else if ( transb ==  't' || transb == 'T' || transb == 'c' || transb == 'C')
    {
        ctransb = CblasTrans;
    }

    cblas_gemm_e4m3e4m3f32(CblasColMajor, ctransa, ctransb,
            _m, _n, _k,
            alpha, (MKL_E4M3*) a, _lda, (MKL_E4M3*) b, _ldb,
            beta, c, _ldc);

}

