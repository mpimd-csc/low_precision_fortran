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
#include <stdint.h>
#include "lpf_internal.h"

#include <stdio.h>
#include <string.h>

void LPF_GLOBAL(s2blacpy,S2BLACPY)(char * uplo, int64_t *m, int64_t *n, float *a, int64_t *lda, lpf_bfloat16_t *b, int64_t *ldb,
                                  lpf_fortran_strlen_t uplo_len)
{
    (void) uplo_len;
    int64_t i,j;
    int64_t LDA  = *lda;
    int64_t LDB  = *ldb;

    if (strncasecmp(uplo, "U", 1) == 0) {
        // UPLO = U
        for ( j = 0; j < *n; j++) {
            float *aa = a + j * LDA;
            lpf_bfloat16_t *bb = b + j * LDB;
            int64_t bound = LPF_MIN(j, *m-1);
            for ( i = 0; i <= bound; i++ ) {
                bb [ i ] = (lpf_bfloat16_t) aa[i];
            }
        }
    } else if (strncasecmp(uplo, "L", 1) == 0) {
        // UPLO = L
        for ( j = 0; j < *n; j++) {
            float *aa = a + j * LDA;
            lpf_bfloat16_t *bb = b + j * LDB;
            int64_t bound = *m-1;
            for ( i = j; i <= bound; i++ ) {
                bb [ i ] = (lpf_bfloat16_t) aa[i];
            }
        }
    } else {
       for ( j = 0; j < *n; j++) {
            float *aa = a + j * LDA;
            lpf_bfloat16_t *bb = b + j * LDB;
            int64_t bound = *m-1;
            for ( i = 0; i <= bound; i++ ) {
                bb [ i ] = (lpf_bfloat16_t) aa[i];
            }
        }
    }
}

#include <ISO_Fortran_binding.h>

void lpf_blas_s2blacpy_fortran_dyn_rank_64(char * uplo, int64_t *m, int64_t *n, float *a, int64_t *lda, CFI_cdesc_t *_b, int64_t *ldb)
{
    lpf_bfloat16_t *b = _b->base_addr;
    LPF_GLOBAL(s2blacpy,S2BLACPY)( uplo, m, n, a, lda, b, ldb, 1);
}

void lpf_blas_s2blacpy_fortran_dyn_rank_32(char * uplo, int32_t *m, int32_t *n, float *a, int32_t *lda, CFI_cdesc_t *_b, int32_t *ldb)
{
    lpf_bfloat16_t *b = _b->base_addr;
    int64_t _m = *m;
    int64_t _n = *n;
    int64_t _lda = *lda;
    int64_t _ldb = *ldb;
    LPF_GLOBAL(s2blacpy,S2BLACPY)( uplo, &_m, &_n, a, &_lda, b, &_ldb, 1);
}

