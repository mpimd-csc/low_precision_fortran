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

void LPF_GLOBAL(blacpy,BLACPY)(const char * uplo, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_bfloat16_t *a, lpf_blas_int_t *lda, lpf_bfloat16_t *b, lpf_blas_int_t *ldb,
                                  lpf_fortran_strlen_t uplo_len)
{
    (void) uplo_len;
    lpf_blas_int_t i,j;
    lpf_blas_int_t LDA  = *lda;
    lpf_blas_int_t LDB  = *ldb;

    if (strncasecmp(uplo, "U", 1) == 0) {
        // UPLO = U
        for ( j = 0; j < *n; j++) {
            lpf_bfloat16_t *aa = a + j * LDA;
            lpf_bfloat16_t *bb = b + j * LDB;
            lpf_blas_int_t bound = LPF_MIN(j, *m-1);
            for ( i = 0; i <= bound; i++ ) {
                bb [ i ] = (lpf_bfloat16_t) aa[i];
            }
        }
    } else if (strncasecmp(uplo, "L", 1) == 0) {
        // UPLO = L
        for ( j = 0; j < *n; j++) {
            lpf_bfloat16_t *aa = a + j * LDA;
            lpf_bfloat16_t *bb = b + j * LDB;
            lpf_blas_int_t bound = *m-1;
            for ( i = j; i <= bound; i++ ) {
                bb [ i ] = (lpf_bfloat16_t) aa[i];
            }
        }
    } else {
       for ( j = 0; j < *n; j++) {
            lpf_bfloat16_t *aa = a + j * LDA;
            lpf_bfloat16_t *bb = b + j * LDB;
            lpf_blas_int_t bound = *m-1;
            for ( i = 0; i <= bound; i++ ) {
                bb [ i ] = (lpf_bfloat16_t) aa[i];
            }
        }
    }
}

void lpf_blas_blacpy_fortran(char * uplo, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_fbfloat16_t *a, lpf_blas_int_t *lda, lpf_fbfloat16_t *b, lpf_blas_int_t *ldb)
{
    LPF_GLOBAL(blacpy,BLACPY)( uplo, m, n, (lpf_bfloat16_t *)a, lda, (lpf_bfloat16_t *)b, ldb, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_blacpy_fortran_dyn_rank(char * uplo, lpf_blas_int_t *m, lpf_blas_int_t *n, CFI_cdesc_t *_a, lpf_blas_int_t *lda, CFI_cdesc_t *_b, lpf_blas_int_t *ldb)
{
    lpf_bfloat16_t *a = _a->base_addr;
    lpf_bfloat16_t *b = _b->base_addr;
    LPF_GLOBAL(blacpy,BLACPY)( uplo, m, n, a, lda, b, ldb, 1);
}
