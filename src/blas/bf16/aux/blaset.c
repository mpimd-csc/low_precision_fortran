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

void LPF_GLOBAL(blaset,BLASET)(const char * uplo, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_bfloat16_t *alpha,
        lpf_bfloat16_t *beta,lpf_bfloat16_t *a , lpf_blas_int_t *lda, lpf_fortran_strlen_t uplo_len)
{
    (void) uplo_len;
    lpf_blas_int_t i,j;
    lpf_blas_int_t LDA  = *lda;
    lpf_bfloat16_t fa = *alpha;
    lpf_bfloat16_t fb = *beta;


    if (strncasecmp(uplo, "U", 1) == 0) {
        // UPLO = U
        for ( j = 0; j < *n; j++) {
            lpf_bfloat16_t *aa = a + j * LDA;
            lpf_blas_int_t bound = LPF_MIN(j, *m-1);
            for ( i = 0; i <= bound; i++ ) {
                aa[i] = fa ;
            }
            aa[bound] = fb;
        }
    } else if (strncasecmp(uplo, "L", 1) == 0) {
        // UPLO = L
        for ( j = 0; j < *n; j++) {
            lpf_bfloat16_t *aa = a + j * LDA;
            lpf_blas_int_t bound = *m-1;
            for ( i = j; i <= bound; i++ ) {
                aa[i] = fa;
            }
            aa[j] = fb;
        }
    } else {
       for ( j = 0; j < *n; j++) {
            lpf_bfloat16_t *aa = a + j * LDA;
            lpf_blas_int_t bound = *m-1;
            for ( i = 0; i <= bound; i++ ) {
                aa[i] = fa;
                if ( i == j ) aa[i] = fb;
            }

        }
    }
}

void hpblas_blaset_fortran(char * uplo, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_fbfloat16_t *alpha, lpf_fbfloat16_t *beta, lpf_fbfloat16_t *a, lpf_blas_int_t *lda)
{
    LPF_GLOBAL(blaset,BLASET)( uplo, m, n, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t*) beta, (lpf_bfloat16_t *)a, lda, 1);
}
