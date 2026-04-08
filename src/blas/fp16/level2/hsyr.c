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

#include <string.h>

/* > \brief \b HSYR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HSYR(UPLO,N,ALPHA,X,INCX,A,LDA) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HSYR   performs the symmetric rank 1 operation */
/* > */
/* >    A := alpha*x*x**T + A, */
/* > */
/* > where alpha is a lpf_float16_t scalar, x is an n element vector and A is an */
/* > n by n symmetric matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array A is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/* >                                  is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. On exit, the */
/* >           upper triangular part of the array A is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           upper triangular part of A is not referenced. On exit, the */
/* >           lower triangular part of the array A is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           LPF_MAX( 1, n ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(hsyr,HSYR)(char *uplo, lpf_blas_int_t *n, lpf_float16_t *alpha, lpf_float16_t *x,
        lpf_blas_int_t *incx, lpf_float16_t *a, lpf_blas_int_t *lda, lpf_fortran_strlen_t uplo_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    lpf_blas_int_t i__, j, ix, jx, kx, info;
    lpf_float16_t temp;

    /*  -- Reference BLAS level2 routine (version 3.4.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    (void) uplo_len;
    /* Function Body */
    info = 0;
    if (! (strncasecmp(uplo, "U", 1) == 0) && ! (strncasecmp(uplo, "L", 1) == 0)) {
            info = 1;
        } else if (*n < 0) {
            info = 2;
        } else if (*incx == 0) {
            info = 5;
        } else if (*lda < LPF_MAX(1,*n)) {
            info = 7;
        }
        if (info != 0) {
            LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HSYR  ", &info, (lpf_fortran_strlen_t)6);
            return;
        }

        /*     Quick return if possible. */

        if (*n == 0 || *alpha == 0.f) {
            return;
        }

        /*     Set the start point in X if the increment is not unity. */

        if (*incx <= 0) {
            kx = 1 - (*n - 1) * *incx;
        } else if (*incx != 1) {
            kx = 1;
        }

        /*     Start the operations. In this version the elements of A are */
        /*     accessed sequentially with one pass through the triangular part */
        /*     of A. */

        if ( strncasecmp(uplo, "U", 1) == 0 ) {

            /*        Form  A  when A is stored in upper triangle. */

            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.f) {
                        temp = *alpha * x[j];
                        i__2 = j;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            a[i__ + j * a_dim1] += x[i__] * temp;
                            /* L10: */
                        }
                    }
                    /* L20: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.f) {
                        temp = *alpha * x[jx];
                        ix = kx;
                        i__2 = j;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            a[i__ + j * a_dim1] += x[ix] * temp;
                            ix += *incx;
                            /* L30: */
                        }
                    }
                    jx += *incx;
                    /* L40: */
                }
            }
        } else {

            /*        Form  A  when A is stored in lower triangle. */

            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.f) {
                        temp = *alpha * x[j];
                        i__2 = *n;
                        for (i__ = j; i__ <= i__2; ++i__) {
                            a[i__ + j * a_dim1] += x[i__] * temp;
                            /* L50: */
                        }
                    }
                    /* L60: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.f) {
                        temp = *alpha * x[jx];
                        ix = jx;
                        i__2 = *n;
                        for (i__ = j; i__ <= i__2; ++i__) {
                            a[i__ + j * a_dim1] += x[ix] * temp;
                            ix += *incx;
                            /* L70: */
                        }
                    }
                    jx += *incx;
                    /* L80: */
                }
            }
        }

        return;

        /*     End of HSYR  . */

    } /* hsyr_ */

void lpf_blas_hsyr_fortran(char *uplo, lpf_blas_int_t *n, lpf_ffloat16_t *alpha, lpf_ffloat16_t *x,
        lpf_blas_int_t *incx, lpf_ffloat16_t *a, lpf_blas_int_t *lda)
{
    LPF_GLOBAL(hsyr,HSYR)(uplo, n, (lpf_float16_t *)alpha, (lpf_float16_t *)x,
        incx, (lpf_float16_t *)a, lda, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hsyr_fortran_dyn_rank(char *uplo, lpf_blas_int_t *n, lpf_ffloat16_t *alpha, CFI_cdesc_t *x,
        lpf_blas_int_t *incx, CFI_cdesc_t *a, lpf_blas_int_t *lda)
{
    lpf_float16_t *x_ptr = x->base_addr;
    lpf_float16_t *a_ptr = a->base_addr;

    LPF_GLOBAL(hsyr,HSYR)(uplo, n, (lpf_float16_t *)alpha, (lpf_float16_t *)x_ptr,
        incx, (lpf_float16_t *)a_ptr, lda, 1);
}
