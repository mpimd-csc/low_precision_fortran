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

/* > \brief \b BSYMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > BSYMV  performs the matrix-vector  operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           upper triangular part of A is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           LPF_MAX( 1, n ). */
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
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. On exit, Y is overwritten by the updated */
/* >           vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
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
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(bsymv,BSYMV)(char *uplo, lpf_blas_int_t *n, lpf_bfloat16_t *alpha, lpf_bfloat16_t *a,
        lpf_blas_int_t *lda, lpf_bfloat16_t *x, lpf_blas_int_t *incx, lpf_bfloat16_t *beta, lpf_bfloat16_t *y, lpf_blas_int_t *
        incy, lpf_fortran_strlen_t uplo_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    lpf_blas_int_t i__, j, ix, iy, jx, jy, kx, ky, info;
    lpf_bfloat16_t temp1, temp2;


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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;
    (void) uplo_len;

    /* Function Body */
    info = 0;
    if (! ( strncasecmp(uplo, "U", 1 ) == 0) && ! (strncasecmp(uplo, "L", 1) == 0 )) {
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*lda < LPF_MAX(1,*n)) {
        info = 5;
    } else if (*incx == 0) {
        info = 7;
    } else if (*incy == 0) {
        info = 10;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSYMV ", &info, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0.f && *beta == 1.f)) {
        return;
    }

    /*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (*n - 1) * *incy;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through the triangular part */
    /*     of A. */

    /*     First form  y := beta*y. */

    if (*beta != 1.f) {
        if (*incy == 1) {
            if (*beta == 0.f) {
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = 0.f;
                    /* L10: */
                }
            } else {
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = *beta * y[i__];
                    /* L20: */
                }
            }
        } else {
            iy = ky;
            if (*beta == 0.f) {
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = 0.f;
                    iy += *incy;
                    /* L30: */
                }
            } else {
                i__1 = *n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = *beta * y[iy];
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if (*alpha == 0.f) {
        return;
    }
    if ( strncasecmp(uplo, "U", 1) == 0 ) {

        /*        Form  y  when A is stored in upper triangle. */

        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    y[i__] += temp1 * a[i__ + j * a_dim1];
                    temp2 += a[i__ + j * a_dim1] * x[i__];
                    /* L50: */
                }
                y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
                /* L60: */
            }
        } else {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[jx];
                temp2 = 0.f;
                ix = kx;
                iy = ky;
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    y[iy] += temp1 * a[i__ + j * a_dim1];
                    temp2 += a[i__ + j * a_dim1] * x[ix];
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
                jx += *incx;
                jy += *incy;
                /* L80: */
            }
        }
    } else {

        /*        Form  y  when A is stored in lower triangle. */

        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                y[j] += temp1 * a[j + j * a_dim1];
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    y[i__] += temp1 * a[i__ + j * a_dim1];
                    temp2 += a[i__ + j * a_dim1] * x[i__];
                    /* L90: */
                }
                y[j] += *alpha * temp2;
                /* L100: */
            }
        } else {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[jx];
                temp2 = 0.f;
                y[jy] += temp1 * a[j + j * a_dim1];
                ix = jx;
                iy = jy;
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    ix += *incx;
                    iy += *incy;
                    y[iy] += temp1 * a[i__ + j * a_dim1];
                    temp2 += a[i__ + j * a_dim1] * x[ix];
                    /* L110: */
                }
                y[jy] += *alpha * temp2;
                jx += *incx;
                jy += *incy;
                /* L120: */
            }
        }
    }

    return;

    /*     End of BSYMV . */

} /* bsymv_ */

void lpf_blas_bsymv_fortran(char *uplo, lpf_blas_int_t *n, lpf_fbfloat16_t *alpha, lpf_fbfloat16_t *a,
        lpf_blas_int_t *lda, lpf_fbfloat16_t *x, lpf_blas_int_t *incx, lpf_fbfloat16_t *beta, lpf_fbfloat16_t *y, lpf_blas_int_t *
        incy)
{
    LPF_GLOBAL(bsymv,BSYMV)(uplo, n, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)a,
        lda, (lpf_bfloat16_t *)x, incx, (lpf_bfloat16_t *)beta, (lpf_bfloat16_t *)y,
        incy, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_bsymv_fortran_dyn_rank(char *uplo, lpf_blas_int_t *n, lpf_fbfloat16_t *alpha, CFI_cdesc_t *_a,
        lpf_blas_int_t *lda, CFI_cdesc_t *_x, lpf_blas_int_t *incx, lpf_fbfloat16_t *beta, CFI_cdesc_t *_y, lpf_blas_int_t *
        incy)
{
    lpf_bfloat16_t *a = _a->base_addr;
    lpf_bfloat16_t *x = _x->base_addr;
    lpf_bfloat16_t *y = _y->base_addr;
    LPF_GLOBAL(bsymv,BSYMV)(uplo, n, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)a,
        lda, (lpf_bfloat16_t *)x, incx, (lpf_bfloat16_t *)beta, (lpf_bfloat16_t *)y,
        incy, 1);
}
