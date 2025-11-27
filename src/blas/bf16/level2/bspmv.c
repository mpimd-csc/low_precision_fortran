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

/* > \brief \b BSPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > BSPMV  performs the matrix-vector operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix, supplied in packed form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the matrix A is supplied in the packed */
/* >           array AP as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   The upper triangular part of A is */
/* >                                  supplied in AP. */
/* > */
/* >              UPLO = 'L' or 'l'   The lower triangular part of A is */
/* >                                  supplied in AP. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. */
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
void LPF_GLOBAL(bspmv,BSPMV)(char *uplo, lpf_blas_int_t *n, lpf_bfloat16_t *alpha, lpf_bfloat16_t *ap,
        lpf_bfloat16_t *x, lpf_blas_int_t *incx, lpf_bfloat16_t *beta, lpf_bfloat16_t *y, lpf_blas_int_t *incy, lpf_fortran_strlen_t
        uplo_len)
{
    /* System generated locals */
    lpf_blas_int_t i__1, i__2;

    /* Local variables */
    lpf_blas_int_t i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
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

    /*     Test the input parameters. */

    /* Parameter adjustments */
    --y;
    --x;
    --ap;
    (void) uplo_len;

    /* Function Body */
    info = 0;
    if (! (strncasecmp(uplo, "U", 1) == 0 ) && ! (strncasecmp(uplo, "L", 1) == 0) ){
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*incx == 0) {
        info = 6;
    } else if (*incy == 0) {
        info = 9;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSPMV ", &info, (lpf_fortran_strlen_t)6);
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

    /*     Start the operations. In this version the elements of the array AP */
    /*     are accessed sequentially with one pass through AP. */

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
    kk = 1;
    if ( strncasecmp(uplo, "U", 1) == 0 ) {

        /*        Form  y  when AP contains the upper triangle. */

        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                k = kk;
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    y[i__] += temp1 * ap[k];
                    temp2 += ap[k] * x[i__];
                    ++k;
                    /* L50: */
                }
                y[j] = y[j] + temp1 * ap[kk + j - 1] + *alpha * temp2;
                kk += j;
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
                i__2 = kk + j - 2;
                for (k = kk; k <= i__2; ++k) {
                    y[iy] += temp1 * ap[k];
                    temp2 += ap[k] * x[ix];
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                y[jy] = y[jy] + temp1 * ap[kk + j - 1] + *alpha * temp2;
                jx += *incx;
                jy += *incy;
                kk += j;
                /* L80: */
            }
        }
    } else {

        /*        Form  y  when AP contains the lower triangle. */

        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[j];
                temp2 = 0.f;
                y[j] += temp1 * ap[kk];
                k = kk + 1;
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    y[i__] += temp1 * ap[k];
                    temp2 += ap[k] * x[i__];
                    ++k;
                    /* L90: */
                }
                y[j] += *alpha * temp2;
                kk += *n - j + 1;
                /* L100: */
            }
        } else {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp1 = *alpha * x[jx];
                temp2 = 0.f;
                y[jy] += temp1 * ap[kk];
                ix = jx;
                iy = jy;
                i__2 = kk + *n - j;
                for (k = kk + 1; k <= i__2; ++k) {
                    ix += *incx;
                    iy += *incy;
                    y[iy] += temp1 * ap[k];
                    temp2 += ap[k] * x[ix];
                    /* L110: */
                }
                y[jy] += *alpha * temp2;
                jx += *incx;
                jy += *incy;
                kk += *n - j + 1;
                /* L120: */
            }
        }
    }

    return;

    /*     End of BSPMV . */

} /* bspmv_ */

void lpf_blas_bspmv_fortran(char *uplo, lpf_blas_int_t *n, lpf_fbfloat16_t *alpha, lpf_fbfloat16_t *ap,
        lpf_fbfloat16_t *x, lpf_blas_int_t *incx, lpf_fbfloat16_t *beta, lpf_fbfloat16_t *y, lpf_blas_int_t *incy)
{
    LPF_GLOBAL(bspmv,BSPMV)(uplo, n, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)ap,
        (lpf_bfloat16_t *)x, incx, (lpf_bfloat16_t *)beta, (lpf_bfloat16_t *)y, incy, 1);
}

