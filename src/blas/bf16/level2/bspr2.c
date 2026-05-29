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

/* > \brief \b BSPR2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
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
/* > BSPR2  performs the symmetric rank 2 operation */
/* > */
/* >    A := alpha*x*y**T + alpha*y*x**T + A, */
/* > */
/* > where alpha is a scalar, x and y are n element vectors and A is an */
/* > n by n symmetric matrix, supplied in packed form. */
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
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the upper triangular part of the */
/* >           updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the lower triangular part of the */
/* >           updated matrix. */
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
void LPF_GLOBAL(bspr2,BSPR2)(char *uplo, int64_t *n, lpf_bfloat16_t *alpha, lpf_bfloat16_t *x,
        int64_t *incx, lpf_bfloat16_t *y, int64_t *incy, lpf_bfloat16_t *ap, lpf_fortran_strlen_t uplo_len)
{
    /* System generated locals */
    int64_t i__1, i__2;

    /* Local variables */
    int64_t i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
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
    --ap;
    --y;
    --x;
    (void) uplo_len;
    /* Function Body */
    info = 0;
    if (! (strncasecmp(uplo, "U", 1) == 0 ) && ! (strncasecmp(uplo, "L", 1) == 0)) {
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*incx == 0) {
        info = 5;
    } else if (*incy == 0) {
        info = 7;
    }
    if (info != 0) {
        int32_t infox = (int32_t)info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSPR2 ", &infox, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.f) {
        return;
    }

    /*     Set up the start points in X and Y if the increments are not both */
    /*     unity. */

    if (*incx != 1 || *incy != 1) {
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
        jx = kx;
        jy = ky;
    }

    /*     Start the operations. In this version the elements of the array AP */
    /*     are accessed sequentially with one pass through AP. */

    kk = 1;
    if (strncasecmp(uplo, "U", 1) == 0 ) {

        /*        Form  A  when upper triangle is stored in AP. */

        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[j] != 0.f || y[j] != 0.f) {
                    temp1 = *alpha * y[j];
                    temp2 = *alpha * x[j];
                    k = kk;
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
                        ++k;
                        /* L10: */
                    }
                }
                kk += j;
                /* L20: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[jx] != 0.f || y[jy] != 0.f) {
                    temp1 = *alpha * y[jy];
                    temp2 = *alpha * x[jx];
                    ix = kx;
                    iy = ky;
                    i__2 = kk + j - 1;
                    for (k = kk; k <= i__2; ++k) {
                        ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
                        ix += *incx;
                        iy += *incy;
                        /* L30: */
                    }
                }
                jx += *incx;
                jy += *incy;
                kk += j;
                /* L40: */
            }
        }
    } else {

        /*        Form  A  when lower triangle is stored in AP. */

        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[j] != 0.f || y[j] != 0.f) {
                    temp1 = *alpha * y[j];
                    temp2 = *alpha * x[j];
                    k = kk;
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
                        ++k;
                        /* L50: */
                    }
                }
                kk = kk + *n - j + 1;
                /* L60: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[jx] != 0.f || y[jy] != 0.f) {
                    temp1 = *alpha * y[jy];
                    temp2 = *alpha * x[jx];
                    ix = jx;
                    iy = jy;
                    i__2 = kk + *n - j;
                    for (k = kk; k <= i__2; ++k) {
                        ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
                        ix += *incx;
                        iy += *incy;
                        /* L70: */
                    }
                }
                jx += *incx;
                jy += *incy;
                kk = kk + *n - j + 1;
                /* L80: */
            }
        }
    }

    return;

    /*     End of BSPR2 . */

} /* bspr2_ */

void lpf_blas_bspr2_fortran_dyn_rank_64(char *uplo, int64_t *n, lpf_fbfloat16_t *alpha, CFI_cdesc_t *_x,
        int64_t *incx, CFI_cdesc_t *_y, int64_t *incy, CFI_cdesc_t *_ap)
{
    lpf_bfloat16_t *x = _x->base_addr;
    lpf_bfloat16_t *y = _y->base_addr;
    lpf_bfloat16_t *ap = _ap->base_addr;
    LPF_GLOBAL(bspr2,BSPR2)(uplo, n, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)x,
        incx, (lpf_bfloat16_t *)y, incy, (lpf_bfloat16_t *)ap, 1);
}

void lpf_blas_bspr2_fortran_dyn_rank_32(char *uplo, int32_t *n, lpf_fbfloat16_t *alpha, CFI_cdesc_t *_x,
        int32_t *incx, CFI_cdesc_t *_y, int32_t *incy, CFI_cdesc_t *_ap)
{
    lpf_bfloat16_t *x = _x->base_addr;
    lpf_bfloat16_t *y = _y->base_addr;
    lpf_bfloat16_t *ap = _ap->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(bspr2,BSPR2)(uplo, &_n, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)x,
        &_incx, (lpf_bfloat16_t *)y, &_incy, (lpf_bfloat16_t *)ap, 1);
}

