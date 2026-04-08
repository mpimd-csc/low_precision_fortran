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

/* > \brief \b HGBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,KL,KU,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HGBMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n band matrix, with kl sub-diagonals and ku super-diagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */
/* > */
/* >              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */
/* > */
/* >              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >           On entry, KL specifies the number of sub-diagonals of the */
/* >           matrix A. KL must satisfy  0 .le. KL. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >           On entry, KU specifies the number of super-diagonals of the */
/* >           matrix A. KU must satisfy  0 .le. KU. */
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
/* >           Before entry, the leading ( kl + ku + 1 ) by n part of the */
/* >           array A must contain the matrix of coefficients, supplied */
/* >           column by column, with the leading diagonal of the matrix in */
/* >           row ( ku + 1 ) of the array, the first super-diagonal */
/* >           starting at position 2 in row ku, the first sub-diagonal */
/* >           starting at position 1 in row ( ku + 2 ), and so on. */
/* >           Elements in the array A that do not correspond to elements */
/* >           in the band matrix (such as the top left ku by ku triangle) */
/* >           are not referenced. */
/* >           The following program segment will transfer a band matrix */
/* >           from conventional full matrix storage to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    K = KU + 1 - J */
/* >                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL ) */
/* >                       A( K + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           ( kl + ku + 1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array of DIMENSION at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
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
/* >          Y is REAL array of DIMENSION at least */
/* >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/* >           Before entry, the incremented array Y must contain the */
/* >           vector y. On exit, Y is overwritten by the updated vector y. */
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

/* > \date November 2015 */

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
void LPF_GLOBAL(hgbmv,HGBMV)(char *trans, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_blas_int_t *kl,
        lpf_blas_int_t *ku, lpf_float16_t *alpha, lpf_float16_t *a, lpf_blas_int_t *lda, lpf_float16_t *x, lpf_blas_int_t *
        incx, lpf_float16_t *beta, lpf_float16_t *y, lpf_blas_int_t *incy, lpf_fortran_strlen_t trans_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    lpf_blas_int_t i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
    lpf_float16_t temp;
    lpf_blas_int_t lenx, leny;

    /*  -- Reference BLAS level2 routine (version 3.6.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2015 */

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

    (void) trans_len;

    /* Function Body */
    info = 0;
    if (! (strncasecmp(trans,"N", 1) == 0) && ! (strncasecmp(trans, "T",1) == 0)
        && ! (strncasecmp(trans, "C", 1) == 0)) {
        info = 1;
    } else if (*m < 0) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*kl < 0) {
        info = 4;
    } else if (*ku < 0) {
        info = 5;
    } else if (*lda < *kl + *ku + 1) {
        info = 8;
    } else if (*incx == 0) {
        info = 10;
    } else if (*incy == 0) {
        info = 13;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HGBMV ", &info, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0.f && *beta == 1.f)) {
        return;
    }

    /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
    /*     up the start points in  X  and  Y. */

    if (strncasecmp(trans, "N", 1) == 0) {
        lenx = *n;
        leny = *m;
    } else {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * *incy;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through the band part of A. */

    /*     First form  y := beta*y. */

    if (*beta != 1.f) {
        if (*incy == 1) {
            if (*beta == 0.f) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = 0.f;
                    /* L10: */
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = *beta * y[i__];
                    /* L20: */
                }
            }
        } else {
            iy = ky;
            if (*beta == 0.f) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = 0.f;
                    iy += *incy;
                    /* L30: */
                }
            } else {
                i__1 = leny;
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
    kup1 = *ku + 1;
    if (strncasecmp(trans, "N", 1) == 0 ) {

        /*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (*incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = *alpha * x[jx];
                k = kup1 - j;
                /* Computing MAX */
                i__2 = 1, i__3 = j - *ku;
                /* Computing MIN */
                i__5 = *m, i__6 = j + *kl;
                i__4 = LPF_MIN(i__5,i__6);
                for (i__ = LPF_MAX(i__2,i__3); i__ <= i__4; ++i__) {
                    y[i__] += temp * a[k + i__ + j * a_dim1];
                    /* L50: */
                }
                jx += *incx;
                /* L60: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = *alpha * x[jx];
                iy = ky;
                k = kup1 - j;
                /* Computing MAX */
                i__4 = 1, i__2 = j - *ku;
                /* Computing MIN */
                i__5 = *m, i__6 = j + *kl;
                i__3 = LPF_MIN(i__5,i__6);
                for (i__ = LPF_MAX(i__4,i__2); i__ <= i__3; ++i__) {
                    y[iy] += temp * a[k + i__ + j * a_dim1];
                    iy += *incy;
                    /* L70: */
                }
                jx += *incx;
                if (j > *ku) {
                    ky += *incy;
                }
                /* L80: */
            }
        }
    } else {

        /*        Form  y := alpha*A**T*x + y. */

        jy = ky;
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.f;
                k = kup1 - j;
                /* Computing MAX */
                i__3 = 1, i__4 = j - *ku;
                /* Computing MIN */
                i__5 = *m, i__6 = j + *kl;
                i__2 = LPF_MIN(i__5,i__6);
                for (i__ = LPF_MAX(i__3,i__4); i__ <= i__2; ++i__) {
                    temp += a[k + i__ + j * a_dim1] * x[i__];
                    /* L90: */
                }
                y[jy] += *alpha * temp;
                jy += *incy;
                /* L100: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.f;
                ix = kx;
                k = kup1 - j;
                /* Computing MAX */
                i__2 = 1, i__3 = j - *ku;
                /* Computing MIN */
                i__5 = *m, i__6 = j + *kl;
                i__4 = LPF_MIN(i__5,i__6);
                for (i__ = LPF_MAX(i__2,i__3); i__ <= i__4; ++i__) {
                    temp += a[k + i__ + j * a_dim1] * x[ix];
                    ix += *incx;
                    /* L110: */
                }
                y[jy] += *alpha * temp;
                jy += *incy;
                if (j > *ku) {
                    kx += *incx;
                }
                /* L120: */
            }
        }
    }

    return;

    /*     End of HGBMV . */

} /* hgbmv_ */

void lpf_blas_hgbmv_fortran(char *trans, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_blas_int_t *kl,
        lpf_blas_int_t *ku, lpf_ffloat16_t *alpha, lpf_ffloat16_t *a, lpf_blas_int_t *lda, lpf_ffloat16_t *x, lpf_blas_int_t *
        incx, lpf_ffloat16_t *beta, lpf_ffloat16_t *y, lpf_blas_int_t *incy)
{
    LPF_GLOBAL(hgbmv,HGBMV)(trans, m, n, kl,
        ku, (lpf_float16_t *)alpha, (lpf_float16_t *)a, lda, (lpf_float16_t *)x,
        incx, (lpf_float16_t *)beta, (lpf_float16_t *)y, incy, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hgbmv_fortran_dyn_rank(char *trans, lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_blas_int_t *kl,
        lpf_blas_int_t *ku, lpf_ffloat16_t *alpha, CFI_cdesc_t *a, lpf_blas_int_t *lda, CFI_cdesc_t *x, lpf_blas_int_t *
        incx, lpf_ffloat16_t *beta, CFI_cdesc_t *y, lpf_blas_int_t *incy)
{
    lpf_float16_t *a_ptr = a->base_addr;
    lpf_float16_t *x_ptr = x->base_addr;
    lpf_float16_t *y_ptr = y->base_addr;

    LPF_GLOBAL(hgbmv,HGBMV)(trans, m, n, kl,
        ku, (lpf_float16_t *)alpha, (lpf_float16_t *)a_ptr, lda, (lpf_float16_t *)x_ptr,
        incx, (lpf_float16_t *)beta, (lpf_float16_t *)y_ptr, incy, 1);
}
