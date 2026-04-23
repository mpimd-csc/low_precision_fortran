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

/* > \brief \b HGER */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HGER   performs the rank 1 operation */
/* > */
/* >    A := alpha*x*y**T + A, */
/* > */
/* > where alpha is a scalar, x is an m element vector, y is an n element */
/* > vector and A is an m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array of dimension at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the m */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. On exit, A is */
/* >           overwritten by the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           LPF_MAX( 1, m ). */
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
void LPF_GLOBAL(hger,HGER)(lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_float16_t *alpha, lpf_float16_t *x,
        lpf_blas_int_t *incx, lpf_float16_t *y, lpf_blas_int_t *incy, lpf_float16_t *a, lpf_blas_int_t *lda)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    lpf_blas_int_t i__, j, ix, jy, kx, info;
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
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (*m < 0) {
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*incx == 0) {
        info = 5;
    } else if (*incy == 0) {
        info = 7;
    } else if (*lda < LPF_MAX(1,*m)) {
        info = 9;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HGER  ", &info, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.f) {
        return;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */

    if (*incy > 0) {
        jy = 1;
    } else {
        jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (y[jy] != 0.f) {
                temp = *alpha * y[jy];
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    a[i__ + j * a_dim1] += x[i__] * temp;
                    /* L10: */
                }
            }
            jy += *incy;
            /* L20: */
        }
    } else {
        if (*incx > 0) {
            kx = 1;
        } else {
            kx = 1 - (*m - 1) * *incx;
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (y[jy] != 0.f) {
                temp = *alpha * y[jy];
                ix = kx;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    a[i__ + j * a_dim1] += x[ix] * temp;
                    ix += *incx;
                    /* L30: */
                }
            }
            jy += *incy;
            /* L40: */
        }
    }

    return;

    /*     End of HGER  . */

} /* hger_ */

void lpf_blas_hger_fortran(lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_ffloat16_t *alpha, lpf_ffloat16_t *x,
        lpf_blas_int_t *incx, lpf_ffloat16_t *y, lpf_blas_int_t *incy, lpf_ffloat16_t *a, lpf_blas_int_t *lda)
{
    LPF_GLOBAL(hger,HGER)(m, n, (lpf_float16_t *)alpha, (lpf_float16_t *)x,
        incx, (lpf_float16_t *)y, incy, (lpf_float16_t *)a, lda);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hger_fortran_dyn_rank(lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_ffloat16_t *alpha, CFI_cdesc_t *x,
        lpf_blas_int_t *incx, CFI_cdesc_t *y, lpf_blas_int_t *incy, CFI_cdesc_t *a, lpf_blas_int_t *lda)
{
    lpf_float16_t *x_ptr = x->base_addr;
    lpf_float16_t *y_ptr = y->base_addr;
    lpf_float16_t *a_ptr = a->base_addr;

    LPF_GLOBAL(hger,HGER)(m, n, (lpf_float16_t *)alpha, (lpf_float16_t *)x_ptr,
        incx, (lpf_float16_t *)y_ptr, incy, (lpf_float16_t *)a_ptr, lda);
}
