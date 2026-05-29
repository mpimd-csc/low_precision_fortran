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

/* > \brief \b HTRMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*) */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HTRMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   x := A*x. */
/* > */
/* >              TRANS = 'T' or 't'   x := A**T*x. */
/* > */
/* >              TRANS = 'C' or 'c'   x := A**T*x. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit */
/* >           triangular as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular matrix and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular matrix and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced either, but are assumed to be unity. */
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
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. On exit, X is overwritten with the */
/* >           tranformed vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
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
void LPF_GLOBAL(htrmv,HTRMV)(char *uplo, char *trans, char *diag, int64_t *n,
        lpf_float16_t *a, int64_t *lda, lpf_float16_t *x, int64_t *incx, lpf_fortran_strlen_t uplo_len,
        lpf_fortran_strlen_t trans_len, lpf_fortran_strlen_t diag_len)
{
    /* System generated locals */
    int64_t a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int64_t i__, j, ix, jx, kx, info;
    lpf_float16_t temp;
    lpf_logical_t nounit;

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
    (void) uplo_len;
    (void) diag_len;
    (void) trans_len;

    /* Function Body */
    info = 0;
    if (! (strncasecmp(uplo, "U", 1) == 0) && ! (strncasecmp(uplo, "L", 1) == 0)) {
        info = 1;
    } else if (! (strncasecmp(trans, "N", 1) == 0) && ! (strncasecmp(trans, "T", 1) == 0 )
                && ! (strncasecmp(trans, "C", 1) == 0)) {
        info = 2;
    } else if (! (strncasecmp(diag, "U", 1) == 0) && ! (strncasecmp(diag, "N", 1) == 0)) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*lda < LPF_MAX(1,*n)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    }
    if (info != 0) {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HTRMV ", &infox, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*n == 0) {
        return;
    }

    nounit = strncasecmp(diag, "N", 1) == 0;

    /*     Set up the start point in X if the increment is not unity. This */
    /*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */

    if (strncasecmp(trans, "N", 1) == 0) {

        /*        Form  x := A*x. */

        if (strncasecmp(uplo, "U", 1) == 0 ) {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.f) {
                        temp = x[j];
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[i__] += temp * a[i__ + j * a_dim1];
                            /* L10: */
                        }
                        if (nounit) {
                            x[j] *= a[j + j * a_dim1];
                        }
                    }
                    /* L20: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.f) {
                        temp = x[jx];
                        ix = kx;
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[ix] += temp * a[i__ + j * a_dim1];
                            ix += *incx;
                            /* L30: */
                        }
                        if (nounit) {
                            x[jx] *= a[j + j * a_dim1];
                        }
                    }
                    jx += *incx;
                    /* L40: */
                }
            }
        } else {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.f) {
                        temp = x[j];
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[i__] += temp * a[i__ + j * a_dim1];
                            /* L50: */
                        }
                        if (nounit) {
                            x[j] *= a[j + j * a_dim1];
                        }
                    }
                    /* L60: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    if (x[jx] != 0.f) {
                        temp = x[jx];
                        ix = kx;
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[ix] += temp * a[i__ + j * a_dim1];
                            ix -= *incx;
                            /* L70: */
                        }
                        if (nounit) {
                            x[jx] *= a[j + j * a_dim1];
                        }
                    }
                    jx -= *incx;
                    /* L80: */
                }
            }
        }
    } else {

        /*        Form  x := A**T*x. */

        if (strncasecmp(uplo, "U", 1) == 0) {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        temp += a[i__ + j * a_dim1] * x[i__];
                        /* L90: */
                    }
                    x[j] = temp;
                    /* L100: */
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        ix -= *incx;
                        temp += a[i__ + j * a_dim1] * x[ix];
                        /* L110: */
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    /* L120: */
                }
            }
        } else {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        temp += a[i__ + j * a_dim1] * x[i__];
                        /* L130: */
                    }
                    x[j] = temp;
                    /* L140: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        ix += *incx;
                        temp += a[i__ + j * a_dim1] * x[ix];
                        /* L150: */
                    }
                    x[jx] = temp;
                    jx += *incx;
                    /* L160: */
                }
            }
        }
    }

    return;

    /*     End of HTRMV . */

} /* htrmv_ */

#include <ISO_Fortran_binding.h>

void lpf_blas_htrmv_fortran_dyn_rank_64(char *uplo, char *trans, char *diag, int64_t *n,
        CFI_cdesc_t *a, int64_t *lda, CFI_cdesc_t *x, int64_t *incx)
{
    lpf_float16_t *a_ptr = a->base_addr;
    lpf_float16_t *x_ptr = x->base_addr;

    LPF_GLOBAL(htrmv,HTRMV)(uplo, trans, diag, n,
        (lpf_float16_t *)a_ptr, lda, (lpf_float16_t *)x_ptr, incx, 1, 1, 1);
}

void lpf_blas_htrmv_fortran_dyn_rank_32(char *uplo, char *trans, char *diag, int32_t *n,
        CFI_cdesc_t *a, int32_t *lda, CFI_cdesc_t *x, int32_t *incx)
{
    lpf_float16_t *a_ptr = a->base_addr;
    lpf_float16_t *x_ptr = x->base_addr;
    int64_t _n = *n;
    int64_t _lda = *lda;
    int64_t _incx = *incx;

    LPF_GLOBAL(htrmv,HTRMV)(uplo, trans, diag, &_n,
        (lpf_float16_t *)a_ptr, &_lda, (lpf_float16_t *)x_ptr, &_incx, 1, 1, 1);
}
