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

/* > \brief \b BTPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL AP(*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > BTPMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix, supplied in packed form. */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular matrix packed sequentially, */
/* >           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/* >           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/* >           respectively, and so on. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular matrix packed sequentially, */
/* >           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/* >           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/* >           respectively, and so on. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced, but are assumed to be unity. */
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
void LPF_GLOBAL(btpmv,BTPMV)(char *uplo, char *trans, char *diag, lpf_blas_int_t *n,
        lpf_bfloat16_t *ap, lpf_bfloat16_t *x, lpf_blas_int_t *incx, lpf_fortran_strlen_t uplo_len, lpf_fortran_strlen_t trans_len,
        lpf_fortran_strlen_t diag_len)
{
    /* System generated locals */
    lpf_blas_int_t i__1, i__2;

    /* Local variables */
    lpf_blas_int_t i__, j, k, kk, ix, jx, kx, info;
    lpf_bfloat16_t temp;
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

    /*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    (void) uplo_len;
    (void) trans_len;
    (void) diag_len;
    /* Function Body */
    info = 0;
    if (! (strncasecmp(uplo, "U", 1) == 0) && ! (strncasecmp(uplo, "L", 1) == 0 )){
        info = 1;
    } else if (! (strncasecmp(trans, "N", 1) == 0) && ! (strncasecmp(trans, "T", 1) == 0 )
                && ! (strncasecmp(trans, "C", 1) == 0)) {
        info = 2;
    } else if (! (strncasecmp(diag, "U", 1) == 0 ) && ! (strncasecmp(diag, "N", 1) == 0 )){
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*incx == 0) {
        info = 7;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BTPMV ", &info, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*n == 0) {
        return;
    }

    nounit = strncasecmp(diag, "N", 1) == 0 ;

    /*     Set up the start point in X if the increment is not unity. This */
    /*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }

    /*     Start the operations. In this version the elements of AP are */
    /*     accessed sequentially with one pass through AP. */

    if ( strncasecmp(trans, "N", 1) == 0) {

        /*        Form  x:= A*x. */

        if ( strncasecmp(uplo, "U", 1 ) == 0 ) {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.f) {
                        temp = x[j];
                        k = kk;
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[i__] += temp * ap[k];
                            ++k;
                            /* L10: */
                        }
                        if (nounit) {
                            x[j] *= ap[kk + j - 1];
                        }
                    }
                    kk += j;
                    /* L20: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.f) {
                        temp = x[jx];
                        ix = kx;
                        i__2 = kk + j - 2;
                        for (k = kk; k <= i__2; ++k) {
                            x[ix] += temp * ap[k];
                            ix += *incx;
                            /* L30: */
                        }
                        if (nounit) {
                            x[jx] *= ap[kk + j - 1];
                        }
                    }
                    jx += *incx;
                    kk += j;
                    /* L40: */
                }
            }
        } else {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.f) {
                        temp = x[j];
                        k = kk;
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[i__] += temp * ap[k];
                            --k;
                            /* L50: */
                        }
                        if (nounit) {
                            x[j] *= ap[kk - *n + j];
                        }
                    }
                    kk -= *n - j + 1;
                    /* L60: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    if (x[jx] != 0.f) {
                        temp = x[jx];
                        ix = kx;
                        i__1 = kk - (*n - (j + 1));
                        for (k = kk; k >= i__1; --k) {
                            x[ix] += temp * ap[k];
                            ix -= *incx;
                            /* L70: */
                        }
                        if (nounit) {
                            x[jx] *= ap[kk - *n + j];
                        }
                    }
                    jx -= *incx;
                    kk -= *n - j + 1;
                    /* L80: */
                }
            }
        }
    } else {

        /*        Form  x := A**T*x. */

        if ( strncasecmp(uplo, "U", 1) == 0) {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= ap[kk];
                    }
                    k = kk - 1;
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        temp += ap[k] * x[i__];
                        --k;
                        /* L90: */
                    }
                    x[j] = temp;
                    kk -= j;
                    /* L100: */
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= ap[kk];
                    }
                    i__1 = kk - j + 1;
                    for (k = kk - 1; k >= i__1; --k) {
                        ix -= *incx;
                        temp += ap[k] * x[ix];
                        /* L110: */
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    kk -= j;
                    /* L120: */
                }
            }
        } else {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= ap[kk];
                    }
                    k = kk + 1;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        temp += ap[k] * x[i__];
                        ++k;
                        /* L130: */
                    }
                    x[j] = temp;
                    kk += *n - j + 1;
                    /* L140: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= ap[kk];
                    }
                    i__2 = kk + *n - j;
                    for (k = kk + 1; k <= i__2; ++k) {
                        ix += *incx;
                        temp += ap[k] * x[ix];
                        /* L150: */
                    }
                    x[jx] = temp;
                    jx += *incx;
                    kk += *n - j + 1;
                    /* L160: */
                }
            }
        }
    }

    return;

    /*     End of BTPMV . */

} /* btpmv_ */

void lpf_blas_btpmv_fortran(char *uplo, char *trans, char *diag, lpf_blas_int_t *n,
        lpf_fbfloat16_t *ap, lpf_fbfloat16_t *x, lpf_blas_int_t *incx)
{
    LPF_GLOBAL(btpmv,BTPMV)(uplo, trans, diag, n,
        (lpf_bfloat16_t *)ap, (lpf_bfloat16_t *)x, incx, 1, 1,1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_btpmv_fortran_dyn_rank(char *uplo, char *trans, char *diag, lpf_blas_int_t *n,
        CFI_cdesc_t *_ap, CFI_cdesc_t *_x, lpf_blas_int_t *incx)
{
    lpf_bfloat16_t *ap = _ap->base_addr;
    lpf_bfloat16_t *x = _x->base_addr;
    LPF_GLOBAL(btpmv,BTPMV)(uplo, trans, diag, n,
        (lpf_bfloat16_t *)ap, (lpf_bfloat16_t *)x, incx, 1, 1, 1);
}
