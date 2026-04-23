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

/* > \brief \b BTBSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,K,LDA,N */
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
/* > BTBSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular band matrix, with ( k + 1 ) */
/* > diagonals. */
/* > */
/* > No test for singularity or near-singularity is included in this */
/* > routine. Such tests must be performed before calling this routine. */
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
/* >           On entry, TRANS specifies the equations to be solved as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   A*x = b. */
/* > */
/* >              TRANS = 'T' or 't'   A**T*x = b. */
/* > */
/* >              TRANS = 'C' or 'c'   A**T*x = b. */
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
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with UPLO = 'U' or 'u', K specifies the number of */
/* >           super-diagonals of the matrix A. */
/* >           On entry with UPLO = 'L' or 'l', K specifies the number of */
/* >           sub-diagonals of the matrix A. */
/* >           K must satisfy  0 .le. K. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
/* >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the upper triangular */
/* >           band part of the matrix of coefficients, supplied column by */
/* >           column, with the leading diagonal of the matrix in row */
/* >           ( k + 1 ) of the array, the first super-diagonal starting at */
/* >           position 2 in row k, and so on. The top left k by k triangle */
/* >           of the array A is not referenced. */
/* >           The following program segment will transfer an upper */
/* >           triangular band matrix from conventional full matrix storage */
/* >           to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = K + 1 - J */
/* >                    DO 10, I = MAX( 1, J - K ), J */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the lower triangular */
/* >           band part of the matrix of coefficients, supplied column by */
/* >           column, with the leading diagonal of the matrix in row 1 of */
/* >           the array, the first sub-diagonal starting at position 1 in */
/* >           row 2, and so on. The bottom right k by k triangle of the */
/* >           array A is not referenced. */
/* >           The following program segment will transfer a lower */
/* >           triangular band matrix from conventional full matrix storage */
/* >           to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = 1 - J */
/* >                    DO 10, I = J, MIN( N, J + K ) */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Note that when DIAG = 'U' or 'u' the elements of the array A */
/* >           corresponding to the diagonal elements of the matrix are not */
/* >           referenced, but are assumed to be unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           ( k + 1 ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is REAL array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element right-hand side vector b. On exit, X is overwritten */
/* >           with the solution vector x. */
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
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(btbsv,BTBSV)(char *uplo, char *trans, char *diag, lpf_blas_int_t *n,
        lpf_blas_int_t *k, lpf_bfloat16_t *a, lpf_blas_int_t *lda, lpf_bfloat16_t *x, lpf_blas_int_t *incx, lpf_fortran_strlen_t
        uplo_len, lpf_fortran_strlen_t trans_len, lpf_fortran_strlen_t diag_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    lpf_blas_int_t i__, j, l, ix, jx, kx, info;
    lpf_bfloat16_t temp;
    lpf_blas_int_t kplus1;
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
    (void) trans_len;
    (void) diag_len;

    /* Function Body */
    info = 0;
    if (! (strncasecmp(uplo, "U", 1) == 0 ) && ! (strncasecmp(uplo, "L", 1) == 0 )) {
        info = 1;
    } else if (! ( strncasecmp(trans, "N", 1) == 0) && ! (strncasecmp(trans, "T", 1) == 0 )
               && ! (strncasecmp(trans, "C", 1) == 0)) {
        info = 2;
    } else if (! (strncasecmp(diag, "U", 1) == 0) && ! (strncasecmp(diag, "N", 1) == 0 )) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*k < 0) {
        info = 5;
    } else if (*lda < *k + 1) {
        info = 7;
    } else if (*incx == 0) {
        info = 9;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BTBSV ", &info, (lpf_fortran_strlen_t)6);
        return ;
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
    /*     accessed by sequentially with one pass through A. */

    if (strncasecmp(trans, "N", 1) == 0 ) {

        /*        Form  x := inv( A )*x. */

        if ( strncasecmp(uplo, "U", 1) == 0) {
            kplus1 = *k + 1;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.f) {
                        l = kplus1 - j;
                        if (nounit) {
                            x[j] /= a[kplus1 + j * a_dim1];
                        }
                        temp = x[j];
                        /* Computing MAX */
                        i__2 = 1, i__3 = j - *k;
                        i__1 = LPF_MAX(i__2,i__3);
                        for (i__ = j - 1; i__ >= i__1; --i__) {
                            x[i__] -= temp * a[l + i__ + j * a_dim1];
                            /* L10: */
                        }
                    }
                    /* L20: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    kx -= *incx;
                    if (x[jx] != 0.f) {
                        ix = kx;
                        l = kplus1 - j;
                        if (nounit) {
                            x[jx] /= a[kplus1 + j * a_dim1];
                        }
                        temp = x[jx];
                        /* Computing MAX */
                        i__2 = 1, i__3 = j - *k;
                        i__1 = LPF_MAX(i__2,i__3);
                        for (i__ = j - 1; i__ >= i__1; --i__) {
                            x[ix] -= temp * a[l + i__ + j * a_dim1];
                            ix -= *incx;
                            /* L30: */
                        }
                    }
                    jx -= *incx;
                    /* L40: */
                }
            }
        } else {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.f) {
                        l = 1 - j;
                        if (nounit) {
                            x[j] /= a[j * a_dim1 + 1];
                        }
                        temp = x[j];
                        /* Computing MIN */
                        i__3 = *n, i__4 = j + *k;
                        i__2 = LPF_MIN(i__3,i__4);
                        for (i__ = j + 1; i__ <= i__2; ++i__) {
                            x[i__] -= temp * a[l + i__ + j * a_dim1];
                            /* L50: */
                        }
                    }
                    /* L60: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    kx += *incx;
                    if (x[jx] != 0.f) {
                        ix = kx;
                        l = 1 - j;
                        if (nounit) {
                            x[jx] /= a[j * a_dim1 + 1];
                        }
                        temp = x[jx];
                        /* Computing MIN */
                        i__3 = *n, i__4 = j + *k;
                        i__2 = LPF_MIN(i__3,i__4);
                        for (i__ = j + 1; i__ <= i__2; ++i__) {
                            x[ix] -= temp * a[l + i__ + j * a_dim1];
                            ix += *incx;
                            /* L70: */
                        }
                    }
                    jx += *incx;
                    /* L80: */
                }
            }
        }
    } else {

        /*        Form  x := inv( A**T)*x. */

        if ( strncasecmp(uplo, "U", 1) == 0 ) {
            kplus1 = *k + 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    l = kplus1 - j;
                    /* Computing MAX */
                    i__2 = 1, i__3 = j - *k;
                    i__4 = j - 1;
                    for (i__ = LPF_MAX(i__2,i__3); i__ <= i__4; ++i__) {
                        temp -= a[l + i__ + j * a_dim1] * x[i__];
                        /* L90: */
                    }
                    if (nounit) {
                        temp /= a[kplus1 + j * a_dim1];
                    }
                    x[j] = temp;
                    /* L100: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = kx;
                    l = kplus1 - j;
                    /* Computing MAX */
                    i__4 = 1, i__2 = j - *k;
                    i__3 = j - 1;
                    for (i__ = LPF_MAX(i__4,i__2); i__ <= i__3; ++i__) {
                        temp -= a[l + i__ + j * a_dim1] * x[ix];
                        ix += *incx;
                        /* L110: */
                    }
                    if (nounit) {
                        temp /= a[kplus1 + j * a_dim1];
                    }
                    x[jx] = temp;
                    jx += *incx;
                    if (j > *k) {
                        kx += *incx;
                    }
                    /* L120: */
                }
            }
        } else {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    l = 1 - j;
                    /* Computing MIN */
                    i__1 = *n, i__3 = j + *k;
                    i__4 = j + 1;
                    for (i__ = LPF_MIN(i__1,i__3); i__ >= i__4; --i__) {
                        temp -= a[l + i__ + j * a_dim1] * x[i__];
                        /* L130: */
                    }
                    if (nounit) {
                        temp /= a[j * a_dim1 + 1];
                    }
                    x[j] = temp;
                    /* L140: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = kx;
                    l = 1 - j;
                    /* Computing MIN */
                    i__4 = *n, i__1 = j + *k;
                    i__3 = j + 1;
                    for (i__ = LPF_MIN(i__4,i__1); i__ >= i__3; --i__) {
                        temp -= a[l + i__ + j * a_dim1] * x[ix];
                        ix -= *incx;
                        /* L150: */
                    }
                    if (nounit) {
                        temp /= a[j * a_dim1 + 1];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    if (*n - j >= *k) {
                        kx -= *incx;
                    }
                    /* L160: */
                }
            }
        }
    }

    return;

    /*     End of BTBSV . */

} /* btbsv_ */

void lpf_blas_btbsv_fortran(char *uplo, char *trans, char *diag, lpf_blas_int_t *n,
        lpf_blas_int_t *k, lpf_fbfloat16_t *a, lpf_blas_int_t *lda, lpf_fbfloat16_t *x, lpf_blas_int_t *incx)
{
    LPF_GLOBAL(btbsv,BTBSV)(uplo, trans, diag, n,
        k, (lpf_bfloat16_t *)a, lda, (lpf_bfloat16_t *)x, incx, 1, 1, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_btbsv_fortran_dyn_rank(char *uplo, char *trans, char *diag, lpf_blas_int_t *n,
        lpf_blas_int_t *k, CFI_cdesc_t *_a, lpf_blas_int_t *lda, CFI_cdesc_t *_x, lpf_blas_int_t *incx)
{
    lpf_bfloat16_t *a = _a->base_addr;
    lpf_bfloat16_t *x = _x->base_addr;
    LPF_GLOBAL(btbsv,BTBSV)(uplo, trans, diag, n,
        k, (lpf_bfloat16_t *)a, lda, (lpf_bfloat16_t *)x, incx, 1, 1, 1);
}
