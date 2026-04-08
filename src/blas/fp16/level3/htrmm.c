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

#ifdef BLAS_IS_MKL
#include <mkl.h>
#endif


#include <string.h>
/* > \brief \b HTRMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER LDA,LDB,M,N */
/*       CHARACTER DIAG,SIDE,TRANSA,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HTRMM  performs one of the matrix-matrix operations */
/* > */
/* >    B := alpha*op( A )*B,   or   B := alpha*B*op( A ), */
/* > */
/* > where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE specifies whether  op( A ) multiplies B from */
/* >           the left or right as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */
/* > */
/* >              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix A is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n'   op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't'   op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c'   op( A ) = A**T. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit triangular */
/* >           as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of B. M must be at */
/* >           least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of B.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, k ), where k is m */
/* >           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
/* >           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
/* >           upper triangular part of the array  A must contain the upper */
/* >           triangular matrix  and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
/* >           lower triangular part of the array  A must contain the lower */
/* >           triangular matrix  and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
/* >           A  are not referenced either,  but are assumed to be  unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/* >           LDA  must be at least  LPF_MAX( 1, m ),  when  SIDE = 'R' or 'r' */
/* >           then LDA must be at least LPF_MAX( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array of DIMENSION ( LDB, n ). */
/* >           Before entry,  the leading  m by n part of the array  B must */
/* >           contain the matrix  B,  and  on exit  is overwritten  by the */
/* >           transformed matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           LPF_MAX( 1, m ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level3 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 3 Blas routine. */
/* > */
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(htrmm,HTRMM)(char *side, char *uplo, char *transa, char *diag,
        lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_float16_t *alpha, lpf_float16_t *a, lpf_blas_int_t *lda, lpf_float16_t *b,
        lpf_blas_int_t *ldb, lpf_fortran_strlen_t side_len, lpf_fortran_strlen_t uplo_len, lpf_fortran_strlen_t transa_len,
        lpf_fortran_strlen_t diag_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    lpf_blas_int_t i__, j, k, info;
    lpf_float16_t temp;
    lpf_logical_t lside;
    lpf_blas_int_t nrowa;
    lpf_logical_t upper;
    lpf_logical_t nounit;


    /*  -- Reference BLAS level3 routine (version 3.4.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Parameters .. */
    /*     .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    (void) diag_len;
    (void) side_len;
    (void) uplo_len;
    (void) transa_len;

    /* Function Body */
    lside = strncasecmp(side, "L", 1) == 0;
    if (lside) {
        nrowa = *m;
    } else {
        nrowa = *n;
    }
    nounit = strncasecmp(diag, "N", 1 ) == 0;
    upper = strncasecmp(uplo, "U", 1 ) == 0;

    info = 0;
    if (! lside && ! (strncasecmp(side, "R", 1) == 0)) {
        info = 1;
    } else if (! upper && ! (strncasecmp(uplo, "L", 1) == 0 )) {
        info = 2;
    } else if (! (strncasecmp(transa, "N", 1) == 0 ) && ! (strncasecmp(transa, "T", 1) == 0) && ! (strncasecmp(transa, "C", 1) == 0 )) {
        info = 3;
    } else if (! (strncasecmp(diag, "U", 1) == 0 ) && ! (strncasecmp(diag, "N", 1) == 0 )) {
        info = 4;
    } else if (*m < 0) {
        info = 5;
    } else if (*n < 0) {
        info = 6;
    } else if (*lda < LPF_MAX(1,nrowa)) {
        info = 9;
    } else if (*ldb < LPF_MAX(1,*m)) {
        info = 11;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HTRMM ", &info, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
        return;
    }

    /*     And when  alpha.eq.zero. */

    if (*alpha == 0.f) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] = 0.f;
                /* L10: */
            }
            /* L20: */
        }
        return;
    }


    /*     Start the operations. */

    if (lside) {
        if (strncasecmp(transa, "N", 1) == 0) {

            /*           Form  B := alpha*A*B. */

            if (upper) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (k = 1; k <= i__2; ++k) {
                        if (b[k + j * b_dim1] != 0.f) {
                            temp = *alpha * b[k + j * b_dim1];
                            i__3 = k - 1;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * a[i__ + k *
                                    a_dim1];
                                /* L30: */
                            }
                            if (nounit) {
                                temp *= a[k + k * a_dim1];
                            }
                            b[k + j * b_dim1] = temp;
                        }
                        /* L40: */
                    }
                    /* L50: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (k = *m; k >= 1; --k) {
                        if (b[k + j * b_dim1] != 0.f) {
                            temp = *alpha * b[k + j * b_dim1];
                            b[k + j * b_dim1] = temp;
                            if (nounit) {
                                b[k + j * b_dim1] *= a[k + k * a_dim1];
                            }
                            i__2 = *m;
                            for (i__ = k + 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * a[i__ + k *
                                    a_dim1];
                                /* L60: */
                            }
                        }
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        } else {

            /*           Form  B := alpha*A**T*B. */

            if (upper) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (i__ = *m; i__ >= 1; --i__) {
                        temp = b[i__ + j * b_dim1];
                        if (nounit) {
                            temp *= a[i__ + i__ * a_dim1];
                        }
                        i__2 = i__ - 1;
                        for (k = 1; k <= i__2; ++k) {
                            temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
                            /* L90: */
                        }
                        b[i__ + j * b_dim1] = *alpha * temp;
                        /* L100: */
                    }
                    /* L110: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        temp = b[i__ + j * b_dim1];
                        if (nounit) {
                            temp *= a[i__ + i__ * a_dim1];
                        }
                        i__3 = *m;
                        for (k = i__ + 1; k <= i__3; ++k) {
                            temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
                            /* L120: */
                        }
                        b[i__ + j * b_dim1] = *alpha * temp;
                        /* L130: */
                    }
                    /* L140: */
                }
            }
        }
    } else {
        if ( strncasecmp(transa, "N", 1) == 0 ) {

            /*           Form  B := alpha*B*A. */

            if (upper) {
                for (j = *n; j >= 1; --j) {
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__1 = *m;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        /* L150: */
                    }
                    i__1 = j - 1;
                    for (k = 1; k <= i__1; ++k) {
                        if (a[k + j * a_dim1] != 0.f) {
                            temp = *alpha * a[k + j * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k *
                                    b_dim1];
                                /* L160: */
                            }
                        }
                        /* L170: */
                    }
                    /* L180: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        /* L190: */
                    }
                    i__2 = *n;
                    for (k = j + 1; k <= i__2; ++k) {
                        if (a[k + j * a_dim1] != 0.f) {
                            temp = *alpha * a[k + j * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k *
                                    b_dim1];
                                /* L200: */
                            }
                        }
                        /* L210: */
                    }
                    /* L220: */
                }
            }
        } else {

            /*           Form  B := alpha*B*A**T. */

            if (upper) {
                i__1 = *n;
                for (k = 1; k <= i__1; ++k) {
                    i__2 = k - 1;
                    for (j = 1; j <= i__2; ++j) {
                        if (a[j + k * a_dim1] != 0.f) {
                            temp = *alpha * a[j + k * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k *
                                    b_dim1];
                                /* L230: */
                            }
                        }
                        /* L240: */
                    }
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[k + k * a_dim1];
                    }
                    if (temp != 1.f) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                            /* L250: */
                        }
                    }
                    /* L260: */
                }
            } else {
                for (k = *n; k >= 1; --k) {
                    i__1 = *n;
                    for (j = k + 1; j <= i__1; ++j) {
                        if (a[j + k * a_dim1] != 0.f) {
                            temp = *alpha * a[j + k * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k *
                                    b_dim1];
                                /* L270: */
                            }
                        }
                        /* L280: */
                    }
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[k + k * a_dim1];
                    }
                    if (temp != 1.f) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                            /* L290: */
                        }
                    }
                    /* L300: */
                }
            }
        }
    }

    return;

    /*     End of HTRMM . */

} /* htrmm_ */

void lpf_blas_htrmm_fortran(char *side, char *uplo, char *transa, char *diag,
        lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_ffloat16_t *alpha, lpf_ffloat16_t *a, lpf_blas_int_t *lda, lpf_ffloat16_t *b,
        lpf_blas_int_t *ldb)
{
    LPF_GLOBAL(htrmm,HTRMM)(side, uplo, transa, diag,
        m, n, (lpf_float16_t *)alpha, (lpf_float16_t *)a, lda, (lpf_float16_t *)b,
        ldb, 1, 1, 1, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_htrmm_fortran_dyn_rank(char *side, char *uplo, char *transa, char *diag,
        lpf_blas_int_t *m, lpf_blas_int_t *n, lpf_ffloat16_t *alpha, CFI_cdesc_t *_a, lpf_blas_int_t *lda, CFI_cdesc_t *_b, lpf_blas_int_t *ldb)
{
    lpf_float16_t *a = _a->base_addr;
    lpf_float16_t *b = _b->base_addr;

    LPF_GLOBAL(htrmm,HTRMM)(side, uplo, transa, diag,
        m, n, (lpf_float16_t *)alpha, (lpf_float16_t *)a, lda, (lpf_float16_t *)b,
        ldb, 1, 1, 1, 1);
}
