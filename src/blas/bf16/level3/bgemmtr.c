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
/* > \brief \b BGEMMTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BGEMMTR(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,BETA, */
/*                         C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,N */
/*       CHARACTER TRANSA,TRANSB, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > BGEMMTR  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*op( A )*op( B ) + beta*C, */
/* > */
/* > where  op( X ) is one of */
/* > */
/* >    op( X ) = X   or   op( X ) = X**T, */
/* > */
/* > alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/* > an n by k matrix,  op( B )  a  k by n matrix and  C an n by n matrix. */
/* > Thereby, the routine only accesses and updates the upper or lower */
/* > triangular part of the result matrix C. This behaviour can be used if */
/* > the resulting matrix C is known to be symmetric. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the lower or the upper */
/* >           triangular part of C is access and updated. */
/* > */
/* >              UPLO = 'L' or 'l', the lower triangular part of C is used. */
/* > */
/* >              UPLO = 'U' or 'u', the upper triangular part of C is used. */
/* > \endverbatim */

/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n',  op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't',  op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c',  op( A ) = A**T. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSB */
/* > \verbatim */
/* >          TRANSB is CHARACTER*1 */
/* >           On entry, TRANSB specifies the form of op( B ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSB = 'N' or 'n',  op( B ) = B. */
/* > */
/* >              TRANSB = 'T' or 't',  op( B ) = B**T. */
/* > */
/* >              TRANSB = 'C' or 'c',  op( B ) = B**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N specifies the number of rows and columns of */
/* >           the matrix C, the number of columns of op(B) and the number */
/* >           of rows of op(A).  N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry,  K  specifies  the number of columns of the matrix */
/* >           op( A ) and the number of rows of the matrix op( B ). K must */
/* >           be at least  zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension ( LDA, ka ), where ka is */
/* >           k  when  TRANSA = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANSA = 'N' or 'n',  the leading  n by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by m  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/* >           LDA must be at least  LPF_MAX( 1, n ), otherwise  LDA must be at */
/* >           least  LPF_MAX( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension ( LDB, kb ), where kb is */
/* >           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/* >           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  n by k  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/* >           LDB must be at least  LPF_MAX( 1, k ), otherwise  LDB must be at */
/* >           least  LPF_MAX( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL. */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension ( LDC, N ) */
/* >           Before entry, the leading  n by n  part of the array  C must */
/* >           contain the matrix  C,  except when  beta  is zero, in which */
/* >           case C need not be set on entry. */
/* >           On exit, the upper or lower triangular part of the matrix */
/* >           C  is overwritten by the n by n matrix */
/* >           ( alpha*op( A )*op( B ) + beta*C ). */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/* >           LPF_MAX( 1, n ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Martin Koehler */

/* > \ingroup gemmtr */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 3 Blas routine. */
/* > */
/* >  -- Written on 19-July-2023. */
/* >     Martin Koehler, MPI Magdeburg */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(bgemmtr,BGEMMTR)(char *uplo, char *transa, char *transb, lpf_blas_int_t
        *n, lpf_blas_int_t *k, lpf_bfloat16_t *alpha, lpf_bfloat16_t *a, lpf_blas_int_t *lda, lpf_bfloat16_t *b, lpf_blas_int_t *
        ldb, lpf_bfloat16_t *beta, lpf_bfloat16_t *c__, lpf_blas_int_t *ldc, lpf_fortran_strlen_t uplo_len, lpf_fortran_strlen_t
        transa_len, lpf_fortran_strlen_t transb_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
                 i__3;

    /* Local variables */
    lpf_blas_int_t i__, j, l, info;
    lpf_logical_t nota, notb;
    lpf_bfloat16_t temp;
    lpf_blas_int_t nrowa, nrowb;
    lpf_logical_t upper;
    lpf_blas_int_t istop;
    lpf_blas_int_t istart;


    /*  -- Reference BLAS level3 routine -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

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

    /*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
    /*     transposed and set  NROWA and NROWB  as the number of rows of  A */
    /*     and  B  respectively. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    (void) uplo_len;
    (void) transa_len;
    (void) transb_len;

    /* Function Body */
    nota = strncasecmp(transa, "N", 1) == 0;
    notb = strncasecmp(transb, "N", 1) == 0;
    if (nota) {
        nrowa = *n;
    } else {
        nrowa = *k;
    }
    if (notb) {
        nrowb = *k;
    } else {
        nrowb = *n;
    }
    upper = strncasecmp(uplo, "U", 1) == 0;

    /*     Test the input parameters. */

    info = 0;
    if (! upper && ! (strncasecmp(uplo, "L", 1) == 0)) {
        info = 1;
    } else if (! nota && ! (strncasecmp(transa, "C", 1) == 0) && ! (strncasecmp(transa, "T", 1) == 0)) {
        info = 2;
    } else if (! notb && ! (strncasecmp(transb, "C", 1) == 0 )&& ! (strncasecmp(transb, "T", 1) == 0 )) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*k < 0) {
        info = 5;
    } else if (*lda < LPF_MAX(1,nrowa)) {
        info = 8;
    } else if (*ldb < LPF_MAX(1,nrowb)) {
        info = 10;
    } else if (*ldc < LPF_MAX(1,*n)) {
        info = 13;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BGEMMTR", &info, (lpf_fortran_strlen_t)7);
        return ;
    }

    /*     Quick return if possible. */

    if (*n == 0) {
        return ;
    }

    /*     And if  alpha.eq.zero. */

    if (*alpha == 0.f) {
        if (*beta == 0.f) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = 0.f;
                    /* L10: */
                }
                /* L20: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    /* L30: */
                }
                /* L40: */
            }
        }
        return;
    }

    /*     Start the operations. */

    if (notb) {
        if (nota) {

            /*           Form  C := alpha*A*B + beta*C. */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = *n;
                }
                if (*beta == 0.f) {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.f;
                        /* L50: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                        /* L60: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    temp = *alpha * b[l + j * b_dim1];
                    i__3 = istop;
                    for (i__ = istart; i__ <= i__3; ++i__) {
                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                        /* L70: */
                    }
                    /* L80: */
                }
                /* L90: */
            }
        } else {

            /*           Form  C := alpha*A**T*B + beta*C */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__) {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        /* L100: */
                    }
                    if (*beta == 0.f) {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    } else {
                        c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
                            i__ + j * c_dim1];
                    }
                    /* L110: */
                }
                /* L120: */
            }
        }
    } else {
        if (nota) {

            /*           Form  C := alpha*A*B**T + beta*C */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = *n;
                }
                if (*beta == 0.f) {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.f;
                        /* L130: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = istop;
                    for (i__ = istart; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                        /* L140: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    temp = *alpha * b[j + l * b_dim1];
                    i__3 = istop;
                    for (i__ = istart; i__ <= i__3; ++i__) {
                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                        /* L150: */
                    }
                    /* L160: */
                }
                /* L170: */
            }
        } else {

            /*           Form  C := alpha*A**T*B**T + beta*C */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = *n;
                }
                i__2 = istop;
                for (i__ = istart; i__ <= i__2; ++i__) {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
                        /* L180: */
                    }
                    if (*beta == 0.f) {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    } else {
                        c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
                            i__ + j * c_dim1];
                    }
                    /* L190: */
                }
                /* L200: */
            }
        }
    }

    return;

    /*     End of BGEMMTR */

} /* bgemmtr_ */

void lpf_blas_bgemmtr_fortran(char *uplo, char *transa, char *transb, lpf_blas_int_t
        *n, lpf_blas_int_t *k, lpf_fbfloat16_t *alpha, lpf_fbfloat16_t *a, lpf_blas_int_t *lda, lpf_fbfloat16_t *b, lpf_blas_int_t *
        ldb, lpf_fbfloat16_t *beta, lpf_fbfloat16_t *c__, lpf_blas_int_t *ldc)
{
    LPF_GLOBAL(bgemmtr,BGEMMTR)(uplo, transa, transb, n, k, (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)a, lda, (lpf_bfloat16_t *)b,
        ldb, (lpf_bfloat16_t *)beta, (lpf_bfloat16_t *)c__, ldc, 1, 1, 1);
}

