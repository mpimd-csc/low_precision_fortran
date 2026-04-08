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

/* > \brief \b BSYMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER LDA,LDB,LDC,M,N */
/*       CHARACTER SIDE,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > BSYMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*A*B + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*B*A + beta*C, */
/* > */
/* > where alpha and beta are scalars,  A is a symmetric matrix and  B and */
/* > C are  m by n matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE  specifies whether  the  symmetric matrix  A */
/* >           appears on the  left or right  in the  operation as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C, */
/* > */
/* >              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C, */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of  the  symmetric  matrix   A  is  to  be */
/* >           referenced as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of the */
/* >                                  symmetric matrix is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of the */
/* >                                  symmetric matrix is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry,  M  specifies the number of rows of the matrix  C. */
/* >           M  must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix C. */
/* >           N  must be at least zero. */
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
/* >          A is REAL array of DIMENSION ( LDA, ka ), where ka is */
/* >           m  when  SIDE = 'L' or 'l'  and is  n otherwise. */
/* >           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of */
/* >           the array  A  must contain the  symmetric matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading m by m upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  symmetric matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  m by m  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  symmetric */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* >           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of */
/* >           the array  A  must contain the  symmetric matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading n by n upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  symmetric matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  n by n  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  symmetric */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/* >           LDA must be at least  LPF_MAX( 1, m ), otherwise  LDA must be at */
/* >           least  LPF_MAX( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array of DIMENSION ( LDB, n ). */
/* >           Before entry, the leading  m by n part of the array  B  must */
/* >           contain the matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           LPF_MAX( 1, m ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array of DIMENSION ( LDC, n ). */
/* >           Before entry, the leading  m by n  part of the array  C must */
/* >           contain the matrix  C,  except when  beta  is zero, in which */
/* >           case C need not be set on entry. */
/* >           On exit, the array  C  is overwritten by the  m by n updated */
/* >           matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
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
void LPF_GLOBAL(bsymm,BSYMM)(char *side, char *uplo, lpf_blas_int_t *m, lpf_blas_int_t *n,
        lpf_bfloat16_t *alpha, lpf_bfloat16_t *a, lpf_blas_int_t *lda, lpf_bfloat16_t *b, lpf_blas_int_t *ldb, lpf_bfloat16_t *beta,
        lpf_bfloat16_t *c__, lpf_blas_int_t *ldc, lpf_fortran_strlen_t side_len, lpf_fortran_strlen_t uplo_len)
{
    /* System generated locals */
    lpf_blas_int_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
                 i__3;

    /* Local variables */
    lpf_blas_int_t i__, j, k, info;
    lpf_bfloat16_t temp1, temp2;
    lpf_blas_int_t nrowa;
    lpf_logical_t upper;


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

    /*     Set NROWA as the number of rows of A. */

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
    (void) side_len;

    /* Function Body */
    if (strncasecmp(side, "L", 1) == 0 ) {
        nrowa = *m;
    } else {
        nrowa = *n;
    }
    upper = strncasecmp(uplo, "U", 1 ) == 0;

    /*     Test the input parameters. */

    info = 0;
    if (! (strncasecmp(side, "L", 1) == 0) && ! (strncasecmp(side, "R", 1) == 0)) {
        info = 1;
    } else if (! upper && ! (strncasecmp(uplo, "L", 1) == 0)) {
        info = 2;
    } else if (*m < 0) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*lda < LPF_MAX(1,nrowa)) {
        info = 7;
    } else if (*ldb < LPF_MAX(1,*m)) {
        info = 9;
    } else if (*ldc < LPF_MAX(1,*m)) {
        info = 12;
    }
    if (info != 0) {
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("BSYMM ", &info, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0.f && *beta == 1.f)) {
        return;
    }

    /*     And when  alpha.eq.zero. */

    if (*alpha == 0.f) {
        if (*beta == 0.f) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = 0.f;
                    /* L10: */
                }
                /* L20: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    /* L30: */
                }
                /* L40: */
            }
        }
        return;
    }

    /*     Start the operations. */

    if ( strncasecmp(side, "L", 1) == 0) {

        /*        Form  C := alpha*A*B + beta*C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp1 = *alpha * b[i__ + j * b_dim1];
                    temp2 = 0.f;
                    i__3 = i__ - 1;
                    for (k = 1; k <= i__3; ++k) {
                        c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
                        temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
                        /* L50: */
                    }
                    if (*beta == 0.f) {
                        c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1]
                            + *alpha * temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1]
                            + temp1 * a[i__ + i__ * a_dim1] + *alpha *
                            temp2;
                    }
                    /* L60: */
                }
                /* L70: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                for (i__ = *m; i__ >= 1; --i__) {
                    temp1 = *alpha * b[i__ + j * b_dim1];
                    temp2 = 0.f;
                    i__2 = *m;
                    for (k = i__ + 1; k <= i__2; ++k) {
                        c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
                        temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
                        /* L80: */
                    }
                    if (*beta == 0.f) {
                        c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1]
                            + *alpha * temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1]
                            + temp1 * a[i__ + i__ * a_dim1] + *alpha *
                            temp2;
                    }
                    /* L90: */
                }
                /* L100: */
            }
        }
    } else {

        /*        Form  C := alpha*B*A + beta*C. */

        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            temp1 = *alpha * a[j + j * a_dim1];
            if (*beta == 0.f) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = temp1 * b[i__ + j * b_dim1];
                    /* L110: */
                }
            } else {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] +
                        temp1 * b[i__ + j * b_dim1];
                    /* L120: */
                }
            }
            i__2 = j - 1;
            for (k = 1; k <= i__2; ++k) {
                if (upper) {
                    temp1 = *alpha * a[k + j * a_dim1];
                } else {
                    temp1 = *alpha * a[j + k * a_dim1];
                }
                i__3 = *m;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
                    /* L130: */
                }
                /* L140: */
            }
            i__2 = *n;
            for (k = j + 1; k <= i__2; ++k) {
                if (upper) {
                    temp1 = *alpha * a[j + k * a_dim1];
                } else {
                    temp1 = *alpha * a[k + j * a_dim1];
                }
                i__3 = *m;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
                    /* L150: */
                }
                /* L160: */
            }
            /* L170: */
        }
    }

    return ;

    /*     End of BSYMM . */

} /* bsymm_ */

void lpf_blas_bsymm_fortran(char *side, char *uplo, lpf_blas_int_t *m, lpf_blas_int_t *n,
        lpf_fbfloat16_t *alpha, lpf_fbfloat16_t *a, lpf_blas_int_t *lda, lpf_fbfloat16_t *b, lpf_blas_int_t *ldb, lpf_fbfloat16_t *beta,
        lpf_fbfloat16_t *c__, lpf_blas_int_t *ldc)
{
    LPF_GLOBAL(bsymm,BSYMM)(side, uplo, m, n,
        (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)a, lda, (lpf_bfloat16_t *)b, ldb, (lpf_bfloat16_t *)beta,
        (lpf_bfloat16_t *)c__, ldc, 1, 1);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_bsymm_fortran_dyn_rank(char *side, char *uplo, lpf_blas_int_t *m, lpf_blas_int_t *n,
        lpf_fbfloat16_t *alpha, CFI_cdesc_t *_a, lpf_blas_int_t *lda, CFI_cdesc_t *_b, lpf_blas_int_t *ldb, lpf_fbfloat16_t *beta, CFI_cdesc_t *_c, lpf_blas_int_t *ldc)
{
    lpf_bfloat16_t *a = _a->base_addr;
    lpf_bfloat16_t *b = _b->base_addr;
    lpf_bfloat16_t *c = _c->base_addr;
    LPF_GLOBAL(bsymm,BSYMM)(side, uplo, m, n,
        (lpf_bfloat16_t *)alpha, (lpf_bfloat16_t *)a, lda, (lpf_bfloat16_t *)b, ldb, (lpf_bfloat16_t *)beta,
        (lpf_bfloat16_t *)c, ldc, 1, 1);
}
