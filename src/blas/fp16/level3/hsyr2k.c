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
/* > \brief \b HSYR2K */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HSYR2K  performs one of the symmetric rank 2k operations */
/* > */
/* >    C := alpha*A*B**T + alpha*B*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*B + alpha*B**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
/* > and  A and B  are  n by k  matrices  in the  first  case  and  k by n */
/* > matrices in the second case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of the  array  C  is to be  referenced  as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/* >                                  is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry,  TRANS  specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   C := alpha*A*B**T + alpha*B*A**T + */
/* >                                        beta*C. */
/* > */
/* >              TRANS = 'T' or 't'   C := alpha*A**T*B + alpha*B**T*A + */
/* >                                        beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'   C := alpha*A**T*B + alpha*B**T*A + */
/* >                                        beta*C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N specifies the order of the matrix C.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/* >           of  columns  of the  matrices  A and B,  and on  entry  with */
/* >           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
/* >           of rows of the matrices  A and B.  K must be at least  zero. */
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
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by n  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDA must be at least  LPF_MAX( 1, n ), otherwise  LDA must */
/* >           be at least  LPF_MAX( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array of DIMENSION ( LDB, kb ), where kb is */
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  k by n  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDB must be at least  LPF_MAX( 1, n ), otherwise  LDB must */
/* >           be at least  LPF_MAX( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array of DIMENSION ( LDC, n ). */
/* >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/* >           upper triangular part of the array C must contain the upper */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           lower triangular part of C is not referenced.  On exit, the */
/* >           upper triangular part of the array  C is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/* >           lower triangular part of the array C must contain the lower */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           upper triangular part of C is not referenced.  On exit, the */
/* >           lower triangular part of the array  C is overwritten by the */
/* >           lower triangular part of the updated matrix. */
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
/* > */
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(hsyr2k,HSYR2K)(char *uplo, char *trans, int64_t *n, int64_t *k,
        lpf_float16_t *alpha, lpf_float16_t *a, int64_t *lda, lpf_float16_t *b, int64_t *ldb, lpf_float16_t *beta,
        lpf_float16_t *c__, int64_t *ldc, lpf_fortran_strlen_t uplo_len, lpf_fortran_strlen_t trans_len)
{
    /* System generated locals */
    int64_t a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
                 i__3;

    /* Local variables */
    int64_t i__, j, l, info;
    lpf_float16_t temp1, temp2;
    int64_t nrowa;
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

    /*     Test the input parameters. */

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
    (void) trans_len;
    /* Function Body */
    if ( strncasecmp(trans, "N", 1) == 0 ) {
        nrowa = *n;
    } else {
        nrowa = *k;
    }
    upper = strncasecmp(uplo, "U", 1) == 0;

    info = 0;
    if (! upper && ! (strncasecmp(uplo, "L", 1) == 0 )) {
        info = 1;
    } else if (! (strncasecmp(trans, "N", 1) == 0 ) && ! (strncasecmp(trans, "T", 1) == 0)  && ! (strncasecmp(trans, "C", 1) == 0)) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*k < 0) {
        info = 4;
    } else if (*lda < LPF_MAX(1,nrowa)) {
        info = 7;
    } else if (*ldb < LPF_MAX(1,nrowa)) {
        info = 9;
    } else if (*ldc < LPF_MAX(1,*n)) {
        info = 12;
    }
    if (info != 0) {
        int32_t infox = info;
        LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)("HSYR2K", &infox, (lpf_fortran_strlen_t)6);
        return;
    }

    /*     Quick return if possible. */

    if (*n == 0 || ((*alpha == 0.f || *k == 0) && *beta == 1.f)) {
        return;
    }

    /*     And when  alpha.eq.zero. */

    if (*alpha == 0.f) {
        if (upper) {
            if (*beta == 0.f) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.f;
                        /* L10: */
                    }
                    /* L20: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                        /* L30: */
                    }
                    /* L40: */
                }
            }
        } else {
            if (*beta == 0.f) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.f;
                        /* L50: */
                    }
                    /* L60: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        }
        return;
    }

    /*     Start the operations. */

    if (strncasecmp(trans, "N", 1) == 0) {

        /*        Form  C := alpha*A*B**T + alpha*B*A**T + C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.f) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.f;
                        /* L90: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                        /* L100: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    if (a[j + l * a_dim1] != 0.f || b[j + l * b_dim1] != 0.f)
                    {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = j;
                        for (i__ = 1; i__ <= i__3; ++i__) {
                            c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
                                i__ + l * a_dim1] * temp1 + b[i__ + l *
                                    b_dim1] * temp2;
                            /* L110: */
                        }
                    }
                    /* L120: */
                }
                /* L130: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.f) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.f;
                        /* L140: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                        /* L150: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    if (a[j + l * a_dim1] != 0.f || b[j + l * b_dim1] != 0.f)
                    {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = *n;
                        for (i__ = j; i__ <= i__3; ++i__) {
                            c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
                                i__ + l * a_dim1] * temp1 + b[i__ + l *
                                    b_dim1] * temp2;
                            /* L160: */
                        }
                    }
                    /* L170: */
                }
                /* L180: */
            }
        }
    } else {

        /*        Form  C := alpha*A**T*B + alpha*B**T*A + C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp1 = 0.f;
                    temp2 = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
                        /* L190: */
                    }
                    if (*beta == 0.f) {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha *
                            temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1]
                            + *alpha * temp1 + *alpha * temp2;
                    }
                    /* L200: */
                }
                /* L210: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n;
                for (i__ = j; i__ <= i__2; ++i__) {
                    temp1 = 0.f;
                    temp2 = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
                        /* L220: */
                    }
                    if (*beta == 0.f) {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha *
                            temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1]
                            + *alpha * temp1 + *alpha * temp2;
                    }
                    /* L230: */
                }
                /* L240: */
            }
        }
    }

    return;

    /*     End of HSYR2K. */

} /* hsyr2k_ */

#include <ISO_Fortran_binding.h>

void lpf_blas_hsyr2k_fortran_dyn_rank_64(char *uplo, char *trans, int64_t *n, int64_t *k,
        lpf_ffloat16_t *alpha, CFI_cdesc_t *_a, int64_t *lda, CFI_cdesc_t *_b, int64_t *ldb,
        lpf_ffloat16_t *beta, CFI_cdesc_t *_c, int64_t *ldc)
{
    lpf_float16_t *a = _a->base_addr;
    lpf_float16_t *b = _b->base_addr;
    lpf_float16_t *c = _c->base_addr;

    LPF_GLOBAL(hsyr2k,HSYR2K)(uplo, trans, n, k,
        (lpf_float16_t *)alpha, (lpf_float16_t *)a, lda, (lpf_float16_t *)b, ldb, (lpf_float16_t *)beta,
        (lpf_float16_t *)c, ldc, 1, 1);
}

void lpf_blas_hsyr2k_fortran_dyn_rank_32(char *uplo, char *trans, int32_t *n, int32_t *k,
        lpf_ffloat16_t *alpha, CFI_cdesc_t *_a, int32_t *lda, CFI_cdesc_t *_b, int32_t *ldb,
        lpf_ffloat16_t *beta, CFI_cdesc_t *_c, int32_t *ldc)
{
    int64_t _n = *n;
    int64_t _k = *k;
    int64_t _lda = *lda;
    int64_t _ldb = *ldb;
    int64_t _ldc = *ldc;

    lpf_float16_t *a = _a->base_addr;
    lpf_float16_t *b = _b->base_addr;
    lpf_float16_t *c = _c->base_addr;

    LPF_GLOBAL(hsyr2k,HSYR2K)(uplo, trans, &_n, &_k,
        (lpf_float16_t *)alpha, (lpf_float16_t *)a, &_lda, (lpf_float16_t *)b, &_ldb, (lpf_float16_t *)beta,
        (lpf_float16_t *)c, &_ldc, 1, 1);
}
