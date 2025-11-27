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

/* > \brief \b BROTM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BROTM(N,SX,INCX,SY,INCY,SPARAM) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SPARAM(5),SX(*),SY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX */
/* > */
/* >    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN */
/* >    (SX**T) */
/* > */
/* >    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE */
/* >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY. */
/* >    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS.. */
/* > */
/* >    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */
/* > */
/* >      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/* >    H=(          )    (          )    (          )    (          ) */
/* >      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */
/* >    SEE  BROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SX */
/* > \verbatim */
/* >          SX is REAL array, dimension N */
/* >         double precision vector with N elements */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of SX */
/* > \endverbatim */
/* > */
/* > \param[in,out] SY */
/* > \verbatim */
/* >          SY is REAL array, dimension N */
/* >         double precision vector with N elements */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of SY */
/* > \endverbatim */
/* > */
/* > \param[in,out] SPARAM */
/* > \verbatim */
/* >          SPARAM is REAL array, dimension 5 */
/* >     SPARAM(1)=SFLAG */
/* >     SPARAM(2)=SH11 */
/* >     SPARAM(3)=SH21 */
/* >     SPARAM(4)=SH12 */
/* >     SPARAM(5)=SH22 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level1 */

/*  ===================================================================== */
void LPF_GLOBAL(brotm,BROTM)(lpf_blas_int_t *n, lpf_bfloat16_t *sx, lpf_blas_int_t *incx, lpf_bfloat16_t *sy,
        lpf_blas_int_t *incy, lpf_bfloat16_t *sparam)
{
    /* Initialized data */

    const lpf_bfloat16_t zero = 0.f;
    const lpf_bfloat16_t two = 2.f;

    /* System generated locals */
    lpf_blas_int_t i__1, i__2;

    /* Local variables */
    lpf_blas_int_t i__;
    lpf_bfloat16_t w, z__;
    lpf_blas_int_t kx, ky;
    lpf_bfloat16_t sh11, sh12, sh21, sh22, sflag;
    lpf_blas_int_t nsteps;


    /*  -- Reference BLAS level1 routine (version 3.4.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Data statements .. */
    /* Parameter adjustments */
    --sparam;
    --sy;
    --sx;

    /* Function Body */
    /*     .. */

    sflag = sparam[1];
    if (*n <= 0 || sflag + two == zero) {
        return;
    }
    if (*incx == *incy && *incx > 0) {

        nsteps = *n * *incx;
        if (sflag < zero) {
            sh11 = sparam[2];
            sh12 = sparam[4];
            sh21 = sparam[3];
            sh22 = sparam[5];
            i__1 = nsteps;
            i__2 = *incx;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                w = sx[i__];
                z__ = sy[i__];
                sx[i__] = w * sh11 + z__ * sh12;
                sy[i__] = w * sh21 + z__ * sh22;
            }
        } else if (sflag == zero) {
            sh12 = sparam[4];
            sh21 = sparam[3];
            i__2 = nsteps;
            i__1 = *incx;
            for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
                w = sx[i__];
                z__ = sy[i__];
                sx[i__] = w + z__ * sh12;
                sy[i__] = w * sh21 + z__;
            }
        } else {
            sh11 = sparam[2];
            sh22 = sparam[5];
            i__1 = nsteps;
            i__2 = *incx;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                w = sx[i__];
                z__ = sy[i__];
                sx[i__] = w * sh11 + z__;
                sy[i__] = -w + sh22 * z__;
            }
        }
    } else {
        kx = 1;
        ky = 1;
        if (*incx < 0) {
            kx = (1 - *n) * *incx + 1;
        }
        if (*incy < 0) {
            ky = (1 - *n) * *incy + 1;
        }

        if (sflag < zero) {
            sh11 = sparam[2];
            sh12 = sparam[4];
            sh21 = sparam[3];
            sh22 = sparam[5];
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                w = sx[kx];
                z__ = sy[ky];
                sx[kx] = w * sh11 + z__ * sh12;
                sy[ky] = w * sh21 + z__ * sh22;
                kx += *incx;
                ky += *incy;
            }
        } else if (sflag == zero) {
            sh12 = sparam[4];
            sh21 = sparam[3];
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                w = sx[kx];
                z__ = sy[ky];
                sx[kx] = w + z__ * sh12;
                sy[ky] = w * sh21 + z__;
                kx += *incx;
                ky += *incy;
            }
        } else {
            sh11 = sparam[2];
            sh22 = sparam[5];
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                w = sx[kx];
                z__ = sy[ky];
                sx[kx] = w * sh11 + z__;
                sy[ky] = -w + sh22 * z__;
                kx += *incx;
                ky += *incy;
            }
        }
    }
    return;
} /* brotm_ */

 void lpf_blas_brotm_fortran(lpf_blas_int_t *n, lpf_fbfloat16_t *sx, lpf_blas_int_t *incx, lpf_fbfloat16_t *sy,
        lpf_blas_int_t *incy, lpf_fbfloat16_t *sparam)
{
    LPF_GLOBAL(brotm,BROTM)(n, (lpf_bfloat16_t *) sx, incx, (lpf_bfloat16_t *) sy, incy, (lpf_bfloat16_t *) sparam);
}
