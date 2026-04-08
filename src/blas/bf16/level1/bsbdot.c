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



/* > \brief \b SDBDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION SDBDOT(N,SB,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL SB */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*),SY(*) */
/*       .. */

/*    PURPOSE */
/*    ======= */

/*    Compute the inner product of two vectors with extended */
/*    precision accumulation. */

/*    Returns S.P. result with dot product accumulated in D.P. */
/*    SDBDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY), */
/*    where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*    defined in a similar way using INCY. */

/*    AUTHOR */
/*    ====== */
/*    Lawson, C. L., (JPL), Hanson, R. J., (SNLA), */
/*    Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL) */

/*    ARGUMENTS */
/*    ========= */

/*    N      (input) INTEGER */
/*           number of elements in input vector(s) */

/*    SB     (input) REAL */
/*           single precision scalar to be added to inner product */

/*    SX     (input) REAL array, dimension (N) */
/*           single precision vector with N elements */

/*    INCX   (input) INTEGER */
/*           storage spacing between elements of SX */

/*    SY     (input) REAL array, dimension (N) */
/*           single precision vector with N elements */

/*    INCY   (input) INTEGER */
/*           storage spacing between elements of SY */

/*    SDBDOT (output) REAL */
/*           single precision dot product (SB if N .LE. 0) */

/*    Further Details */
/*    =============== */

/*    REFERENCES */

/*    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*    Krogh, Basic linear algebra subprograms for Fortran */
/*    usage, Algorithm No. 539, Transactions on Mathematical */
/*    Software 5, 3 (September 1979), pp. 308-323. */

/*    REVISION HISTORY  (YYMMDD) */

/*    791001  DATE WRITTEN */
/*    890531  Changed all specific intrinsics to generic.  (WRB) */
/*    890831  Modified array declarations.  (WRB) */
/*    890831  REVISION DATE from Version 3.2 */
/*    891214  Prologue converted to Version 4.0 format.  (BAB) */
/*    920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*    920501  Reformatted the REFERENCES section.  (WRB) */
/*    070118  Reformat to LAPACK coding style */

/*    ===================================================================== */

/*       .. Local Scalars .. */
/*       DOUBLE PRECISION DBDOT */
/*       INTEGER I,KX,KY,NS */
/*       .. */
/*       .. Intrinsic Functions .. */
/*       INTRINSIC DBLE */
/*       .. */
/*       DBDOT = SB */
/*       IF (N.LE.0) THEN */
/*          SDBDOT = DBDOT */
/*          RETURN */
/*       END IF */
/*       IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN */

/*       Code for equal and positive increments. */

/*          NS = N*INCX */
/*          DO I = 1,NS,INCX */
/*             DBDOT = DBDOT + DBLE(SX(I))*DBLE(SY(I)) */
/*          END DO */
/*       ELSE */

/*       Code for unequal or nonpositive increments. */

/*          KX = 1 */
/*          KY = 1 */
/*          IF (INCX.LT.0) KX = 1 + (1-N)*INCX */
/*          IF (INCY.LT.0) KY = 1 + (1-N)*INCY */
/*          DO I = 1,N */
/*             DBDOT = DBDOT + DBLE(SX(KX))*DBLE(SY(KY)) */
/*             KX = KX + INCX */
/*             KY = KY + INCY */
/*          END DO */
/*       END IF */
/*       SDBDOT = DBDOT */
/*       RETURN */
/*       END */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
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
int16_t LPF_GLOBAL(bsbdot,BSBDOT)(lpf_blas_int_t *n, lpf_bfloat16_t *sb, lpf_bfloat16_t *sx, lpf_blas_int_t *incx, lpf_bfloat16_t *sy,
        lpf_blas_int_t *incy)
{
    /* System generated locals */
    lpf_blas_int_t i__1, i__2;
    lpf_bfloat16_t ret_val;

    /* Local variables */
    lpf_blas_int_t i__, ns, kx, ky;
    float dbdot;


    /*  -- Reference BLAS level1 routine (version 3.4.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  PURPOSE */
    /*  ======= */

    /*  Compute the inner product of two vectors with extended */
    /*  precision accumulation. */

    /*  Returns S.P. result with dot product accumulated in D.P. */
    /*  SDBDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY), */
    /*  where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
    /*  defined in a similar way using INCY. */

    /*  AUTHOR */
    /*  ====== */
    /*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA), */
    /*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL) */

    /*  ARGUMENTS */
    /*  ========= */

    /*  N      (input) INTEGER */
    /*         number of elements in input vector(s) */

    /*  SB     (input) REAL */
    /*         single precision scalar to be added to inner product */

    /*  SX     (input) REAL array, dimension (N) */
    /*         single precision vector with N elements */

    /*  INCX   (input) INTEGER */
    /*         storage spacing between elements of SX */

    /*  SY     (input) REAL array, dimension (N) */
    /*         single precision vector with N elements */

    /*  INCY   (input) INTEGER */
    /*         storage spacing between elements of SY */

    /*  SDBDOT (output) REAL */
    /*         single precision dot product (SB if N .LE. 0) */

    /*  Further Details */
    /*  =============== */

    /*  REFERENCES */

    /*  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
    /*  Krogh, Basic linear algebra subprograms for Fortran */
    /*  usage, Algorithm No. 539, Transactions on Mathematical */
    /*  Software 5, 3 (September 1979), pp. 308-323. */

    /*  REVISION HISTORY  (YYMMDD) */

    /*  791001  DATE WRITTEN */
    /*  890531  Changed all specific intrinsics to generic.  (WRB) */
    /*  890831  Modified array declarations.  (WRB) */
    /*  890831  REVISION DATE from Version 3.2 */
    /*  891214  Prologue converted to Version 4.0 format.  (BAB) */
    /*  920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
    /*  920501  Reformatted the REFERENCES section.  (WRB) */
    /*  070118  Reformat to LAPACK coding style */

    /*  ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    dbdot = *sb;
    if (*n <= 0) {
        ret_val = dbdot;
        RETURN_BF16(ret_val);
    }
    if (*incx == *incy && *incx > 0) {

        /*     Code for equal and positive increments. */

        ns = *n * *incx;
        i__1 = ns;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dbdot += (float) sx[i__] * (float) sy[i__];
        }
    } else {

        /*     Code for unequal or nonpositive increments. */

        kx = 1;
        ky = 1;
        if (*incx < 0) {
            kx = (1 - *n) * *incx + 1;
        }
        if (*incy < 0) {
            ky = (1 - *n) * *incy + 1;
        }
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            dbdot += (float) sx[kx] * (float) sy[ky];
            kx += *incx;
            ky += *incy;
        }
    }
    ret_val = dbdot;
    /** return (lpf_bfloat16_t) ret_val; */
    RETURN_BF16(ret_val);
} /* sdbdot_ */

 lpf_fbfloat16_t lpf_blas_bsbdot_fortran(lpf_blas_int_t *n, lpf_fbfloat16_t *sb, lpf_fbfloat16_t *sx, lpf_blas_int_t *incx, lpf_fbfloat16_t *sy,
        lpf_blas_int_t *incy)
{
    lpf_fbfloat16_t r;
    r.value = LPF_GLOBAL(bsbdot,BSBDOT)(n,(lpf_bfloat16_t *) sb, (lpf_bfloat16_t *) sx, incx, (lpf_bfloat16_t *) sy, incy);
    return r;
}

#include <ISO_Fortran_binding.h>

lpf_fbfloat16_t lpf_blas_bsbdot_fortran_dyn_rank(lpf_blas_int_t *n, lpf_fbfloat16_t *sb, CFI_cdesc_t *_sx, lpf_blas_int_t *incx, CFI_cdesc_t *_sy, lpf_blas_int_t *incy)
{
    lpf_bfloat16_t *sx = _sx->base_addr;
    lpf_bfloat16_t *sy = _sy->base_addr;
    return lpf_blas_bsbdot_fortran(n, sb, (lpf_fbfloat16_t *)sx, incx, (lpf_fbfloat16_t *)sy, incy);
}

