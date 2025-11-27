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

/* > \brief \b BSWAP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BSWAP(N,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*),SY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    interchanges two vectors. */
/* >    uses unrolled loops for increments equal to 1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup single_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(bswap,BSWAP)(lpf_blas_int_t *n, lpf_bfloat16_t *sx, lpf_blas_int_t *incx, lpf_bfloat16_t *sy,
        lpf_blas_int_t *incy)
{
    /* System generated locals */
    lpf_blas_int_t i__1;

    /* Local variables */
    lpf_blas_int_t i__, m, ix, iy, mp1;
    lpf_bfloat16_t stemp;


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
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
        return;
    }
    if (*incx == 1 && *incy == 1) {

        /*       code for both increments equal to 1 */


        /*       clean-up loop */

        m = *n % 3;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                stemp = sx[i__];
                sx[i__] = sy[i__];
                sy[i__] = stemp;
            }
            if (*n < 3) {
                return;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 3) {
            stemp = sx[i__];
            sx[i__] = sy[i__];
            sy[i__] = stemp;
            stemp = sx[i__ + 1];
            sx[i__ + 1] = sy[i__ + 1];
            sy[i__ + 1] = stemp;
            stemp = sx[i__ + 2];
            sx[i__ + 2] = sy[i__ + 2];
            sy[i__ + 2] = stemp;
        }
    } else {

        /*       code for unequal increments or equal increments not equal */
        /*         to 1 */

        ix = 1;
        iy = 1;
        if (*incx < 0) {
            ix = (-(*n) + 1) * *incx + 1;
        }
        if (*incy < 0) {
            iy = (-(*n) + 1) * *incy + 1;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            stemp = sx[ix];
            sx[ix] = sy[iy];
            sy[iy] = stemp;
            ix += *incx;
            iy += *incy;
        }
    }
    return;
} /* bswap_ */

 void lpf_blas_bswap_fortran(lpf_blas_int_t * n, lpf_fbfloat16_t *sx, lpf_blas_int_t *incx, lpf_fbfloat16_t *sy,
        lpf_blas_int_t *incy)
{
    LPF_GLOBAL(bswap,BSWAP)(n, (lpf_bfloat16_t *) sx, incx, (lpf_bfloat16_t *) sy, incy);
}

