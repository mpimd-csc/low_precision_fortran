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

/* > \brief \b SAXPBY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SAXPBY(N,SA,SX,INCX,SY,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL SA, SB */
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
/* >    SAXPBY constant times a vector plus a constant times a vector. */
/* >    uses unrolled loops for increments equal to one. */
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
void LPF_GLOBAL(haxpby,HAXPBY)(lpf_blas_int_t *n, lpf_float16_t *sa, lpf_float16_t *sx, lpf_blas_int_t *incx,
        lpf_float16_t *sb, lpf_float16_t *sy, lpf_blas_int_t *incy)
{
    /* System generated locals */
    lpf_blas_int_t i__1;

    /* Local variables */
    lpf_blas_int_t i__, m, ix, iy, mp1;


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
    if (*sa == 0.f && *sb == 1.f) {
        return;
    }
    if (*incx == 1 && *incy == 1) {

        /*        code for both increments equal to 1 */


        /*        clean-up loop */

        m = *n % 4;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                sy[i__] = *sb * sy[i__] + *sa * sx[i__];
            }
        }
        if (*n < 4) {
            return;
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 4) {
            sy[i__] = *sb * sy[i__] + *sa * sx[i__];
            sy[i__ + 1] = *sb * sy[i__ + 1] + *sa * sx[i__ + 1];
            sy[i__ + 2] = *sb * sy[i__ + 2] + *sa * sx[i__ + 2];
            sy[i__ + 3] = *sb * sy[i__ + 3] + *sa * sx[i__ + 3];
        }
    } else {

        /*        code for unequal increments or equal increments */
        /*          not equal to 1 */

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
            sy[iy] = *sb * sy[iy] + *sa * sx[ix];
            ix += *incx;
            iy += *incy;
        }
    }
    return;
} /* saxpby_ */

 void lpf_blas_haxpby_fortran( lpf_blas_int_t *n, lpf_ffloat16_t *sa, lpf_ffloat16_t *sx, lpf_blas_int_t *incx, lpf_ffloat16_t * sb, lpf_ffloat16_t *sy, lpf_blas_int_t * incy )
{
    LPF_GLOBAL(haxpby,HAXPBY)( n, (lpf_float16_t * ) sa, (lpf_float16_t * ) sx, incx, (lpf_float16_t * ) sb, (lpf_float16_t*) sy, incy);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_haxpby_fortran_dyn_rank(lpf_blas_int_t *n, lpf_ffloat16_t *sa, CFI_cdesc_t *_sx, lpf_blas_int_t *incx, lpf_ffloat16_t *sb, CFI_cdesc_t *_sy, lpf_blas_int_t *incy)
{
    lpf_float16_t *sx = _sx->base_addr;
    lpf_float16_t *sy = _sy->base_addr;
    lpf_blas_haxpby_fortran(n, sa, (lpf_ffloat16_t *)sx, incx, sb, (lpf_ffloat16_t *)sy, incy);
}

