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

/* > \brief \b HROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE HROT(N,SX,INCX,SY,INCY,C,S) */

/*       .. Scalar Arguments .. */
/*       REAL C,S */
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
/* >    applies a plane rotation. */
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
void LPF_GLOBAL(hrot,HROT)(lpf_blas_int_t *n, lpf_float16_t *sx, lpf_blas_int_t *incx, lpf_float16_t *sy,
        lpf_blas_int_t *incy, lpf_float16_t *c__, lpf_float16_t *s)
{
    /* System generated locals */
    lpf_blas_int_t i__1;

    /* Local variables */
    lpf_blas_int_t i__, ix, iy;
    lpf_float16_t stemp;


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
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
        return;
    }
    if (*incx == 1 && *incy == 1) {

        /*       code for both increments equal to 1 */

        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            stemp = *c__ * sx[i__] + *s * sy[i__];
            sy[i__] = *c__ * sy[i__] - *s * sx[i__];
            sx[i__] = stemp;
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
            stemp = *c__ * sx[ix] + *s * sy[iy];
            sy[iy] = *c__ * sy[iy] - *s * sx[ix];
            sx[ix] = stemp;
            ix += *incx;
            iy += *incy;
        }
    }
    return;
} /* hrot_ */

 void lpf_blas_hrot_fortran(lpf_blas_int_t *n, lpf_ffloat16_t *sx, lpf_blas_int_t *incx, lpf_ffloat16_t *sy,
        lpf_blas_int_t *incy, lpf_ffloat16_t *c__, lpf_ffloat16_t *s)
{
    LPF_GLOBAL(hrot,HROT)(n, (lpf_float16_t *) sx, incx, (lpf_float16_t *) sy,
                incy, (lpf_float16_t *) c__, (lpf_float16_t *) s);
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hrot_fortran_dyn_rank(lpf_blas_int_t *n, CFI_cdesc_t *_sx, lpf_blas_int_t *incx, CFI_cdesc_t *_sy,
        lpf_blas_int_t *incy, lpf_ffloat16_t *c__, lpf_ffloat16_t *s)
{
    lpf_float16_t *sx = _sx->base_addr;
    lpf_float16_t *sy = _sy->base_addr;
    lpf_blas_hrot_fortran(n, (lpf_ffloat16_t *)sx, incx, (lpf_ffloat16_t *)sy,
                incy, c__, s);
}
