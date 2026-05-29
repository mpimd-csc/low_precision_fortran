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

/* > \brief \b BROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BROT(N,SX,INCX,SY,INCY,C,S) */

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
void LPF_GLOBAL(brot,BROT)(int64_t *n, lpf_bfloat16_t *sx, int64_t *incx, lpf_bfloat16_t *sy,
        int64_t *incy, lpf_bfloat16_t *c__, lpf_bfloat16_t *s)
{
    /* System generated locals */
    int64_t i__1;

    /* Local variables */
    int64_t i__, ix, iy;
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
} /* brot_ */

#include <ISO_Fortran_binding.h>

void lpf_blas_brot_fortran_dyn_rank_64(int64_t *n, CFI_cdesc_t *_sx, int64_t *incx, CFI_cdesc_t *_sy,
        int64_t *incy, lpf_fbfloat16_t *c__, lpf_fbfloat16_t *s)
{
    lpf_bfloat16_t *sx = _sx->base_addr;
    lpf_bfloat16_t *sy = _sy->base_addr;
    LPF_GLOBAL(brot,BROT)(n, (lpf_bfloat16_t *)sx, incx, (lpf_bfloat16_t *)sy,
                incy, (lpf_bfloat16_t *)c__, (lpf_bfloat16_t *)s);
}

void lpf_blas_brot_fortran_dyn_rank_32(int32_t *n, CFI_cdesc_t *_sx, int32_t *incx, CFI_cdesc_t *_sy,
        int32_t *incy, lpf_fbfloat16_t *c__, lpf_fbfloat16_t *s)
{
    lpf_bfloat16_t *sx = _sx->base_addr;
    lpf_bfloat16_t *sy = _sy->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    LPF_GLOBAL(brot,BROT)(&_n, (lpf_bfloat16_t *)sx, &_incx, (lpf_bfloat16_t *)sy,
                &_incy, (lpf_bfloat16_t *)c__, (lpf_bfloat16_t *)s);
}
