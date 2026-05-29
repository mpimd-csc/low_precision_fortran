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
#include <stdio.h>
#include "lpf_internal.h"

/* > \brief \b BDOT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION BDOT(N,SX,INCX,SY,INCY) */

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
/* >    BDOT forms the dot product of two vectors. */
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
int16_t LPF_GLOBAL(bdot,BDOT)(int64_t *n, lpf_bfloat16_t *sx, int64_t *incx, lpf_bfloat16_t *sy, int64_t *incy)
{
    /* System generated locals */
    int64_t i__1;
    lpf_bfloat16_t ret_val;

    /* Local variables */
    int64_t i__, m, ix, iy, mp1;
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
    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0) {
        RETURN_BF16(ret_val);
    }
    if (*incx == 1 && *incy == 1) {

        /*        code for both increments equal to 1 */


        /*        clean-up loop */

        m = *n % 5;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                stemp += sx[i__] * sy[i__];
            }
            if (*n < 5) {
                ret_val = stemp;
                RETURN_BF16(ret_val);
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 5) {
            stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] +
                sx[i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] +
                sx[i__ + 4] * sy[i__ + 4];
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
            stemp += sx[ix] * sy[iy];
            ix += *incx;
            iy += *incy;
        }
    }
    ret_val = stemp;

    /** lpf_blas_bf16_i16_t r; */
    /** r.f16 = ret_val; */
    /** printf("%d\n",(int) r.i16); */
    /** return r.i16; */
    RETURN_BF16(ret_val);
} /* bdot_ */


#include <ISO_Fortran_binding.h>

lpf_fbfloat16_t lpf_blas_bdot_fortran_dyn_rank_64(int64_t *n, CFI_cdesc_t *_sx, int64_t *incx, CFI_cdesc_t *_sy, int64_t *incy)
{
    lpf_bfloat16_t *sx = _sx->base_addr;
    lpf_bfloat16_t *sy = _sy->base_addr;
    lpf_fbfloat16_t r;
    r.value = LPF_GLOBAL(bdot,BDOT)(n, sx, incx, sy, incy);
    return r;
}

lpf_fbfloat16_t lpf_blas_bdot_fortran_dyn_rank_32(int32_t *n, CFI_cdesc_t *_sx, int32_t *incx, CFI_cdesc_t *_sy, int32_t *incy)
{
    lpf_bfloat16_t *sx = _sx->base_addr;
    lpf_bfloat16_t *sy = _sy->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    int64_t _incy = *incy;
    lpf_fbfloat16_t r;
    r.value = LPF_GLOBAL(bdot,BDOT)(&_n, sx, &_incx, sy, &_incy);
    return r;
}
