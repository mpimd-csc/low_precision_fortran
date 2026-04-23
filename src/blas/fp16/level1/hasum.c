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

#include <stdint.h>
#include "lpf_internal.h"
#include "fp16_helper.h"




/* > \brief \b HASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION HASUM(N,SX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL SX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    HASUM takes the sum of the absolute values. */
/* >    uses unrolled loops for increment equal to one. */
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
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */


int16_t LPF_GLOBAL(hasum,HASUM)(lpf_blas_int_t *n, lpf_float16_t *sx, lpf_blas_int_t *incx)
{
    /* System generated locals */
    lpf_blas_int_t i__1, i__2;
    lpf_float16_t ret_val;

    /* Local variables */
    lpf_blas_int_t i__, m, mp1, nincx;
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
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0) {
        RETURN_FP16(ret_val);
    }
    if (*incx == 1) {
        /*        code for increment equal to 1 */


        /*        clean-up loop */

        m = *n % 6;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                stemp += abs_fp16(sx[i__]);
            }
            if (*n < 6) {
                ret_val = stemp;
                RETURN_FP16(ret_val);
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 6) {
            stemp = stemp
                + abs_fp16(sx[i__])
                + abs_fp16(sx[i__ + 1])
                + abs_fp16(sx[i__ + 2])
                + abs_fp16(sx[i__ + 3])
                + abs_fp16(sx[i__ + 4])
                + abs_fp16(sx[i__ + 5]);
        }
    } else {

        /*        code for increment not equal to 1 */

        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            stemp += abs_fp16(sx[i__]);
        }
    }
    ret_val = stemp;
    RETURN_FP16(ret_val);
} /* hasum_ */

lpf_ffloat16_t lpf_blas_hasum_fortran(lpf_blas_int_t *n, lpf_ffloat16_t *sx, lpf_blas_int_t *incx)
{
    lpf_ffloat16_t r;
    r.value = LPF_GLOBAL(hasum,HASUM)(n, (lpf_float16_t*) sx, incx);
    return r;
}

#include <ISO_Fortran_binding.h>

lpf_ffloat16_t lpf_blas_hasum_fortran_dyn_rank(lpf_blas_int_t *n, CFI_cdesc_t *_sx, lpf_blas_int_t *incx)
{
    lpf_float16_t *sx = _sx->base_addr;
    return lpf_blas_hasum_fortran(n, (lpf_ffloat16_t *)sx, incx);
}

