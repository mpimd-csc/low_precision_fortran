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

/* > \brief \b IHAMAX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IHAMAX(N,SX,INCX) */

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
/* >    IHAMAX finds the index of the first element having maximum absolute value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup aux_blas */

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
lpf_blas_int_t LPF_GLOBAL(ihamax,IHAMAX)(lpf_blas_int_t *n, lpf_float16_t *sx, lpf_blas_int_t *incx)
{
    /* System generated locals */
    lpf_blas_int_t ret_val, i__1;
    lpf_float16_t r__1;

    /* Local variables */
    lpf_blas_int_t i__, ix;
    lpf_float16_t smax;


    /*  -- Reference BLAS level1 routine (version 3.6.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2015 */

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
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
        return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
        return ret_val;
    }
    if (*incx == 1) {

        /*        code for increment equal to 1 */

        smax = abs_fp16(sx[1]);
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            r__1 = sx[i__];
            if (abs_fp16(r__1) > smax) {
                ret_val = i__;
                smax = abs_fp16(r__1);
            }
        }
    } else {

        /*        code for increment not equal to 1 */

        ix = 1;
        smax = abs_fp16(sx[1]);
        ix += *incx;
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            r__1 = abs_fp16(sx[ix]);
            if ( r__1 > smax) {
                ret_val = i__;
                smax = r__1;
            }
            ix += *incx;
        }
    }
    return ret_val;
} /* ihamax_ */

lpf_blas_int_t lpf_blas_ihamax_fortran(lpf_blas_int_t *n, lpf_ffloat16_t *sx, lpf_blas_int_t *incx)
{
    return LPF_GLOBAL(ihamax,IHAMAX)(n, (lpf_float16_t *) sx, incx);
}
