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
#include <math.h>
#include "lpf_internal.h"
#include "fp16_helper.h"

/* > \brief \b HNRM2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       REAL FUNCTION HNRM2(N,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > HNRM2 returns the euclidean norm of a vector via the function */
/* > name, so that */
/* > */
/* >    HNRM2 := sqrt( x'*x ). */
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
/* >  -- This version written on 25-October-1982. */
/* >     Modified on 14-October-1993 to inline the call to SLASSQ. */
/* >     Sven Hammarling, Nag Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
int16_t LPF_GLOBAL(hnrm2_fp32,HNRM2_FP32)(lpf_blas_int_t *n, lpf_float16_t *x, lpf_blas_int_t *incx)
{
    /* System generated locals */
    lpf_blas_int_t i__1, i__2;
    lpf_float16_t ret_val;
    float r__1;

    /* Builtin functions */
    /* Local variables */
    lpf_blas_int_t ix;
    float ssq, norm, scale, absxi;


    /*  -- Reference BLAS level1 routine (version 3.4.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1 || *incx < 1) {
        norm = 0.f;
    } else if (*n == 1) {
        norm = fabsf((float)x[1]);
    } else {
        scale = 0.f;
        ssq = 1.f;
        /*        The following loop is equivalent to this call to the LAPACK */
        /*        auxiliary routine: */
        /*        CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
            if (x[ix] != 0.f) {
                r__1 = (float) x[ix];
                absxi = fabsf(r__1);
                if (scale < absxi) {
                    /* Computing 2nd power */
                    r__1 = scale / absxi;
                    ssq = ssq * (r__1 * r__1) + 1.f;
                    scale = absxi;
                } else {
                    /* Computing 2nd power */
                    r__1 = absxi / scale;
                    ssq += r__1 * r__1;
                }
            }
            /* L10: */
        }
        norm = scale * sqrtf(ssq);
    }

    ret_val = norm;
    /** return ret_val; */
    RETURN_FP16(ret_val);
    /*     End of HNRM2. */

} /* hnrm2_ */


 lpf_ffloat16_t lpf_blas_hnrm2_fp32_fortran(lpf_blas_int_t *n, lpf_ffloat16_t *x, lpf_blas_int_t *incx)
{
    lpf_ffloat16_t r;
    r.value = LPF_GLOBAL(hnrm2_fp32,HNRM2_FP32)(n, (lpf_float16_t *)x, incx);
    return r;
}

#include <ISO_Fortran_binding.h>

lpf_ffloat16_t lpf_blas_hnrm2_fp32_fortran_dyn_rank(lpf_blas_int_t *n, CFI_cdesc_t *_x, lpf_blas_int_t *incx)
{
    lpf_float16_t *x = _x->base_addr;
    return lpf_blas_hnrm2_fp32_fortran(n, (lpf_ffloat16_t *)x, incx);
}

