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
#include "bf16_helper.h"

/* > \brief \b BROTG */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE BROTG(SA,SB,C,S) */

/*       .. Scalar Arguments .. */
/*       REAL C,S,SA,SB */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    BROTG construct givens plane rotation. */
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
/* > \endverbatim */
/* > */
/*  ===================================================================== */
void LPF_GLOBAL(brotg,BROTG)(lpf_bfloat16_t *sa, lpf_bfloat16_t *sb, lpf_bfloat16_t *c__, lpf_bfloat16_t *s)
{
    /* System generated locals */
    lpf_bfloat16_t r__1, r__2;

    /* Local variables */
    lpf_bfloat16_t r__, z__, roe, scale;


    const lpf_bfloat16_t c_b2 = 1.f;
    /*  -- Reference BLAS level1 routine (version 3.4.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    roe = *sb;
    if ( abs_bf16(*sa) > abs_bf16(*sb)) {
        roe = *sa;
    }
    scale = abs_bf16(*sa) + abs_bf16(*sb);
    if (scale == 0.f) {
        *c__ = 1.f;
        *s = 0.f;
        r__ = 0.f;
        z__ = 0.f;
    } else {
        /* Computing 2nd power */
        r__1 = *sa / scale;
        /* Computing 2nd power */
        r__2 = *sb / scale;
        r__ = scale * sqrtf((float)(r__1 * r__1 + r__2 * r__2));
        r__ = sign_bf16(c_b2, roe) * r__;
        *c__ = *sa / r__;
        *s = *sb / r__;
        z__ = 1.f;
        if ( abs_bf16(*sa) > abs_bf16(*sb)) {
            z__ = *s;
        }
        if ( abs_bf16(*sb) >= abs_bf16(*sa) && *c__ != 0.f) {
            z__ = 1.f / *c__;
        }
    }
    *sa = r__;
    *sb = z__;
    return ;
} /* brotg_ */


void lpf_blas_brotg_fortran(lpf_fbfloat16_t *sa, lpf_fbfloat16_t *sb, lpf_fbfloat16_t *c__, lpf_fbfloat16_t *s)
{
    LPF_GLOBAL(brotg,BROTG)((lpf_bfloat16_t *) sa, (lpf_bfloat16_t *) sb, (lpf_bfloat16_t *) c__, (lpf_bfloat16_t *) s);
}
