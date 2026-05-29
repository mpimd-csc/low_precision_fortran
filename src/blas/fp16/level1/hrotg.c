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
#include "fp16_helper.h"
#include "lpf_internal.h"
#include <math.h>
#include <stdint.h>

void LPF_GLOBAL(hrotg, HROTG)(lpf_float16_t* sa, lpf_float16_t* sb,
                              lpf_float16_t* c__, lpf_float16_t* s)
{

    lpf_float16_t r__1, r__2;

    lpf_float16_t r__, z__, roe, scale;

    const lpf_float16_t c_b2 = 1.f;

    roe = *sb;
    if (abs_fp16(*sa) > abs_fp16(*sb))
    {
        roe = *sa;
    }
    scale = abs_fp16(*sa) + abs_fp16(*sb);
    if (scale == 0.f)
    {
        *c__ = 1.f;
        *s = 0.f;
        r__ = 0.f;
        z__ = 0.f;
    }
    else
    {

        r__1 = *sa / scale;

        r__2 = *sb / scale;
        r__ = scale * sqrtf((float)(r__1 * r__1 + r__2 * r__2));
        r__ = sign_fp16(c_b2, roe) * r__;
        *c__ = *sa / r__;
        *s = *sb / r__;
        z__ = 1.f;
        if (abs_fp16(*sa) > abs_fp16(*sb))
        {
            z__ = *s;
        }
        if (abs_fp16(*sb) >= abs_fp16(*sa) && *c__ != 0.f)
        {
            z__ = 1.f / *c__;
        }
    }
    *sa = r__;
    *sb = z__;
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_hrotg_fortran(lpf_ffloat16_t* sa, lpf_ffloat16_t* sb,
                            lpf_ffloat16_t* c__, lpf_ffloat16_t* s)
{
    LPF_GLOBAL(hrotg, HROTG)((lpf_float16_t*)sa, (lpf_float16_t*)sb,
                             (lpf_float16_t*)c__, (lpf_float16_t*)s);
}
