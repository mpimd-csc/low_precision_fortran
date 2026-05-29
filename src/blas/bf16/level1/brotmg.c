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
#include "bf16_helper.h"
#include "lpf_internal.h"
#include <stdint.h>

void LPF_GLOBAL(brotmg, BROTMG)(lpf_bfloat16_t* sd1, lpf_bfloat16_t* sd2,
                                lpf_bfloat16_t* sx1, lpf_bfloat16_t* sy1,
                                lpf_bfloat16_t* sparam)
{

    const lpf_bfloat16_t zero = 0.f;
    const lpf_bfloat16_t one = 1.f;
    const lpf_bfloat16_t two = 2.f;
    const lpf_bfloat16_t gam = 32.f;
    const lpf_bfloat16_t gamsq = 1024.f;
    const lpf_bfloat16_t rgamsq = 9.76562e-4f;

    lpf_bfloat16_t r__1;

    lpf_bfloat16_t su, sp1, sp2, sq1, sq2, sh11, sh12, sh21, sh22, sflag, stemp;

    --sparam;

    if (*sd1 < zero)
    {

        sflag = -one;
        sh11 = zero;
        sh12 = zero;
        sh21 = zero;
        sh22 = zero;

        *sd1 = zero;
        *sd2 = zero;
        *sx1 = zero;
    }
    else
    {

        sp2 = *sd2 * *sy1;
        if (sp2 == zero)
        {
            sflag = -two;
            sparam[1] = sflag;
            return;
        }

        sp1 = *sd1 * *sx1;
        sq2 = sp2 * *sy1;
        sq1 = sp1 * *sx1;

        if (abs_bf16(sq1) > abs_bf16(sq2))
        {
            sh21 = -(*sy1) / *sx1;
            sh12 = sp2 / sp1;

            su = one - sh12 * sh21;

            if (su > zero)
            {
                sflag = zero;
                *sd1 /= su;
                *sd2 /= su;
                *sx1 *= su;
            }
        }
        else
        {
            if (sq2 < zero)
            {

                sflag = -one;
                sh11 = zero;
                sh12 = zero;
                sh21 = zero;
                sh22 = zero;

                *sd1 = zero;
                *sd2 = zero;
                *sx1 = zero;
            }
            else
            {
                sflag = one;
                sh11 = sp1 / sp2;
                sh22 = *sx1 / *sy1;
                su = one + sh11 * sh22;
                stemp = *sd2 / su;
                *sd2 = *sd1 / su;
                *sd1 = stemp;
                *sx1 = *sy1 * su;
            }
        }

        if (*sd1 != zero)
        {
            while (*sd1 <= rgamsq || *sd1 >= gamsq)
            {
                if (sflag == zero)
                {
                    sh11 = one;
                    sh22 = one;
                    sflag = -one;
                }
                else
                {
                    sh21 = -one;
                    sh12 = one;
                    sflag = -one;
                }
                if (*sd1 <= rgamsq)
                {

                    r__1 = gam;
                    *sd1 *= r__1 * r__1;
                    *sx1 /= gam;
                    sh11 /= gam;
                    sh12 /= gam;
                }
                else
                {

                    r__1 = gam;
                    *sd1 /= r__1 * r__1;
                    *sx1 *= gam;
                    sh11 *= gam;
                    sh12 *= gam;
                }
            }
        }
        if (*sd2 != zero)
        {
            while (abs_bf16(*sd2) <= rgamsq || abs_bf16(*sd2) >= gamsq)
            {
                if (sflag == zero)
                {
                    sh11 = one;
                    sh22 = one;
                    sflag = -one;
                }
                else
                {
                    sh21 = -one;
                    sh12 = one;
                    sflag = -one;
                }
                if (abs_bf16(*sd2) <= rgamsq)
                {

                    r__1 = gam;
                    *sd2 *= r__1 * r__1;
                    sh21 /= gam;
                    sh22 /= gam;
                }
                else
                {

                    r__1 = gam;
                    *sd2 /= r__1 * r__1;
                    sh21 *= gam;
                    sh22 *= gam;
                }
            }
        }
    }
    if (sflag < zero)
    {
        sparam[2] = sh11;
        sparam[3] = sh21;
        sparam[4] = sh12;
        sparam[5] = sh22;
    }
    else if (sflag == zero)
    {
        sparam[3] = sh21;
        sparam[4] = sh12;
    }
    else
    {
        sparam[2] = sh11;
        sparam[5] = sh22;
    }
    sparam[1] = sflag;
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_brotmg_fortran_dyn_rank(lpf_fbfloat16_t* sd1,
                                      lpf_fbfloat16_t* sd2,
                                      lpf_fbfloat16_t* sx1,
                                      lpf_fbfloat16_t* sy1,
                                      CFI_cdesc_t* _sparam)
{
    lpf_bfloat16_t* sparam = _sparam->base_addr;
    LPF_GLOBAL(brotmg, BROTMG)((lpf_bfloat16_t*)sd1, (lpf_bfloat16_t*)sd2,
                               (lpf_bfloat16_t*)sx1, (lpf_bfloat16_t*)sy1,
                               (lpf_bfloat16_t*)sparam);
}
