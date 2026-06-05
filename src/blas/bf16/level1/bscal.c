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
#include "lpf_internal.h"
#include <math.h>
#include <stdint.h>

/**
 * @brief Bfloat16 Scaling (BSCAL).
 *
 * Performs the operation: x := alpha * x
 *
 * @param[in] n Number of elements in vector x.
 * @param[in] sa Scalar multiplier alpha.
 * @param[in,out] sx Vector x.
 * @param[in] incx Increment for the elements of x.
 */
void LPF_GLOBAL(bscal, BSCAL)(int64_t* n, lpf_bfloat16_t* sa,
                              lpf_bfloat16_t* sx, int64_t* incx)
{

    int64_t i__1, i__2;

    int64_t i__, m, mp1, nincx;

    --sx;

    if (*n <= 0 || *incx <= 0)
    {
        return;
    }
    if (*incx == 1)
    {

        m = *n % 5;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                sx[i__] = *sa * sx[i__];
            }
            if (*n < 5)
            {
                return;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 5)
        {
            sx[i__] = *sa * sx[i__];
            sx[i__ + 1] = *sa * sx[i__ + 1];
            sx[i__ + 2] = *sa * sx[i__ + 2];
            sx[i__ + 3] = *sa * sx[i__ + 3];
            sx[i__ + 4] = *sa * sx[i__ + 4];
        }
    }
    else
    {

        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            sx[i__] = *sa * sx[i__];
        }
    }
    return;
}

#include <ISO_Fortran_binding.h>

void lpf_blas_bscal_fortran_dyn_rank_64(int64_t* n, lpf_fbfloat16_t* sa,
                                        CFI_cdesc_t* _sx, int64_t* incx)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    LPF_GLOBAL(bscal, BSCAL)(n, (lpf_bfloat16_t*)sa, sx, incx);
}

void lpf_blas_bscal_fortran_dyn_rank_32(int32_t* n, lpf_fbfloat16_t* sa,
                                        CFI_cdesc_t* _sx, int32_t* incx)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    LPF_GLOBAL(bscal, BSCAL)(&_n, (lpf_bfloat16_t*)sa, sx, &_incx);
}
