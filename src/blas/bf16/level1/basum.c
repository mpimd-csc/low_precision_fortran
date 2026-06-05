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

/**
 * @brief Bfloat16 Sum of Absolute Values (BASUM).
 *
 * Performs the operation: result := sum(|x|)
 *
 * @param[in] n Number of elements in vector x.
 * @param[in] sx Vector x.
 * @param[in] incx Increment for the elements of x.
 * @return The sum of absolute values of vector x.
 */
int16_t LPF_GLOBAL(basum, BASUM)(int64_t* n, lpf_bfloat16_t* sx, int64_t* incx)
{

    int64_t i__1, i__2;
    lpf_bfloat16_t ret_val;

    int64_t i__, m, mp1, nincx;
    lpf_bfloat16_t stemp;

    --sx;

    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0 || *incx <= 0)
    {
        RETURN_BF16(ret_val);
    }
    if (*incx == 1)
    {

        m = *n % 6;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                stemp += abs_bf16(sx[i__]);
            }
            if (*n < 6)
            {
                ret_val = stemp;
                RETURN_BF16(ret_val);
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 6)
        {
            stemp = stemp + abs_bf16(sx[i__]) + abs_bf16(sx[i__ + 1]) +
                    abs_bf16(sx[i__ + 2]) + abs_bf16(sx[i__ + 3]) +
                    abs_bf16(sx[i__ + 4]) + abs_bf16(sx[i__ + 5]);
        }
    }
    else
    {

        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            stemp += abs_bf16(sx[i__]);
        }
    }
    ret_val = stemp;
    RETURN_BF16(ret_val);
}

#include <ISO_Fortran_binding.h>

lpf_fbfloat16_t lpf_blas_basum_fortran_dyn_rank_64(int64_t* n, CFI_cdesc_t* _sx,
                                                   int64_t* incx)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    lpf_fbfloat16_t r;
    r.value = LPF_GLOBAL(basum, BASUM)(n, sx, incx);
    return r;
}

lpf_fbfloat16_t lpf_blas_basum_fortran_dyn_rank_32(int32_t* n, CFI_cdesc_t* _sx,
                                                   int32_t* incx)
{
    lpf_bfloat16_t* sx = _sx->base_addr;
    int64_t _n = *n;
    int64_t _incx = *incx;
    lpf_fbfloat16_t r;
    r.value = LPF_GLOBAL(basum, BASUM)(&_n, sx, &_incx);
    return r;
}
