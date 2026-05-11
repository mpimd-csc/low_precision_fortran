// SPDX-License-Identifier LGPL-3.0-or-later
/*
 * Copyright (C)  2025 by Martin Koehler
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include "fp8_e5m2.h"
#include "fp8_e5m2_bits.h"

/**
 * fp8_e5m2_isnan - Check if an FP8 E5M2 value is Not-a-Number.
 *
 * In E5M2, NaN is encoded as exponent=31 and a non-zero mantissa.
 */
int fp8_e5m2_isnan(fp8_e5m2_t fp8)
{
    return (FP8_E5M2_GET_EXP(fp8) == FP8_E5M2_E_MAX) && (FP8_E5M2_GET_MANT(fp8) != 0);
}

/**
 * fp8_e5m2_isinf - Check if an FP8 E5M2 value represents infinity.
 *
 * In E5M2, infinity is encoded as exponent=31 and a zero mantissa.
 */
int fp8_e5m2_isinf(fp8_e5m2_t fp8)
{
    return (FP8_E5M2_GET_EXP(fp8) == FP8_E5M2_E_MAX) && (FP8_E5M2_GET_MANT(fp8) == 0);
}

/**
 * fp8_e5m2_iszero - Check if an FP8 E5M2 value is zero (positive or negative).
 *
 * E5M2 has two zero representations: +0 (0x00) and -0 (0x80),
 * analogous to IEEE 754 signed zeros.
 */
int fp8_e5m2_iszero(fp8_e5m2_t fp8)
{
    return (fp8 == FP8_E5M2_ZERO1) || (fp8 == FP8_E5M2_ZERO2);
}
