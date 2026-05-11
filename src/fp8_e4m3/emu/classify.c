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

#include "fp8_e4m3.h"
#include "fp8_e4m3_bits.h"

/**
 * fp8_e4m3_isnan - Check if an FP8 E4M3 value is Not-a-Number.
 *
 * In E4M3, NaN is encoded as exponent=15, mantissa=7.
 * Both 0x7F (+NaN) and 0xFF (-NaN) match this pattern.
 * Note: E4M3 has no dedicated infinity; INF_* aliases map to NaN.
 */
int fp8_e4m3_isnan(fp8_e4m3_t fp8)
{
    return (fp8 == FP8_E4M3_NAN_POS) || (fp8 == FP8_E4M3_NAN_NEG);
}

/**
 * fp8_e4m3_isinf - Check if an FP8 E4M3 value represents infinity.
 *
 * Since E4M3 has no dedicated infinity encoding,
 * FP8_E4M3_INF_POS and FP8_E4M3_INF_NEG are aliased to the NaN encodings,
 * so this function is equivalent to isnan for the bit-patterns that would
 * be infinity in IEEE 754.
 */
int fp8_e4m3_isinf(fp8_e4m3_t fp8)
{
    return (fp8 == FP8_E4M3_INF_POS) || (fp8 == FP8_E4M3_INF_NEG);
}

/**
 * fp8_e4m3_iszero - Check if an FP8 E4M3 value is zero (positive or negative).
 *
 * E4M3 has two zero representations: +0 (0x00) and -0 (0x80),
 * analogous to IEEE 754 signed zeros.
 */
int fp8_e4m3_iszero(fp8_e4m3_t fp8)
{
    return (fp8 == FP8_E4M3_ZERO1) || (fp8 == FP8_E4M3_ZERO2);
}


