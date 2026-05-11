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
#include "fp32.h"

/**
 * Construct an FP8 E4M3 value from separate sign, exponent, and mantissa fields.
 *
 * Packs the three components into the 8-bit FP8 E4M3 encoding using the format
 * [S][EEEE][MMM]. If the resulting encoding would be NaN (exponent field = 15
 * and mantissa = 7), the NaN representation is returned instead.
 *
 * @param sign  Sign bit: 0 for positive, 1 for negative.
 * @param expi  Unbiased exponent in range [-6, 8]. Added to FP8_E4M3_E_BIAS (7) to
 *              produce the biased exponent field.
 * @param mant  Mantissa bits (0..7). Only the low 3 bits are used.
 * @return     Packed fp8_e4m3_t value, or NaN if the encoding is all-ones.
 */
fp8_e4m3_t fp8_e4m3_from_components(uint8_t sign, int8_t expi, uint8_t mant)
{
    uint8_t expb = expi + FP8_E4M3_E_BIAS;
    if ( expb >= FP8_E4M3_E_MAX && mant >= 0x07) return (sign) ? (FP8_E4M3_NAN_NEG) : (FP8_E4M3_NAN_POS);
    return (sign << FP8_E4M3_S_SHIFT) | ( expb << FP8_E4M3_E_SHIFT) | (mant & FP8_E4M3_M_MASK);
}


/**
 * Convert an FP8 E4M3 value to a 32-bit IEEE float.
 *
 * Handles all special cases: positive/negative zero maps to +/-0.0f, and
 * NaN encodings map to NAN. Normal numbers reconstruct the implicit leading 1
 * and apply the biased exponent. Subnormal numbers (exponent field = 0) have
 * no implicit bit and use the fixed exponent (1 - FP8_E4M3_E_BIAS).
 *
 * @param fp8  Raw FP8 E4M3 value.
 * @return     Equivalent float value.
 */
float fp8_e4m3_to_float(fp8_e4m3_t fp8) {
    float result = 0.0f;
    int8_t expo = FP8_E4M3_GET_EXP(fp8);
    uint8_t sign = FP8_E4M3_GET_SIGN(fp8);
    uint8_t mantinssa = FP8_E4M3_GET_MANT(fp8);

    // Handle special cases (zero, infinity, NaN)
    if (fp8 == FP8_E4M3_ZERO1) {
        return 0.0f;
    } else if (fp8 == FP8_E4M3_ZERO2 ) {
        return -0.0f;
    } else if (fp8 == FP8_E4M3_NAN_POS || fp8 == FP8_E4M3_NAN_NEG) {
        return NAN;
    }


    if ( expo > 0  ) {
        /* We have an implicit 1 in front of the mantissa */
        expo -= FP8_E4M3_E_BIAS;
        float m = ( 1.0f + (FP8_E4M3_FP8_E4M3_MANTISSA_MULTIPLIER  * mantinssa ));
        float v = ldexpf(m, expo);
        return (sign)?(-v):v;

    } else { // Subnormal numbers, no implicit digit
        float m = ( 0.0f + (FP8_E4M3_FP8_E4M3_MANTISSA_MULTIPLIER  * mantinssa ));
        float v = ldexpf(m, 1-FP8_E4M3_E_BIAS);
        return (sign)?(-v):v;
    }

}

/**
 * Convert a 32-bit IEEE float to an FP8 E4M3 value with round-to-nearest-even.
 *
 * Decomposes the float into sign, exponent, and mantissa fields, then rounds
 * the wider mantissa to fit 3 bits using a guard/round/sticky scheme. Three
 * cases are handled:
 *
 *   1. Overflow  (exponent > max normal): returns NaN (E4M3 has no infinity).
 *   2. Subnormal (exponent below min normal): shifts mantissa into subnormal
 *      range with rounding; flushes to signed zero if too small.
 *   3. Normal    (in representable range): rounds the FP32 mantissa to 3 bits,
 *      re-normalizes if rounding overflowed, and saturates to NaN if the
 *      result would exceed the largest normal value.
 *
 * @param x  Input float value.
 * @return   Nearest FP8 E4M3 encoding, or NaN on overflow.
 */
fp8_e4m3_t fp8_e4m3_from_float(float x){
    union {
        float f32;
        uint32_t i32;
    } fp32 = {.f32 = x};

    if (x == 0.0f) return FP8_E4M3_ZERO1;
    if (x == -0.0f) return FP8_E4M3_ZERO2;
    if (x == INFINITY) return FP8_E4M3_INF_POS;
    if (x == -INFINITY) return FP8_E4M3_INF_NEG;
    if (x == NAN) return FP8_E4M3_NAN_POS;
    if (x == -NAN) return FP8_E4M3_NAN_NEG;

    uint32_t sign = FP32_GET_SIGN(fp32.i32);
    uint32_t  expo = FP32_GET_EXP(fp32.i32);
    uint32_t mant = FP32_GET_MANT(fp32.i32);
    int32_t expu = expo - FP32_E_BIAS;

    /* Out of range, absolute value too large */
    if ( expu > FP8_E4M3_E_MAX_U ) {
        return (sign)?(FP8_E4M3_INF_NEG):(FP8_E4M3_INF_POS);
    }
    /* Eventually a subnormal number */
    if ( expu < FP8_E4M3_E_MIN_U-FP8_E4M3_M_BITS-1) {
         return (sign << FP8_E4M3_S_SHIFT);
    }
    if ( expu < FP8_E4M3_E_MIN_U) {
        mant = (1 << FP32_M_BITS) | mant;   // Add the implicit 1
        int shift = (FP32_M_BITS - FP8_E4M3_M_BITS + (FP8_E4M3_E_MIN_U - expu));
        uint8_t sticky = 0;
        for ( int i = 0; i < shift-2; i++) {
            sticky |= ( mant >> i ) & 0x01;
        }
        uint8_t guard = (mant >> (shift-1)) & 1;
        uint8_t round = (mant >> (shift-2)) & 1;
        uint8_t increment = 0;
        if ( guard && ( round || sticky || (mant >> shift) & 0x01)) increment = 1;

        mant = (mant >> shift) + increment;

        expu = 0;
        if (mant >= ( 1 << FP8_E4M3_M_BITS )){
            expu = FP8_E4M3_E_MIN_U + FP8_E4M3_E_BIAS;
        }
        return (sign << FP8_E4M3_S_SHIFT) | (expu << FP8_E4M3_M_BITS) | ((mant) & FP8_E4M3_M_MASK);
    }

    /* in the range */
    /* add the implicit 1 */
    mant |= (1<<FP32_M_BITS);

    int shift = (FP32_M_BITS - FP8_E4M3_M_BITS);
    uint8_t sticky = 0;
    for (int i = 0; i < shift - 2; i++) {
        sticky |= (mant >> i) & 0x01;
    }
    uint8_t round = (mant >> (shift-2)) & 0x01;
    uint8_t guard = (mant >> (shift-1)) & 0x01;

    uint32_t increment = 0;
    if ( guard && ( round || sticky || (mant >> shift) & 0x01)) increment = 1;

    mant = (mant >> shift) + increment;

    /* normalize again */
    if (mant >= (1<<(FP8_E4M3_M_BITS+1))){
        expu++;
        mant >>= 1;
    }

    /* Overflow after rounding */
    if ( expu >= FP8_E4M3_E_MAX_U && (mant & FP8_E4M3_M_MASK)  == 0x07) {
        return (sign)?(FP8_E4M3_INF_NEG):(FP8_E4M3_INF_POS);
    }

    expu = expu + FP8_E4M3_E_BIAS;
    uint8_t res = (sign << FP8_E4M3_S_SHIFT) | (expu << FP8_E4M3_E_SHIFT) |((mant) & FP8_E4M3_M_MASK);
    return res;
}


