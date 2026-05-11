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
 * @file base_arithmetic.c
 * @brief Software-emulated FP8 E4M3 arithmetic: add, subtract, multiply, divide.
 *
 * Implements four basic arithmetic operations on the OCP FP8 E4M3 format
 * using integer-only arithmetic with guard/round/sticky (GRS) rounding.
 *
 * Special-value propagation rules:
 *   - 1/+0 = +Inf,  1/-0 = -Inf
 *   - Inf + Inf = Inf,  Inf - Inf = NaN
 *   - 0/0 = NaN,  0 * Inf = NaN
 *   - x * Inf = signed Inf (sign = XOR of operand signs)
 *   - Any NaN operand propagates NaN
 */

/**
 * Right-shift a mantissa with GRS rounding for exponent alignment.
 *
 * Used when aligning mantissae of different exponents before addition.
 * Bits shifted out beyond the round position are OR'd into a sticky bit;
 * if the round or sticky bits are nonzero, the LSB of the shifted result
 * is set (round-up).
 *
 * If @p shift is so large that no significant bits remain (>= 2*(FP8_E4M3_M_BITS+1)),
 * returns zero.
 *
 * @param mant  Mantissa with GRS extension bits (guard, round, sticky).
 * @param shift Number of positions to right-shift.
 * @return      Shifted mantissa with rounding applied.
 */
static inline fp_mant_t add_shift_helper(fp_mant_t mant, int shift){
    if (shift >= 2*(FP8_E4M3_M_BITS+1)) {
        return 0;
    } else {
        uint32_t sticky = 0;
        for (int i = 0; i < shift-1 ; ++i) {
            sticky |= (mant >> i) & 1;
        }
        uint8_t round  = (mant >> (shift-1 )) & 0x03;
        mant = mant >> shift;
        mant |= (round | sticky ) ? 0x01 : 0x00;
    }
    return mant;
}

/** @macro NORMALIZE_WITH_GRS
 * Normalize a mantissa with 3 GRS extension bits to the range [2^(FP8_E4M3_M_BITS+3), 2^(FP8_E4M3_M_BITS+4)).
 *
 * Shifts the mantissa right (and increments the exponent) if it overflowed,
 * or left (decrements the exponent) if it is too small, so the result occupies
 * the bit positions [FP8_E4M3_M_BITS+3 .. FP8_E4M3_M_BITS+1 .. 0] = 1.xxx GRS.
 *
 * @param mant  Mantissa value (modified in place).
 * @param expo  Exponent value (modified in place).
 */
#define NORMALIZE_WITH_GRS(mant, expo) do { \
    while (mant >= (1 << (FP8_E4M3_M_BITS + 4))) { \
        mant >>= 1; \
        expo++; \
    } \
    while (mant < (1 << (FP8_E4M3_M_BITS + 3))) { \
        mant <<= 1; \
        expo--; \
    }  \
} while(0)

/** @macro NORMALIZE_WITH_GRS_STICKY
 * Like NORMALIZE_WITH_GRS, but accumulates shifted-out bits into a sticky flag.
 *
 * When right-shifting the mantissa during normalization, the LSB is OR'd
 * into @p sticky so that no rounding information is lost.  Left-shifting
 * does not affect sticky (the bits shifted in are zeros).
 *
 * @param mant   Mantissa value (modified in place).
 * @param expo   Exponent value (modified in place).
 * @param sticky Sticky bit accumulator (modified in place).
 */
#define NORMALIZE_WITH_GRS_STICKY(mant, expo, sticky) do { \
    while (mant >= (1 << (FP8_E4M3_M_BITS + 4))) { \
        sticky |= (mant & 0x1); \
        mant >>= 1; \
        expo++; \
    } \
    while (mant < (1 << (FP8_E4M3_M_BITS + 3))) { \
        mant <<= 1; \
        expo--; \
    }  \
} while(0)



/**
 * Add two FP8 E4M3 values (a + b) with GRS rounding.
 *
 * Handles NaN/Inf propagation, sign-magnitude addition/subtraction when
 * signs differ, exponent alignment via add_shift_helper(), normalization,
 * and rounding.  Returns signed zero when the result is zero.
 *
 * @param a  First operand.
 * @param b  Second operand.
 * @return   FP8 E4M3 result of a + b.
 */
fp8_e4m3_t fp8_e4m3_add(fp8_e4m3_t a, fp8_e4m3_t b)
{
    uint8_t sa = FP8_E4M3_GET_SIGN(a);
    uint8_t sb = FP8_E4M3_GET_SIGN(b);

    if (fp8_e4m3_isnan(a)) return a;
    if (fp8_e4m3_isnan(b)) return b;
    if (fp8_e4m3_isinf(a)) return a;
    if (fp8_e4m3_isinf(b)) return b;

    int8_t expa = FP8_E4M3_GET_EXP(a);
    int8_t expb = FP8_E4M3_GET_EXP(b);
    fp_mant_t manta = FP8_E4M3_GET_MANT(a);
    fp_mant_t mantb = FP8_E4M3_GET_MANT(b);

    /* Add hidden 1 if necessary */
    if ( expa != 0 ) {
        manta = (1 << FP8_E4M3_M_BITS ) | manta;
        expa  -= FP8_E4M3_E_BIAS;
    } else {
        expa = 1 - FP8_E4M3_E_BIAS;
    }
    if ( expb != 0 ) {
        mantb = (1 << FP8_E4M3_M_BITS ) | mantb;
        expb  -= FP8_E4M3_E_BIAS;
    } else {
        expb = 1 - FP8_E4M3_E_BIAS;
    }

    /* add GRS bits */
    manta <<= 3;
    mantb <<= 3;

    /* Adjust the Mantissa */
    if (expa > expb) {
        int shift = expa - expb;
        mantb = add_shift_helper(mantb, shift);
        expb = expa;
    } else if (expb > expa) {
        int shift = expb - expa;
        manta = add_shift_helper(manta, shift);
        expa = expb;
    }

    uint8_t result_sign = 0;
    fp_mant_t result_mant;
    if (sa == sb) {
        result_mant = manta + mantb;
        result_sign = sa;
    } else {
        if (manta > mantb) {
            result_mant = manta - mantb;
            result_sign = sa;
        } else {
            result_mant = mantb - manta;
            result_sign = sb;
        }
    }
    int8_t result_exp = expa;

    if ( result_mant == 0x00 ) {
        return FP8_E4M3_ZERO1;
    }
    NORMALIZE_WITH_GRS(result_mant, result_exp);

    /* Round */
    uint8_t guard = (result_mant >> 2) & 1;
    uint8_t round = (result_mant >> 1) & 1;
    uint8_t sticky = result_mant & 0x1;

    uint8_t increment = 0;
    if (guard && (round || sticky || ((result_mant >> 3) & 1))) {
        increment = 1;
    }
    result_mant = (result_mant >> 3) + increment;

    // Check overflow in mantisse
    if (result_mant >= (1 << (FP8_E4M3_M_BITS + 1))) {
        result_mant >>= 1;
        result_exp++;
    }

    /* We end in the subnormal area */
    if (result_exp < FP8_E4M3_E_MIN_U) {
        int shift = FP8_E4M3_E_MIN_U - result_exp;
        result_mant >>= shift;
        result_exp = - FP8_E4M3_E_BIAS;
    }

    result_mant = result_mant & FP8_E4M3_M_MASK;
    result_exp  = result_exp + FP8_E4M3_E_BIAS;

    if (result_exp >= FP8_E4M3_E_MAX && result_mant >= 0x07) return (result_sign) ? FP8_E4M3_NAN_NEG : FP8_E4M3_NAN_POS;
    uint8_t res =  (result_sign << 7) | ((result_exp & 0xF) << FP8_E4M3_M_BITS) | result_mant;
    return res;

}

/**
 * Subtract two FP8 E4M3 values (a - b).
 *
 * Implemented by flipping the sign bit of @p b and delegating to
 * fp8_e4m3_add().  The only special case handled directly is
 * Inf - Inf, which returns NaN (per the rules at the top of this file).
 *
 * @param a  Minuend.
 * @param b  Subtrahend.
 * @return   FP8 E4M3 result of a - b, or NaN for Inf - Inf.
 */
fp8_e4m3_t fp8_e4m3_sub(fp8_e4m3_t a, fp8_e4m3_t b)
{
    if ( (a == FP8_E4M3_INF_POS && b == FP8_E4M3_INF_NEG)
         || ( a == FP8_E4M3_INF_NEG && b == FP8_E4M3_INF_POS) )
        return FP8_E4M3_NAN_POS;
    return fp8_e4m3_add(a, b ^ (1<<FP8_E4M3_S_SHIFT));
}


/**
 * Multiply two FP8 E4M3 values (a * b) with GRS rounding.
 *
 * Decomposes both operands into sign, unbiased exponent, and mantissa,
 * multiplies the mantissae, adds exponents, and normalizes the product
 * using NORMALIZE_WITH_GRS_STICKY.  Subnormal inputs are normalized
 * before multiplication; subnormal results are handled with sticky-bit
 * accumulation.  Overflow saturates to NaN (E4M3 has no infinity).
 * Zero * Inf returns NaN.
 *
 * @param a  Multiplicand.
 * @param b  Multiplier.
 * @return   FP8 E4M3 result of a * b.
 */
fp8_e4m3_t fp8_e4m3_mul(fp8_e4m3_t a, fp8_e4m3_t b) {
    uint8_t sa = FP8_E4M3_GET_SIGN(a);
    uint8_t sb = FP8_E4M3_GET_SIGN(b);

    if (fp8_e4m3_isnan(a)) return FP8_E4M3_NAN_POS;
    if (fp8_e4m3_isnan(b)) return FP8_E4M3_NAN_POS;
    if (fp8_e4m3_isinf(a)) return (sa^sb) ? FP8_E4M3_INF_NEG: FP8_E4M3_INF_POS;
    if (fp8_e4m3_isinf(b)) return (sa^sb) ? FP8_E4M3_INF_NEG: FP8_E4M3_INF_POS;
    if (fp8_e4m3_isinf(a) && fp8_e4m3_iszero(b)) return FP8_E4M3_NAN_POS;
    if (fp8_e4m3_isinf(b) && fp8_e4m3_iszero(a)) return FP8_E4M3_NAN_POS;

      int8_t expa = FP8_E4M3_GET_EXP(a);
    int8_t expb = FP8_E4M3_GET_EXP(b);
    fp_mant_t manta = FP8_E4M3_GET_MANT(a);
    fp_mant_t mantb = FP8_E4M3_GET_MANT(b);

    if ((expa == 0 && manta == 0) || (expb == 0 && mantb == 0)) return 0;

    /* Add hidden 1 if necessary */
    if ( expa != 0 ) {
        manta = (1 << FP8_E4M3_M_BITS ) | manta;
        expa  -= FP8_E4M3_E_BIAS;
    } else {
        /* a is subnormal*/
        expa = 1 - FP8_E4M3_E_BIAS;
        while ( manta < (1<<FP8_E4M3_M_BITS)) {
            manta <<= 1;
            expa--;
        }
    }
    if ( expb != 0 ) {
        mantb = (1 << FP8_E4M3_M_BITS ) | mantb;
        expb  -= FP8_E4M3_E_BIAS;
    } else {
        expb = 1 - FP8_E4M3_E_BIAS;
        while ( mantb < (1<<FP8_E4M3_M_BITS)) {
            mantb <<= 1;
            expb--;
        }
    }

    int exp_res = expa + expb;

    uint32_t mant_res = manta * mantb;

    /* Add GRS */
    mant_res <<= 3;

    uint8_t sticky = (mant_res & FP8_E4M3_M_MASK) ? 1:0;
    mant_res >>= FP8_E4M3_M_BITS;

    NORMALIZE_WITH_GRS_STICKY(mant_res, exp_res, sticky);

    /* Result is subnormal */
    if ( exp_res < FP8_E4M3_E_MIN_U ) {
        int shift = FP8_E4M3_E_MIN_U-exp_res;
        for ( int i = 0; i < shift; i++) {
            sticky |= (mant_res >> i) & 0x1;
        }
        mant_res >>= shift;
    }

    /* Result is numerically zero */
    if ( exp_res < FP8_E4M3_E_MIN_U - FP8_E4M3_M_BITS -1) {
        return (sa^sb)<<FP8_E4M3_S_SHIFT;
    }

    /* We generate a zero */
    if (mant_res == 0x00) {
        return ((sa^sb)<<FP8_E4M3_S_SHIFT);
    }

    /* rounding */
    sticky |= mant_res & 1;
    uint8_t guard = (mant_res >> 2) & 1;
    uint8_t round = (mant_res >> 1) & 1;
    uint8_t last  = (mant_res >> 3) & 1;
    uint8_t increment = 0;

    if ( guard && ( round || sticky || last) ) increment = 1;

    /* Remove GRS and increment if necessary */
    mant_res = (mant_res >> 3) + increment;

    if ( exp_res < FP8_E4M3_E_MIN_U ) {
        exp_res = 0x00;
        if ( mant_res >= (1<<FP8_E4M3_M_BITS)) {
            exp_res = FP8_E4M3_E_MIN_U + FP8_E4M3_E_BIAS;
            mant_res = 0x00;
        }
        return ((sa^sb)<< FP8_E4M3_S_SHIFT) | (exp_res << FP8_E4M3_E_SHIFT) | (mant_res & FP8_E4M3_M_MASK);
    } else {
        if ( mant_res >= (1<<(FP8_E4M3_M_BITS+1))) {
            mant_res >>= 1;
            exp_res++;
        }
        exp_res += FP8_E4M3_E_BIAS;
        if (exp_res > FP8_E4M3_E_MAX ) return (sa ^ sb) ? FP8_E4M3_INF_NEG : FP8_E4M3_INF_POS;
        return ((sa ^ sb) << FP8_E4M3_S_SHIFT) | (exp_res << FP8_E4M3_M_BITS) | (mant_res & FP8_E4M3_M_MASK);
    }
}

/**
 * Divide two FP8 E4M3 values (a / b) with GRS rounding.
 *
 * Decomposes operands into sign, exponent, and mantissa, computes the
 * quotient via integer division with extended precision (the dividend is
 * left-shifted by FP8_E4M3_M_BITS + 3 + COMPUTFP8_E4M3_E_BITS to produce fractional bits),
 * then normalizes and rounds.  Division by zero returns signed NaN; 0/0
 * returns NaN.  Subnormal inputs are normalized before division; subnormal
 * results are flushed or shifted with sticky-bit accumulation.
 *
 * @param a  Dividend.
 * @param b  Divisor (must not be zero).
 * @return   FP8 E4M3 result of a / b, or NaN on 0/0 or div-by-zero.
 */
fp8_e4m3_t fp8_e4m3_div(fp8_e4m3_t a, fp8_e4m3_t b) {
    uint8_t sa = FP8_E4M3_GET_SIGN(a);
    uint8_t sb = FP8_E4M3_GET_SIGN(b);

    if (fp8_e4m3_isnan(a)) return a;
    if (fp8_e4m3_isinf(a)) return a;
    if (fp8_e4m3_isnan(b)) return b;
    if (fp8_e4m3_isinf(b)) return b;

    if (fp8_e4m3_iszero(a) && fp8_e4m3_iszero(b)) return FP8_E4M3_NAN_POS;
    if (fp8_e4m3_iszero(b)) return (sa^sb) ? FP8_E4M3_INF_NEG: FP8_E4M3_INF_POS;

    int8_t expa = FP8_E4M3_GET_EXP(a);
    int8_t expb = FP8_E4M3_GET_EXP(b);
    fp_mant_t manta = FP8_E4M3_GET_MANT(a);
    fp_mant_t mantb = FP8_E4M3_GET_MANT(b);

    if ((expa == 0 && manta == 0) || (expb == 0 && mantb == 0)) return 0;

    /* Add hidden 1 if necessary */
    if ( expa != 0 ) {
        manta = (1 << FP8_E4M3_M_BITS ) | manta;
        expa  -= FP8_E4M3_E_BIAS;
    } else {
        /* a is subnormal*/
        expa = 1 - FP8_E4M3_E_BIAS;
        while ( manta < (1<<FP8_E4M3_M_BITS)) {
            manta <<= 1;
            expa--;
        }
    }
    if ( expb != 0 ) {
        mantb = (1 << FP8_E4M3_M_BITS ) | mantb;
        expb  -= FP8_E4M3_E_BIAS;
    } else {
        expb = 1 - FP8_E4M3_E_BIAS;
        while ( mantb < (1<<FP8_E4M3_M_BITS)) {
            mantb <<= 1;
            expb--;
        }
    }

#define COMPUTFP8_E4M3_E_BITS FP8_E4M3_M_BITS
#define COMPUTFP8_E4M3_E_BITFP8_E4M3_S_MASK FP8_E4M3_M_MASK
    int exp_res = expa - expb;

    uint32_t mant_res = (manta<< ( FP8_E4M3_M_BITS + 3 + COMPUTFP8_E4M3_E_BITS) ) / mantb;

    uint8_t sticky = (mant_res & COMPUTFP8_E4M3_E_BITFP8_E4M3_S_MASK) ? 1:0;
    mant_res >>= COMPUTFP8_E4M3_E_BITS;
#undef COMPUTFP8_E4M3_E_BITS
#undef COMPUTFP8_E4M3_E_BITFP8_E4M3_S_MASK

    NORMALIZE_WITH_GRS_STICKY(mant_res, exp_res, sticky);

    /* Result is subnormal */
    if ( exp_res < FP8_E4M3_E_MIN_U ) {
        int shift = FP8_E4M3_E_MIN_U-exp_res;
        for ( int i = 0; i < shift; i++) {
            sticky |= (mant_res >> i) & 0x1;
        }
        mant_res >>= shift;
    }

    /* Result is numerically zero */
    if ( exp_res < FP8_E4M3_E_MIN_U - FP8_E4M3_M_BITS -1) {
        return (sa^sb)<<FP8_E4M3_S_SHIFT;
    }

    /* We generate a zero */
    if (mant_res == 0x00) {
        return ((sa^sb)<<FP8_E4M3_S_SHIFT);
    }

    /* rounding */
    sticky |= mant_res & 1;
    uint8_t guard = (mant_res >> 2) & 1;
    uint8_t round = (mant_res >> 1) & 1;
    uint8_t last  = (mant_res >> 3) & 1;
    uint8_t increment = 0;

    if ( guard && ( round || sticky || last) ) increment = 1;

    /* Remove GRS and increment if necessary */
    mant_res = (mant_res >> 3) + increment;

    if ( exp_res < FP8_E4M3_E_MIN_U ) {
        exp_res = 0x00;
        if ( mant_res >= (1<<FP8_E4M3_M_BITS)) {
            exp_res = FP8_E4M3_E_MIN_U + FP8_E4M3_E_BIAS;
            mant_res = 0x00;
        }
        return ((sa^sb)<< FP8_E4M3_S_SHIFT) | (exp_res << FP8_E4M3_E_SHIFT) | (mant_res & FP8_E4M3_M_MASK);
    } else {
        if ( mant_res >= (1<<(FP8_E4M3_M_BITS+1))) {
            mant_res >>= 1;
            exp_res++;
        }
        exp_res += FP8_E4M3_E_BIAS;
        if (exp_res > FP8_E4M3_E_MAX ) return (sa ^ sb) ? FP8_E4M3_INF_NEG : FP8_E4M3_INF_POS;
        return ((sa ^ sb) << FP8_E4M3_S_SHIFT) | (exp_res << FP8_E4M3_M_BITS) | (mant_res & FP8_E4M3_M_MASK);
    }
}

/**
 * Negate an FP8 E4M3 value (unary minus).
 *
 * Flips the sign bit (bit 7) of @p x via XOR.  This is valid for all
 * FP8 E4M3 encodings: normal numbers, subnormals, signed zeros, and NaNs
 * all have their sign at bit 7, so a single XOR with (1 << FP8_E4M3_S_SHIFT)
 * produces the correctly-negated encoding.  No exponent or mantissa
 * bits are affected.
 *
 * @param x  FP8 E4M3 value to negate.
 * @return   -x in FP8 E4M3 encoding.
 */
fp8_e4m3_t fp8_e4m3_unitary_minus(fp8_e4m3_t x)
{
    return x ^ (1 << FP8_E4M3_S_SHIFT);
}
