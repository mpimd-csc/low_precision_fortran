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

#include "fp8_e4m3/emu/fp8_e4m3_bits.h"
#include "lpf_fp8_e5m2_emu.h"
#include "fp8_e5m2_bits.h"

/**
 * @file base_arithmetic.c
 * @brief Software-emulated FP8 E5M2 arithmetic: add, subtract, multiply, divide.
 */

/**
 * Right-shift a mantissa with GRS rounding for exponent alignment.
 */
static inline fp_mant_t add_shift_helper(fp_mant_t mant, int shift){
    if (shift >= 2*(FP8_E5M2_M_BITS+1)) {
        return 0;
    } else {
        uint32_t sticky = 0;
        for (int i = 0; i < shift-1 ; ++i) {
            sticky |= (mant >> i) & 1;
        }
        uint8_t round  = (mant >> (shift-1 )) & 1;
        mant = mant >> shift;
        mant |= (round | sticky ) ? 0x01 : 0x00;
    }
    return mant;
}

/** @macro NORMALIZE_WITH_GRS
 * Normalize a mantissa to the range [2^(FP8_E5M2_M_BITS+3), 2^(FP8_E5M2_M_BITS+4)).
 */
#define NORMALIZE_WITH_GRS(mant, expo) do { \
    while (mant >= (1 << (FP8_E5M2_M_BITS + 4))) { \
        mant >>= 1; \
        expo++; \
    } \
    while (mant < (1 << (FP8_E5M2_M_BITS + 3))) { \
        mant <<= 1; \
        expo--; \
    }  \
} while(0)

/** @macro NORMALIZE_WITH_GRS_STICKY
 * Like NORMALIZE_WITH_GRS, but accumulates shifted-out bits into a sticky flag.
 */
#define NORMALIZE_WITH_GRS_STICKY(mant, expo, sticky) do { \
    while (mant >= (1 << (FP8_E5M2_M_BITS + 4))) { \
        sticky |= (mant & 0x1); \
        mant >>= 1; \
        expo++; \
    } \
    while (mant < (1 << (FP8_E5M2_M_BITS + 3))) { \
        mant <<= 1; \
        expo--; \
    }  \
} while(0)

fp8_e5m2_t fp8_e5m2_add(fp8_e5m2_t a, fp8_e5m2_t b)
{
    uint8_t sa = FP8_E5M2_GET_SIGN(a);
    uint8_t sb = FP8_E5M2_GET_SIGN(b);

    if (fp8_e5m2_isnan(a)) return a;
    if (fp8_e5m2_isnan(b)) return b;
    if (fp8_e5m2_isinf(a)) return a;
    if (fp8_e5m2_isinf(b)) return b;

    int8_t expa = FP8_E5M2_GET_EXP(a);
    int8_t expb = FP8_E5M2_GET_EXP(b);
    fp_mant_t manta = FP8_E5M2_GET_MANT(a);
    fp_mant_t mantb = FP8_E5M2_GET_MANT(b);

    if ( expa != 0 ) {
        manta = (1 << FP8_E5M2_M_BITS ) | manta;
        expa  -= FP8_E5M2_E_BIAS;
    } else {
        expa = 1 - FP8_E5M2_E_BIAS;
    }
    if ( expb != 0 ) {
        mantb = (1 << FP8_E5M2_M_BITS ) | mantb;
        expb  -= FP8_E5M2_E_BIAS;
    } else {
        expb = 1 - FP8_E5M2_E_BIAS;
    }

    manta <<= 3;
    mantb <<= 3;

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
        return FP8_E5M2_ZERO1;
    }
    NORMALIZE_WITH_GRS(result_mant, result_exp);

    uint8_t guard = (result_mant >> 2) & 1;
    uint8_t round = (result_mant >> 1) & 1;
    uint8_t sticky = result_mant & 0x1;

    uint8_t increment = 0;
    if (guard && (round || sticky || ((result_mant >> 3) & 1))) {
        increment = 1;
    }
    result_mant = (result_mant >> 3) + increment;

    if (result_mant >= (1 << (FP8_E5M2_M_BITS + 1))) {
        result_mant >>= 1;
        result_exp++;
    }

    if (result_exp < FP8_E5M2_E_MIN_U) {
        int shift = FP8_E5M2_E_MIN_U - result_exp;
        result_mant >>= shift;
        result_exp = - FP8_E5M2_E_BIAS;
    }

    result_mant = result_mant & FP8_E5M2_M_MASK;
    result_exp  = result_exp + FP8_E5M2_E_BIAS;

    if (result_exp >= FP8_E5M2_E_MAX) return (result_sign) ? FP8_E5M2_INF_NEG : FP8_E5M2_INF_POS;
    uint8_t res =  (result_sign << 7) | ((result_exp & 0x1F) << FP8_E5M2_M_BITS) | result_mant;
    return res;
}

fp8_e5m2_t fp8_e5m2_sub(fp8_e5m2_t a, fp8_e5m2_t b)
{
    if ( (a == FP8_E5M2_INF_POS && b == FP8_E5M2_INF_POS)
         || ( a == FP8_E5M2_INF_NEG && b == FP8_E5M2_INF_NEG) )
        return FP8_E5M2_NAN_POS1;
    return fp8_e5m2_add(a, b ^ (1<<FP8_E5M2_S_SHIFT));
}

fp8_e5m2_t fp8_e5m2_mul(fp8_e5m2_t a, fp8_e5m2_t b) {
    uint8_t sa = FP8_E5M2_GET_SIGN(a);
    uint8_t sb = FP8_E5M2_GET_SIGN(b);

    if (fp8_e5m2_isnan(a)) return FP8_E5M2_NAN_POS1;
    if (fp8_e5m2_isnan(b)) return FP8_E5M2_NAN_POS1;
    if (fp8_e5m2_isinf(a) && fp8_e5m2_iszero(b)) return FP8_E5M2_NAN_POS1;
    if (fp8_e5m2_isinf(b) && fp8_e5m2_iszero(a)) return FP8_E5M2_NAN_POS1;
    if (fp8_e5m2_isinf(a)) return (sa^sb) ? FP8_E5M2_INF_NEG : FP8_E5M2_INF_POS;
    if (fp8_e5m2_isinf(b)) return (sa^sb) ? FP8_E5M2_INF_NEG : FP8_E5M2_INF_POS;

    int8_t expa = FP8_E5M2_GET_EXP(a);
    int8_t expb = FP8_E5M2_GET_EXP(b);
    fp_mant_t manta = FP8_E5M2_GET_MANT(a);
    fp_mant_t mantb = FP8_E5M2_GET_MANT(b);

    if ((expa == 0 && manta == 0) || (expb == 0 && mantb == 0)) return 0;

    if ( expa != 0 ) {
        manta = (1 << FP8_E5M2_M_BITS ) | manta;
        expa  -= FP8_E5M2_E_BIAS;
    } else {
        expa = 1 - FP8_E5M2_E_BIAS;
        while ( manta < (1<<FP8_E5M2_M_BITS)) {
            manta <<= 1;
            expa--;
        }
    }
    if ( expb != 0 ) {
        mantb = (1 << FP8_E5M2_M_BITS ) | mantb;
        expb  -= FP8_E5M2_E_BIAS;
    } else {
        expb = 1 - FP8_E5M2_E_BIAS;
        while ( mantb < (1<<FP8_E5M2_M_BITS)) {
            mantb <<= 1;
            expb--;
        }
    }

    int exp_res = expa + expb;
    uint32_t mant_res = manta * mantb;
    mant_res <<= 3;

    uint8_t sticky = (mant_res & FP8_E5M2_M_MASK) ? 1:0;
    mant_res >>= FP8_E5M2_M_BITS;

    NORMALIZE_WITH_GRS_STICKY(mant_res, exp_res, sticky);

    if ( exp_res < FP8_E5M2_E_MIN_U ) {
        int shift = FP8_E5M2_E_MIN_U-exp_res;
        for ( int i = 0; i < shift; i++) {
            sticky |= (mant_res >> i) & 0x1;
        }
        mant_res >>= shift;
    }

    if ( exp_res < FP8_E5M2_E_MIN_U - FP8_E5M2_M_BITS -1) {
        return (sa^sb)<<FP8_E5M2_S_SHIFT;
    }

    if (mant_res == 0x00) {
        return ((sa^sb)<<FP8_E5M2_S_SHIFT);
    }

    sticky |= mant_res & 1;
    uint8_t guard = (mant_res >> 2) & 1;
    uint8_t round = (mant_res >> 1) & 1;
    uint8_t last  = (mant_res >> 3) & 1;
    uint8_t increment = 0;

    if ( guard && ( round || sticky || last) ) increment = 1;

    mant_res = (mant_res >> 3) + increment;

    if ( exp_res < FP8_E5M2_E_MIN_U ) {
        exp_res = 0x00;
        if ( mant_res >= (1<<FP8_E5M2_M_BITS)) {
            exp_res = FP8_E5M2_E_MIN_U + FP8_E5M2_E_BIAS;
            mant_res = 0x00;
        }
        return ((sa^sb)<< FP8_E5M2_S_SHIFT) | (exp_res << FP8_E5M2_E_SHIFT) | (mant_res & FP8_E5M2_M_MASK);
    } else {
        if ( mant_res >= (1<<(FP8_E5M2_M_BITS+1))) {
            mant_res >>= 1;
            exp_res++;
        }
        exp_res += FP8_E5M2_E_BIAS;
        if (exp_res > FP8_E5M2_E_MAX ) return (sa ^ sb) ? FP8_E5M2_INF_NEG : FP8_E5M2_INF_POS;
        return ((sa ^ sb) << FP8_E5M2_S_SHIFT) | (exp_res << FP8_E5M2_E_SHIFT) | (mant_res & FP8_E5M2_M_MASK);
    }
}

fp8_e5m2_t fp8_e5m2_div(fp8_e5m2_t a, fp8_e5m2_t b) {
    uint8_t sa = FP8_E5M2_GET_SIGN(a);
    uint8_t sb = FP8_E5M2_GET_SIGN(b);

    if (fp8_e5m2_isnan(a)) return a;

    if (fp8_e5m2_isinf(a)) {
        if ( fp8_e5m2_isinf(b)) return FP8_E5M2_NAN_POS1;
        return (sa == sb) ? FP8_E4M3_INF_POS : FP8_E4M3_INF_NEG;
    }
    if (fp8_e5m2_isnan(b)) return b;
    if (fp8_e5m2_isinf(b)) {
        if (fp8_e5m2_iszero(a)) return FP8_E5M2_ZERO1;
        return (sa^sb) ? FP8_E5M2_ZERO2 : FP8_E5M2_ZERO1;
    }

    if (fp8_e5m2_iszero(a) && fp8_e5m2_iszero(b)) return FP8_E5M2_NAN_POS1;
    if (fp8_e5m2_iszero(b)) return (sa^sb) ? FP8_E5M2_INF_NEG: FP8_E5M2_INF_POS;

    int8_t expa = FP8_E5M2_GET_EXP(a);
    int8_t expb = FP8_E5M2_GET_EXP(b);
    fp_mant_t manta = FP8_E5M2_GET_MANT(a);
    fp_mant_t mantb = FP8_E5M2_GET_MANT(b);

    if ((expa == 0 && manta == 0) || (expb == 0 && mantb == 0)) return 0;

    if ( expa != 0 ) {
        manta = (1 << FP8_E5M2_M_BITS ) | manta;
        expa  -= FP8_E5M2_E_BIAS;
    } else {
        expa = 1 - FP8_E5M2_E_BIAS;
        while ( manta < (1<<FP8_E5M2_M_BITS)) {
            manta <<= 1;
            expa--;
        }
    }
    if ( expb != 0 ) {
        mantb = (1 << FP8_E5M2_M_BITS ) | mantb;
        expb  -= FP8_E5M2_E_BIAS;
    } else {
        expb = 1 - FP8_E5M2_E_BIAS;
        while ( mantb < (1<<FP8_E5M2_M_BITS)) {
            mantb <<= 1;
            expb--;
        }
    }

    int exp_res = expa - expb;
    uint32_t mant_res = (manta << (FP8_E5M2_M_BITS + 3 + FP8_E5M2_M_BITS)) / mantb;
    uint8_t sticky = (mant_res & FP8_E5M2_M_MASK) ? 1:0;
    mant_res >>= FP8_E5M2_M_BITS;

    NORMALIZE_WITH_GRS_STICKY(mant_res, exp_res, sticky);

    if ( exp_res < FP8_E5M2_E_MIN_U ) {
        int shift = FP8_E5M2_E_MIN_U-exp_res;
        for ( int i = 0; i < shift; i++) {
            sticky |= (mant_res >> i) & 0x1;
        }
        mant_res >>= shift;
    }

    if ( exp_res < FP8_E5M2_E_MIN_U - FP8_E5M2_M_BITS -1) {
        return (sa^sb)<<FP8_E5M2_S_SHIFT;
    }

    if (mant_res == 0x00) {
        return ((sa^sb)<<FP8_E5M2_S_SHIFT);
    }

    sticky |= mant_res & 1;
    uint8_t guard = (mant_res >> 2) & 1;
    uint8_t round = (mant_res >> 1) & 1;
    uint8_t last  = (mant_res >> 3) & 1;
    uint8_t increment = 0;

    if ( guard && ( round || sticky || last) ) increment = 1;

    mant_res = (mant_res >> 3) + increment;

    if ( exp_res < FP8_E5M2_E_MIN_U ) {
        exp_res = 0x00;
        if ( mant_res >= (1<<FP8_E5M2_M_BITS)) {
            exp_res = FP8_E5M2_E_MIN_U + FP8_E5M2_E_BIAS;
            mant_res = 0x00;
        }
        return ((sa^sb)<< FP8_E5M2_S_SHIFT) | (exp_res << FP8_E5M2_E_SHIFT) | (mant_res & FP8_E5M2_M_MASK);
    } else {
        if ( mant_res >= (1<<(FP8_E5M2_M_BITS+1))) {
            mant_res >>= 1;
            exp_res++;
        }
        exp_res += FP8_E5M2_E_BIAS;
        if (exp_res > FP8_E5M2_E_MAX ) return (sa ^ sb) ? FP8_E5M2_INF_NEG : FP8_E5M2_INF_POS;
        return ((sa ^ sb) << FP8_E5M2_S_SHIFT) | (exp_res << FP8_E5M2_E_SHIFT) | (mant_res & FP8_E5M2_M_MASK);
    }
}

fp8_e5m2_t fp8_e5m2_unitary_minus(fp8_e5m2_t x)
{
    return x ^ (1 << FP8_E5M2_S_SHIFT);
}


fp8_e5m2_t fp8_e5m2_pow(fp8_e5m2_t a, fp8_e5m2_t b)
{
    float _a = fp8_e5m2_to_float(a);
    float _b = fp8_e5m2_to_float(b);
    float _c = powf(_a, _b);
    return fp8_e5m2_from_float(_c);
}

