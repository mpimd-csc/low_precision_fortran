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
#ifndef LPF_FP8_E5M2_BITS_H
#define LPF_FP8_E5M2_BITS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <math.h>

    typedef uint16_t fp_mant_t;  /* Wide type for mantissa arithmetic (avoids overflow during intermediate ops) */

    /* Bit-field masks, shifts, and widths
     *
     *   FP8_E5M2_S_SHIFT=7  FP8_E5M2_E_SHIFT=2        FP8_E5M2_M_SHIFT=0 (implicit)
     *      |       |                   |
     *      v       v                   v
     *   [S] [E E E E E] [M M]
     *    7   6 5 4 3 2   1 0
     */

#define FP8_E5M2_M_MASK 0x03  /* Mantissa bit mask: bits [1:0] */
#define FP8_E5M2_M_BITS 2     /* Mantissa width in bits */
#define FP8_E5M2_E_MASK 0x7C  /* Exponent bit mask: bits [6:2] */
#define FP8_E5M2_E_SHIFT FP8_E5M2_M_BITS /* Exponent field position (shift amount from bit 0) */
#define FP8_E5M2_E_BITS 5     /* Exponent width in bits */
#define FP8_E5M2_S_MASK 0x80  /* Sign bit mask: bit [7] */
#define FP8_E5M2_S_SHIFT (FP8_E5M2_E_BITS+FP8_E5M2_M_BITS) /* Sign field position (shift amount from bit 0) */

    /* Exponent parameters
     *
     * The 5-bit exponent field ranges from 0 to 31.
     *   - Field = 0  : subnormal numbers (no implicit leading 1; effective exponent = 1 - BIAS = -14)
     *   - Field = 1..30: normal numbers  (implicit leading 1; exponent = field - BIAS)
     *   - Field = 31 : Special values (Infinity and NaN)
     *
     * FP8_E5M2_E_MAX_U / FP8_E5M2_E_MIN_U give the unbiased (actual) exponent range for normal numbers.
     */

#define FP8_E5M2_E_BIAS 15    /* Exponent bias: actual exponent = field - FP8_E5M2_E_BIAS */
#define FP8_E5M2_E_MAX  31    /* Maximum exponent field value (all 1s) */
#define FP8_E5M2_E_MAX_U 15   /* Maximum unbiased normal exponent (field=30 minus FP8_E5M2_E_BIAS) */
#define FP8_E5M2_E_MIN_U -14  /* Minimum unbiased normal exponent (field=1 minus FP8_E5M2_E_BIAS) */

#define FP8_E5M2_MANTISSA_MULTIPLIER  0.25  /* 2^(-FP8_E5M2_M_BITS): scales raw mantissa bits to fractional value */

    /*
     * Special values and boundary encodings
     *
     * Following IEEE 754 conventions for E5M2.
     *
     * Zeros:     Two representations (+0 and -0).
     * Infinity:  Exponent=31, Mantissa=0. Both +Inf and -Inf are defined.
     * NaN:       Exponent=31, Mantissa > 0. Multiple NaN encodings exist.
     * Max norm:  Exponent=30, Mantissa=3  =>  +/- 2^15 * 1.75 = +/- 57344
     * Min norm:  Exponent=1,  Mantissa=0  =>  +/- 2^(-14)
     * Max sub:   Exponent=0,  Mantissa=3  =>  +/- 0.75 * 2^(-14)
     * Min sub:   Exponent=0,  Mantissa=1  =>  +/- 2^(-16)
     */

#define FP8_E5M2_ZERO1    0x00  /* Positive zero: sign=0, exp=0, mant=0 */
#define FP8_E5M2_ZERO2    0x80  /* Negative zero: sign=1, exp=0, mant=0 */
#define FP8_E5M2_NAN_POS1  0x7D  /* Positive NaN: sign=0, exp=31, mant=1 */
#define FP8_E5M2_NAN_POS2  0x7E  /* Positive NaN: sign=0, exp=31, mant=2 */
#define FP8_E5M2_NAN_POS3  0x7F  /* Positive NaN: sign=0, exp=31, mant=3 */
#define FP8_E5M2_NAN_NEG1  0xFD  /* Negative NaN: sign=1, exp=31, mant=1 */
#define FP8_E5M2_NAN_NEG2  0xFE  /* Negative NaN: sign=1, exp=31, mant=2 */
#define FP8_E5M2_NAN_NEG3  0xFF  /* Negative NaN: sign=1, exp=31, mant=3 */
#define FP8_E5M2_INF_POS   0x7C  /* Positive Infinity: sign=0, exp=31, mant=0 */
#define FP8_E5M2_INF_NEG   0xFC  /* Negative Infinity: sign=1, exp=31, mant=0 */
#define FP8_E5M2_MAX_NUM_POS 0x7B   /* Largest positive normal: 57344 */
#define FP8_E5M2_MAX_NUM_NEG 0xFB   /* Largest negative normal: -57344 */
#define FP8_E5M2_MIN_NUM_POS 0x04   /* Smallest positive normal: 2^-14 */
#define FP8_E5M2_MIN_NUM_NEG 0x84   /* Smallest negative normal: -2^-14 */
#define FP8_E5M2_MAX_SUB_NUM_POS 0x03   /* Largest positive subnormal: 0.75 * 2^-14 */
#define FP8_E5M2_MAX_SUB_NUM_NEG 0x83   /* Largest negative subnormal: -0.75 * 2^-14 */
#define FP8_E5M2_MIN_SUB_NUM_POS 0x01   /* Smallest positive subnormal: 2^-16 */
#define FP8_E5M2_MIN_SUB_NUM_NEG 0x81   /* Smallest negative subnormal: -2^-16 */
#define FP8_E5M2_EPSILON 0x34 /* Epsilon 0.25 = 2 ^ -2 */

    /* Field extraction helpers
     *
     * FP8_E5M2_GET_SIGN  -> 0 or 1 (7)
     * FP8_E5M2_GET_EXP   -> raw exponent field (2..6)
     * FP8_E5M2_GET_MANT  -> raw mantissa field (0..1)
     *
     * To reconstruct the sign-adjusted value:
     *   sign  = FP8_E5M2_GET_SIGN(x)  ?  -1.0 : 1.0;
     *   exp   = FP8_E5M2_GET_EXP(x);
     *   mant  = FP8_E5M2_GET_MANT(x);
     *   if (exp == 0)  -> subnormal:  sign * (mant * FP8_E5M2_MANTISSA_MULTIPLIER) * 2^(1 - FP8_E5M2_E_BIAS)
     *   else if (exp == 31) -> special:  Infinity or NaN
     *   else           -> normal:     sign * (1.0 + mant * FP8_E5M2_MANTISSA_MULTIPLIER) * 2^(exp - FP8_E5M2_E_BIAS)
     */

#define FP8_E5M2_GET_SIGN(x) (((x) & FP8_E5M2_S_MASK) >> (FP8_E5M2_S_SHIFT))
#define FP8_E5M2_GET_EXP(x)  (((x) & FP8_E5M2_E_MASK) >> (FP8_E5M2_E_SHIFT))
#define FP8_E5M2_GET_MANT(x) (((x) & FP8_E5M2_M_MASK))

#ifdef __cplusplus
}
#endif

#endif /* LPF_FP8_E5M2_BITS_H */
