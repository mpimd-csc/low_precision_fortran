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
#ifndef LPF_FP8_E4M3_BITS_H
#define LPF_FP8_E4M3_BITS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <math.h>

    typedef uint16_t fp_mant_t;  /* Wide type for mantissa arithmetic (avoids overflow during intermediate ops) */

    /* Bit-field masks, shifts, and widths
     *
     *   FP8_E4M3_S_SHIFT=7  FP8_E4M3_E_SHIFT=3        M_SHIFT=0 (implicit)
     *      |       |                   |
     *      v       v                   v
     *   [S] [E E E E] [M M M]
     *    7   6 5 4 3   2 1 0
     */

#define FP8_E4M3_M_MASK 0x07  /* Mantissa bit mask: bits [2:0] */
#define FP8_E4M3_M_BITS 3     /* Mantissa width in bits */
#define FP8_E4M3_E_MASK 0x78  /* Exponent bit mask: bits [6:3] */
#define FP8_E4M3_E_SHIFT FP8_E4M3_M_BITS /* Exponent field position (shift amount from bit 0) */
#define FP8_E4M3_E_BITS 4     /* Exponent width in bits */
#define FP8_E4M3_S_MASK 0x80  /* Sign bit mask: bit [7] */
#define FP8_E4M3_S_SHIFT (FP8_E4M3_E_BITS+FP8_E4M3_M_BITS) /* Sign field position (shift amount from bit 0) */

    /* Exponent parameters
     *
     * The 4-bit exponent field ranges from 0 to 15.
     *   - Field = 0  : subnormal numbers (no implicit leading 1; effective exponent = 1 - BIAS = -6)
     *   - Field = 1..15: normal numbers  (implicit leading 1; exponent = field - BIAS)
     *   - Field = 15, mantissa = 7 (all 1s): NaN
     *   - Field = 15, mantissa = 0..6: largest valid normal numbers (no infinity)
     *
     * FP8_E4M3_E_MAX_U / FP8_E4M3_E_MIN_U give the unbiased (actual) exponent range for normal numbers.
     */

#define FP8_E4M3_E_BIAS 7    /* Exponent bias: actual exponent = field - FP8_E4M3_E_BIAS */
#define FP8_E4M3_E_MAX  15   /* Maximum exponent field value (all 1s) */
#define FP8_E4M3_E_MAX_U 8   /* Maximum unbiased normal exponent (FP8_E4M3_E_MAX - FP8_E4M3_E_BIAS) */
#define FP8_E4M3_E_MIN_U -6  /* Minimum unbiased normal exponent (field=1 minus FP8_E4M3_E_BIAS) */

#define FP8_E4M3_FP8_E4M3_MANTISSA_MULTIPLIER  0.125  /* 2^(-FP8_E4M3_M_BITS): scales raw mantissa bits to fractional value */

    /*
     * Special values and boundary encodings
     *
     * Unlike IEEE 754, E4M3 has no dedicated infinity bit pattern.
     * Infinity maps to NaN instead (FP8_E4M3_INF_* == FP8_E4M3_NAN_*).
     *
     * Zeros:     Two representations (+0 and -0), like IEEE 754.
     * NaN:       Only when exponent=15 AND mantissa=7 (all 1s).
     *            Both +NaN (0x7F) and -NaN (0xFF) are defined.
     * Max norm:  exponent=15, mantissa=6  =>  +/- 2^8 * 1.75 = +/- 448
     * Min norm:  exponent=1,  mantissa=0  =>  +/- 2^(-6)
     * Max sub:   exponent=0,  mantissa=7  =>  +/- 0.875 * 2^(-6)
     * Min sub:   exponent=0,  mantissa=1  =>  +/- 2^(-9)
     */

#define FP8_E4M3_ZERO1    0x00  /* Positive zero: sign=0, exp=0, mant=0 */
#define FP8_E4M3_ZERO2    0x80  /* Negative zero: sign=1, exp=0, mant=0 */
#define FP8_E4M3_NAN_POS  0x7F  /* Positive NaN: sign=0, exp=15, mant=7 */
#define FP8_E4M3_INF_POS  FP8_E4M3_NAN_POS   /* No infinity in E4M3; maps to NaN */
#define FP8_E4M3_NAN_NEG  0xFF  /* Negative NaN: sign=1, exp=15, mant=7 */
#define FP8_E4M3_INF_NEG  FP8_E4M3_NAN_NEG   /* No infinity in E4M3; maps to NaN */
#define FP8_E4M3_MAX_NUM_POS 0x7E   /* Largest positive normal: 448 */
#define FP8_E4M3_MAX_NUM_NEG 0xFE   /* Largest negative normal: -448 */
#define FP8_E4M3_MIN_NUM_POS 0x08   /* Smallest positive normal: 2^-6 */
#define FP8_E4M3_MIN_NUM_NEG 0x88   /* Smallest negative normal: -2^-6 */
#define FP8_E4M3_MAX_SUB_NUM_POS 0x07   /* Largest positive subnormal: 0.875 * 2^-6 */
#define FP8_E4M3_MAX_SUB_NUM_NEG 0x87   /* Largest negative subnormal: -0.875 * 2^-6 */
#define FP8_E4M3_MIN_SUB_NUM_POS 0x01   /* Smallest positive subnormal: 2^-9 */
#define FP8_E4M3_MIN_SUB_NUM_NEG 0x81   /* Smallest negative subnormal: -2^-9 */
#define FP8_E4M3_EPSILON 0x20 /* Epsilon 0.125 = 2 ^ -3 */

    /* Field extraction helpers
     *
     * FP8_E4M3_GET_SIGN  -> 0 or 1 (7)
     * FP8_E4M3_GET_EXP   -> raw exponent field (3..6)
     * FP8_E4M3_GET_MANT  -> raw mantissa field (0..2)
     *
     * To reconstruct the sign-adjusted value:
     *   sign  = FP8_E4M3_GET_SIGN(x)  ?  -1.0 : 1.0;
     *   exp   = FP8_E4M3_GET_EXP(x);
     *   mant  = FP8_E4M3_GET_MANT(x);
     *   if (exp == 0)  -> subnormal:  sign * (mant * FP8_E4M3_FP8_E4M3_MANTISSA_MULTIPLIER) * 2^(1 - FP8_E4M3_E_BIAS)
     *   else           -> normal:     sign * (1.0 + mant * FP8_E4M3_FP8_E4M3_MANTISSA_MULTIPLIER) * 2^(exp - FP8_E4M3_E_BIAS)
     */

#define FP8_E4M3_GET_SIGN(x) (((x) & FP8_E4M3_S_MASK) >> (FP8_E4M3_S_SHIFT))
#define FP8_E4M3_GET_EXP(x)  (((x) & FP8_E4M3_E_MASK) >> (FP8_E4M3_E_SHIFT))
#define FP8_E4M3_GET_MANT(x)  ((x) & FP8_E4M3_M_MASK)



#ifdef __cplusplus
}
#endif

#endif /* LPF_FP8_E4M3_BITS_H */
