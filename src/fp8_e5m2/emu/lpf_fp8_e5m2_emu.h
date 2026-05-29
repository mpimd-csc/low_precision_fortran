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
#ifndef LPF_FP8_E5M2_EMU_H
#define LPF_FP8_E5M2_EMU_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <math.h>
#include "fp8_e5m2_bits.h"

    /*
     * FP8 E5M2 floating-point format (OCP standard)
     * =============================================
     *
     * Bit layout (8 bits total):
     *   [7]     Sign      (1 bit)
     *   [6:2]   Exponent  (5 bits, biased by 15)
     *   [1:0]   Mantissa  (2 bits, implicit leading 1 for normals)
     *
     * Normal number range:   [2^-14, 57344]  and  [-57344, -2^-14]
     * Subnormal number range: [2^-16, 0.75*2^-14]  and  [-0.75*2^-14, -2^-16]
     *
     * References:
     *  OCP 8-bit Floating Point Specification (OFP8)
     */
    typedef uint8_t fp8_e5m2_t;  /* Raw FP8 E5M2 value stored in a byte */

    /* Classification */
    int fp8_e5m2_isnan(fp8_e5m2_t fp8);
    int fp8_e5m2_isinf(fp8_e5m2_t fp8);
    int fp8_e5m2_iszero(fp8_e5m2_t fp8);

    /* Conversion */
    fp8_e5m2_t fp8_e5m2_from_components(uint8_t sign, int8_t expi, uint8_t mant);
    float fp8_e5m2_to_float(fp8_e5m2_t fp8);
    fp8_e5m2_t fp8_e5m2_from_float(float x);

    /* Arithmetic */
    fp8_e5m2_t fp8_e5m2_add(fp8_e5m2_t a, fp8_e5m2_t b);
    fp8_e5m2_t fp8_e5m2_sub(fp8_e5m2_t a, fp8_e5m2_t b);
    fp8_e5m2_t fp8_e5m2_mul(fp8_e5m2_t a, fp8_e5m2_t b);
    fp8_e5m2_t fp8_e5m2_div(fp8_e5m2_t a, fp8_e5m2_t b);
    fp8_e5m2_t fp8_e5m2_unitary_minus(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_pow(fp8_e5m2_t a, fp8_e5m2_t b);

    /* functions (non lookup) */
    fp8_e5m2_t fp8_e5m2_abs(fp8_e5m2_t x);


    /* Math functions (lookup-table based) */
    fp8_e5m2_t fp8_e5m2_acos(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_acosh(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_asin(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_asinh(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_atan(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_atanh(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_ceiling(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_cos(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_cosh(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_erf(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_erfc(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_exp(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_floor(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_gamma(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_log(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_log2(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_log10(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_log_gamma(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_sin(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_sinh(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_sqrt(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_tan(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_cotan(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_tanh(fp8_e5m2_t x);

    fp8_e5m2_t fp8_e5m2_bessel_j0(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_bessel_j1(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_bessel_y0(fp8_e5m2_t x);
    fp8_e5m2_t fp8_e5m2_bessel_y1(fp8_e5m2_t x);

    fp8_e5m2_t fp8_e5m2_atan2(fp8_e5m2_t x, fp8_e5m2_t y);
    fp8_e5m2_t fp8_e5m2_hypot(fp8_e5m2_t x, fp8_e5m2_t y);

    /* Debug */
    void fp8_e5m2_dump_all(void);

    /* Defines */
#define __FP8_E5M2_MIN_EXP__ (-14)
#define __FP8_E5M2_MAX_EXP__ (15)


#ifdef __cplusplus
}
#endif

#endif /* LPF_FP8_E5M2_H */
