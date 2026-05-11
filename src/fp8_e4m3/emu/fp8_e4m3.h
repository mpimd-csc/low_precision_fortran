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
#ifndef LPF_FP8_E4M3_H
#define LPF_FP8_E4M3_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <math.h>
#include "fp8_e4m3_bits.h"

    /*
     * FP8 E4M3 floating-point format (OCP standard)
     * =============================================
     *
     * Bit layout (8 bits total):
     *   [7]     Sign      (1 bit)
     *   [6:3]   Exponent  (4 bits, biased by 7)
     *   [2:0]   Mantissa  (3 bits, implicit leading 1 for normals)
     *
     * Key difference from IEEE 754:
     *   There is NO infinity representation. When the exponent field is all 1s
     *   and the mantissa field is all 1s, the encoding is NaN (not infinity).
     *   All other exponent-all-1s patterns are valid normal numbers.
     *   This extends the dynamic range compared to IEEE-style FP8 (E5M2).
     *
     * Normal number range:   [2^-6, 448]  and  [-448, -2^-6]
     * Subnormal number range: [2^-9, 0.875*2^-6]  and  [-0.875*2^-6, -2^-9]
     *
     * References:
     *  OCP 8-bit Floating Point Specification (OFP8)
     */
    typedef uint8_t fp8_e4m3_t;  /* Raw FP8 E4M3 value stored in a byte */

    /* Classification */
    int fp8_e4m3_isnan(fp8_e4m3_t fp8);
    int fp8_e4m3_isinf(fp8_e4m3_t fp8);
    int fp8_e4m3_iszero(fp8_e4m3_t fp8);

    /* Conversion */
    fp8_e4m3_t fp8_e4m3_from_components(uint8_t sign, int8_t expi, uint8_t mant);
    float fp8_e4m3_to_float(fp8_e4m3_t fp8);
    fp8_e4m3_t fp8_e4m3_from_float(float x);

    /* Arithmetic */
    fp8_e4m3_t fp8_e4m3_add(fp8_e4m3_t a, fp8_e4m3_t b);
    fp8_e4m3_t fp8_e4m3_sub(fp8_e4m3_t a, fp8_e4m3_t b);
    fp8_e4m3_t fp8_e4m3_mul(fp8_e4m3_t a, fp8_e4m3_t b);
    fp8_e4m3_t fp8_e4m3_div(fp8_e4m3_t a, fp8_e4m3_t b);
    fp8_e4m3_t fp8_e4m3_unitary_minus(fp8_e4m3_t x);

    /* functions (non lookup) */
    fp8_e4m3_t fp8_e4m3_abs(fp8_e4m3_t x);


    /* Math functions (lookup-table based) */
    fp8_e4m3_t fp8_e4m3_acos(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_acosh(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_asin(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_asinh(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_atan(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_atanh(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_ceiling(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_cos(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_cosh(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_erf(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_erfc(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_exp(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_floor(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_gamma(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_log(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_log10(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_log_gamma(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_sin(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_sinh(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_sqrt(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_tan(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_cotan(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_tanh(fp8_e4m3_t x);

    fp8_e4m3_t fp8_e4m3_bessel_j0(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_bessel_j1(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_bessel_y0(fp8_e4m3_t x);
    fp8_e4m3_t fp8_e4m3_bessel_y1(fp8_e4m3_t x);

    fp8_e4m3_t fp8_e4m3_atan2(fp8_e4m3_t x, fp8_e4m3_t y);
    fp8_e4m3_t fp8_e4m3_hypot(fp8_e4m3_t x, fp8_e4m3_t y);

    /* Debug */
    void fp8_e4m3_dump_all(void);


#ifdef __cplusplus
}
#endif

#endif /* LPF_FP8_E4M3_H */
