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

#include <stddef.h>
#include "lpf_fp8_e5m2_emu.h"
#include "fp8_e5m2_lut.h"

fp8_e5m2_t fp8_e5m2_acos(fp8_e5m2_t x)      { return fp8_e5m2_acos_table[x]; }
fp8_e5m2_t fp8_e5m2_acosh(fp8_e5m2_t x)     { return fp8_e5m2_acosh_table[x]; }
fp8_e5m2_t fp8_e5m2_asin(fp8_e5m2_t x)      { return fp8_e5m2_asin_table[x]; }
fp8_e5m2_t fp8_e5m2_asinh(fp8_e5m2_t x)     { return fp8_e5m2_asinh_table[x]; }
fp8_e5m2_t fp8_e5m2_atan(fp8_e5m2_t x)      { return fp8_e5m2_atan_table[x]; }
fp8_e5m2_t fp8_e5m2_atanh(fp8_e5m2_t x)     { return fp8_e5m2_atanh_table[x]; }
fp8_e5m2_t fp8_e5m2_ceiling(fp8_e5m2_t x)   { return fp8_e5m2_ceiling_table[x]; }
fp8_e5m2_t fp8_e5m2_cos(fp8_e5m2_t x)       { return fp8_e5m2_cos_table[x]; }
fp8_e5m2_t fp8_e5m2_cosh(fp8_e5m2_t x)      { return fp8_e5m2_cosh_table[x]; }
fp8_e5m2_t fp8_e5m2_erf(fp8_e5m2_t x)       { return fp8_e5m2_erf_table[x]; }
fp8_e5m2_t fp8_e5m2_erfc(fp8_e5m2_t x)      { return fp8_e5m2_erfc_table[x]; }
fp8_e5m2_t fp8_e5m2_exp(fp8_e5m2_t x)       { return fp8_e5m2_exp_table[x]; }
fp8_e5m2_t fp8_e5m2_floor(fp8_e5m2_t x)     { return fp8_e5m2_floor_table[x]; }
fp8_e5m2_t fp8_e5m2_gamma(fp8_e5m2_t x)     { return fp8_e5m2_gamma_table[x]; }
fp8_e5m2_t fp8_e5m2_log(fp8_e5m2_t x)       { return fp8_e5m2_log_table[x]; }
fp8_e5m2_t fp8_e5m2_log2(fp8_e5m2_t x)       { return fp8_e5m2_log2_table[x]; }
fp8_e5m2_t fp8_e5m2_log10(fp8_e5m2_t x)     { return fp8_e5m2_log10_table[x]; }
fp8_e5m2_t fp8_e5m2_log_gamma(fp8_e5m2_t x) { return fp8_e5m2_log_gamma_table[x]; }
fp8_e5m2_t fp8_e5m2_sin(fp8_e5m2_t x)       { return fp8_e5m2_sin_table[x]; }
fp8_e5m2_t fp8_e5m2_sinh(fp8_e5m2_t x)      { return fp8_e5m2_sinh_table[x]; }
fp8_e5m2_t fp8_e5m2_sqrt(fp8_e5m2_t x)      { return fp8_e5m2_sqrt_table[x]; }
fp8_e5m2_t fp8_e5m2_tan(fp8_e5m2_t x)       { return fp8_e5m2_tan_table[x]; }
fp8_e5m2_t fp8_e5m2_cotan(fp8_e5m2_t x)     { return fp8_e5m2_cotan_table[x]; }
fp8_e5m2_t fp8_e5m2_tanh(fp8_e5m2_t x)      { return fp8_e5m2_tanh_table[x]; }
fp8_e5m2_t fp8_e5m2_bessel_j0(fp8_e5m2_t x) { return fp8_e5m2_bessel_j0_table[x]; }
fp8_e5m2_t fp8_e5m2_bessel_j1(fp8_e5m2_t x) { return fp8_e5m2_bessel_j1_table[x]; }
fp8_e5m2_t fp8_e5m2_bessel_y0(fp8_e5m2_t x) { return fp8_e5m2_bessel_y0_table[x]; }
fp8_e5m2_t fp8_e5m2_bessel_y1(fp8_e5m2_t x) { return fp8_e5m2_bessel_y1_table[x]; }
fp8_e5m2_t fp8_e5m2_atan2(fp8_e5m2_t x, fp8_e5m2_t y)
{
    return fp8_e5m2_atan2_table[(size_t) x * 256 +(size_t)y];
}

fp8_e5m2_t fp8_e5m2_hypot(fp8_e5m2_t x, fp8_e5m2_t y)
{
    return fp8_e5m2_hypot_table[(size_t) x * 256 +(size_t)y];
}
