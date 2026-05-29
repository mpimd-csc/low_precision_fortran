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
#include "lpf_fp8_e4m3_emu.h"
#include "fp8_e4m3_lut.h"

fp8_e4m3_t fp8_e4m3_acos(fp8_e4m3_t x)      { return fp8_e4m3_acos_table[x]; }
fp8_e4m3_t fp8_e4m3_acosh(fp8_e4m3_t x)     { return fp8_e4m3_acosh_table[x]; }
fp8_e4m3_t fp8_e4m3_asin(fp8_e4m3_t x)      { return fp8_e4m3_asin_table[x]; }
fp8_e4m3_t fp8_e4m3_asinh(fp8_e4m3_t x)     { return fp8_e4m3_asinh_table[x]; }
fp8_e4m3_t fp8_e4m3_atan(fp8_e4m3_t x)      { return fp8_e4m3_atan_table[x]; }
fp8_e4m3_t fp8_e4m3_atanh(fp8_e4m3_t x)     { return fp8_e4m3_atanh_table[x]; }
fp8_e4m3_t fp8_e4m3_ceiling(fp8_e4m3_t x)   { return fp8_e4m3_ceiling_table[x]; }
fp8_e4m3_t fp8_e4m3_cos(fp8_e4m3_t x)       { return fp8_e4m3_cos_table[x]; }
fp8_e4m3_t fp8_e4m3_cosh(fp8_e4m3_t x)      { return fp8_e4m3_cosh_table[x]; }
fp8_e4m3_t fp8_e4m3_erf(fp8_e4m3_t x)       { return fp8_e4m3_erf_table[x]; }
fp8_e4m3_t fp8_e4m3_erfc(fp8_e4m3_t x)      { return fp8_e4m3_erfc_table[x]; }
fp8_e4m3_t fp8_e4m3_exp(fp8_e4m3_t x)       { return fp8_e4m3_exp_table[x]; }
fp8_e4m3_t fp8_e4m3_floor(fp8_e4m3_t x)     { return fp8_e4m3_floor_table[x]; }
fp8_e4m3_t fp8_e4m3_gamma(fp8_e4m3_t x)     { return fp8_e4m3_gamma_table[x]; }
fp8_e4m3_t fp8_e4m3_log(fp8_e4m3_t x)       { return fp8_e4m3_log_table[x]; }
fp8_e4m3_t fp8_e4m3_log2(fp8_e4m3_t x)       { return fp8_e4m3_log2_table[x]; }
fp8_e4m3_t fp8_e4m3_log10(fp8_e4m3_t x)     { return fp8_e4m3_log10_table[x]; }
fp8_e4m3_t fp8_e4m3_log_gamma(fp8_e4m3_t x) { return fp8_e4m3_log_gamma_table[x]; }
fp8_e4m3_t fp8_e4m3_sin(fp8_e4m3_t x)       { return fp8_e4m3_sin_table[x]; }
fp8_e4m3_t fp8_e4m3_sinh(fp8_e4m3_t x)      { return fp8_e4m3_sinh_table[x]; }
fp8_e4m3_t fp8_e4m3_sqrt(fp8_e4m3_t x)      { return fp8_e4m3_sqrt_table[x]; }
fp8_e4m3_t fp8_e4m3_tan(fp8_e4m3_t x)       { return fp8_e4m3_tan_table[x]; }
fp8_e4m3_t fp8_e4m3_cotan(fp8_e4m3_t x)       { return fp8_e4m3_cotan_table[x]; }
fp8_e4m3_t fp8_e4m3_tanh(fp8_e4m3_t x)      { return fp8_e4m3_tanh_table[x]; }
fp8_e4m3_t fp8_e4m3_bessel_j0(fp8_e4m3_t x) { return fp8_e4m3_bessel_j0_table[x]; }
fp8_e4m3_t fp8_e4m3_bessel_j1(fp8_e4m3_t x) { return fp8_e4m3_bessel_j1_table[x]; }
fp8_e4m3_t fp8_e4m3_bessel_y0(fp8_e4m3_t x) { return fp8_e4m3_bessel_y0_table[x]; }
fp8_e4m3_t fp8_e4m3_bessel_y1(fp8_e4m3_t x) { return fp8_e4m3_bessel_y1_table[x]; }
fp8_e4m3_t fp8_e4m3_atan2(fp8_e4m3_t x, fp8_e4m3_t y)
{

    return fp8_e4m3_atan2_table[(size_t) x * 256 +(size_t)y];
}

fp8_e4m3_t fp8_e4m3_hypot(fp8_e4m3_t x, fp8_e4m3_t y)
{

    return fp8_e4m3_hypot_table[(size_t) x * 256 +(size_t)y];
}
