// SPDX-License-Identifier LGPL-3.0-or-later
/*
   This file is part of LPF, a Low Precision helper for Fortran
   Copyright (C) 2025 Martin Koehler

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
   */


#ifndef FP16_HELPER_H_P7LSSW6X
#define FP16_HELPER_H_P7LSSW6X

#include <stdint.h>
#include "lpf_internal.h"
#include <math.h>
#include "constants.h"

#define HIDDEN __attribute__((visibility("hidden")))

typedef union {
    _Float16 f16;
    int16_t  i16;
    uint16_t u16;
} fp16_handler_t;

/* Exports for C */

lpf_float16_t abs_fp16(lpf_float16_t x);
lpf_float16_t sqrt_fp16(lpf_float16_t x);
lpf_float16_t sign_fp16(lpf_float16_t a, lpf_float16_t b);
lpf_float16_t diff_fp16(lpf_float16_t a, lpf_float16_t b);

/* Operations with real64 */
void __fp16_helper_add_fp16_real64(int16_t *out, int16_t a, double b);
void __fp16_helper_add_real64_fp16(int16_t *out, double a, int16_t b);
void __fp16_helper_sub_fp16_real64(int16_t *out, int16_t a, double b);
void __fp16_helper_sub_real64_fp16(int16_t *out, double a, int16_t b);
void __fp16_helper_mul_fp16_real64(int16_t *out, int16_t a, double b);
void __fp16_helper_mul_real64_fp16(int16_t *out, double a, int16_t b);
void __fp16_helper_div_fp16_real64(int16_t *out, int16_t a, double b);
void __fp16_helper_div_real64_fp16(int16_t *out, double a, int16_t b);
#endif /* end of include guard: FP16_HELPER_H_P7LSSW6X */
