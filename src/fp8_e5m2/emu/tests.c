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

/**
 * @file tests.c
 * @brief Correctness tests for FP8 E5M2 conversions and arithmetic.
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "lpf_fp8_e5m2_emu.h"
#include "fp8_e5m2_bits.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

/**
 * Test that converting an FP8 E5M2 value to float and back results in the same encoding.
 * This is a basic consistency check.
 */
void test_roundtrip_exhaustive() {
    int err = 0;
    printf("Running exhaustive roundtrip tests (FP8 -> Float -> FP8)...\n");

    for (uint16_t x = 0; x < 256; x++) {
        fp8_e5m2_t a = (fp8_e5m2_t)x;
        float f = fp8_e5m2_to_float(a);
        fp8_e5m2_t b = fp8_e5m2_from_float(f);

        if (fp8_e5m2_isnan(a)) {
            if (!fp8_e5m2_isnan(b)) {
                printf("NaN roundtrip failed: 0x%02x -> %f -> 0x%02x\n", a, f, b);
                err++;
            }
        } else if (a != b) {
            printf("Roundtrip failed: 0x%02x -> %20.16f -> 0x%02x\n", a, f, b);
            err++;
        }
        /** printf("0x%02x -> %20.16f -> 0x%02x %s\n", a, f, b, ( (fp8_e5m2_isnan(a) && fp8_e5m2_isnan(b)) || a == b) ? "OK" : "ERROR"); */
    }

    if (err) {
        printf("Roundtrip tests failed with %d errors.\n", err);
        abort();
    } else {
        printf("Roundtrip tests successful.\n");
    }
}

/**
 * Test specific boundary values and special encodings.
 */
void test_special_values() {
    printf("Testing special values...\n");

    // Zeros
    assert(fp8_e5m2_to_float(FP8_E5M2_ZERO1) == 0.0f);
    assert(fp8_e5m2_to_float(FP8_E5M2_ZERO2) == -0.0f);
    assert(fp8_e5m2_from_float(0.0f) == FP8_E5M2_ZERO1);
    assert(fp8_e5m2_from_float(-0.0f) == FP8_E5M2_ZERO2);

    // Infinities
    assert(fp8_e5m2_to_float(FP8_E5M2_INF_POS) == INFINITY);
    assert(fp8_e5m2_to_float(FP8_E5M2_INF_NEG) == -INFINITY);
    assert(fp8_e5m2_from_float(INFINITY) == FP8_E5M2_INF_POS);
    assert(fp8_e5m2_from_float(-INFINITY) == FP8_E5M2_INF_NEG);

    // NaNs
    assert(isnan(fp8_e5m2_to_float(FP8_E5M2_NAN_POS1)));
    assert(isnan(fp8_e5m2_to_float(FP8_E5M2_NAN_NEG1)));
    assert(isnan(fp8_e5m2_to_float(FP8_E5M2_NAN_POS3)));
    assert(fp8_e5m2_isnan(fp8_e5m2_from_float(NAN)));

    // Boundaries
    printf("Max Pos: %f (expected 57344.0)\n", fp8_e5m2_to_float(FP8_E5M2_MAX_NUM_POS));
    assert(fp8_e5m2_to_float(FP8_E5M2_MAX_NUM_POS) == 57344.0f);

    printf("Min Pos: %f (expected 2^-14)\n", fp8_e5m2_to_float(FP8_E5M2_MIN_NUM_POS));
    assert(fp8_e5m2_to_float(FP8_E5M2_MIN_NUM_POS) == powf(2.0f, -14.0f));

    printf("Min Sub Pos: %f (expected 2^-16)\n", fp8_e5m2_to_float(FP8_E5M2_MIN_SUB_NUM_POS));
    assert(fp8_e5m2_to_float(FP8_E5M2_MIN_SUB_NUM_POS) == powf(2.0f, -16.0f));

    printf("Special values tests successful.\n");
}

void test_add() {
    int err = 0;
    printf("Running exhaustive add tests...\n");
    for (uint16_t x = 0; x < 256; x++) {
        for (uint16_t y = 0; y < 256; y++) {
            fp8_e5m2_t a = (fp8_e5m2_t)x;
            fp8_e5m2_t b = (fp8_e5m2_t)y;
            fp8_e5m2_t c = fp8_e5m2_add(a, b);
            float aa = fp8_e5m2_to_float(a);
            float bb = fp8_e5m2_to_float(b);
            fp8_e5m2_t cc = fp8_e5m2_from_float(aa + bb);
            float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
            int is_ok = (d == 0.0f || isnanf(d));
            /** printf("%20.16f (0x%02x) + %20.16f (0x%02x) = %20.16f (0x%02x) \t ref: %20.16f (0x%02x) \t diff = %20.16f %s\n", */
                   /** aa, a, bb, b, fp8_e5m2_to_float(c), c, fp8_e5m2_to_float(cc), cc, d, is_ok ? "OK" : "ERROR"); */
            if (!is_ok) err++;
        }
    }
    if (err) { printf("Add tests failed.\n"); abort(); }
    printf("Add tests successful.\n");
}

void test_sub() {
    int err = 0;
    printf("Running exhaustive sub tests...\n");
    for (uint16_t x = 0; x < 256; x++) {
        for (uint16_t y = 0; y < 256; y++) {
            fp8_e5m2_t a = (fp8_e5m2_t)x;
            fp8_e5m2_t b = (fp8_e5m2_t)y;
            fp8_e5m2_t c = fp8_e5m2_sub(a, b);
            float aa = fp8_e5m2_to_float(a);
            float bb = fp8_e5m2_to_float(b);
            fp8_e5m2_t cc = fp8_e5m2_from_float(aa - bb);
            float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
            int is_ok = (d == 0.0f || isnanf(d));
            /** printf("%20.16f (0x%02x) - %20.16f (0x%02x) = %20.16f (0x%02x) \t ref: %20.16f (0x%02x) \t diff = %20.16f %s\n", */
                   /** aa, a, bb, b, fp8_e5m2_to_float(c), c, fp8_e5m2_to_float(cc), cc, d, is_ok ? "OK" : "ERROR"); */
            if (!is_ok) err++;
        }
    }
    if (err) { printf("Sub tests failed.\n"); abort(); }
    printf("Sub tests successful.\n");
}

void test_mul() {
    int err = 0;
    printf("Running exhaustive mul tests...\n");
    for (uint16_t x = 0; x < 256; x++) {
        for (uint16_t y = 0; y < 256; y++) {
            fp8_e5m2_t a = (fp8_e5m2_t)x;
            fp8_e5m2_t b = (fp8_e5m2_t)y;
            fp8_e5m2_t c = fp8_e5m2_mul(a, b);
            float aa = fp8_e5m2_to_float(a);
            float bb = fp8_e5m2_to_float(b);
            fp8_e5m2_t cc = fp8_e5m2_from_float(aa * bb);
            float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
            int is_ok = (d == 0.0f || isnanf(d));
            /** printf("%20.16f (0x%02x) * %20.16f (0x%02x) = %20.16f (0x%02x) \t ref: %20.16f (0x%02x) \t diff = %20.16f %s\n", */
                   /** aa, a, bb, b, fp8_e5m2_to_float(c), c, fp8_e5m2_to_float(cc), cc, d, is_ok ? "OK" : "ERROR"); */
            if (!is_ok) err++;
        }
    }
    if (err) { printf("Mul tests failed.\n"); abort(); }
    printf("Mul tests successful.\n");
}

void test_div() {
    int err = 0;
    printf("Running exhaustive div tests...\n");
    for (uint16_t x = 0; x < 256; x++) {
        for (uint16_t y = 0; y < 256; y++) {
            fp8_e5m2_t a = (fp8_e5m2_t)x;
            fp8_e5m2_t b = (fp8_e5m2_t)y;
            fp8_e5m2_t c = fp8_e5m2_div(a, b);
            float aa = fp8_e5m2_to_float(a);
            float bb = fp8_e5m2_to_float(b);
            fp8_e5m2_t cc = fp8_e5m2_from_float(aa / bb);
            float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
            int is_ok = (d == 0.0f || isnanf(d));
            /** printf("%20.16f (0x%02x) / %20.16f (0x%02x) = %20.16f (0x%02x) \t ref: %20.16f (0x%02x) \t diff = %20.16f %s\n", */
                   /** aa, a, bb, b, fp8_e5m2_to_float(c), c, fp8_e5m2_to_float(cc), cc, d, is_ok ? "OK" : "ERROR"); */
            if (!is_ok) err++;
        }
    }
    if (err) { printf("Div tests failed.\n"); abort(); }
    printf("Div tests successful.\n");
}


void test_abs() {
    int err = 0;
    printf("Running exhaustive abs tests...\n");
    for (uint16_t x = 0; x < 256; x++) {
        fp8_e5m2_t a = (fp8_e5m2_t)x;
        fp8_e5m2_t c = fp8_e5m2_abs(a);
        float aa = fp8_e5m2_to_float(a);
        fp8_e5m2_t cc = fp8_e5m2_from_float(fabsf(aa));
        float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
        int is_ok = (d == 0.0f || isnanf(d));
        /** printf("%20.16f (0x%02x) abs = %20.16f (0x%02x) \t ref: %20.16f (0x%02x) \t diff = %20.16f %s\n", */
               /** aa, a, fp8_e5m2_to_float(c), c, fp8_e5m2_to_float(cc), cc, d, is_ok ? "OK" : "ERROR"); */
        if (!is_ok) err++;
    }
    if (err) { printf("Abs tests failed.\n"); abort(); }
    printf("Abs tests successful.\n");
}


typedef fp8_e5m2_t (*math_func_t)(fp8_e5m2_t);

void run_unary_test(const char* name, math_func_t func, float (*ref_func)(float)) {
    int err = 0;
    printf("Testing %s...\n", name);
    for (uint16_t x = 0; x < 256; x++) {
        fp8_e5m2_t a = (fp8_e5m2_t)x;
        fp8_e5m2_t c = func(a);
        float aa = fp8_e5m2_to_float(a);
        float ref_val = ref_func(aa);
        fp8_e5m2_t cc = fp8_e5m2_from_float(ref_val);
        float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
        if (!(d == 0.0f || isnanf(d))) {
            printf("%s failed: 0x%02x -> %20.16f vs ref %20.16f (diff %20.16f)\n",
                   name, a, fp8_e5m2_to_float(c), fp8_e5m2_to_float(cc), d);
            err++;
        }
    }
    if (err) {
        printf("%s tests failed with %d errors.\n", name, err);
        abort();
    }
    printf("%s tests successful.\n", name);
}

float cotan_ref(float x) {
    return -tanf(M_PI/2.0f + x);
}

void test_math_lut() {
    printf("Running math LUT tests...\n");
    run_unary_test("acos", fp8_e5m2_acos, acosf);
    run_unary_test("acosh", fp8_e5m2_acosh, acoshf);
    run_unary_test("asin", fp8_e5m2_asin, asinf);
    run_unary_test("asinh", fp8_e5m2_asinh, asinhf);
    run_unary_test("atan", fp8_e5m2_atan, atanf);
    run_unary_test("atanh", fp8_e5m2_atanh, atanhf);
    run_unary_test("ceiling", fp8_e5m2_ceiling, ceilf);
    run_unary_test("cos", fp8_e5m2_cos, cosf);
    run_unary_test("cosh", fp8_e5m2_cosh, coshf);
    run_unary_test("erf", fp8_e5m2_erf, erff);
    run_unary_test("erfc", fp8_e5m2_erfc, erfcf);
    run_unary_test("exp", fp8_e5m2_exp, expf);
    run_unary_test("floor", fp8_e5m2_floor, floorf);
    run_unary_test("gamma", fp8_e5m2_gamma, tgammaf);
    run_unary_test("log", fp8_e5m2_log, logf);
    run_unary_test("log10", fp8_e5m2_log10, log10f);
    run_unary_test("log_gamma", fp8_e5m2_log_gamma, lgammaf);
    run_unary_test("sin", fp8_e5m2_sin, sinf);
    run_unary_test("sinh", fp8_e5m2_sinh, sinhf);
    run_unary_test("sqrt", fp8_e5m2_sqrt, sqrtf);
    run_unary_test("tan", fp8_e5m2_tan, tanf);
    run_unary_test("cotan", fp8_e5m2_cotan, cotan_ref);
    run_unary_test("tanh", fp8_e5m2_tanh, tanhf);
    run_unary_test("bessel_j0", fp8_e5m2_bessel_j0, j0f);
    run_unary_test("bessel_j1", fp8_e5m2_bessel_j1, j1f);
    run_unary_test("bessel_y0", fp8_e5m2_bessel_y0, y0f);
    run_unary_test("bessel_y1", fp8_e5m2_bessel_y1, y1f);

    // Binary functions
    printf("Testing atan2...\n");
    int err_atan2 = 0;
    for (uint16_t x = 0; x < 256; x++) {
        for (uint16_t y = 0; y < 256; y++) {
            fp8_e5m2_t a = (fp8_e5m2_t)x;
            fp8_e5m2_t b = (fp8_e5m2_t)y;
            fp8_e5m2_t c = fp8_e5m2_atan2(a, b);
            float aa = fp8_e5m2_to_float(a);
            float bb = fp8_e5m2_to_float(b);
            fp8_e5m2_t cc = fp8_e5m2_from_float(atan2f(aa, bb));
            float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
            if (!(d == 0.0f || isnanf(d))) {
                printf("atan2 failed: 0x%02x, 0x%02x -> %20.16f vs ref %20.16f (diff %20.16f)\n",
                       a, b, fp8_e5m2_to_float(c), fp8_e5m2_to_float(cc), d);
                err_atan2++;
            }
        }
    }
    if (err_atan2) { printf("atan2 tests failed with %d errors.\n", err_atan2); abort(); }
    printf("atan2 tests successful.\n");

    printf("Testing hypot...\n");
    int err_hypot = 0;
    for (uint16_t x = 0; x < 256; x++) {
        for (uint16_t y = 0; y < 256; y++) {
            fp8_e5m2_t a = (fp8_e5m2_t)x;
            fp8_e5m2_t b = (fp8_e5m2_t)y;
            fp8_e5m2_t c = fp8_e5m2_hypot(a, b);
            float aa = fp8_e5m2_to_float(a);
            float bb = fp8_e5m2_to_float(b);
            fp8_e5m2_t cc = fp8_e5m2_from_float(hypotf(aa, bb));
            float d = fabsf(fp8_e5m2_to_float(c) - fp8_e5m2_to_float(cc));
            if (!(d == 0.0f || isnanf(d))) {
                printf("hypot failed: 0x%02x, 0x%02x -> %20.16f vs ref %20.16f (diff %20.16f)\n",
                       a, b, fp8_e5m2_to_float(c), fp8_e5m2_to_float(cc), d);
                err_hypot++;
            }
        }
    }
    if (err_hypot) { printf("hypot tests failed with %d errors.\n", err_hypot); abort(); }
    printf("hypot tests successful.\n");

    printf("Math LUT tests successful.\n");
}

int main(int argc, char **argv) {
    (void)argc; (void)argv;

    test_special_values();
    test_abs();
    test_math_lut();
    test_roundtrip_exhaustive();
    test_add();
    test_sub();
    test_mul();
    test_div();

    printf("\n==== Full Encoding Dump ====\n");
    /** fp8_e5m2_dump_all(); */

    printf("\nAll conversion and arithmetic tests passed!\n");
    return 0;
}
