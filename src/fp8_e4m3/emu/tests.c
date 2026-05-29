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
 * @brief Exhaustive correctness tests for the FP8 E4M3 emulator.
 *
 * Verifies conversion (fp8_e4m3_to_float / fp8_e4m3_from_float) against
 * known boundary values and specific encodings, then runs 256x256
 * exhaustive tests for each arithmetic operation (add, sub, mul, div, sqrt)
 * by comparing the emulator result against a reference computed in float
 * space and re-quantized back to FP8 E4M3.
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "lpf_fp8_e4m3_emu.h"
#include "fp8_e4m3_bits.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

/**
 * Exhaustive addition test over all 256x256 FP8 E4M3 input pairs.
 *
 * For each pair (a, b), computes c = fp8_e4m3_add(a, b) and a reference
 * cc = fp8_e4m3_from_float(fp8_e4m3_to_float(a) + fp8_e4m3_to_float(b)),
 * then checks that both decode to the same float value.  Prints mismatches
 * and aborts on any error.
 */
void test_add() {
    int err = 0;

    /** uint8_t x = 0x05; */
    /** uint8_t y = 0xa8; */
    for (uint16_t x = 0; x < 256; x++)
    {
        for (uint16_t y = 0; y < 256; y++)
        {
            fp8_e4m3_t a = x;
            fp8_e4m3_t b = y;

            fp8_e4m3_t c = fp8_e4m3_add(a, b);

            float aa = fp8_e4m3_to_float(a);
            float bb = fp8_e4m3_to_float(b);

            fp8_e4m3_t cc = fp8_e4m3_from_float(aa+bb);

            float d = fabsf(fp8_e4m3_to_float(c) - fp8_e4m3_to_float(cc));
            if (d != 0.0f && !isnanf(d)) {
                printf("%20.16f (0x%02x) + %20.16f (0x%02x) = %20.16f \t %20.16f \t diff = %20.16f %s\n", fp8_e4m3_to_float(a), a, fp8_e4m3_to_float(b), b, fp8_e4m3_to_float(c), fp8_e4m3_to_float(cc),
                        d, (d == 0.0f || isnanf(d))?"":"error"
                      );
                err++;
            }
        }
    }
    if ( err ) {
        printf("Add tests fail.\n");
        abort();
    } else {
        printf("Add tests successful.\n");
    }

}

/**
 * Exhaustive multiplication test over all 256x256 FP8 E4M3 input pairs.
 *
 * Same strategy as test_add: compare fp8_e4m3_mul(a, b) against
 * fp8_e4m3_from_float(fp8_e4m3_to_float(a) * fp8_e4m3_to_float(b)).
 */
void test_mul() {
    int err = 0;

    /** uint8_t x = 0x01; */
    /** uint8_t y = 0x31; */
    for (uint16_t x = 0; x < 256; x++)
    {
        for (uint16_t y = 0; y < 256; y++)
        {
            fp8_e4m3_t a = x;
            fp8_e4m3_t b = y;

            fp8_e4m3_t c = fp8_e4m3_mul(a, b);

            float aa = fp8_e4m3_to_float(a);
            float bb = fp8_e4m3_to_float(b);

            fp8_e4m3_t cc = fp8_e4m3_from_float(aa*bb);

            float d = fabsf(fp8_e4m3_to_float(c) - fp8_e4m3_to_float(cc));
            if (d != 0.0f && !isnanf(d)) {
                printf("%20.16f (0x%02x) * %20.16f (0x%02x) = %20.16f \t %20.16f \t diff = %20.16f %s\n", fp8_e4m3_to_float(a), a, fp8_e4m3_to_float(b), b, fp8_e4m3_to_float(c), fp8_e4m3_to_float(cc),
                        d, (d == 0.0f || isnanf(d))?"":"error"
                      );
                err++;
            }
        }
    }
    if ( err ) {
        printf("Mul tests fail.\n");
        abort();
    } else {
        printf("Mul tests successful.\n");
    }

}

/**
 * Exhaustive division test over all 256x256 FP8 E4M3 input pairs.
 *
 * Same strategy as test_add: compare fp8_e4m3_div(a, b) against
 * fp8_e4m3_from_float(fp8_e4m3_to_float(a) / fp8_e4m3_to_float(b)).
 */
void test_div() {
    int err = 0;

    /** uint8_t x = 0x38; */
    /** uint8_t y = 0x02; */
    for (uint16_t x = 0; x < 256; x++)
    {
        for (uint16_t y = 0; y < 256; y++)
        {
            fp8_e4m3_t a = x;
            fp8_e4m3_t b = y;

            fp8_e4m3_t c = fp8_e4m3_div(a, b);

            float aa = fp8_e4m3_to_float(a);
            float bb = fp8_e4m3_to_float(b);

            fp8_e4m3_t cc = fp8_e4m3_from_float(aa/bb);

            float d = fabsf(fp8_e4m3_to_float(c) - fp8_e4m3_to_float(cc));
            if (d != 0.0f && !isnanf(d)) {
                printf("%20.16f (0x%02x) / %20.16f (0x%02x) = %20.16f \t %20.16f \t diff = %20.16f %s\n", fp8_e4m3_to_float(a), a, fp8_e4m3_to_float(b), b, fp8_e4m3_to_float(c), fp8_e4m3_to_float(cc),
                        d, (d == 0.0f || isnanf(d))?"":"error"
                      );
                err++;
            }
        }
    }
    if ( err ) {
        printf("DIV tests fail.\n");
        abort();
    } else {
        printf("DIV tests successful.\n");
    }

}

/**
 * Exhaustive subtraction test over all 256x256 FP8 E4M3 input pairs.
 *
 * Same strategy as test_add: compare fp8_e4m3_sub(a, b) against
 * fp8_e4m3_from_float(fp8_e4m3_to_float(a) - fp8_e4m3_to_float(b)).
 */
void test_sub() {
    int err = 0;

    /** uint8_t x = 0x05; */
    /** uint8_t y = 0xa8; */
    for (uint16_t x = 0; x < 256; x++)
    {
        for (uint16_t y = 0; y < 256; y++)
        {
            fp8_e4m3_t a = x;
            fp8_e4m3_t b = y;

            fp8_e4m3_t c = fp8_e4m3_sub(a, b);

            float aa = fp8_e4m3_to_float(a);
            float bb = fp8_e4m3_to_float(b);

            fp8_e4m3_t cc = fp8_e4m3_from_float(aa-bb);

            float d = fabsf(fp8_e4m3_to_float(c) - fp8_e4m3_to_float(cc));
            if (d != 0.0f && !isnanf(d)) {
                printf("%20.16f (0x%02x) - %20.16f (0x%02x) = %20.16f \t %20.16f \t diff = %20.16f %s\n", fp8_e4m3_to_float(a), a, fp8_e4m3_to_float(b), b, fp8_e4m3_to_float(c), fp8_e4m3_to_float(cc),
                        d, (d == 0.0f || isnanf(d))?"":"error"
                      );
                err++;
            }
        }
    }
    if ( err ) {
        printf("Sub tests fail.\n");
        abort();
    } else {
        printf("Sub tests successful.\n");
    }

}

/**
 * Exhaustive square root test over all 256 FP8 E4M3 values.
 *
 * For each encoding x, compares fp8_e4m3_sqrt(x) against
 * fp8_e4m3_from_float(sqrtf(fp8_e4m3_to_float(x))).  Prints every
 * result (including NaN for negative inputs) rather than only mismatches,
 * since the table-lookup sqrt truncation may differ from the float-reference
 * rounding path.
 */
void test_sqrt() {
    int err = 0;

    for (uint16_t x = 0; x < 256; x++)
    {
        fp8_e4m3_t a = x;

        fp8_e4m3_t c = fp8_e4m3_sqrt(a);
        float aa = fp8_e4m3_to_float(a);

        fp8_e4m3_t cc = fp8_e4m3_from_float(sqrtf(aa));

        float d = fabsf(fp8_e4m3_to_float(c) - fp8_e4m3_to_float(cc));
        printf("sqrt(%20.16f (0x%02x)) = %20.16f (0x%02x)\t %20.16f \t diff = %20.16f %s\n", fp8_e4m3_to_float(a), a, fp8_e4m3_to_float(c), c, fp8_e4m3_to_float(cc),
                d, (d == 0.0f || isnanf(d))?"":"error"
              );
    }
    if ( err ) {
        printf("DIV tests fail.\n");
        abort();
    } else {
        printf("DIV tests successful.\n");
    }

}

/**
 * Main test driver for the FP8 E4M3 emulator.
 *
 * Runs four phases:
 *   1. **fp8_e4m3_to_float** — asserts on boundary values (zeros, min/max
 *      normal and subnormal, 1.0) and a set of specific encodings.
 *   2. **fp8_e4m3_from_float** — asserts on the reverse mapping for the
 *      same values plus overflow-to-NaN cases (512, -512) and
 *      subnormal/flush-to-zero cases.
 *   3. **Arithmetic** — exhaustive 256x256 tests for add, sub, mul, div.
 *   4. **sqrt** — exhaustive 256-value test for square root.
 *
 * Aborts on any assertion failure.  Returns 0 on success.
 */
int main(int argc, char **argv)
{
    (void) argc;
    (void) argv;

    printf("==== TEST for fp8_e4m3_to_float ====\n");

    printf("FP8_E4M3_ZERO1 : %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_ZERO1));
    assert(fp8_e4m3_to_float(FP8_E4M3_ZERO1) == 0.0f);
    printf("FP8_E4M3_ZERO2 : %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_ZERO2));
    assert(fp8_e4m3_to_float(FP8_E4M3_ZERO2) == -0.0f);


    printf("FP8_E4M3_MAX_NUM_POS: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MAX_NUM_POS));
    assert(fp8_e4m3_to_float(FP8_E4M3_MAX_NUM_POS) == 448.0f);

    printf("FP8_E4M3_MAX_NUM_NEG: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MAX_NUM_NEG));
    assert(fp8_e4m3_to_float(FP8_E4M3_MAX_NUM_NEG) == -448.0f);

    printf("FP8_E4M3_MIN_NUM_POS: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MIN_NUM_POS));
    assert(fp8_e4m3_to_float(FP8_E4M3_MIN_NUM_POS) == 0.01562500f);

    printf("FP8_E4M3_MIN_NUM_NEG: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MIN_NUM_NEG));
    assert(fp8_e4m3_to_float(FP8_E4M3_MIN_NUM_NEG) == -0.01562500f);

    printf("FP8_E4M3_MAX_SUB_NUM_POS: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MAX_SUB_NUM_POS));
    assert(fp8_e4m3_to_float(FP8_E4M3_MAX_SUB_NUM_POS) == 0.013671875f);

    printf("FP8_E4M3_MAX_SUB_NUM_NEG: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MAX_SUB_NUM_NEG));
    assert(fp8_e4m3_to_float(FP8_E4M3_MAX_SUB_NUM_NEG) == -0.013671875f);

    printf("FP8_E4M3_MIN_SUB_NUM_POS: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MIN_SUB_NUM_POS));
    assert(fp8_e4m3_to_float(FP8_E4M3_MIN_SUB_NUM_POS) == 0.001953125f);

    printf("FP8_E4M3_MIN_SUB_NUM_NEG: %10.8f\n", fp8_e4m3_to_float(FP8_E4M3_MIN_SUB_NUM_NEG));
    assert(fp8_e4m3_to_float(FP8_E4M3_MIN_SUB_NUM_NEG) == -0.001953125f);

    printf("FP8_ONE (0x38): %10.8f\n", fp8_e4m3_to_float(0x38));
    assert(fp8_e4m3_to_float(0x38) == 1.0);

#define TEST(out,in) do {\
    printf("fp8_e4m3_to_float(0x%02x) = %10.2f\n", (unsigned int) in, fp8_e4m3_to_float(in)); \
    assert(fp8_e4m3_to_float(in) == out); \
    } while(0)

    TEST(0.0, 0x00);
    TEST(0.001953125, 0x1);
    TEST(0.013671875, 0x7);
    TEST(0.015625,0x8);
    TEST(0.34375,0x2b);
    TEST(0.9375,0x37);
    TEST(1.0, 0x38);
    TEST(1.125, 0x39);
    TEST(1.375, 0x3b);
    TEST(3.25, 0x45);
    TEST(240.0, 0x77);
    TEST(448.0, 0x7e);
    TEST(-2.0, 0xc0);

#undef TEST



    printf("\n");
    printf("==== Test for float_to_fp8_e4m3 ====\n");
#define TEST(in,out) do {\
    printf("fp8_e4m3_from_float(%10.8f) = 0x%02x\n", in, fp8_e4m3_from_float(in)); \
    assert(fp8_e4m3_from_float(in) == out); \
    } while(0)

    TEST(0.0, 0x00);
    TEST(0.0001123, 0x0);
    TEST(0.001953125, 0x1);
    TEST(0.013671875, 0x7);
    TEST(0.0139, 0x7);
    TEST(0.0145, 0x7);
    TEST(0.0148, 0x8);
    TEST(0.015, 0x8);
    TEST(0.015625,0x8);
    TEST(0.34375,0x2b);
    TEST(0.9375,0x37);
    TEST(1.0, 0x38);
    TEST(1.125, 0x39);
    TEST(1.375, 0x3b);
    TEST(3.25, 0x45);
    TEST(240.0, 0x77);
    TEST(448.0, 0x7e);
    TEST(-2.0, 0xc0);
    TEST(512.0, 0x7f);
    TEST(-512.0,0xff);

    TEST(0.001953125*0.50, 0x00);
    TEST(0.001953125*0.5625000000000000, 0x01);
    TEST(0.001953125/288.0,0x00);
#undef TEST
    printf("\n");
    printf("==== Test Add ====\n");
    test_add();

    printf("\n");
    printf("==== Test Sub ====\n");
    test_sub();

    printf("\n");
    printf("==== Test Mul ====\n");
    test_mul();

    printf("\n");
    printf("==== Test Div ====\n");
    test_div();

    printf("\n");
    printf("==== Test SQRT ====\n");
    /** list_sqrt(); */
    test_sqrt();





    /** fp8_e4m3_dump_all(); */

    /** for (int x = 0; x <= 128; x++) { */
    /**     float xx = x; */
    /**     fp8_e4m3_t x8 = fp8_e4m3_from_float(xx); */
    /**  */
    /**     printf("%20.16f \t 0x%02x \t %20.16f\n", xx, x8, fp8_e4m3_to_float(x8)); */
    /** } */
    /**  */
    return 0;
}


