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

/*
 * gen_lookup_tables.c - Lookup table generator for fp8_e4m3 math functions.
 *
 * Generates a C header file with precomputed 256-entry lookup tables for
 * unary math functions (sin, cos, sqrt, etc.). Each table entry is computed
 * by the float round-trip:
 *   fp8_e4m3_from_float(func(fp8_e4m3_to_float(i)))
 *
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "fp8_e4m3.h"
#include "fp8_e4m3_bits.h"

typedef float (*math_func_t)(float);
typedef float (*math2_func_t)(float, float);

typedef struct {
    const char *name;
    const char *cname;
    math_func_t func;
} func_entry_t;

typedef struct {
    const char *name;
    const char *cname;
    math2_func_t func;
} func2_entry_t;

float cotanf(float x)
{
    float arg = M_PI/2.0f + x;
    return -tanf(arg);

}

static func_entry_t supported_functions[] = {
    {"acos",       "acosf",    acosf},
    {"acosh",      "acoshf",   acoshf},
    {"asin",       "asinf",    asinf},
    {"asinh",      "asinhf",   asinhf},
    {"atan",       "atanf",    atanf},
    {"atanh",      "atanhf",   atanhf},
    {"ceiling",    "ceilf",    ceilf},
    {"cos",        "cosf",     cosf},
    {"cosh",       "coshf",    coshf},
    {"erf",        "erff",     erff},
    {"erfc",       "erfcf",    erfcf},
    {"exp",        "expf",     expf},
    {"floor",      "floorf",   floorf},
    {"gamma",      "tgammaf",  tgammaf},
    {"log",        "logf",     logf},
    {"log10",      "log10f",   log10f},
    {"log_gamma",  "lgammaf",  lgammaf},
    {"sin",        "sinf",     sinf},
    {"sinh",       "sinhf",    sinhf},
    {"sqrt",       "sqrtf",    sqrtf},
    {"tan",        "tanf",     tanf},
    {"cotan",        "cotanf",     cotanf},
    {"tanh",       "tanhf",    tanhf},
    {"bessel_j0",  "j0f",      j0f},
    {"bessel_j1",  "j1f",      j1f},
    {"bessel_y0",  "y0f",      y0f},
    {"bessel_y1",  "y1f",      y1f},
};

static func2_entry_t supported_functions2[] = {
    {"atan2",       "atan2f",    atan2f},
    {"hypot",      "hypof",    hypotf}
};




#define NUM_FUNCTIONS ((int)(sizeof(supported_functions) / sizeof(supported_functions[0])))
#define NUM_FUNCTIONS2 ((int)(sizeof(supported_functions2) / sizeof(supported_functions2[0])))

static void print_list(void)
{
    for (int i = 0; i < NUM_FUNCTIONS; i++) {
        printf("%-12s (%s)\n", supported_functions[i].name, supported_functions[i].cname);
    }
}

static void generate_table(FILE *out, const func_entry_t *entry)
{
    fprintf(out, "/* Lookup table for fp8_e4m3 %s: table[i] = "
            "fp8_e4m3_from_float(%s(fp8_e4m3_to_float(i))) */\n",
            entry->name, entry->cname);
    fprintf(out, "static const uint8_t fp8_e4m3_%s_table[256] = {\n",
            entry->name);

    for (int i = 0; i < 256; i++) {
        fp8_e4m3_t val = (uint8_t)i;
        uint8_t result;

        if (fp8_e4m3_isnan(val)) {
            result = FP8_E4M3_GET_SIGN(val) ? FP8_E4M3_NAN_NEG : FP8_E4M3_NAN_POS;
        } else {
            float f = fp8_e4m3_to_float(val);
            float r = entry->func(f);
            result = fp8_e4m3_from_float(r);
        }

        fprintf(out, "    0x%02X%s", result, (i < 255) ? "," : "");
        if ((i + 1) % 8 == 0) fprintf(out, "\n");
    }
    fprintf(out, "};\n\n");
}

static void generate_table2(FILE *out, const func2_entry_t *entry)
{
    fprintf(out, "/* Lookup table for fp8_e4m3 %s: table[i] = "
            "fp8_e4m3_from_float(%s(fp8_e4m3_to_float(i), fp8_e4m3_to_float(i))) */\n",
            entry->name, entry->cname);
    fprintf(out, "static const uint8_t fp8_e4m3_%s_table[256*256] = {\n",
            entry->name);

    int cnt = 0;
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            fp8_e4m3_t val1 = (uint8_t)i;
            fp8_e4m3_t val2 = (uint8_t)j;
            uint8_t result;

            if (fp8_e4m3_isnan(val1)) {
                result = FP8_E4M3_GET_SIGN(val1) ? FP8_E4M3_NAN_NEG : FP8_E4M3_NAN_POS;
            } else if (fp8_e4m3_isnan(val2)) {
                result = FP8_E4M3_GET_SIGN(val2) ? FP8_E4M3_NAN_NEG : FP8_E4M3_NAN_POS;
            }  else {
                float f1 = fp8_e4m3_to_float(val1);
                float f2 = fp8_e4m3_to_float(val2);
                float r = entry->func(f1, f2);
                result = fp8_e4m3_from_float(r);
            }

            fprintf(out, "    0x%02X%s", result, ( cnt+1 < 256*256 ) ? "," : "");
            if ((cnt + 1) % 8 == 0) fprintf(out, "\n");
            cnt++;
        }
    }
    fprintf(out, "};\n\n");
}


static void write_header(FILE *out)
{
    fprintf(out,
        "/* AUTO-GENERATED by gen_lookup_tables -- DO NOT EDIT */\n"
        "/* SPDX-License-Identifier LGPL-3.0-or-later */\n"
        "/*\n"
        " * Copyright (C)  2025 by Martin Koehler\n"
        " *\n"
        " * This program is free software; you can redistribute it and/or\n"
        " * modify it under the terms of the GNU Lesser General Public\n"
        " * License as published by the Free Software Foundation; either\n"
        " * version 3 of the License, or (at your option) any later version.\n"
        " *\n"
        " * This program is distributed in the hope that it will be useful,\n"
        " * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        " * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"
        " * Lesser General Public License for more details.\n"
        " */\n"
        "#ifndef LPF_FP8_E4M3_LUT_H\n"
        "#define LPF_FP8_E4M3_LUT_H\n\n"
        "#include <stdint.h>\n\n"
        "#ifdef __cplusplus\n"
        "extern \"C\" {\n"
        "#endif\n\n"
    );

    for (int i = 0; i < NUM_FUNCTIONS; i++) {
        generate_table(out, &supported_functions[i]);
    }
    for (int i = 0; i < NUM_FUNCTIONS2; i++) {
        generate_table2(out, &supported_functions2[i]);
    }


    fprintf(out,
        "#ifdef __cplusplus\n"
        "}\n"
        "#endif\n\n"
        "#endif /* LPF_FP8_E4M3_LUT_H */\n"
    );
}

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;

    write_header(stdout);
    return 0;
}
