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


#include "lpf_fp8_e5m2_emu.h"
#include "fp8_e5m2_bits.h"
#include <stdio.h>

#if 0
/**
 * Format a 32-bit value as a space-separated binary string (groups of 4 bits, MSB first).
 *
 * Produces a string like "0000 0000 0000 0000 0000 0000 0000 0000".
 * Uses a thread-unsafe static buffer — not suitable for concurrent use.
 *
 * @param x  32-bit unsigned integer to format.
 * @return   Pointer to a static buffer containing the formatted string.
 */
static const char * str_bit_32(uint32_t x) {
    static char ret[50];
    int pos = 0;
    int c = 0;
    int i;
    for ( i = 31; i >= 0; i--) {
        ret[pos] = ((x>>i)&0x01) ? '1':'0';
        pos ++;
        c ++;
        if (c == 4 && i > 0) {
            ret[pos] = ' ';
            pos++;
            c = 0;
        }
    }
    ret[pos] = 0;
    return ret;
}

/**
 * Format a 16-bit value as a space-separated binary string (groups of 4 bits, MSB first).
 *
 * Produces a string like "0000 0000 0000 0000".
 * Uses a thread-unsafe static buffer — not suitable for concurrent use.
 *
 * @param x  16-bit unsigned integer to format.
 * @return   Pointer to a static buffer containing the formatted string.
 */
static const char * str_bit_16(uint16_t x) {
    static char ret[50];
    int pos = 0;
    int c = 0;
    int i;
    for ( i = 15; i >= 0; i--) {
        ret[pos] = ((x>>i)&0x01) ? '1':'0';
        pos ++;
        c ++;
        if (c == 4 && i > 0) {
            ret[pos] = ' ';
            pos++;
            c = 0;
        }
    }
    ret[pos] = 0;
    return ret;
}
#endif

/**
 * Format an 8-bit value as a space-separated binary string (groups of 4 bits, MSB first).
 *
 * Produces a string like "0000 0000", aligning with the FP8 E5M2 bit-field
 * layout: [S] [EEEEE] [MM].
 * Uses a thread-unsafe static buffer — not suitable for concurrent use.
 *
 * @param x  8-bit unsigned integer to format.
 * @return   Pointer to a static buffer containing the formatted string.
 */
static const char * str_bit_8(uint8_t x) {
    static char ret[50];
    int pos = 0;
    int c = 0;
    int i;
    for ( i = 7; i >= 0; i--) {
        ret[pos] = ((x>>i)&0x01) ? '1':'0';
        pos ++;
        c ++;
        if (c == 4 && i > 0) {
            ret[pos] = ' ';
            pos++;
            c = 0;
        }
    }
    ret[pos] = 0;
    return ret;
}


/**
 * Print a table of all 256 FP8 E5M2 encodings to stdout.
 *
 * Each line shows the hex code, the 8-bit binary representation, and the float value
 * obtained via fp8_e5m2_to_float().  Useful for verification and
 * visual inspection of the full encoding space.
 */
void fp8_e5m2_dump_all()
{
    uint16_t i;
    for (i = 0; i < 256; ++i) {
        printf("FP8_E5M2 (0x%02x, %8s) = %10.8f\n", i, str_bit_8(i), fp8_e5m2_to_float(i));
    }
}
