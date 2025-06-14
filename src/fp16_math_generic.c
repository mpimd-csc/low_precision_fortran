//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of FP16_SUPPORTS, an FP16 Helper for Fortran
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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "fp16_helper.h"

/*
 * ABS
 */
HIDDEN void __fp16_helper_abs(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    if ( _a.f16 < 0 ) {
        r->f16 = - _a.f16;
    } else {
        r->f16 =  _a.f16;
    }
}

/*
 * EPSILON
 */
HIDDEN void __fp16_helper_epsilon(int16_t *out)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    r->f16 = __FLT16_EPSILON__;
}

/*
 * Exponent
 */
HIDDEN void __fp16_helper_exponent(int *out, int16_t in)
{
    /** uint16_t number = * (uint16_t*) &in; */
    /** uint16_t exponent = (number >> 10) & 0x1F; */
    /** uint16_t fraction = number & 0x3FF; */
    /** int ex = exponent - 15; */
    /** if ( fraction == 0 && exponent == 0 ) { */
    /**     *out = 0; */
    /**     return; */
    /** } */
    /** *out = (int) ex; */
    fp16_handler_t r = {.i16 = in};
    float x = (float) r.f16;
    frexpf(x, out);
}

/*
 * Fraction
 */
HIDDEN void __fp16_helper_fraction(int16_t *out, int16_t in)
{
    fp16_handler_t *rout = (fp16_handler_t * ) out;
    fp16_handler_t r = {.i16 = in};
    float x = (float) r.f16;
    int dummy = 0;
    float frac = frexpf(x, &dummy);
    rout->f16 = (_Float16) frac;
}

