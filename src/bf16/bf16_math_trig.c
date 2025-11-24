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
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "bf16_helper.h"

static float rad2deg(float x)
{
    return (x * 180.0f)/M_PI;
}

static float deg2rad(float x)
{
    return (x * M_PI)/180.0f;
}


/*
 * ACOSD
 */
HIDDEN void __bf16_helper_acosd(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = (float) acosf(a.bf16);
    r->bf16 = rad2deg(x);
}

/*
 * ASIND
 */
HIDDEN void __bf16_helper_asind(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = (float) asinf(a.bf16);
    r->bf16 = rad2deg(x);
}

/*
 * ATAND
 */
HIDDEN void __bf16_helper_atand(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = (float) atanf(a.bf16);
    r->bf16 = rad2deg(x);
}

/*
 * ATAN2D
 */
HIDDEN void __bf16_helper_atan2d(int16_t *out, int16_t in1, int16_t in2)
{
    bf16_handler_t a = {.i16 = in1};
    bf16_handler_t b = {.i16 = in2};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = (float) atan2f(a.bf16, b.bf16);
    r->bf16 = rad2deg(x);
}

/*
 * COSD
 */
HIDDEN void __bf16_helper_cosd(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = deg2rad(a.bf16);
    r->bf16 = (__bf16) cosf(x);
}

/*
 * SIND
 */
HIDDEN void __bf16_helper_sind(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = deg2rad(a.bf16);
    r->bf16 = (__bf16) sinf(x);
}

/*
 * TAND
 */
HIDDEN void __bf16_helper_tand(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float x = deg2rad(a.bf16);
    r->bf16 = (__bf16) tanf(x);
}

/*
 * COTAN(X) is translated into -TAN(X+PI/2)
 */
HIDDEN void __bf16_helper_cotan(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float arg = M_PI/2.0f + (float) a.bf16;
    r->bf16 = (__bf16) -tanf(arg);
}

/*
 * COTAND(X) is translated into -TAND(X+90) for REAL argument.
 */
HIDDEN void __bf16_helper_cotand(int16_t *out, int16_t in)
{
    bf16_handler_t a = {.i16 = in};
    bf16_handler_t *r = (bf16_handler_t *) out;

    float arg = (float) a.bf16 + 90.0f;
    arg = deg2rad(arg);
    r->bf16 = (__bf16) -tanf(arg);
}





