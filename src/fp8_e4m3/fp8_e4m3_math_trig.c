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

#include "fp8_e4m3/emu/fp8_e4m3.h"
#include "fp8_e4m3_helper.h"

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
HIDDEN void __fp8_e4m3_helper_acosd(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float x = acosf(y);
    r->fp8_e4m3 = fp8_e4m3_from_float(rad2deg(x));
}

/*
 * ASIND
 */
HIDDEN void __fp8_e4m3_helper_asind(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float x = asinf(y);
    r->fp8_e4m3 = fp8_e4m3_from_float(rad2deg(x));
}

/*
 * ATAND
 */
HIDDEN void __fp8_e4m3_helper_atand(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float x = atanf(y);
    r->fp8_e4m3 = fp8_e4m3_from_float(rad2deg(x));
}

/*
 * ATAN2D
 */
HIDDEN void __fp8_e4m3_helper_atan2d(int8_t *out, int8_t in1, int8_t in2)
{
    fp8_e4m3_handler_t a = {.i8 = in1};
    fp8_e4m3_handler_t b = {.i8 = in2};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y1 = fp8_e4m3_to_float(a.fp8_e4m3);
    float y2 = fp8_e4m3_to_float(b.fp8_e4m3);
    float x = atan2f(y1, y2);
    r->fp8_e4m3 = fp8_e4m3_from_float(rad2deg(x));
}

/*
 * COSD
 */
HIDDEN void __fp8_e4m3_helper_cosd(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float x = deg2rad(y);
    r->fp8_e4m3 = fp8_e4m3_cos(fp8_e4m3_from_float(x));
}

/*
 * SIND
 */
HIDDEN void __fp8_e4m3_helper_sind(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float x = deg2rad(y);
    r->fp8_e4m3 = fp8_e4m3_sin(fp8_e4m3_from_float(x));
}

/*
 * TAND
 */
HIDDEN void __fp8_e4m3_helper_tand(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float x = deg2rad(y);
    r->fp8_e4m3 = fp8_e4m3_tan(fp8_e4m3_from_float(x));
}

/*
 * COTAN(X) is translated into -TAN(X+PI/2)
 */
HIDDEN void __fp8_e4m3_helper_cotan(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    r->fp8_e4m3 = fp8_e4m3_cotan(a.fp8_e4m3);
}

/*
 * COTAND(X) is translated into -TAND(X+90) for REAL argument.
 */
HIDDEN void __fp8_e4m3_helper_cotand(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t a = {.i8 = in};
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t *) out;

    float y = fp8_e4m3_to_float(a.fp8_e4m3);
    float arg = deg2rad(y);
    r->fp8_e4m3 = fp8_e4m3_cotan(fp8_e4m3_from_float((arg)));
}





