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

#include "fp8_e4m3_helper.h"
#include "emu/fp8_e4m3.h"
#include "lpf_internal.h"

/*
 * Function: acos ( C: acosf )
 */
HIDDEN void __fp8_e4m3_helper_acos(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_acos(_a.fp8_e4m3);
}

/*
 * Function: acosh ( C: acoshf )
 */
HIDDEN void __fp8_e4m3_helper_acosh(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_acosh(_a.fp8_e4m3);
}

/*
 * Function: asin ( C: asinf )
 */
HIDDEN void __fp8_e4m3_helper_asin(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_asin(_a.fp8_e4m3);
}

/*
 * Function: asinh ( C: asinhf )
 */
HIDDEN void __fp8_e4m3_helper_asinh(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_asinh(_a.fp8_e4m3);
}

/*
 * Function: atan ( C: atanf )
 */
HIDDEN void __fp8_e4m3_helper_atan(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_atan(_a.fp8_e4m3);
}

/*
 * Function: atanh ( C: atanhf )
 */
HIDDEN void __fp8_e4m3_helper_atanh(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_atanh(_a.fp8_e4m3);
}

/*
 * Function: bessel_j0 ( C: j0f )
 */
HIDDEN void __fp8_e4m3_helper_bessel_j0(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_bessel_j0(_a.fp8_e4m3);
}

/*
 * Function: bessel_j1 ( C: j1f )
 */
HIDDEN void __fp8_e4m3_helper_bessel_j1(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_bessel_j1(_a.fp8_e4m3);
}

/*
 * Function: bessel_y0 ( C: y0f )
 */
HIDDEN void __fp8_e4m3_helper_bessel_y0(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_bessel_y0(_a.fp8_e4m3);
}

/*
 * Function: bessel_y1 ( C: y1f )
 */
HIDDEN void __fp8_e4m3_helper_bessel_y1(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_bessel_y1(_a.fp8_e4m3);
}

/*
 * Function: ceiling ( C: ceilf )
 */
HIDDEN void __fp8_e4m3_helper_ceiling(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_ceiling(_a.fp8_e4m3);
}

/*
 * Function: cos ( C: cosf )
 */
HIDDEN void __fp8_e4m3_helper_cos(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_cos(_a.fp8_e4m3);
}

/*
 * Function: cosh ( C: coshf )
 */
HIDDEN void __fp8_e4m3_helper_cosh(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_cosh(_a.fp8_e4m3);
}

/*
 * Function: erf ( C: erff )
 */
HIDDEN void __fp8_e4m3_helper_erf(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_erf(_a.fp8_e4m3);
}

/*
 * Function: erfc ( C: erfcf )
 */
HIDDEN void __fp8_e4m3_helper_erfc(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_erfc(_a.fp8_e4m3);
}

/*
 * Function: exp ( C: expf )
 */
HIDDEN void __fp8_e4m3_helper_exp(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_exp(_a.fp8_e4m3);
}

/*
 * Function: floor ( C: floorf )
 */
HIDDEN void __fp8_e4m3_helper_floor(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_floor(_a.fp8_e4m3);
}

/*
 * Function: gamma ( C: tgammaf )
 */
HIDDEN void __fp8_e4m3_helper_gamma(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_gamma(_a.fp8_e4m3);
}

/*
 * Function: log ( C: logf )
 */
HIDDEN void __fp8_e4m3_helper_log(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_log(_a.fp8_e4m3);
}

/*
 * Function: log10 ( C: log10f )
 */
HIDDEN void __fp8_e4m3_helper_log10(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_log10(_a.fp8_e4m3);
}

/*
 * Function: log_gamma ( C: lgammaf )
 */
HIDDEN void __fp8_e4m3_helper_log_gamma(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_log_gamma(_a.fp8_e4m3);
}

/*
 * Function: sin ( C: sinf )
 */
HIDDEN void __fp8_e4m3_helper_sin(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_sin(_a.fp8_e4m3);
}

/*
 * Function: sinh ( C: sinhf )
 */
HIDDEN void __fp8_e4m3_helper_sinh(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_sinh(_a.fp8_e4m3);
}

/*
 * Function: sqrt ( C: sqrtf )
 */
HIDDEN void __fp8_e4m3_helper_sqrt(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_sqrt(_a.fp8_e4m3);
}

/*
 * Function: tan ( C: tanf )
 */
HIDDEN void __fp8_e4m3_helper_tan(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_tan(_a.fp8_e4m3);
}

/*
 * Function: tanh ( C: tanhf )
 */
HIDDEN void __fp8_e4m3_helper_tanh(int8_t *out, int8_t in)
{
    fp8_e4m3_handler_t *r = (fp8_e4m3_handler_t * ) out;
    fp8_e4m3_handler_t _a = { .i8 = in};
    r->fp8_e4m3 = fp8_e4m3_tanh(_a.fp8_e4m3);
}

