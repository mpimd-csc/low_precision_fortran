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

#include "fp16_helper.h"

/*
 * Function: acos ( C: acosf )
 */
HIDDEN void __fp16_helper_acos(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) acosf ((float) _a.f16);
}

/*
 * Function: acosh ( C: acoshf )
 */
HIDDEN void __fp16_helper_acosh(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) acoshf ((float) _a.f16);
}

/*
 * Function: asin ( C: asinf )
 */
HIDDEN void __fp16_helper_asin(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) asinf ((float) _a.f16);
}

/*
 * Function: asinh ( C: asinhf )
 */
HIDDEN void __fp16_helper_asinh(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) asinhf ((float) _a.f16);
}

/*
 * Function: atan ( C: atanf )
 */
HIDDEN void __fp16_helper_atan(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) atanf ((float) _a.f16);
}

/*
 * Function: atanh ( C: atanhf )
 */
HIDDEN void __fp16_helper_atanh(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) atanhf ((float) _a.f16);
}

/*
 * Function: bessel_j0 ( C: j0f )
 */
HIDDEN void __fp16_helper_bessel_j0(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) j0f ((float) _a.f16);
}

/*
 * Function: bessel_j1 ( C: j1f )
 */
HIDDEN void __fp16_helper_bessel_j1(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) j1f ((float) _a.f16);
}

/*
 * Function: bessel_y0 ( C: y0f )
 */
HIDDEN void __fp16_helper_bessel_y0(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) y0f ((float) _a.f16);
}

/*
 * Function: bessel_y1 ( C: y1f )
 */
HIDDEN void __fp16_helper_bessel_y1(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) y1f ((float) _a.f16);
}

/*
 * Function: ceiling ( C: ceilf )
 */
HIDDEN void __fp16_helper_ceiling(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) ceilf ((float) _a.f16);
}

/*
 * Function: cos ( C: cosf )
 */
HIDDEN void __fp16_helper_cos(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) cosf ((float) _a.f16);
}

/*
 * Function: cosh ( C: coshf )
 */
HIDDEN void __fp16_helper_cosh(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) coshf ((float) _a.f16);
}

/*
 * Function: erf ( C: erff )
 */
HIDDEN void __fp16_helper_erf(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) erff ((float) _a.f16);
}

/*
 * Function: erfc ( C: erfcf )
 */
HIDDEN void __fp16_helper_erfc(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) erfcf ((float) _a.f16);
}

/*
 * Function: exp ( C: expf )
 */
HIDDEN void __fp16_helper_exp(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) expf ((float) _a.f16);
}

/*
 * Function: floor ( C: floorf )
 */
HIDDEN void __fp16_helper_floor(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) floorf ((float) _a.f16);
}

/*
 * Function: gamma ( C: tgammaf )
 */
HIDDEN void __fp16_helper_gamma(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) tgammaf ((float) _a.f16);
}

/*
 * Function: log ( C: logf )
 */
HIDDEN void __fp16_helper_log(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) logf ((float) _a.f16);
}

/*
 * Function: log10 ( C: log10f )
 */
HIDDEN void __fp16_helper_log10(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) log10f ((float) _a.f16);
}

/*
 * Function: log_gamma ( C: lgammaf )
 */
HIDDEN void __fp16_helper_log_gamma(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) lgammaf ((float) _a.f16);
}

/*
 * Function: sin ( C: sinf )
 */
HIDDEN void __fp16_helper_sin(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) sinf ((float) _a.f16);
}

/*
 * Function: sinh ( C: sinhf )
 */
HIDDEN void __fp16_helper_sinh(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) sinhf ((float) _a.f16);
}

/*
 * Function: sqrt ( C: sqrtf )
 */
HIDDEN void __fp16_helper_sqrt(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) sqrtf ((float) _a.f16);
}

/*
 * Function: tan ( C: tanf )
 */
HIDDEN void __fp16_helper_tan(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) tanf ((float) _a.f16);
}

/*
 * Function: tanh ( C: tanhf )
 */
HIDDEN void __fp16_helper_tanh(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) tanhf ((float) _a.f16);
}

