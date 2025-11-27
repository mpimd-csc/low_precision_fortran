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
#include <stdint.h>
#include <math.h>

#include "fp16_helper.h"
#include "lpf_internal.h"

/*
 *
 */
lpf_float16_t diff_fp16(lpf_float16_t a, lpf_float16_t b)
{
    if ( a < b )
        return b - a;
    else
        return a - b;

}

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

lpf_float16_t abs_fp16(lpf_float16_t x)
{
    if ( x < 0)
        return -x;
    else
        return x;
}

/*
 * DIGITS
 */
HIDDEN void __fp16_helper_digits(int *out)
{
    *out =__FLT16_MANT_DIG__;
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
 * Huge
 */
HIDDEN void __fp16_helper_huge(int16_t *out)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    r->f16 = __FLT16_MAX__;
}

/*
 * Tiny
 */
HIDDEN void __fp16_helper_tiny(int16_t *out)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    r->f16 = __FLT16_MIN__;
}


/*
 * Minexponent
 */
HIDDEN void __fp16_helper_minexponent(int * out)
{
    *out = __FLT16_MIN_EXP__;
}

/*
 * Maxexponent
 */
HIDDEN void __fp16_helper_maxexponent(int * out)
{
    *out = __FLT16_MAX_EXP__;
}



/*
 * Exponent
 */
HIDDEN void __fp16_helper_exponent(int *out, int16_t in)
{
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

/*
 * ERFC_SCALED
 * We use, like gfortran
 *
 *   ERFC_SCALED(x) = ERFC(x) * EXP(X**2)
 * in higher precision to compute it. In this case intermediate is double!
 */
HIDDEN void __fp16_helper_erfc_scaled(int16_t *out, int16_t in)
{
    fp16_handler_t *rout = (fp16_handler_t * ) out;
    fp16_handler_t r = {.i16 = in};
    double x = (double) r.f16;
    double p1 = erfc(x);
    double p2 = exp(x * x);
    x = p1 * p2;
    rout->f16 = (_Float16) x;
}



/*
 * MOD
 */
HIDDEN void __fp16_helper_mod(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t a = {.i16 = in1};
    fp16_handler_t b = {.i16 = in2};
    fp16_handler_t *rout = (fp16_handler_t * ) out;
    rout -> f16 = (_Float16) fmodf((float) a.f16, (float) b.f16);
}

/*
 * MODULO
 */
HIDDEN void __fp16_helper_modulo(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t a = {.i16 = in1};
    fp16_handler_t b = {.i16 = in2};
    fp16_handler_t *rout = (fp16_handler_t * ) out;

    float z = fmod((float) a.f16, (float) b.f16);
    if (z < 0)
        z += (float) b.f16;
    rout->f16 = (_Float16) z;

}

/*
 * Nearest
 */
HIDDEN void __fp16_helper_nearest(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t a = {.i16 = in1};
    fp16_handler_t b = {.i16 = in2};
    fp16_handler_t *rout = (fp16_handler_t * ) out;

    if ( a.u16  == 0 && b.f16 < 0) {
        rout -> f16 = a.f16;
        return;
    }
    if ( a.u16 == __UINT16_MAX__ && b.f16 > 0) {
        rout -> f16 = a.f16;
        return;
    }

    if ( b.f16 < 0 ) {
        rout->u16 = a.u16 - 1;
    } else {
        rout->u16 = a.u16 + 1;
    }

}

/*
 * NINT
 */
HIDDEN void __fp16_helper_nint(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16  = rintf((float) _a.f16);
}


/*
 * PRECISION
 */
HIDDEN void __fp16_helper_precision(int *out)
{
    *out  = __FLT16_DIG__;
}

/*
 * RANGE
 */
HIDDEN void __fp16_helper_range(int *out)
{
    if ( - __FLT16_MIN_10_EXP__ < __FLT16_MAX_10_EXP__)
        *out = - __FLT16_MIN_10_EXP__;
    else
        *out = __FLT16_MAX_10_EXP__;
}

/*
 * SCALE
 */
HIDDEN void __fp16_helper_scale(int16_t *out, int16_t in, int s)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};

    r->f16 = (_Float16) ( (float)_a.f16 * powf(2.0,(float)s) );

}

/*
 * SIGN
 */
HIDDEN void __fp16_helper_sign(int16_t *out, int16_t in1, int16_t  in2)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in1};
    fp16_handler_t _b = { .i16 = in2};

    _Float16 aa = _a.f16;
    if ( aa < 0 )
        aa = -aa;

    if (_b.f16 < 0)
        r->f16 = -aa;
    else
        r->f16 = aa;
}

lpf_float16_t sign_fp16(lpf_float16_t a, lpf_float16_t b)
{
    if (a < 0 )
    {
        a = -a;
    }
    if ( b < 0)
        return -a;
    else
        return a;
}

/*
 * IS INF
 */
HIDDEN void __fp16_helper_isinf(int *out, int16_t in)
{
    fp16_handler_t _a = { .i16 = in};
    *out = isinf((float) _a.f16);
}

/*
 * IS NAN
 */
HIDDEN void __fp16_helper_isnan(int *out, int16_t in)
{
    fp16_handler_t _a = { .i16 = in};
    *out = isnan((float) _a.f16);
}


