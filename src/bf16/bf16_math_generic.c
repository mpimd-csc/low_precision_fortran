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
#include "lpf_internal.h"


lpf_bfloat16_t diff_bf16(lpf_bfloat16_t a, lpf_bfloat16_t b)
{
    if ( a < b )
        return b - a ;
    else
        return a - b;
}
/*
 * ABS
 */
HIDDEN void __bf16_helper_abs(int16_t *out, int16_t in)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = in};
    if ( _a.bf16 < 0 ) {
        r->bf16 = - _a.bf16;
    } else {
        r->bf16 =  _a.bf16;
    }
}

lpf_bfloat16_t abs_bf16(lpf_bfloat16_t a)
{
    if ( a < 0)
        return -a;
    else
        return a;
}


/*
 * DIGITS
 */
HIDDEN void __bf16_helper_digits(int *out)
{
    *out =__BFLT16_MANT_DIG__;
}


/*
 * EPSILON
 */
HIDDEN void __bf16_helper_epsilon(int16_t *out)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    r->bf16 = __BFLT16_EPSILON__;
}

/*
 * Huge
 */
HIDDEN void __bf16_helper_huge(int16_t *out)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    r->bf16 = __BFLT16_MAX__;
}

/*
 * Tiny
 */
HIDDEN void __bf16_helper_tiny(int16_t *out)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    r->bf16 = __BFLT16_MIN__;
}


/*
 * Minexponent
 */
HIDDEN void __bf16_helper_minexponent(int * out)
{
    *out = __BFLT16_MIN_EXP__;
}

/*
 * Maxexponent
 */
HIDDEN void __bf16_helper_maxexponent(int * out)
{
    *out = __BFLT16_MAX_EXP__;
}



/*
 * Exponent
 */
HIDDEN void __bf16_helper_exponent(int *out, int16_t in)
{
    bf16_handler_t r = {.i16 = in};
    float x = (float) r.bf16;
    frexpf(x, out);
}

/*
 * Fraction
 */
HIDDEN void __bf16_helper_fraction(int16_t *out, int16_t in)
{
    bf16_handler_t *rout = (bf16_handler_t * ) out;
    bf16_handler_t r = {.i16 = in};
    float x = (float) r.bf16;
    int dummy = 0;
    float frac = frexpf(x, &dummy);
    rout->bf16 = (__bf16) frac;
}

/*
 * ERFC_SCALED
 * We use, like gfortran
 *
 *   ERFC_SCALED(x) = ERFC(x) * EXP(X**2)
 * in higher precision to compute it. In this case intermediate is double!
 */
HIDDEN void __bf16_helper_erfc_scaled(int16_t *out, int16_t in)
{
    bf16_handler_t *rout = (bf16_handler_t * ) out;
    bf16_handler_t r = {.i16 = in};
    double x = (double) r.bf16;
    double p1 = erfc(x);
    double p2 = exp(x * x);
    x = p1 * p2;
    rout->bf16 = (__bf16) x;
}



/*
 * MOD
 */
HIDDEN void __bf16_helper_mod(int16_t *out, int16_t in1, int16_t in2)
{
    bf16_handler_t a = {.i16 = in1};
    bf16_handler_t b = {.i16 = in2};
    bf16_handler_t *rout = (bf16_handler_t * ) out;
    rout -> bf16 = (__bf16) fmodf((float) a.bf16, (float) b.bf16);
}

/*
 * MODULO
 */
HIDDEN void __bf16_helper_modulo(int16_t *out, int16_t in1, int16_t in2)
{
    bf16_handler_t a = {.i16 = in1};
    bf16_handler_t b = {.i16 = in2};
    bf16_handler_t *rout = (bf16_handler_t * ) out;

    float z = fmod((float) a.bf16, (float) b.bf16);
    if (z < 0)
        z += (float) b.bf16;
    rout->bf16 = (__bf16) z;

}

/*
 * Nearest
 */
HIDDEN void __bf16_helper_nearest(int16_t *out, int16_t in1, int16_t in2)
{
    bf16_handler_t a = {.i16 = in1};
    bf16_handler_t b = {.i16 = in2};
    bf16_handler_t *rout = (bf16_handler_t * ) out;

    if ( a.u16  == 0 && b.bf16 < 0) {
        rout -> bf16 = a.bf16;
        return;
    }
    if ( a.u16 == __UINT16_MAX__ && b.bf16 > 0) {
        rout -> bf16 = a.bf16;
        return;
    }

    if ( b.bf16 < 0 ) {
        rout->u16 = a.u16 - 1;
    } else {
        rout->u16 = a.u16 + 1;
    }

}

/*
 * NINT
 */
HIDDEN void __bf16_helper_nint(int16_t *out, int16_t in)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = in};
    r->bf16  = rintf((float) _a.bf16);
}


/*
 * PRECISION
 */
HIDDEN void __bf16_helper_precision(int *out)
{
    *out  = __BFLT16_DIG__;
}

/*
 * RANGE
 */
HIDDEN void __bf16_helper_range(int *out)
{
    if ( - __BFLT16_MIN_10_EXP__ < __BFLT16_MAX_10_EXP__)
        *out = - __BFLT16_MIN_10_EXP__;
    else
        *out = __BFLT16_MAX_10_EXP__;
}

/*
 * SCALE
 */
HIDDEN void __bf16_helper_scale(int16_t *out, int16_t in, int s)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = in};

    r->bf16 = (__bf16) ( (float)_a.bf16 * powf(2.0,(float)s) );

}

/*
 * SIGN
 */
HIDDEN void __bf16_helper_sign(int16_t *out, int16_t in1, int16_t  in2)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = in1};
    bf16_handler_t _b = { .i16 = in2};

    __bf16 aa = _a.bf16;
    if ( aa < 0 )
        aa = -aa;

    if (_b.bf16 < 0)
        r->bf16 = -aa;
    else
        r->bf16 = aa;
}

lpf_bfloat16_t sign_bf16(lpf_bfloat16_t a, lpf_bfloat16_t b)
{
    if ( a < 0 )
        a = -a;
    if ( b < 0 )
        return -a ;
    else
        return a;
}

/*
 * IS INF
 */
HIDDEN void __bf16_helper_isinf(int *out, int16_t in)
{
    bf16_handler_t _a = { .i16 = in};
    *out = isinf((float) _a.bf16);
}

/*
 * IS NAN
 */
HIDDEN void __bf16_helper_isnan(int *out, int16_t in)
{
    bf16_handler_t _a = { .i16 = in};
    *out = isnan((float) _a.bf16);
}


