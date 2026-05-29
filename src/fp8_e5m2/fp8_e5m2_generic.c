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

#include "emu/lpf_fp8_e5m2_emu.h"
#include "emu/fp8_e5m2_bits.h"
#include "fp8_e5m2_helper.h"
#include "lpf_internal.h"

/*
 * ABS
 */
HIDDEN void __fp8_e5m2_helper_abs(int8_t *out, int8_t in)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = in};
    r->fp8_e5m2 = fp8_e5m2_abs(_a.fp8_e5m2);
}

/*
 * DIGITS
 */
HIDDEN void __fp8_e5m2_helper_digits(int *out)
{
    *out = FP8_E5M2_M_BITS + 1;
}


/*
 * EPSILON
 */
HIDDEN void __fp8_e5m2_helper_epsilon(int8_t *out)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    r->fp8_e5m2 = FP8_E5M2_EPSILON;
}

/*
 * Huge
 */
HIDDEN void __fp8_e5m2_helper_huge(int8_t *out)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    r->fp8_e5m2 = FP8_E5M2_MAX_NUM_POS;
}

/*
 * Tiny
 */
HIDDEN void __fp8_e5m2_helper_tiny(int8_t *out)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    r->fp8_e5m2 = FP8_E5M2_MIN_NUM_POS;
}


/*
 * Minexponent
 */
HIDDEN void __fp8_e5m2_helper_minexponent(int * out)
{
    *out = FP8_E5M2_E_MIN_U ;
}

/*
 * Maxexponent
 */
HIDDEN void __fp8_e5m2_helper_maxexponent(int * out)
{
    *out = FP8_E5M2_E_MAX_U;
}



/*
 * Exponent
 */
HIDDEN void __fp8_e5m2_helper_exponent(int *out, int8_t in)
{
    fp8_e5m2_handler_t r = {.i8 = in};
    if (r.fp8_e5m2 == FP8_E5M2_ZERO1 || r.fp8_e5m2 == FP8_E5M2_ZERO2) {
        *out = 0;
        return;
    }
    *out = FP8_E5M2_GET_EXP(r.fp8_e5m2) - FP8_E5M2_E_BIAS + 1;
}

/*
 * Fraction
 */
HIDDEN void __fp8_e5m2_helper_fraction(int8_t *out, int8_t in)
{
    fp8_e5m2_handler_t *rout = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t r = {.i8 = in};
    if (r.fp8_e5m2 == FP8_E5M2_ZERO1 || r.fp8_e5m2 == FP8_E5M2_ZERO2) {
        *out = 0;
        return;
    }

    uint8_t s = FP8_E5M2_GET_SIGN(r.fp8_e5m2);
    uint8_t m = FP8_E5M2_GET_MANT(r.fp8_e5m2);
    rout->fp8_e5m2 = fp8_e5m2_from_components(s, -1, m);
}

/*
 * ERFC_SCALED
 * We use, like gfortran
 *
 *   ERFC_SCALED(x) = ERFC(x) * EXP(X**2)
 * in higher precision to compute it. In this case intermediate is double!
 */
HIDDEN void __fp8_e5m2_helper_erfc_scaled(int8_t *out, int8_t in)
{
    fp8_e5m2_handler_t *rout = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t r = {.i8 = in};
    double x = (double) fp8_e5m2_to_float(r.fp8_e5m2);
    double p1 = erfc(x);
    double p2 = exp(x * x);
    x = p1 * p2;
    rout->fp8_e5m2 = fp8_e5m2_from_float(x);
}



/*
 * MOD
 */
HIDDEN void __fp8_e5m2_helper_mod(int8_t *out, int8_t in1, int8_t in2)
{
    fp8_e5m2_handler_t a = {.i8 = in1};
    fp8_e5m2_handler_t b = {.i8 = in2};
    fp8_e5m2_handler_t *rout = (fp8_e5m2_handler_t * ) out;
    float _a = fp8_e5m2_to_float(a.fp8_e5m2);
    float _b = fp8_e5m2_to_float(b.fp8_e5m2);

    rout -> fp8_e5m2 = fp8_e5m2_from_float(fmodf(_a, _b));
}

/*
 * MODULO
 */
HIDDEN void __fp8_e5m2_helper_modulo(int8_t *out, int8_t in1, int8_t in2)
{
    fp8_e5m2_handler_t a = {.i8 = in1};
    fp8_e5m2_handler_t b = {.i8 = in2};
    fp8_e5m2_handler_t *rout = (fp8_e5m2_handler_t * ) out;

    float _a = fp8_e5m2_to_float(a.fp8_e5m2);
    float _b = fp8_e5m2_to_float(b.fp8_e5m2);

    float z = fmodf(_a, _b);
    if (z < 0)
        z += (float) _b;
    rout->fp8_e5m2 = fp8_e5m2_from_float(z);
}

/*
 * Nearest
 */
HIDDEN void __fp8_e5m2_helper_nearest(int8_t *out, int8_t in1, int8_t in2)
{
    fp8_e5m2_handler_t a = {.i8 = in1};
    fp8_e5m2_handler_t b = {.i8 = in2};
    fp8_e5m2_handler_t *rout = (fp8_e5m2_handler_t * ) out;
    uint8_t bsign = FP8_E5M2_GET_SIGN(b.u8);

    if ( a.u8  == 0 && bsign > 0) {
        rout -> fp8_e5m2 = a.fp8_e5m2;
        return;
    }
    if ( a.u8 ==  __UINT8_MAX__ && bsign == 0) {
        rout -> fp8_e5m2 = a.fp8_e5m2;
        return;
    }

    if ( bsign ) {
        uint8_t s = FP8_E5M2_GET_SIGN(a.u8);
        int8_t e = FP8_E5M2_GET_EXP(a.u8)-FP8_E5M2_E_BIAS;
        uint8_t m = FP8_E5M2_GET_MANT(a.u8);
        if (m == 0x00)  {
            m = FP8_E5M2_M_MASK;
            e = e - 1;
        }
        rout->fp8_e5m2 = fp8_e5m2_from_components(s, e, m);
    } else {
        uint8_t s = FP8_E5M2_GET_SIGN(a.u8);
        int8_t e = FP8_E5M2_GET_EXP(a.u8)-FP8_E5M2_E_BIAS;
        uint8_t m = FP8_E5M2_GET_MANT(a.u8);
        m = m + 1;
        if (m >  FP8_E5M2_M_MASK)  {
            m = 0x00;
            e = e + 1;
        }
        rout->fp8_e5m2 = fp8_e5m2_from_components(s, e, m);
    }

}

/*
 * NINT
 */
HIDDEN void __fp8_e5m2_helper_nint(int8_t *out, int8_t in)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = in};
    float a = fp8_e5m2_to_float(_a.fp8_e5m2);

    r->fp8_e5m2  = fp8_e5m2_from_float(rintf(a));
}


/*
 * PRECISION
 */
HIDDEN void __fp8_e5m2_helper_precision(int *out)
{
    /* The result has the value INT((DIGITS( x) - 1) * LOG10(RADIX( x))). If RADIX( x) is an integral power of 10, 1 is added */
    *out  = 0;
}

/*
 * RANGE
 */
HIDDEN void __fp8_e5m2_helper_range(int *out)
{
    /* INT(LOG10( HUGE( x) )) */
    *out = 4;
}

/*
 * SCALE
 */
HIDDEN void __fp8_e5m2_helper_scale(int8_t *out, int8_t in, int scl)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = in};

    uint8_t s = FP8_E5M2_GET_SIGN(_a.u8);
    int8_t e = FP8_E5M2_GET_EXP(_a.u8)-FP8_E5M2_E_BIAS;
    uint8_t m = FP8_E5M2_GET_MANT(_a.u8);

    e += scl;
    if ( e > FP8_E5M2_E_MAX_U) {
        r->fp8_e5m2 = s ? FP8_E5M2_INF_NEG : FP8_E5M2_INF_POS;
    } else if ( e < FP8_E5M2_E_MIN_U ) {
        r->fp8_e5m2 = s ? FP8_E5M2_ZERO2 : FP8_E5M2_ZERO1;
    } else {
        r->fp8_e5m2 = fp8_e5m2_from_components(s, e, m);
    }
}

/*
 * SIGN
 */
HIDDEN void __fp8_e5m2_helper_sign(int8_t *out, int8_t in1, int8_t  in2)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = in1};
    fp8_e5m2_handler_t _b = { .i8 = in2};
    uint8_t bsign = FP8_E5M2_GET_SIGN(_b.u8);
    int8_t e = FP8_E5M2_GET_EXP(_a.u8)-FP8_E5M2_E_BIAS;
    uint8_t m = FP8_E5M2_GET_MANT(_a.u8);
    r->fp8_e5m2 = fp8_e5m2_from_components(bsign, e, m);
}

/*
 * IS INF
 */
HIDDEN void __fp8_e5m2_helper_isinf(int *out, int8_t in)
{
    fp8_e5m2_handler_t _a = { .i8 = in};
    *out = fp8_e5m2_isinf(_a.fp8_e5m2);
}

/*
 * IS NAN
 */
HIDDEN void __fp8_e5m2_helper_isnan(int *out, int8_t in)
{
    fp8_e5m2_handler_t _a = { .i8 = in};
    *out = fp8_e5m2_isnan(_a.fp8_e5m2);
}


