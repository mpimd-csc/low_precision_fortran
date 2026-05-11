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

#include "fp8_e5m2_helper.h"
#include "fp8_e5m2/emu/fp8_e5m2.h"
/*
 * Assignment and Constructors
 */
HIDDEN void __fp8_e5m2_helper_set_from_int(int8_t *out, int in)
{
    fp8_e5m2_handler_t * x = (fp8_e5m2_handler_t *) out;
    x->fp8_e5m2 = fp8_e5m2_from_float((float) in);
}

HIDDEN void __fp8_e5m2_helper_set_from_float(int8_t *out, float in)
{
    fp8_e5m2_handler_t * x = (fp8_e5m2_handler_t *) out;
    x->fp8_e5m2 = fp8_e5m2_from_float(in);
}

HIDDEN void __fp8_e5m2_helper_set_from_double(int8_t *out, double in)
{
    fp8_e5m2_handler_t * x = (fp8_e5m2_handler_t *) out;
    x->fp8_e5m2 = fp8_e5m2_from_float(in);
}

HIDDEN float __fp8_e5m2_helper_get_float(int8_t in)
{
    fp8_e5m2_handler_t x;
    x.i8 = in;
    return (float)  fp8_e5m2_to_float(x.fp8_e5m2);
}

/*
 * Operator(+)
 */
HIDDEN void __fp8_e5m2_helper_add_fp8_e5m2_fp8_e5m2(int8_t *out, int8_t a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_handler_t _b = { .i8 = b};

    r->fp8_e5m2 = fp8_e5m2_add(_a.fp8_e5m2, _b.fp8_e5m2);
}

HIDDEN void __fp8_e5m2_helper_add_fp8_e5m2_real(int8_t *out, int8_t a, float b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_t _b = fp8_e5m2_from_float(b);

    r->fp8_e5m2 = fp8_e5m2_add(_a.fp8_e5m2, _b);
}

/*
 * Operator(-)
 */
HIDDEN void __fp8_e5m2_helper_sub_fp8_e5m2_fp8_e5m2(int8_t *out, int8_t a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_handler_t _b = { .i8 = b};

    r->fp8_e5m2 = fp8_e5m2_sub(_a.fp8_e5m2, _b.fp8_e5m2);
}

HIDDEN void __fp8_e5m2_helper_sub_fp8_e5m2_real(int8_t *out, int8_t a, float b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_t _b = fp8_e5m2_from_float(b);
    r->fp8_e5m2 = fp8_e5m2_sub(_a.fp8_e5m2, _b);
}

HIDDEN void __fp8_e5m2_helper_sub_real_fp8_e5m2(int8_t *out, float a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _b = { .i8 = b};
    fp8_e5m2_t _a = fp8_e5m2_from_float(a);
    r->fp8_e5m2 = fp8_e5m2_sub(_a, _b.fp8_e5m2);
}

HIDDEN void __fp8_e5m2_helper_unitary_minus(int8_t *out, int8_t a)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a };
    r->fp8_e5m2 = fp8_e5m2_unitary_minus(_a.fp8_e5m2) ;
}


/*
 * Operator(*)
 */
HIDDEN void __fp8_e5m2_helper_mul_fp8_e5m2_fp8_e5m2(int8_t *out, int8_t a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_handler_t _b = { .i8 = b};

    r->fp8_e5m2 = fp8_e5m2_mul(_a.fp8_e5m2,_b.fp8_e5m2);
}

HIDDEN void __fp8_e5m2_helper_mul_fp8_e5m2_real(int8_t *out, int8_t a, float b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_t _b = fp8_e5m2_from_float(b);

    r->fp8_e5m2 = fp8_e5m2_mul(_a.fp8_e5m2, _b) ;
}

/*
 * Operator(/)
 */
HIDDEN void __fp8_e5m2_helper_div_fp8_e5m2_fp8_e5m2(int8_t *out, int8_t a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_handler_t _b = { .i8 = b};

    r->fp8_e5m2 = fp8_e5m2_div(_a.fp8_e5m2,  _b.fp8_e5m2);
}

HIDDEN void __fp8_e5m2_helper_div_fp8_e5m2_real(int8_t *out, int8_t a, float b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = a};
    fp8_e5m2_t _b = fp8_e5m2_from_float(b);

    r->fp8_e5m2 = fp8_e5m2_div(_a.fp8_e5m2, _b);
}

HIDDEN void __fp8_e5m2_helper_div_real_fp8_e5m2(int8_t *out, float a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _b = { .i8 = b};
    fp8_e5m2_t _a = fp8_e5m2_from_float(a);

    r->fp8_e5m2 = fp8_e5m2_div(_a, _b.fp8_e5m2);
}


/*
 * Operator (**)
 */
HIDDEN void __fp8_e5m2_helper_pow_fp8_e5m2_fp8_e5m2(int8_t *out, int8_t a, int8_t b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;

    float _a = fp8_e5m2_to_float(a);
    float _b = fp8_e5m2_to_float(b);
    float _r = powf(_a, _b);

    r->fp8_e5m2 = fp8_e5m2_from_float(_r) ;
}

HIDDEN void __fp8_e5m2_helper_pow_fp8_e5m2_real(int8_t *out, int8_t a, float b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    float _a = fp8_e5m2_to_float(a);
    float _r = powf(_a, b);
    r->fp8_e5m2 = fp8_e5m2_from_float(_r) ;
}

HIDDEN void __fp8_e5m2_helper_pow_fp8_e5m2_int(int8_t *out, int8_t a, int b)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    float _a = fp8_e5m2_to_float(a);
    float _r = powf(_a, b);
    r->fp8_e5m2 = fp8_e5m2_from_float(_r) ;
}


