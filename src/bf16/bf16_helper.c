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
/*
 * Assignment and Constructors
 */
HIDDEN void __bf16_helper_set_from_int(int16_t *out, int in)
{
    bf16_handler_t * x = (bf16_handler_t *) out;
    x->bf16 = (__bf16) in;
}

HIDDEN void __bf16_helper_set_from_float(int16_t *out, float in)
{
    bf16_handler_t * x = (bf16_handler_t *) out;
    x->bf16 = (__bf16) in;
}

HIDDEN void __bf16_helper_set_from_double(int16_t *out, double in)
{
    bf16_handler_t * x = (bf16_handler_t *) out;
    x->bf16 = (__bf16) in;
}

HIDDEN float __bf16_helper_get_float(int16_t in)
{
    bf16_handler_t x;
    x.i16 = in;
    return (float) x.bf16;
}

/*
 * Operator(+)
 */
HIDDEN void __bf16_helper_add_bf16_bf16(int16_t *out, int16_t a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = _a.bf16 + _b.bf16;
}

HIDDEN void __bf16_helper_add_bf16_real(int16_t *out, int16_t a, float b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = _a.bf16 + (__bf16) b;
}

/*
 * Operator(-)
 */
HIDDEN void __bf16_helper_sub_bf16_bf16(int16_t *out, int16_t a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = _a.bf16 - _b.bf16;
}

HIDDEN void __bf16_helper_sub_bf16_real(int16_t *out, int16_t a, float b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = _a.bf16 - (__bf16) b;
}

HIDDEN void __bf16_helper_sub_real_bf16(int16_t *out, float a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = (__bf16 ) a - _b.bf16;
}

HIDDEN void __bf16_helper_unitary_minus(int16_t *out, int16_t a)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = (__bf16 ) - _a.bf16;
}




/*
 * Operator(*)
 */
HIDDEN void __bf16_helper_mul_bf16_bf16(int16_t *out, int16_t a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = _a.bf16 * _b.bf16;
}

HIDDEN void __bf16_helper_mul_bf16_real(int16_t *out, int16_t a, float b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = _a.bf16 * (__bf16) b;
}

/*
 * Operator(/)
 */
HIDDEN void __bf16_helper_div_bf16_bf16(int16_t *out, int16_t a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = _a.bf16 / _b.bf16;
}

HIDDEN void __bf16_helper_div_bf16_real(int16_t *out, int16_t a, float b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = _a.bf16 / (__bf16) b;
}

HIDDEN void __bf16_helper_div_real_bf16(int16_t *out, float a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = (__bf16 ) a / _b.bf16;
}


/*
 * Operator (**)
 */
HIDDEN void __bf16_helper_pow_bf16_bf16(int16_t *out, int16_t a, int16_t b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};
    bf16_handler_t _b = { .i16 = b};

    r->bf16 = (__bf16) powf((float) _a.bf16 ,(float)_b.bf16);
}

HIDDEN void __bf16_helper_pow_bf16_real(int16_t *out, int16_t a, float b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = (__bf16) powf((float) _a.bf16 , b);
}

HIDDEN void __bf16_helper_pow_bf16_int(int16_t *out, int16_t a, int b)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = a};

    r->bf16 = (__bf16) powf((float) _a.bf16 ,(float) b);
}


