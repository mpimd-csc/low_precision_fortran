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
 * Assignment and Constructors
 */
HIDDEN void __fp16_helper_set_from_int(int16_t *out, int in)
{
    fp16_handler_t * x = (fp16_handler_t *) out;
    x->f16 = (_Float16) in;
}

HIDDEN void __fp16_helper_set_from_float(int16_t *out, float in)
{
    fp16_handler_t * x = (fp16_handler_t *) out;
    x->f16 = (_Float16) in;
}

HIDDEN void __fp16_helper_set_from_double(int16_t *out, double in)
{
    fp16_handler_t * x = (fp16_handler_t *) out;
    x->f16 = (_Float16) in;
}

HIDDEN float __fp16_helper_get_float(int16_t in)
{
    fp16_handler_t x;
    x.i16 = in;
    return (float) x.f16;
}

/*
 * Operator(+)
 */
HIDDEN void __fp16_helper_add_fp16_fp16(int16_t *out, int16_t a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};
    fp16_handler_t _b = { .i16 = b};

    r->f16 = _a.f16 + _b.f16;
}

HIDDEN void __fp16_helper_add_fp16_real(int16_t *out, int16_t a, float b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};

    r->f16 = _a.f16 + (_Float16) b;
}

/*
 * Operator(-)
 */
HIDDEN void __fp16_helper_sub_fp16_fp16(int16_t *out, int16_t a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};
    fp16_handler_t _b = { .i16 = b};

    r->f16 = _a.f16 - _b.f16;
}

HIDDEN void __fp16_helper_sub_fp16_real(int16_t *out, int16_t a, float b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};

    r->f16 = _a.f16 - (_Float16) b;
}

HIDDEN void __fp16_helper_sub_real_fp16(int16_t *out, float a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _b = { .i16 = b};

    r->f16 = (_Float16 ) a + _b.f16;
}



/*
 * Operator(*)
 */
HIDDEN void __fp16_helper_mul_fp16_fp16(int16_t *out, int16_t a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};
    fp16_handler_t _b = { .i16 = b};

    r->f16 = _a.f16 * _b.f16;
}

HIDDEN void __fp16_helper_mul_fp16_real(int16_t *out, int16_t a, float b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};

    r->f16 = _a.f16 * (_Float16) b;
}

/*
 * Operator(/)
 */
HIDDEN void __fp16_helper_div_fp16_fp16(int16_t *out, int16_t a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};
    fp16_handler_t _b = { .i16 = b};

    r->f16 = _a.f16 / _b.f16;
}

HIDDEN void __fp16_helper_div_fp16_real(int16_t *out, int16_t a, float b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};

    r->f16 = _a.f16 - (_Float16) b;
}

HIDDEN void __fp16_helper_div_real_fp16(int16_t *out, float a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _b = { .i16 = b};

    r->f16 = (_Float16 ) a / _b.f16;
}


/*
 * Operator (**)
 */
HIDDEN void __fp16_helper_pow_fp16_fp16(int16_t *out, int16_t a, int16_t b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};
    fp16_handler_t _b = { .i16 = b};

    r->f16 = (_Float16) powf((float) _a.f16 ,(float)_b.f16);
}

HIDDEN void __fp16_helper_pow_fp16_real(int16_t *out, int16_t a, float b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};

    r->f16 = (_Float16) powf((float) _a.f16 , b);
}

HIDDEN void __fp16_helper_pow_fp16_int(int16_t *out, int16_t a, int b)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = a};

    r->f16 = (_Float16) powf((float) _a.f16 ,(float) b);
}

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

