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
 * Function: atan2 ( C: atan2f )
 */
HIDDEN void __fp16_helper_atan2(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in1};
    fp16_handler_t _b = { .i16 = in2};
    r->f16 = (_Float16) atan2f ((float) _a.f16, (float) _b.f16);
}

/*
 * Function: bessel_jn ( C: jnf )
 */
HIDDEN void __fp16_helper_bessel_jn(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in1};
    fp16_handler_t _b = { .i16 = in2};
    r->f16 = (_Float16) jnf ((float) _a.f16, (float) _b.f16);
}

/*
 * Function: bessel_yn ( C: ynf )
 */
HIDDEN void __fp16_helper_bessel_yn(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in1};
    fp16_handler_t _b = { .i16 = in2};
    r->f16 = (_Float16) ynf ((float) _a.f16, (float) _b.f16);
}

/*
 * Function: hypot ( C: hypotf )
 */
HIDDEN void __fp16_helper_hypot(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in1};
    fp16_handler_t _b = { .i16 = in2};
    r->f16 = (_Float16) hypotf ((float) _a.f16, (float) _b.f16);
}

