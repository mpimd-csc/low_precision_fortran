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
 * Function: atan2 ( C: atan2f )
 */
HIDDEN void __bf16_helper_atan2(int16_t *out, int16_t in1, int16_t in2)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = in1};
    bf16_handler_t _b = { .i16 = in2};
    r->bf16 = (__bf16) atan2f ((float) _a.bf16, (float) _b.bf16);
}

/*
 * Function: bessel_jn ( C: jnf )
 */
HIDDEN void __bf16_helper_bessel_jn(int16_t *out, int in1, int16_t in2)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _b = { .i16 = in2};
    r->bf16 = (__bf16) jnf(in1, (float) _b.bf16);
}

/*
 * Function: bessel_yn ( C: ynf )
 */
HIDDEN void __bf16_helper_bessel_yn(int16_t *out, int in1, int16_t in2)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _b = { .i16 = in2};
    r->bf16 = (__bf16) ynf (in1, (float) _b.bf16);
}

/*
 * Function: hypot ( C: hypotf )
 */
HIDDEN void __bf16_helper_hypot(int16_t *out, int16_t in1, int16_t in2)
{
    bf16_handler_t *r = (bf16_handler_t * ) out;
    bf16_handler_t _a = { .i16 = in1};
    bf16_handler_t _b = { .i16 = in2};
    r->bf16 = (__bf16) hypotf ((float) _a.bf16, (float) _b.bf16);
}

