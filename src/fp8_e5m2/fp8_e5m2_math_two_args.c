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

#include "fp8_e5m2/emu/lpf_fp8_e5m2_emu.h"
#include "fp8_e5m2_helper.h"

/*
 * Function: atan2 ( C: atan2f )
 */
HIDDEN void __fp8_e5m2_helper_atan2(int8_t *out, int8_t in1, int8_t in2)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = in1};
    fp8_e5m2_handler_t _b = { .i8 = in2};
    r->fp8_e5m2 = fp8_e5m2_atan2(_a.fp8_e5m2, _b.fp8_e5m2);
}

/*
 * Function: bessel_jn ( C: jnf )
 */
HIDDEN void __fp8_e5m2_helper_bessel_jn(int8_t *out, int in1, int8_t in2)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _b = { .i8 = in2};
    r->fp8_e5m2 = fp8_e5m2_from_float(jnf (in1, fp8_e5m2_to_float( _b.fp8_e5m2)));
}

/*
 * Function: bessel_yn ( C: ynf )
 */
HIDDEN void __fp8_e5m2_helper_bessel_yn(int8_t *out, int in1, int8_t in2)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _b = { .i8 = in2};
    r->fp8_e5m2 = fp8_e5m2_from_float(ynf (in1, fp8_e5m2_to_float( _b.fp8_e5m2)));
}

/*
 * Function: hypot ( C: hypotf )
 */
HIDDEN void __fp8_e5m2_helper_hypot(int8_t *out, int8_t in1, int8_t in2)
{
    fp8_e5m2_handler_t *r = (fp8_e5m2_handler_t * ) out;
    fp8_e5m2_handler_t _a = { .i8 = in1};
    fp8_e5m2_handler_t _b = { .i8 = in2};
    r->fp8_e5m2 = fp8_e5m2_hypot(_a.fp8_e5m2, _b.fp8_e5m2);
}

