//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of LPF-BLAS, a BLAS Implementation for Half-Precision
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
#include <stdint.h>

#include "lpf_internal.h"

__attribute__((constructor))
void lpf_init() {
#ifdef LPF_HAVE_FP16
    if ( ! lpf_cpu_has_basic_fp16() ) {
        fprintf(stderr, "The CPU does not support fp16 numbers. Abort.\n");
        abort();
    }
#endif
#ifdef LPF_HAVE_BF16
    if ( ! lpf_cpu_has_basic_bf16() ) {
        fprintf(stderr, "The CPU does not support bf16 numbers. Abort.\n");
        abort();
    }
#endif



}
