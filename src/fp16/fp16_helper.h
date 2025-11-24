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


#ifndef FP16_HELPER_H_P7LSSW6X
#define FP16_HELPER_H_P7LSSW6X

#include <stdint.h>
#include <math.h>
#include "constants.h"

#define HIDDEN __attribute__((visibility("hidden")))

typedef union {
    _Float16 f16;
    int16_t  i16;
    uint16_t u16;
} fp16_handler_t;


#endif /* end of include guard: FP16_HELPER_H_P7LSSW6X */
