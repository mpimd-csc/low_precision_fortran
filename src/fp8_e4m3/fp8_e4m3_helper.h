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


#ifndef BF16_HELPER_H_P7LSSW6X
#define BF16_HELPER_H_P7LSSW6X

#include <stdint.h>
#include "emu/fp8_e4m3.h"
#include "lpf_internal.h"
#include <math.h>
#include "constants.h"

#define HIDDEN __attribute__((visibility("hidden")))

typedef union {
    fp8_e4m3_t fp8_e4m3;
    int8_t  i8;
    uint8_t u8;
} fp8_e4m3_handler_t;

#endif /* end of include guard: BF16_HELPER_H_P7LSSW6X */
