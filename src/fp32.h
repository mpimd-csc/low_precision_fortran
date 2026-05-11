/// SPDX-License-Identifier LGPL-3.0-or-later
/*
 * Copyright (C)  2025 by Martin Koehler
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
#ifndef LPF_FP32_H
#define LPF_FP32_H

/* FP32 */
#define FP32_M_MASK 0x07FFFFF   // 0b0000 0000 0111 1111 1111 1111 1111 1111
#define FP32_M_BITS 23
#define FP32_E_MASK 0x7F800000  // 0b0111 1111 1000 0000 0000 0000 0000 0000
#define FP32_E_SHIFT FP32_M_BITS
#define FP32_E_BITS 8
#define FP32_S_MASK 0x80000000  // 0b1000 0000 0000 0000 0000 0000 0000 0000
#define FP32_S_SHIFT (FP32_E_BITS+FP32_M_BITS)
#define FP32_E_BIAS 127

#define FP32_GET_SIGN(x) (((x) & FP32_S_MASK) >> (FP32_S_SHIFT))
#define FP32_GET_EXP(x) (((x) & FP32_E_MASK) >> (FP32_E_SHIFT))
#define FP32_GET_MANT(x) (((x) & FP32_M_MASK))

#endif
