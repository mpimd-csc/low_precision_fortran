//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of HPBLAS, a BLAS Implementation for Half-Precision
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

// Function to get CPUID information
static void cpuid(uint32_t eax_in, uint32_t ecx_in, uint32_t *eax, uint32_t *ebx, uint32_t *ecx, uint32_t *edx) {
    asm volatile (
        "cpuid"
        : "=a" (*eax), "=b" (*ebx), "=c" (*ecx), "=d" (*edx)
        : "a" (eax_in), "c" (ecx_in)
    );
    /** printf("eax_in: %lx, ecx_in: %lx, eax = %lx, ebx = %lx, ecx = %lx, edx = %lx\n", eax_in, ecx_in, *eax, *ebx, *ecx, *edx ); */
}

// Function to check for AVX512FP16 flag
static int cpu_has_avx512fp16() {
    uint32_t eax, ebx, ecx, edx;

    // CPUID with EAX = 7, ECX = 0
    cpuid(7, 0, &eax, &ebx, &ecx, &edx);

    // Check bit 23 of ECX (AVX512FP16)
    return (ecx & (1 << 23)) != 0;
}

// Function to check for F16C flag
static int cpu_has_f16c() {
    uint32_t eax, ebx, ecx, edx;

    // CPUID with EAX = 1
    cpuid(1, 0, &eax, &ebx, &ecx, &edx);

    // Check bit 29 of ECX (F16C)
    int ret = (ecx & (1 << 29)) != 0;
    return ret;
}

// Function to check for AVX512F flag
static int cpu_has_avx512f() {
    uint32_t eax, ebx, ecx, edx;

    // CPUID with EAX = 7, ECX = 0
    cpuid(7, 0, &eax, &ebx, &ecx, &edx);

    // Check bit 16 of EBX (AVX512F)
    return (ebx & (1 << 16)) != 0;
}

// Function to check for AVX512BF16 flag
static int cpu_has_avx512bf16() {
    uint32_t eax, ebx, ecx, edx;

    // CPUID with EAX = 7, ECX = 0
    cpuid(7, 0, &eax, &ebx, &ecx, &edx);

    // Check bit 5 of ECX (AVX512BF16)
    return (ecx & (1 << 5)) != 0;
}

int lpf_cpu_has_basic_fp16()
{
	return cpu_has_f16c();
}

int lpf_cpu_has_full_fp16() {
	return cpu_has_avx512f() && cpu_has_avx512fp16();
}

int lpf_cpu_has_basic_bf16() {
	return 1;
}

int lpf_cpu_has_full_bf16() {
	return cpu_has_avx512bf16();
}


