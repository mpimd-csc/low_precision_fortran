// SPDX-License-Identifier LGPL-3.0-or-later
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
#ifndef LPF_INTERNAL_H_JPVCTVA6
#define LPF_INTERNAL_H_JPVCTVA6

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
    /* Integer for the hidden value */
#if defined(__INTEL_LLVM_COMPILER) || defined(__ICC)
    /* Intel Compiler (oneAPI and classic) */
    typedef size_t lpf_fortran_strlen_t;
#elif defined (__PGI) || defined(__NVCOMPILER)
    typedef int lpf_blas_fortran_strlen_t;
#elif defined(__aocc__)
    /* CLANG/FLANG as AMD AOCC */
    typedef size_t lpf_fortran_strlen_t;
#elif defined(__clang__)
    /* CLANG/FLANG */
    typedef size_t lpf_fortran_strlen_t;
#elif __GNUC__ > 7
    /* GNU 8.x and newer */
    typedef size_t lpf_fortran_strlen_t;
#else
    /* GNU 4.x - 7.x */
    typedef int lpf_fortran_strlen_t;
#endif


    /* Logical */
#if __GNUC__ >= 5 && !defined (__clang__)
#define lpf_logical_t int_least32_t
#else
#define lpf_logical_t int
#endif

#ifndef lpf_float16_t
#define lpf_float16_t _Float16
#endif

#ifndef lpf_bfloat16_t
#define lpf_bfloat16_t __bf16
#endif


#define LPF_MIN(A,B) (((A)<(B))?(A):(B))
#define LPF_MAX(A,B) (((A)>(B))?(A):(B))
#define LPF_ABS(X) (((X)<(0))?-(X):(X))

#include "lpf_mangle.h"
#include <ISO_Fortran_binding.h>

    typedef void (*lpf_blas_xerbla_func_c_t)(char *str, int32_t *info, lpf_fortran_strlen_t len);
    typedef void (*lpf_blas_xerbla_func_t)(CFI_cdesc_t *msg_desc, int32_t *info);

    void LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)(char *str, int32_t *info, lpf_fortran_strlen_t len);
    void LPF_GLOBAL(lpf_blas_xerbla_set, LPF_BLAS_XERBLA_SET)(void * func, int32_t *style);


    int lpf_cpu_has_basic_fp16(void);
    int lpf_cpu_has_full_fp16(void);
    int lpf_cpu_has_basic_bf16(void);
    int lpf_cpu_has_full_bf16(void);

    double lpf_get_wtime(void);

    /* Bit Mask representation tweak for fp16 */
    typedef union {
        lpf_float16_t f16;
        uint16_t u16;
    } lpf_fp16_u16_t;

    /* Bit Mask representation tweak for fp16 */
    typedef union {
        lpf_bfloat16_t f16;
        uint16_t u16;
    } lpf_bf16_u16_t;


    /* Structures to handle the ISO_C_Binding correctly. */
    typedef struct _lpf_ffloat16_t {
        int16_t value;
    } lpf_ffloat16_t;
    typedef struct _lpf_fbfloat16_t {
        int16_t value;
    } lpf_fbfloat16_t;



#define RETURN_FP16(X) do { lpf_fp16_u16_t ret; ret.f16=(X) ; return ret.u16; } while(0)
#define MAKE_FP16(X) { lpf_fp16_u16_t ret, ret.f16 = (X), ret.u16; }

#define RETURN_BF16(X) do { lpf_bf16_u16_t ret; ret.f16=(X) ; return ret.u16; } while(0)
#define MAKE_BF16(X) { lpf_bf16_u16_t ret, ret.f16 = (X), ret.u16; }




#ifdef __cplusplus
}
#endif
#endif /* end of include guard: LPF_INTERNAL_H_JPVCTVA6 */
