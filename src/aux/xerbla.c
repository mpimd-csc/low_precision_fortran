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

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <ISO_Fortran_binding.h>
#include "lpf_internal.h"

static lpf_blas_xerbla_func_t xerbla_function = NULL;
static int xerbla_function_set_from_fortran = 0;

void LPF_GLOBAL(lpf_blas_xerbla, LPF_BLAS_XERBLA)(char *str, lpf_blas_int_t *info, lpf_fortran_strlen_t len)
{
    char tmp[32];

    if ( xerbla_function != NULL && xerbla_function_set_from_fortran ) {
        CFI_index_t ext[1];
        CFI_CDESC_T(1) msg_desc_r1;
        CFI_cdesc_t *msg_desc = (CFI_cdesc_t * ) &msg_desc_r1;
        ext[0] =1;
        CFI_establish(msg_desc, (void *) str, CFI_attribute_pointer, CFI_type_char, len, 1, ext );
        xerbla_function(msg_desc, info);

        return;
    } else if ( xerbla_function != NULL && xerbla_function_set_from_fortran == 0) {
        void (*xerbla_function_c)(char *, lpf_blas_int_t*, lpf_fortran_strlen_t);
        *((void **)&xerbla_function_c) = *((void**)&xerbla_function);
        xerbla_function_c(str, info, len);
        return;
    }

    tmp[31] = 0;
    strncpy(tmp, str, LPF_MAX(len,31));
    tmp[len] = 0;

    fprintf(stderr, "** On entry to %s parameter number %d had an illegal value\n", tmp, (int) *info);

    return;

}

void LPF_GLOBAL(lpf_blas_xerbla_set,LPF_BLAS_XERBLA_SET)(void* func)
{
    xerbla_function = func;
    xerbla_function_set_from_fortran = 0;
    return;
}

void lpf_blas_xerbla_set_function_fortran(void * funcptr)
{
    xerbla_function  = funcptr;
    /** xerbla_function_set_from_fortran = 1; */
    return;
}
