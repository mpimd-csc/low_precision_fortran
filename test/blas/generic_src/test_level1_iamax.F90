!  SPDX-License-Identifier: LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program; if not, write to the Free Software Foundation,
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!


#ifdef LPF_FP8_E4M3
#define DT fp8_e4m3
#define BLASMOD lpf_blas_fp8_e4m3
#define TYPEMOD lpf_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
#define DT fp8_e5m2
#define BLASMOD lpf_blas_fp8_e5m2
#define TYPEMOD lpf_fp8_e5m2
#endif
#ifdef LPF_FP16
#define DT fp16
#define BLASMOD lpf_blas_fp16
#define TYPEMOD lpf_fp16
#endif
#ifdef LPF_BF16
#define DT bf16
#define BLASMOD lpf_blas_bf16
#define TYPEMOD lpf_bf16
#endif

program test_level1_iamax
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_n0()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical()
        integer(int64) :: n = 10
        integer(int64) :: incx = 1
        type(DT), dimension(10) :: sx
        integer(int64) :: result

        sx = [1.0, 2.0, 3.0, 10.0, 5.0, 6.0, 7.0, 8.0, 9.0, 4.0]

        result = iamax(n, sx, incx)
        call check_integer("iamax_typical", result, 4_int64)
    end subroutine

    subroutine test_n0()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1
        type(DT), dimension(0) :: sx
        integer(int64) :: result

        result = iamax(n, sx, incx)
        ! iamax(0) behavior is implementation defined, usually 0 or 1.
        ! We just check it doesn't crash.
    end subroutine

    subroutine test_stride()
        integer(int64) :: n = 5
        integer(int64) :: incx = 2
        type(DT), dimension(9) :: sx
        integer(int64) :: result

        sx = [1.0, 0.0, 2.0, 0.0, 10.0, 0.0, 4.0, 0.0, 5.0]

        result = iamax(n, sx, incx)
        call check_integer("iamax_stride", result, 3_int64)
    end subroutine

end program test_level1_iamax
