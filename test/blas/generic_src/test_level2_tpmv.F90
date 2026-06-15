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

program test_level2_tpmv
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical('U', 'N', 'N')
    call test_typical('L', 'N', 'N')
    call test_typical('U', 'T', 'N')
    call test_typical('L', 'T', 'N')

    call test_typical('U', 'N', 'U')
    call test_typical('L', 'N', 'U')
    call test_typical('U', 'T', 'U')
    call test_typical('L', 'T', 'U')

    call test_stride('U', 'N', 'N')
    call test_stride('L', 'N', 'N')
    call test_stride('U', 'T', 'N')
    call test_stride('L', 'T', 'N')

    call test_stride('U', 'N', 'U')
    call test_stride('L', 'N', 'U')
    call test_stride('U', 'T', 'U')
    call test_stride('L', 'T', 'U')


    call test_edge()

    call test_summary()

contains

    subroutine test_typical(uplo, trans, diag)
        character(len=1) :: uplo, trans, diag
        integer(int64) :: n = 3
        integer(int64) :: incx = 1
        type(DT), dimension(6) :: ap
        type(DT), dimension(3,3) :: a, b
        type(DT), dimension(3) :: x, x_old
        integer :: i

        ap = 0.0
        ! For ap = 0.0, A is zero matrix. Let's use 1.0 for identity.
        ap = 1.0
        x = [1.0, 2.0, 3.0]

        a = reshape ( (/ 1, 0, 0, 1, 1, 0, 1, 1, 1 /), (/3,3/))

        if (uplo.eq.'L') then
            b = TRANSPOSE(a)
            a = b
        end if

        call gemv(trans, n, n, DT(1.0), a, n, x, 1_int64, DT(0.0), x_old, 1_int64)

        call tpmv(uplo, trans, diag, n, ap, x, incx)

        do i = 1, 3
            call check_dt_real64("tpmv_typical_"//uplo//trans//diag, int(i,int64), x(i), dble(x_old(i)), GENERIC_TOL)
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1
        type(DT), dimension(0) :: ap
        type(DT), dimension(0) :: x

        call tpmv('U', 'N', 'N', n, ap, x, incx)
        call tpmv('L', 'N', 'N', n, ap, x, incx)
    end subroutine

    subroutine test_stride(uplo, trans, diag)
        character(len=1) :: uplo, trans, diag
        integer(int64) :: n = 2
        integer(int64) :: incx = 2
        type(DT), dimension(3) :: ap
        type(DT), dimension(3) :: x, x_old
        integer(int64) :: i
        type(DT), dimension(2,2) :: a, b

        ap = 1.0
        x = [1.0, 0.0, 2.0]
        x_old = x

        a = reshape ( (/ 1, 0,  1, 1 /), (/2,2/))

        if (uplo.eq.'L') then
            b = TRANSPOSE(a)
            a = b
        end if

        call gemv(trans, n, n, DT(1.0), a, n, x, incx, DT(0.0), x_old, incx)

        call tpmv(uplo, trans, diag, n, ap, x, incx)

        do i = 1, 2
            call check_dt_real64("tpmv_stride_"//uplo//trans//diag, i, x(2*i-1), dble(x_old(2*i-1)), GENERIC_TOL)
        end do
    end subroutine

end program test_level2_tpmv
