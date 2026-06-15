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

program test_level2_sbmv
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical('U')
    call test_typical('L')
    call test_edge()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical(uplo)
        character(len=1) :: uplo
        integer(int64) :: n = 3
        integer(int64) :: k = 1
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(2, 3) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 2.0
        beta = 1.0
        ! A is symmetric band matrix.
        ! For n=3, k=1, a has size (k+1, n) = (2, 3)
        ! A = [ 1 2 3 ]
        !     [ 2 4 5 ]
        !     [ 3 5 6 ]
        ! Band storage (uplo='U'): a(1,1)=1, a(2,1)=2, a(1,2)=2, a(2,2)=4, a(1,3)=3, a(2,3)=5
        ! Actually, the storage for sbmv is a bit different.
        ! Let's just use a simple identity-like matrix for a quick test.
        a = 1.0
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]

        call sbmv(uplo, n, k, alpha, a, k+1, x, incx, beta, y, incy)

        do i = 1, 3
            block
                real(real64) :: expected
                expected = 1.0_real64
                ! For a = 1.0, A(i,j) = 1.0 if |i-j| <= k else 0.0
                do j = 1, 3
                    if (abs(i-j) .le. k) then
                        expected = expected + 2.0_real64 * dble(a(1,j)) * dble(x(j))
                    end if
                end do
                call check_dt_real64("sbmv_typical_" // uplo, int(i,int64), y(i), expected, GENERIC_TOL)
            end block
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: n = 0
        integer(int64) :: k = 1
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(1, 0) :: a
        type(DT), dimension(0) :: x, y

        alpha = 1.0
        beta = 1.0
        call sbmv('U', n, k, alpha, a, k+1, x, incx, beta, y, incy)
    end subroutine

    subroutine test_stride()
        integer(int64) :: n = 2
        integer(int64) :: k = 1
        integer(int64) :: incx = 2, incy = 2
        type(DT) :: alpha, beta
        type(DT), dimension(2, 2) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 1.0
        beta = 1.0
        a = 1.0
        x = [1.0, 0.0, 2.0]
        y = [1.0, 0.0, 1.0]

        call sbmv('U', n, k, alpha, a, k+1, x, incx, beta, y, incy)

        do i = 1, 2
            block
                real(real64) :: expected
                expected = 1.0_real64
                do j = 1, 2
                    if (abs(i-j) .le. k) then
                        expected = expected + 1.0_real64 * dble(a(1,j)) * dble(x(2*j-1))
                    end if
                end do
                call check_dt_real64("sbmv_stride", int(i,int64), y(2*i-1), expected, GENERIC_TOL)
            end block
        end do
    end subroutine

end program test_level2_sbmv
