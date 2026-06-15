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

program test_level2_ger
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_edge()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical()
        integer(int64) :: m = 3, n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha
        type(DT), dimension(3) :: x, y
        type(DT), dimension(3, 3) :: a, a_old
        integer :: i, j

        alpha = 2.0
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]
        a = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        a_old = a

        call ger(m, n, alpha, x, incx, y, incy, a, m)

        do i = 1, 3
            do j = 1, 3
                call check_dt_real64("ger_typical", int(j*10+i,int64), a(j,i), &
                    & dble(a_old(j,i)) + 2.0_real64 * dble(x(j)) * dble(y(i)), GENERIC_TOL)
            end do
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: m = 0, n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha
        type(DT), dimension(0) :: x
        type(DT), dimension(3) :: y
        type(DT), dimension(0, 3) :: a

        alpha = 1.0
        call ger(m, n, alpha, x, incx, y, incy, a, int(1, int64))
    end subroutine

    subroutine test_stride()
        integer(int64) :: m = 2, n = 2
        integer(int64) :: incx = 2, incy = 2
        type(DT) :: alpha
        type(DT), dimension(3) :: x, y
        type(DT), dimension(2, 2) :: a, a_old
        integer :: i, j

        alpha = 1.0
        x = [1.0, 0.0, 2.0]
        y = [1.0, 0.0, 1.0]
        a = reshape([1.0, 2.0, 3.0, 4.0], [2, 2])
        a_old = a

        call ger(m, n, alpha, x, incx, y, incy, a, m)

        do i = 1, 2
            do j = 1, 2
                call check_dt_real64("ger_stride", int(j*10+i,int64), a(j,i), &
                    & dble(a_old(j,i)) + 1.0_real64 * dble(x(2*j-1)) * dble(y(2*i-1)), GENERIC_TOL)
            end do
        end do
    end subroutine

end program test_level2_ger
