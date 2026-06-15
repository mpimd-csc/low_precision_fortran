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

program test_level2_spr2
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical('U')
    call test_typical('L')
    call test_stride('U')
    call test_stride('L')



    call test_edge()

    call test_summary()

contains

    subroutine test_typical(uplo)
        character(len=1) :: uplo
        integer(int64) :: n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha
        type(DT), dimension(3) :: x, y
        type(DT), dimension(6) :: ap, ap_old
        integer(int64) :: i, j, je, js

        alpha = 2.0
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]
        ap = 1.0
        ap_old = ap

        call spr2(uplo, n, alpha, x, incx, y, incy, ap)

        do i = 1, 3
            if ( uplo .eq. 'U') then
                js = 1
                je = i
            else
                js = i
                je = n
            end if
            do j = js, je
                block
                    real(real64) :: expected
                    integer(int64) :: pos, i2
                    if ( uplo .eq. 'U') then
                        pos = ((i * ( i - 1))/2 ) + j
                    else
                        i2 = i - 1
                        pos = i2 * n + ( i2 - i2*i2) / 2 + (j-i+1)
                    end if
                    expected = dble(ap_old(pos)) + dble(alpha) *  dble(x(i)) * dble(y(j)) &
                        & +  dble(alpha) *  dble(x(j)) * dble(y(i))
                    call check_dt_real64("spr2_typical_" // uplo, int(pos,int64), ap(pos), expected, GENERIC_TOL)
                end block
            end do
        end do


    end subroutine

    subroutine test_edge()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha
        type(DT), dimension(0) :: x, y
        type(DT), dimension(0) :: ap

        alpha = 1.0
        call spr2('U', n, alpha, x, incx, y, incy, ap)
    end subroutine

    subroutine test_stride(uplo)
        character(len=1) :: uplo
        integer(int64) :: n = 3
        integer(int64) :: incx = 2, incy = 2
        type(DT) :: alpha
        type(DT), dimension(5) :: x, y
        type(DT), dimension(6) :: ap, ap_old
        integer(int64) :: i, j, js, je

        alpha = 1.0
        x = [1.0, 0.0, 2.0, 0.0, 3.0 ]
        y = [1.0, 0.0, 1.0, 0.0, 1.0 ]
        ap = 1.0
        ap_old = ap

        call spr2(uplo, n, alpha, x, incx, y, incy, ap)

        do i = 1, 3
            if ( uplo .eq. 'U') then
                js = 1
                je = i
            else
                js = i
                je = n
            end if
            do j = js, je
                block
                    real(real64) :: expected
                    integer(int64) :: pos, i2
                    if ( uplo .eq. 'U') then
                        pos = ((i * ( i - 1))/2 ) + j
                    else
                        i2 = i - 1
                        pos = i2 * n + ( i2 - i2*i2) / 2 + (j-i+1)
                    end if
                    expected = dble(ap_old(pos)) + dble(alpha) *  dble(x(incx*i-1)) * dble(y(incy*j-1)) &
                        & +  dble(alpha) *  dble(x(incx*j-1)) * dble(y(incy*i-1))
                    call check_dt_real64("spr2_stride_" // uplo, int(pos,int64), ap(pos), expected, GENERIC_TOL)
                end block
            end do
        end do



    end subroutine

end program test_level2_spr2
