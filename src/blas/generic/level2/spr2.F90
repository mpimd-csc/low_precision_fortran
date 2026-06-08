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

#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_spr2_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_spr2_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    subroutine spr2_impl(uplo, n, alpha, x, incx, y, incy, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: x(*)
        type(DT), intent(in) :: y(*)
        type(DT), intent(inout) :: ap(*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, iy, j, jx, jy, k, kk, kx, ky

        zero = 0.0
        if ((n.eq.0) .or. (alpha.eq.zero)) return

        if (incx.le.0) then
            kx = 1 - (n-1)*incx
        else if (incx.ne.1) then
            kx = 1
        end if

        if (incy.le.0) then
            ky = 1 - (n-1)*incy
        else if (incy.ne.1) then
            ky = 1
        end if

        kk = 1
        if (uplo == 'U' .or. uplo == 'u') then
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    k = kk;
                    do i = 1, j
                        temp = alpha*(x(i)*y(j) + y(i)*x(j))
                        ap(k) = ap(k) + temp
                        k = k + 1
                    end do
                    kk = kk + j
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    k = kk
                    ix = kx
                    iy = ky
                    do i = 1, j
                        ap(k) = ap(k) + alpha*(x(ix)*y(jy) + y(iy)*x(jx))
                        ix = ix + incx
                        k = k + 1
                    end do
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                end do
            end if
        else
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    k = kk
                    do i = j, n
                        ap(k) = ap(k) + alpha*(x(i)*y(j) + y(i)*x(j))
                        k = k + 1
                    end do
                    kk = kk + n - j + 1
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    k = kk
                    ix = jx
                    iy = jy
                    do i = j, n
                        ap(k) = ap(k) + alpha*(x(ix)*y(jy) + y(iy)*x(jx))
                        ix = ix + incx
                        k = k + 1
                    end do
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + n - j + 1
                end do
            end if
        end if
    end subroutine

    module subroutine spr2_64(uplo, n, alpha, x, incx, y, incy, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: ap(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pap(:)
        integer(int64) :: lenx, leny, lenap
        integer(int64) :: lx(1), ly(1), lap(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        lenap = (n * (n + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        call spr2_impl(uplo, n, alpha, px, incx, py, incy, pap)
    end subroutine

    module subroutine spr2_32(uplo, n, alpha, x, incx, y, incy, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: ap(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pap(:)
        integer(int64) :: lenx, leny, lenap
        integer(int64) :: lx(1), ly(1), lap(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (int(n, int64) - 1) * abs(int(incy, int64)))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        lenap = (int(n, int64) * (int(n, int64) + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        call spr2_impl(uplo, int(n, int64), alpha, px, int(incx, int64), py, int(incy, int64), pap)
    end subroutine

end submodule
