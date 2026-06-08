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
submodule (lpf_blas_fp8_e5m2) lpf_blas_symv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_symv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine symv_impl(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, lda, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(in) :: x(*)
        type(DT), intent(inout) :: y(*)

        type(DT) :: temp1, temp2, one, zero
        integer(int64) :: i, ix, iy, j, jx, jy, kx, ky

        one = 1.0
        zero = 0.0

        if ((n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return

        if (incx.gt.0) then
            kx = 1
        else
            kx = 1 - (n-1)*incx
        end if

        if (incy.gt.0) then
            ky = 1
        else
            ky = 1 - (n-1)*incy
        end if

        if (beta.ne.one) then
            if (incy.eq.1) then
                if (beta.eq.zero) then
                    do i = 1, n
                        y(i) = zero
                    end do
                else
                    do i = 1, n
                        y(i) = beta*y(i)
                    end do
                end if
            else
                iy = ky
                if (beta.eq.zero) then
                    do i = 1, n
                        y(iy) = zero
                        iy = iy + incy
                    end do
                else
                    do i = 1, n
                        y(iy) = beta*y(iy)
                        iy = iy + incy
                    end do
                end if
            end if
        end if

        if (alpha.eq.zero) return

        if (uplo == 'U' .or. uplo == 'u') then
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    temp1 = alpha*x(j)
                    temp2 = zero
                    do i = 1, j - 1
                        y(i) = y(i) + temp1*a(i,j)
                        temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + temp1*a(j,j) + alpha*temp2
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    temp1 = alpha*x(jx)
                    temp2 = zero
                    ix = kx
                    iy = ky
                    do i = 1, j - 1
                        y(iy) = y(iy) + temp1*a(i,j)
                        temp2 = temp2 + a(i,j)*x(ix)
                        ix = ix + incx
                        iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                end do
            end if
        else
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    temp1 = alpha*x(j)
                    temp2 = zero
                    y(j) = y(j) + temp1*a(j,j)
                    do i = j + 1, n
                        y(i) = y(i) + temp1*a(i,j)
                        temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + alpha*temp2
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    temp1 = alpha*x(jx)
                    temp2 = zero
                    y(jy) = y(jy) + temp1*a(j,j)
                    ix = jx
                    iy = jy
                    do i = j + 1, n
                        ix = ix + incx
                        iy = iy + incy
                        y(iy) = y(iy) + temp1*a(i,j)
                        temp2 = temp2 + a(i,j)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                end do
            end if
        end if
    end subroutine

    module subroutine symv_64(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, lda, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        call symv_impl(uplo, n, alpha, pa, lda, px, incx, beta, py, incy)
    end subroutine

    module subroutine symv_32(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, incy, lda, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (int(n, int64) - 1) * abs(int(incy, int64)))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        call symv_impl(uplo, int(n, int64), alpha, pa, int(lda, int64), px, int(incx, int64), beta, py, int(incy, int64))
    end subroutine

end submodule
