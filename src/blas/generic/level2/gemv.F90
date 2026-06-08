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
submodule (lpf_blas_fp8_e5m2) lpf_blas_gemv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_gemv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine gemv_impl(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: trans
        integer(int64), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(in) :: x(*)
        type(DT), intent(inout) :: y(*)

        type(DT) :: temp, one, zero
        integer(int64) :: i, ix, iy, j, jx, jy, kx, ky, lenx, leny

        one = 1.0
        zero = 0.0

        if ((m.eq.0) .or. (n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return

        if (trans == 'N' .or. trans == 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        if (incx.gt.0) then
            kx = 1
        else
            kx = 1 - (lenx-1)*incx
        end if

        if (incy.gt.0) then
            ky = 1
        else
            ky = 1 - (leny-1)*incy
        end if

        if (beta.ne.one) then
            if (incy.eq.1) then
                if (beta.eq.zero) then
                    do i = 1, leny
                        y(i) = zero
                    end do
                else
                    do i = 1, leny
                        y(i) = beta*y(i)
                    end do
                end if
            else
                iy = ky
                if (beta.eq.zero) then
                    do i = 1, leny
                        y(iy) = zero
                        iy = iy + incy
                    end do
                else
                    do i = 1, leny
                        y(iy) = beta*y(iy)
                        iy = iy + incy
                    end do
                end if
            end if
        end if

        if (alpha.eq.zero) return

        if (trans == 'N' .or. trans == 'n') then
            jx = kx
            if (incy.eq.1) then
                do j = 1, n
                    temp = alpha*x(jx)
                    do i = 1, m
                        y(i) = y(i) + temp*a(i,j)
                    end do
                    jx = jx + incx
                end do
            else
                do j = 1, n
                    temp = alpha*x(jx)
                    iy = ky
                    do i = 1, m
                        y(iy) = y(iy) + temp*a(i,j)
                        iy = iy + incy
                    end do
                    jx = jx + incx
                end do
            end if
        else
            jy = ky
            if (incx.eq.1) then
                do j = 1, n
                    temp = zero
                    do i = 1, m
                        temp = temp + a(i,j)*x(i)
                    end do
                    y(jy) = y(jy) + alpha*temp
                    jy = jy + incy
                end do
            else
                do j = 1, n
                    temp = zero
                    ix = kx
                    do i = 1, m
                        temp = temp + a(i,j)*x(ix)
                        ix = ix + incx
                    end do
                    y(jy) = y(jy) + alpha*temp
                    jy = jy + incy
                end do
            end if
        end if
    end subroutine

    module subroutine gemv_64(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: trans
        integer(int64), intent(in) :: incx, incy, lda, m, n
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

        if (trans == 'N' .or. trans == 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        lenx = 1 + ( (lenx - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        leny = 1 + ( (leny - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        call gemv_impl(trans, m, n, alpha, pa, lda, px, incx, beta, py, incy)
    end subroutine

    module subroutine gemv_32(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: trans
        integer(int32), intent(in) :: incx, incy, lda, m, n
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

        if (trans == 'N' .or. trans == 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        lenx = 1 + ( (lenx - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        leny = 1 + ( (leny - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        call gemv_impl(trans, int(m, int64), int(n, int64), &
            & alpha, pa, int(lda, int64), px, int(incx, int64), beta, py, int(incy, int64))
    end subroutine

end submodule
