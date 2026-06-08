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
submodule (lpf_blas_fp8_e5m2) lpf_blas_ger_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_ger_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine ger_impl(m,n,alpha,x,incx,y,incy,a,lda)
        implicit none
        integer(int64), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: x(*)
        type(DT), intent(in) :: y(*)
        type(DT), intent(inout) :: a(lda,*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, j, jy, kx

        zero = 0.0
        if ((m.eq.0) .or. (n.eq.0) .or. (alpha.eq.zero)) return

        if (incy.gt.0) then
            jy = 1
        else
            jy = 1 - (n-1)*incy
        end if

        if (incx.eq.1) then
            do j = 1, n
                if (y(jy).ne.zero) then
                    temp = alpha*y(jy)
                    do i = 1, m
                        a(i,j) = a(i,j) + x(i)*temp
                    end do
                end if
                jy = jy + incy
            end do
        else
            if (incx.gt.0) then
                kx = 1
            else
                kx = 1 - (m-1)*incx
            end if
            do j = 1, n
                if (y(jy).ne.zero) then
                    temp = alpha*y(jy)
                    ix = kx
                    do i = 1, m
                        a(i,j) = a(i,j) + x(ix)*temp
                        ix = ix + incx
                    end do
                end if
                jy = jy + incy
            end do
        end if
    end subroutine

    module subroutine ger_64(m,n,alpha,x,incx,y,incy,a,lda)
        implicit none
        integer(int64), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: a(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (m - 1) * abs(incx))
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

        call ger_impl(m, n, alpha, px, incx, py, incy, pa, lda)
    end subroutine

    module subroutine ger_32(m,n,alpha,x,incx,y,incy,a,lda)
        implicit none
        integer(int32), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: a(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(m, int64) - 1) * abs(int(incx, int64)))
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

        call ger_impl(int(m, int64), int(n, int64), alpha, px, int(incx, int64), py, int(incy, int64), pa, int(lda, int64))
    end subroutine

end submodule
