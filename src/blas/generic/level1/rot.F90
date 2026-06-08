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
submodule (lpf_blas_fp8_e5m2) lpf_blas_rot_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rot_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rot_impl(n,sx,incx,sy,incy,c,s)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(inout) :: sx(*)
        type(DT), intent(inout) :: sy(*)
        type(DT), intent(in) :: c, s

        type(DT) :: stemp
        integer(int64) :: i, ix, iy

        if (n.le.0) return

        if (incx.eq.1 .and. incy.eq.1) then
            do i = 1, n
                stemp = c*sx(i) + s*sy(i)
                sy(i) = c*sy(i) - s*sx(i)
                sx(i) = stemp
            end do
        else
            ix = 1
            iy = 1
            if (incx.lt.0) ix = (-n+1)*incx + 1
            if (incy.lt.0) iy = (-n+1)*incy + 1
            do i = 1, n
                stemp = c*sx(ix) + s*sy(iy)
                sy(iy) = c*sy(iy) - s*sx(ix)
                sx(ix) = stemp
                ix = ix + incx
                iy = iy + incy
            end do
        end if
    end subroutine

    module subroutine rot_64(n,sx,incx,sy,incy,c,s)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)
        type(DT), intent(in) :: c, s

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        call rot_impl(n, psx, incx, psy, incy, c, s)
    end subroutine

    module subroutine rot_32(n,sx,incx,sy,incy,c,s)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)
        type(DT), intent(in) :: c, s

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        integer(int32) :: lenx, leny
        integer(int32) :: lx(1), ly(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        call rot_impl(int(n, int64), psx, int(incx, int64), psy, int(incy, int64), c, s)
    end subroutine

end submodule
