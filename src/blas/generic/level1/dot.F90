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
submodule (lpf_blas_fp8_e5m2) lpf_blas_dot_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_dot_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function dot_impl(n,sx,incx,sy,incy) result(out)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sx(*)
        type(DT), intent(in) :: sy(*)
        type(DT) :: out

        type(DT) :: stemp
        integer(int64) :: i, ix, iy, m, mp1

        stemp = 0.0
        out = 0.0
        if (n.le.0) return

        if (incx.eq.1 .and. incy.eq.1) then
            m = mod(n, 5)
            if (m.ne.0) then
                do i = 1, m
                    stemp = stemp + sx(i)*sy(i)
                end do
                if (n.lt.5) then
                    out = stemp
                    return
                end if
            end if
            mp1 = m + 1
            do i = mp1, n, 5
                stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) + &
                        sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
            end do
        else
            ix = 1
            iy = 1
            if (incx.lt.0) ix = (-n+1)*incx + 1
            if (incy.lt.0) iy = (-n+1)*incy + 1
            do i = 1, n
                stemp = stemp + sx(ix)*sy(iy)
                ix = ix + incx
                iy = iy + incy
            end do
        end if
        out = stemp
    end function

    module function dot_64(n,sx,incx,sy,incy) result(out)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(in) :: sx(..)
        type(DT), target, intent(in) :: sy(..)
        type(DT) :: out

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

        out = dot_impl(n, psx, incx, psy, incy)
    end function

    module function dot_32(n,sx,incx,sy,incy) result(out)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), target, intent(in) :: sx(..)
        type(DT), target, intent(in) :: sy(..)
        type(DT) :: out

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

        out = dot_impl(int(n, int64), psx, int(incx, int64), psy, int(incy, int64))
    end function

end submodule
