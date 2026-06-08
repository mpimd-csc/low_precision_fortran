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
submodule (lpf_blas_fp8_e5m2) lpf_blas_sbdot_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_sbdot_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function sbdot_impl(n,sb,sx,incx,sy,incy) result(out)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sb
        type(DT), intent(in) :: sx(*)
        type(DT), intent(in) :: sy(*)
        type(DT) :: out

        real(real32) :: dsdot
        integer(int64) :: i, kx, ky, ns

        dsdot = sb
        if (n.le.0) then
            out = dsdot
            return
        end if

        if (incx.eq.incy .and. incx.gt.0) then
            ns = n*incx
            do i = 1, ns, incx
                dsdot = dsdot + real(sx(i)*sy(i))
            end do
        else
            kx = 1
            ky = 1
            if (incx.lt.0) kx = 1 + (1-n)*incx
            if (incy.lt.0) ky = 1 + (1-n)*incy
            do i = 1, n
                dsdot = dsdot + real(sx(kx)*sy(ky))
                kx = kx + incx
                ky = ky + incy
            end do
        end if
        out = dsdot
    end function

    module function sbdot_64(n,sb,sx,incx,sy,incy) result(out)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sb
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

        out = sbdot_impl(n, sb, psx, incx, psy, incy)
    end function

    module function sbdot_32(n,sb,sx,incx,sy,incy) result(out)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sb
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

        out = sbdot_impl(int(n, int64), sb, psx, int(incx, int64), psy, int(incy, int64))
    end function

end submodule
