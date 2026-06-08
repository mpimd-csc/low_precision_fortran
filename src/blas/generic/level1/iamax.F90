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
submodule (lpf_blas_fp8_e5m2) lpf_blas_samax_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_samax_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function iamax_impl(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: sx(*)
        integer(int64) :: out

        type(DT) :: smax
        integer(int64) :: i, ix

        out = 0
        if (n.lt.1 .or. incx.le.0) return
        out = 1
        if (n.eq.1) return

        if (incx.eq.1) then
            smax = abs(sx(1))
            do i = 2, n
                if (abs(sx(i)).gt.smax) then
                    out = i
                    smax = abs(sx(i))
                end if
            end do
        else
            ix = 1
            smax = abs(sx(1))
            ix = ix + incx
            do i = 2, n
                if (abs(sx(ix)).gt.smax) then
                    out = i
                    smax = abs(sx(ix))
                end if
                ix = ix + incx
            end do
        end if
    end function

    module function iamax_64(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        integer(int64) :: out

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int64) :: len
        integer(int64) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        out = iamax_impl(n, sxi, incx)
    end function

    module function iamax_32(n,sx,incx) result(out)
        implicit none
        integer(int32), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        integer(int32) :: out

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int32) :: len
        integer(int32) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        out = int(iamax_impl(int(n, int64), sxi, int(incx, int64)), int32)
    end function

end submodule
