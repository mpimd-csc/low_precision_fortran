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
submodule (lpf_blas_fp8_e5m2) lpf_blas_scal_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_scal_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine scal_impl(n,sa,sx,incx)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: sa
        type(DT), intent(inout) :: sx(*)

        integer(int64) :: i, m, mp1, nincx
        type(DT) :: one

        one = 1.0
        if (n.le.0 .or. incx.le.0 .or. sa.eq.one) return

        if (incx.eq.1) then
            m = mod(n, 5)
            if (m.ne.0) then
                do i = 1, m
                    sx(i) = sa*sx(i)
                end do
                if (n.lt.5) return
            end if
            mp1 = m + 1
            do i = mp1, n, 5
                sx(i) = sa*sx(i)
                sx(i+1) = sa*sx(i+1)
                sx(i+2) = sa*sx(i+2)
                sx(i+3) = sa*sx(i+3)
                sx(i+4) = sa*sx(i+4)
            end do
        else
            nincx = n*incx
            do i = 1, nincx, incx
                sx(i) = sa*sx(i)
            end do
        end if
    end subroutine

    module subroutine scal_64(n,sa,sx,incx)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: sa
        type(DT), target, intent(inout) :: sx(..)

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int64) :: len
        integer(int64) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        call scal_impl(n, sa, sxi, incx)
    end subroutine

    module subroutine scal_32(n,sa,sx,incx)
        implicit none
        integer(int32), intent(in) :: incx, n
        type(DT), intent(in) :: sa
        type(DT), target, intent(inout) :: sx(..)

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int32) :: len
        integer(int32) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        call scal_impl(int(n, int64), sa, sxi, int(incx, int64))
    end subroutine

end submodule
