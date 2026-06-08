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
submodule (lpf_blas_fp8_e5m2) lpf_blas_rotm_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rotm_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rotm_impl(n,sx,incx,sy,incy,sparams)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(*)
        type(DT), target, intent(inout) :: sy(*)
        type(DT), intent(in) :: sparams(*)

        type(DT) :: sflag, sh11, sh12, sh21, sh22, two, w, z, zero
        integer(int64) :: i, kx, ky, nsteps

        zero = 0.0
        two = 2.0
        sflag = sparams(1)
        if (n.le.0 .or. (sflag + two .eq. zero)) return

        if (incx.eq.incy .and. incx.gt.0) then
            nsteps = n*incx
            if (sflag.lt.zero) then
                sh11 = sparams(2)
                sh12 = sparams(4)
                sh21 = sparams(3)
                sh22 = sparams(5)
                do i = 1, nsteps, incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z*sh12
                    sy(i) = w*sh21 + z*sh22
                end do
            else if (sflag.eq.zero) then
                sh12 = sparams(4)
                sh21 = sparams(3)
                do i = 1, nsteps, incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w + z*sh12
                    sy(i) = w*sh21 + z
                end do
            else
                sh11 = sparams(2)
                sh22 = sparams(5)
                do i = 1, nsteps, incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z
                    sy(i) = -w + sh22*z
                end do
            end if
        else
            kx = 1
            ky = 1
            if (incx.lt.0) kx = 1 + (1-n)*incx
            if (incy.lt.0) ky = 1 + (1-n)*incy

            if (sflag.lt.zero) then
                sh11 = sparams(2)
                sh12 = sparams(4)
                sh21 = sparams(3)
                sh22 = sparams(5)
                do i = 1, n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w*sh11 + z*sh12
                    sy(ky) = w*sh21 + z*sh22
                    kx = kx + incx
                    ky = ky + incy
                end do
            else if (sflag.eq.zero) then
                sh12 = sparams(4)
                sh21 = sparams(3)
                do i = 1, n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w + z*sh12
                    sy(ky) = w*sh21 + z
                    kx = kx + incx
                    ky = ky + incy
                end do
            else
                sh11 = sparams(2)
                sh22 = sparams(5)
                do i = 1, n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w*sh11 + z
                    sy(ky) = -w + sh22*z
                    kx = kx + incx
                    ky = ky + incy
                end do
            end if
        end if
    end subroutine

    module subroutine rotm_64(n,sx,incx,sy,incy,sparam)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)
        type(DT), target, intent(in) :: sparam(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        type(DT), CONTIGUOUS, pointer :: psparams(:)
        integer(int64) :: lenx, leny, lenp
        integer(int64) :: lx(1), ly(1), lp(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        lenp = 5
        lp(1) = lenp
        ptr = c_loc(sparam)
        call c_f_pointer(ptr, psparams, lp)

        call rotm_impl(n, psx, incx, psy, incy, psparams)
    end subroutine

    module subroutine rotm_32(n,sx,incx,sy,incy,sparam)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)
        type(DT), target, intent(in) :: sparam(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        type(DT), CONTIGUOUS, pointer :: psparams(:)
        integer(int32) :: lenx, leny, lenp
        integer(int32) :: lx(1), ly(1), lp(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        lenp = 5
        lp(1) = lenp
        ptr = c_loc(sparam)
        call c_f_pointer(ptr, psparams, lp)

        call rotm_impl(int(n, int64), psx, int(incx, int64), psy, int(incy, int64), psparams)
    end subroutine

end submodule
