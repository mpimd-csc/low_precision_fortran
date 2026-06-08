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
submodule (lpf_blas_fp8_e5m2) lpf_blas_axpy_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
 submodule (lpf_blas_fp8_e4m3) lpf_blas_axpy_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding


    contains

        pure subroutine axpy_impl(n,sa,sx,incx,sy,incy)
            implicit none

            type(DT) , intent(in) :: sa
            integer(int64), intent(in) :: incx, incy, n
            type(DT), intent(in) :: sx(*)
            type(DT), intent(inout) :: sy(*)

            !     .. local scalars ..
            integer(int64) i,ix,iy
            !     ..
            if (n.le.0) return
            if (sa.eq.0.0) return
            if (incx.eq.1 .and. incy.eq.1) then
                do i = 1,n
                    sy(i) = sy(i) + sa*sx(i)
                end do
            else
                !        code for unequal increments or equal increments
                !          not equal to 1
                !
                ix = 1
                iy = 1
                if (incx.lt.0) ix = (-n+1)*incx + 1
                if (incy.lt.0) iy = (-n+1)*incy + 1
                do i = 1,n
                    sy(iy) = sy(iy) + sa*sx(ix)
                    ix = ix + incx
                    iy = iy + incy
                end do
            end if
            return
        end subroutine

        module subroutine axpy_64(n,sa,sx,incx,sy,incy)
            implicit none

            type(DT) , intent(in) :: sa
            integer(int64), intent(in) :: incx, incy, n
            type(DT), target, intent(in) :: sx(..)
            type(DT), target, intent(inout) :: sy(..)

            type(DT), CONTIGUOUS, pointer, dimension(:) :: psx
            type(DT), CONTIGUOUS, pointer, dimension(:) :: psy

            type(c_ptr) :: ptr
            integer(int64) :: ls(1)
            integer(int64) :: total_size

            ! Flatten the array sx
            total_size = size(sx)
            ls(1) = total_size
            ptr = c_loc(sx)
            call c_f_pointer(ptr, psx, ls)

            total_size = size(sy)
            ls(1) = total_size
            ptr = c_loc(sy)
            call c_f_pointer(ptr, psy, ls)


            call axpy_impl(n, sa, psx, incx, psy, incy)

        end subroutine

        module subroutine axpy_32(n,sa,sx,incx,sy,incy)
            implicit none

            type(DT) , intent(in) :: sa
            integer(int32), intent(in) :: incx, incy, n
            type(DT), target, intent(in) :: sx(..)
            type(DT), target, intent(inout) :: sy(..)

            type(DT), CONTIGUOUS, pointer, DIMENSION(:) :: psx
            type(DT), CONTIGUOUS, pointer, DIMENSION(:) :: psy

            integer(int64) :: cn, cincx, cincy
            integer(int64) :: total_size
            type(c_ptr) :: ptr
            integer(int64) :: ls(1)

            cn = int(n, int64)
            cincx = int(incx, int64)
            cincy = int(incy, int64)

            ! Flatten the array sx
            total_size = size(sx)
            ls(1) = total_size
            ptr = c_loc(sx)
            call c_f_pointer(ptr, psx, ls)

            total_size = size(sy)
            ls(1) = total_size
            ptr = c_loc(sy)
            call c_f_pointer(ptr, psy, ls)


            call axpy_impl(cn, sa, psx, cincx, psy, cincy)

        end subroutine


    end submodule
