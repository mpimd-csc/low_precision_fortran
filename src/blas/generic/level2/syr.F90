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
submodule (lpf_blas_fp8_e5m2) lpf_blas_syr_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_syr_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding
    implicit none

    contains

        subroutine syr_impl(uplo,n,alpha,x,incx,a,lda)
            implicit none

            character, intent(in) :: uplo
            integer(int64), intent(in) :: incx, lda, n
            type(DT), intent(in) :: alpha
            type(DT), intent(in) :: x(*)
            type(DT), intent(inout) :: a(lda, *)

            type(DT) :: temp, zero
            integer(int64) :: i, ix, j, jx, kx

            zero = 0.0

            if ((n.eq.0) .or. (alpha.eq.zero)) return

            if (incx.le.0) then
                kx = 1 - (n-1)*incx
            else if (incx.ne.1) then
                kx = 1
            end if

            if (uplo == 'U' .or. uplo == 'u') then
                if (incx.eq.1) then
                    do j = 1,n
                        if (x(j).ne.zero) then
                            temp = alpha*x(j)
                            do i = 1,j
                                a(i,j) = a(i,j) + x(i)*temp
                            end do
                        end if
                    end do
                else
                    jx = kx
                    do j = 1,n
                        if (x(jx).ne.zero) then
                            temp = alpha*x(jx)
                            ix = kx
                            do i = 1,j
                                a(i,j) = a(i,j) + x(ix)*temp
                                ix = ix + incx
                            end do
                        end if
                        jx = jx + incx
                    end do
                end if
            else
                if (incx.eq.1) then
                    do j = 1,n
                        if (x(j).ne.zero) then
                            temp = alpha*x(j)
                            do i = j,n
                                a(i,j) = a(i,j) + x(i)*temp
                            end do
                        end if
                    end do
                else
                    jx = kx
                    do j = 1,n
                        if (x(jx).ne.zero) then
                            temp = alpha*x(jx)
                            ix = jx
                            do i = j,n
                                a(i,j) = a(i,j) + x(ix)*temp
                                ix = ix + incx
                            end do
                        end if
                        jx = jx + incx
                    end do
                end if
            end if
            return
        end subroutine

        module subroutine syr_64(uplo, n, alpha, x, incx, a, lda)
            implicit none
            character, intent(in) :: uplo
            integer(int64), intent(in) :: incx, n, lda
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..)
            type(DT), target, intent(inout) :: a(..)

            type(DT), CONTIGUOUS, pointer :: px(:)
            type(DT), CONTIGUOUS, pointer :: pa(:,:)
            integer(int64) :: lenx
            integer(int64) :: lx(1)
            type(c_ptr) :: ptr

            lenx = 1 + ( (n - 1) * abs(incx))
            lx(1) = lenx
            ptr = c_loc(x)
            call c_f_pointer(ptr, px, lx)

            ptr = c_loc(a)
            call c_f_pointer(ptr, pa,  [lda, n])

            call syr_impl(uplo, n, alpha, px, incx, pa, lda)
        end subroutine

        module subroutine syr_32(uplo, n, alpha, x, incx, a, lda)
            implicit none
            character, intent(in) :: uplo
            integer(int32), intent(in) :: incx, n, lda
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..)
            type(DT), target, intent(inout) :: a(..)

            type(DT), CONTIGUOUS, pointer :: px(:)
            type(DT), CONTIGUOUS, pointer :: pap(:,:)
            integer(int64) :: lenx
            integer(int64) :: lx(1)
            type(c_ptr) :: ptr

            lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
            lx(1) = lenx
            ptr = c_loc(x)
            call c_f_pointer(ptr, px, lx)

            ptr = c_loc(a)
            call c_f_pointer(ptr, pap, [lda, n])

            call syr_impl(uplo, int(n, int64), alpha, px, int(incx, int64), pap, int(lda, int64))
        end subroutine

    end submodule

