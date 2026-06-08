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
submodule (lpf_blas_fp8_e5m2) lpf_blas_tpmv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_tpmv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine tpmv_impl(uplo, trans, diag, n, ap, x, incx)
        implicit none
        character, intent(in) :: uplo, trans, diag
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: ap(*)
        type(DT), intent(inout) :: x(*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, j, jx, k, kk, kx
        logical :: nounit

        zero = 0.0
        if (n == 0) return

        nounit = (diag == 'N' .or. diag == 'n')

        if (incx <= 0) then
            kx = 1 - (n-1)*incx
        else if (incx /= 1) then
            kx = 1
        end if

        if (trans == 'N' .or. trans == 'n') then
            if (uplo == 'U' .or. uplo == 'u') then
                kk = 1
                if (incx == 1) then
                    do j = 1, n
                        if (x(j) /= zero) then
                            temp = x(j)
                            k = kk
                            do i = 1, j - 1
                                x(i) = x(i) + temp*ap(k)
                                k = k + 1
                            end do
                            if (nounit) x(j) = x(j)*ap(kk+j-1)
                        end if
                        kk = kk + j
                    end do
                else
                    jx = kx
                    do j = 1, n
                        if (x(jx) /= zero) then
                            temp = x(jx)
                            ix = kx
                            do k = kk, kk + j - 2
                                x(ix) = x(ix) + temp*ap(k)
                                ix = ix + incx
                            end do
                            if (nounit) x(jx) = x(jx)*ap(kk+j-1)
                        end if
                        jx = jx + incx
                        kk = kk + j
                    end do
                end if
            else
                kk = (n * (n + 1)) / 2
                if (incx == 1) then
                    do j = n, 1, -1
                        if (x(j) /= zero) then
                            temp = x(j)
                            k = kk
                            do i = n, j + 1, -1
                                x(i) = x(i) + temp*ap(k)
                                k = k - 1
                            end do
                            if (nounit) x(j) = x(j)*ap(kk-n+j)
                        end if
                        kk = kk - (n-j+1)
                    end do
                else
                    kx = kx + (n-1)*incx
                    jx = kx
                    do j = n, 1, -1
                        if (x(jx) /= zero) then
                            temp = x(jx)
                            ix = kx
                            do k = kk, kk - (n - (j + 1)), -1
                                x(ix) = x(ix) + temp*ap(k)
                                ix = ix - incx
                            end do
                            if (nounit) x(jx) = x(jx)*ap(kk-n+j)
                        end if
                        jx = jx - incx
                        kk = kk - (n-j+1)
                    end do
                end if
            end if
        else
            if (uplo == 'U' .or. uplo == 'u') then
                kk = (n * (n + 1)) / 2
                if (incx == 1) then
                    do j = n, 1, -1
                        temp = x(j)
                        if (nounit) temp = temp*ap(kk)
                        k = kk - 1
                        do i = j - 1, 1, -1
                            temp = temp + ap(k)*x(i)
                            k = k - 1
                        end do
                        x(j) = temp
                        kk = kk - j
                    end do
                else
                    jx = kx + (n-1)*incx
                    do j = n, 1, -1
                        temp = x(jx)
                        ix = jx
                        if (nounit) temp = temp*ap(kk)
                        do k = kk - 1, kk - j + 1, -1
                            ix = ix - incx
                            temp = temp + ap(k)*x(ix)
                        end do
                        x(jx) = temp
                        jx = jx - incx
                        kk = kk - j
                    end do
                end if
            else
                kk = 1
                if (incx == 1) then
                    do j = 1, n
                        temp = x(j)
                        if (nounit) temp = temp*ap(kk)
                        k = kk + 1
                        do i = j + 1, n
                            temp = temp + ap(k)*x(i)
                            k = k + 1
                        end do
                        x(j) = temp
                        kk = kk + (n-j+1)
                    end do
                else
                    jx = kx
                    do j = 1, n
                        temp = x(jx)
                        ix = jx
                        if (nounit) temp = temp*ap(kk)
                        do k = kk + 1, kk + n - j
                            ix = ix + incx
                            temp = temp + ap(k)*x(ix)
                        end do
                        x(jx) = temp
                        jx = jx + incx
                        kk = kk + (n-j+1)
                    end do
                end if
            end if
        end if
    end subroutine

    module subroutine tpmv_64(uplo, trans, diag, n, ap, x, incx)
        implicit none
        character, intent(in) :: uplo, trans, diag
        integer(int64), intent(in) :: incx, n
        type(DT), target, intent(in) :: ap(..)
        type(DT), target, intent(inout) :: x(..)

        type(DT), CONTIGUOUS, pointer :: pap(:)
        type(DT), CONTIGUOUS, pointer :: px(:)
        integer(int64) :: lenx, lenap
        integer(int64) :: lx(1), lap(1)
        type(c_ptr) :: ptr

        lenap = (n * (n + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        call tpmv_impl(uplo, trans, diag, n, pap, px, incx)
    end subroutine

    module subroutine tpmv_32(uplo, trans, diag, n, ap, x, incx)
        implicit none
        character, intent(in) :: uplo, trans, diag
        integer(int32), intent(in) :: incx, n
        type(DT), target, intent(in) :: ap(..)
        type(DT), target, intent(inout) :: x(..)

        type(DT), CONTIGUOUS, pointer :: pap(:)
        type(DT), CONTIGUOUS, pointer :: px(:)
        integer(int64) :: lenx, lenap
        integer(int64) :: lx(1), lap(1)
        type(c_ptr) :: ptr

        lenap = (int(n, int64) * (int(n, int64) + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        call tpmv_impl(uplo, trans, diag, int(n, int64), pap, px, int(incx, int64))
    end subroutine

end submodule
