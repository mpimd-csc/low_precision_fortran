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
submodule (lpf_blas_fp8_e5m2) lpf_blas_trmm_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_trmm_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine trmm_impl(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        implicit none
        character, intent(in) :: side, uplo, transa, diag
        integer(int64), intent(in) :: m, n, lda, ldb
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(inout) :: b(ldb,*)

        type(DT) :: temp, zero, one
        integer(int64) :: i, j, k, nrowa
        logical :: lside, nounit, upper

        zero = 0.0
        one = 1.0
        lside = (side == 'L' .or. side == 'l')
        nounit = (diag == 'N' .or. diag == 'n')
        upper = (uplo == 'U' .or. uplo == 'u')

        if (lside) then
            nrowa = m
        else
            nrowa = n
        end if

        if ((m == 0) .or. (n == 0)) return

        if (alpha == zero) then
            do j = 1, n
                do i = 1, m
                    b(i,j) = zero
                end do
            end do
            return
        end if

        if (lside) then
            if (transa == 'N' .or. transa == 'n') then
                if (upper) then
                    do j = 1, n
                        do k = 1, m
                            if (b(k,j) /= zero) then
                                temp = alpha*b(k,j)
                                do i = 1, k - 1
                                    b(i,j) = b(i,j) + temp*a(i,k)
                                end do
                                if (nounit) temp = temp*a(k,k)
                                b(k,j) = temp
                            end if
                        end do
                    end do
                else
                    do j = 1, n
                        do k = m, 1, -1
                            if (b(k,j) /= zero) then
                                temp = alpha*b(k,j)
                                b(k,j) = temp
                                if (nounit) b(k,j) = b(k,j)*a(k,k)
                                do i = k + 1, m
                                    b(i,j) = b(i,j) + temp*a(i,k)
                                end do
                            end if
                        end do
                    end do
                end if
            else
                if (upper) then
                    do j = 1, n
                        do i = m, 1, -1
                            temp = b(i,j)
                            if (nounit) temp = temp*a(i,i)
                            do k = 1, i - 1
                                temp = temp + a(k,i)*b(k,j)
                            end do
                            b(i,j) = alpha*temp
                        end do
                    end do
                else
                    do j = 1, n
                        do i = 1, m
                            temp = b(i,j)
                            if (nounit) temp = temp*a(i,i)
                            do k = i + 1, m
                                temp = temp + a(k,i)*b(k,j)
                            end do
                            b(i,j) = alpha*temp
                        end do
                    end do
                end if
            end if
        else
            if (transa == 'N' .or. transa == 'n') then
                if (upper) then
                    do j = n, 1, -1
                        temp = alpha
                        if (nounit) temp = temp*a(j,j)
                        do i = 1, m
                            b(i,j) = temp*b(i,j)
                        end do
                        do k = 1, j - 1
                            if (a(k,j) /= zero) then
                                temp = alpha*a(k,j)
                                do i = 1, m
                                    b(i,j) = b(i,j) + temp*b(i,k)
                                end do
                            end if
                        end do
                    end do
                else
                    do j = 1, n
                        temp = alpha
                        if (nounit) temp = temp*a(j,j)
                        do i = 1, m
                            b(i,j) = temp*b(i,j)
                        end do
                        do k = j + 1, n
                            if (a(k,j) /= zero) then
                                temp = alpha*a(k,j)
                                do i = 1, m
                                    b(i,j) = b(i,j) + temp*b(i,k)
                                end do
                            end if
                        end do
                    end do
                end if
            else
                if (upper) then
                    do k = 1, n
                        do j = 1, k - 1
                            if (a(j,k) /= zero) then
                                temp = alpha*a(j,k)
                                do i = 1, m
                                    b(i,j) = b(i,j) + temp*b(i,k)
                                end do
                            end if
                        end do
                        temp = alpha
                        if (nounit) temp = temp*a(k,k)
                        if (temp /= one) then
                            do i = 1, m
                                b(i,k) = temp*b(i,k)
                            end do
                        end if
                    end do
                else
                    do k = n, 1, -1
                        do j = k + 1, n
                            if (a(j,k) /= zero) then
                                temp = alpha*a(j,k)
                                do i = 1, m
                                    b(i,j) = b(i,j) + temp*b(i,k)
                                end do
                            end if
                        end do
                        temp = alpha
                        if (nounit) temp = temp*a(k,k)
                        if (temp /= one) then
                            do i = 1, m
                                b(i,k) = temp*b(i,k)
                            end do
                        end if
                    end do
                end if
            end if
        end if
    end subroutine

    module subroutine trmm_64(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        implicit none
        character, intent(in) :: side, uplo, transa, diag
        integer(int64), intent(in) :: m, n, lda, ldb
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(inout) :: b(..)

        type(DT), CONTIGUOUS, pointer :: pa(:)
        type(DT), CONTIGUOUS, pointer :: pb(:)
        integer(int64) :: la(1), lb(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        total_size = size(b)
        lb(1) = total_size
        ptr = c_loc(b)
        call c_f_pointer(ptr, pb, lb)

        call trmm_impl(side, uplo, transa, diag, m, n, alpha, pa, lda, pb, ldb)
    end subroutine

    module subroutine trmm_32(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        implicit none
        character, intent(in) :: side, uplo, transa, diag
        integer(int32), intent(in) :: m, n, lda, ldb
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(inout) :: b(..)

        type(DT), CONTIGUOUS, pointer :: pa(:)
        type(DT), CONTIGUOUS, pointer :: pb(:)
        integer(int64) :: la(1), lb(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        total_size = size(b)
        lb(1) = total_size
        ptr = c_loc(b)
        call c_f_pointer(ptr, pb, lb)

        call trmm_impl(side, uplo, transa, diag, int(m, int64), int(n, int64), alpha, pa, int(lda, int64), pb, int(ldb, int64))
    end subroutine

end submodule
