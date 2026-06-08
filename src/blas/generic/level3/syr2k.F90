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
submodule (lpf_blas_fp8_e5m2) lpf_blas_syr2k_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_syr2k_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine syr2k_impl(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        implicit none
        character, intent(in) :: trans, uplo
        integer(int64), intent(in) :: k, lda, ldb, ldc, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(in) :: b(ldb,*)
        type(DT), intent(inout) :: c(ldc,*)

        type(DT) :: temp1, temp2, zero
        integer(int64) :: i, j, l, nrowa
        logical :: upper

        zero = 0.0
        upper = (uplo == 'U' .or. uplo == 'u')

        if (trans == 'N' .or. trans == 'n') then
            nrowa = n
        else
            nrowa = k
        end if

        if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == 1.0))) return

        if (alpha == zero) then
            if (upper) then
                if (beta == zero) then
                    do j = 1, n
                        do i = 1, j
                            c(i,j) = zero
                        end do
                    end do
                else
                    do j = 1, n
                        do i = 1, j
                            c(i,j) = beta*c(i,j)
                        end do
                    end do
                end if
            else
                if (beta == zero) then
                    do j = 1, n
                        do i = j, n
                            c(i,j) = zero
                        end do
                    end do
                else
                    do j = 1, n
                        do i = j, n
                            c(i,j) = beta*c(i,j)
                        end do
                    end do
                end if
            end if
            return
        end if

        if (trans == 'N' .or. trans == 'n') then
            if (upper) then
                do j = 1, n
                    if (beta == zero) then
                        do i = 1, j
                            c(i,j) = zero
                        end do
                    else if (beta /= 1.0) then
                        do i = 1, j
                            c(i,j) = beta*c(i,j)
                        end do
                    end if
                    do l = 1, k
                        if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                            temp1 = alpha*b(j,l)
                            temp2 = alpha*a(j,l)
                            do i = 1, j
                                c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                            end do
                        end if
                    end do
                end do
            else
                do j = 1, n
                    if (beta == zero) then
                        do i = j, n
                            c(i,j) = zero
                        end do
                    else if (beta /= 1.0) then
                        do i = j, n
                            c(i,j) = beta*c(i,j)
                        end do
                    end if
                    do l = 1, k
                        if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                            temp1 = alpha*b(j,l)
                            temp2 = alpha*a(j,l)
                            do i = j, n
                                c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                            end do
                        end if
                    end do
                end do
            end if
        else
            if (upper) then
                do j = 1, n
                    do i = 1, j
                        temp1 = zero
                        temp2 = zero
                        do l = 1, k
                            temp1 = temp1 + a(l,i)*b(l,j)
                            temp2 = temp2 + b(l,i)*a(l,j)
                        end do
                        if (beta == zero) then
                            c(i,j) = alpha*temp1 + alpha*temp2
                        else
                            c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                        end if
                    end do
                end do
            else
                do j = 1, n
                    do i = j, n
                        temp1 = zero
                        temp2 = zero
                        do l = 1, k
                            temp1 = temp1 + a(l,i)*b(l,j)
                            temp2 = temp2 + b(l,i)*a(l,j)
                        end do
                        if (beta == zero) then
                            c(i,j) = alpha*temp1 + alpha*temp2
                        else
                            c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                        end if
                    end do
                end do
            end if
        end if
    end subroutine

    module subroutine syr2k_64(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        implicit none
        character, intent(in) :: trans, uplo
        integer(int64), intent(in) :: k, lda, ldb, ldc, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: b(..)
        type(DT), target, intent(inout) :: c(..)

        type(DT), CONTIGUOUS, pointer :: pa(:)
        type(DT), CONTIGUOUS, pointer :: pb(:)
        type(DT), CONTIGUOUS, pointer :: pc(:)
        integer(int64) :: la(1), lb(1), lc(1)
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

        total_size = size(c)
        lc(1) = total_size
        ptr = c_loc(c)
        call c_f_pointer(ptr, pc, lc)

        call syr2k_impl(uplo, trans, n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc)
    end subroutine

    module subroutine syr2k_32(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        implicit none
        character, intent(in) :: trans, uplo
        integer(int32), intent(in) :: k, lda, ldb, ldc, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: b(..)
        type(DT), target, intent(inout) :: c(..)

        type(DT), CONTIGUOUS, pointer :: pa(:)
        type(DT), CONTIGUOUS, pointer :: pb(:)
        type(DT), CONTIGUOUS, pointer :: pc(:)
        integer(int64) :: la(1), lb(1), lc(1)
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

        total_size = size(c)
        lc(1) = total_size
        ptr = c_loc(c)
        call c_f_pointer(ptr, pc, lc)

        call syr2k_impl(uplo, trans, int(n, int64), int(k, int64), &
            & alpha, pa, int(lda, int64), pb, int(ldb, int64), beta, pc, int(ldc, int64))
    end subroutine

end submodule
