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

module lpf_blas_scale_diag_right
    use iso_fortran_env, only: real32, real64
    use lpf_types
    implicit none

    private

    public :: scale_diag_right

    interface scale_diag_right
        procedure :: scale_diag_right_fp32
        procedure :: scale_diag_right_fp64
    end interface

contains

    subroutine scale_diag_right_fp32(m, n, a, lda, dr, info)
        use lpf_xerbla
        integer(lpf_default_int_kind), intent(in) :: m, n, lda
        integer(lpf_default_int_kind), intent(inout) :: info
        real(real32), intent(inout), dimension(lda, *) :: a
        real(real32), intent(out), dimension(*) :: dr

        integer(lpf_default_int_kind) :: k

        external :: lpf_blas_xerbla

        ! Check arguments
        info = 0
        if ( m .lt. 0 ) then
            info = -1
        else if ( n .lt. 0 ) then
            info = -2
        else if ( lda .lt. max(1, m)) then
            info = -4
        end if

        if ( info .ne. 0) then
            call lpf_blas_xerbla("scale_diag_right / scale_diag_right_fp32", info)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if

        do k = 1, n
            dr(k) = 1.0/maxval(abs(a(1:m,k)))
            a(1:m,k) = a(1:m,k) * dr(k)
        end do

    end subroutine

    subroutine scale_diag_right_fp64(m, n, a, lda, dr, info)
        integer(lpf_default_int_kind), intent(in) :: m, n, lda
        integer(lpf_default_int_kind), intent(inout) :: info
        real(real64), intent(inout), dimension(lda, *) :: a
        real(real64), intent(out), dimension(*) :: dr

        integer(lpf_default_int_kind) :: k

        external :: lpf_blas_xerbla

        ! Check arguments
        info = 0
        if ( m .lt. 0 ) then
            info = -1
        else if ( n .lt. 0 ) then
            info = -2
        else if ( lda .lt. max(1, m)) then
            info = -4
        end if

        if ( info .ne. 0) then
            call lpf_blas_xerbla("scale_diag_right / scale_diag_right_fp64", info)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if

        do k = 1, n
            dr(k) = 1.0D0/maxval(abs(a(1:m,k)))
            a(1:m,k) = a(1:m,k) * dr(k)
        end do


    end subroutine

end module
