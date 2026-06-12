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

submodule (lpf_blas_fp16) lpf_blas_fp16_scale_diag_right
    use lpf_fp16
    use iso_fortran_env, only: real32, real64, int32, int64
    use lpf_types
    implicit none

contains

    module subroutine scale_diag_right_fp16_64(m, n, a, lda, dr, info)
        integer(int64), intent(in) :: m, n, lda
        integer(int64), intent(inout) :: info
        type(fp16), intent(inout), dimension(lda, *) :: a
        type(fp16), intent(out), dimension(*) :: dr

        integer(int64) :: k, l
        integer(int32) :: iinfo
        type(fp16) :: max_value, value

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
            iinfo = int(info, int32)
            call lpf_blas_xerbla("scale_diag_right / scale_diag_right_fp16", iinfo)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if

        do k = 1, n
            max_value = abs(a(1, k))
            do l = 1, m
                value = abs(a(l, k))
                if (value .gt. max_value) then
                    max_value = value
                end if
            end do
            dr(k) = FP16(1.0)/max_value
            a(1:m,k) = a(1:m,k) * dr(k)
        end do
    end subroutine

    module subroutine scale_diag_right_fp16_32(m, n, a, lda, dr, info)
        integer(int32), intent(in) :: m, n, lda
        integer(int32), intent(inout) :: info
        type(fp16), intent(inout), dimension(lda, *) :: a
        type(fp16), intent(out), dimension(*) :: dr

        integer(int64) :: iinfo

        call scale_diag_right_fp16_64(int(m, int64), int(n, int64), a, int(lda, int64), dr, iinfo)

        info = int(iinfo, int32)

    end subroutine
end submodule

