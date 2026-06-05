!  SPDX-License-Identifier LGPL-3.0-or-later
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
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
!
submodule (lpf_fp8_e4m3) lpf_fp8_e4m3_math_norm2
    use iso_fortran_env
    implicit none

contains
    module pure function fp8_e4m3_norm2_r1(vec) result(out)
        type(fp8_e4m3), intent(in) :: vec(:)
        type(fp8_e4m3) :: out

        integer(int64) :: i, n
        type(fp8_e4m3) :: scales, ssq, absxi, one, zero

        out = 0.0

        one  = 1.0
        zero = 0.0

        n = size(vec)

        if (n .lt. 1) then
            out = 0.0
            return
        else if ( n .eq. 1 ) then
            out = abs(vec(1))
        else
            scales = 0.0
            ssq = 1.0

            do i = 1, n
                if ( vec(i).ne.zero) then
                    absxi = abs(vec(i))
                    if (scales.lt.absxi) then
                        ssq = one + ssq* (scales/absxi)**2.0
                        scales = absxi
                    else
                        ssq = ssq + (absxi/scales)**2.0
                    end if
                end if
            end do
            out = scales*sqrt(ssq)
        end if
        return
    end function

    module pure function fp8_e4m3_norm2_r2(vec) result(out)
        type(fp8_e4m3), intent(in) :: vec(:, :)
        type(fp8_e4m3) :: out

        integer(int64) :: i, j, n, m
        type(fp8_e4m3) :: scales, ssq, absxi, one, zero

        out = 0.0

        one  = 1.0
        zero = 0.0

        n = size(vec, 1)
        m = size(vec, 2)

        if (n .lt. 1 .or. m .lt.1 ) then
            out = 0.0
            return
        else if ( n .eq. 1 .and. m .eq. 1 ) then
            out = abs(vec(1 , 1 ))
        else
            scales = 0.0
            ssq = 1.0

            do j = 1, m
                do i = 1, n
                    if ( vec(i, j ).ne.zero) then
                        absxi = abs(vec(i, j))
                        if (scales.lt.absxi) then
                            ssq = one + ssq* (scales/absxi)**2.0
                            scales = absxi
                        else
                            ssq = ssq + (absxi/scales)**2.0
                        end if
                    end if
                end do
            end do
            out = scales*sqrt(ssq)
        end if
        return
    end function

    module pure function fp8_e4m3_norm2_r3(vec) result(out)
        type(fp8_e4m3), intent(in) :: vec(:, :, :)
        type(fp8_e4m3) :: out

        integer(int64) :: i, j, k, n, m, p
        type(fp8_e4m3) :: scales, ssq, absxi, one, zero

        out = 0.0

        one  = 1.0
        zero = 0.0

        n = size(vec, 1)
        m = size(vec, 2)
        p = size(vec, 3)

        if (n .lt. 1 .or. m .lt. 1 .or. p .lt. 1 ) then
            out = 0.0
            return
        else if ( n .eq. 1 .and. m .eq. 1 .and. p .eq. 1) then
            out = abs(vec(1,1,1))
        else
            scales = 0.0
            ssq = 1.0

            do k  = 1, p
                do j = 1, m
                    do i = 1, n
                        if ( vec(i, j , k ).ne.zero) then
                            absxi = abs(vec(i, j, k ))
                            if (scales.lt.absxi) then
                                ssq = one + ssq* (scales/absxi)**2.0
                                scales = absxi
                            else
                                ssq = ssq + (absxi/scales)**2.0
                            end if
                        end if
                    end do
                end do
            end do
            out = scales*sqrt(ssq)
        end if
        return
    end function

    module pure function fp8_e4m3_norm2_r1_dim(vec, dim) result(out)
        type(fp8_e4m3), intent(in) :: vec(:)
        integer, intent(in) :: dim
        type(fp8_e4m3) :: out

        out = 0.0
        if ( dim .eq. 1 ) then
            out = fp8_e4m3_norm2_r1(vec)
        end if
        return
    end function

    module pure function fp8_e4m3_norm2_r2_dim(vec, dim) result(out)
        type(fp8_e4m3), intent(in) :: vec(:, :)
        integer, intent(in) :: dim
        type(fp8_e4m3), dimension(size(vec,merge(2, 1, dim == 1))) :: out

        integer(int64) :: j, m, n

        out = 0.0

        n = size(vec, 1)
        m = size(vec, 2)

        if ( dim .eq. 1 ) then
            do j = 1, m
                out(j) = fp8_e4m3_norm2_r1(vec(:,j))
            end do
        else if ( dim .eq. 2 ) then
            do j = 1, n
                out(j) = fp8_e4m3_norm2_r1(vec(j,:))
            end do
        end if
        return
    end function

    module pure function fp8_e4m3_norm2_r3_dim(vec, dim) result(out)
        type(fp8_e4m3), intent(in) :: vec(:, :, :)
        integer, intent(in) :: dim
        type(fp8_e4m3), dimension( size(vec, merge(2, 1, dim == 1)), &
            & size(vec, merge(2, 3, dim == 3))) :: out

        integer(int64) :: i, j, m, n, p

        out = 0.0

        n = size(vec, 1)
        m = size(vec, 2)
        p = size(vec, 3)

        if ( dim .eq. 1 ) then
            do i = 1, m
                do j = 1, p
                    out(i, j) = fp8_e4m3_norm2_r1(vec(:, i , j))
                end do
            end do
        else if ( dim .eq. 2 ) then
            do i = 1, n
                do j = 1, p
                    out(i, j) = fp8_e4m3_norm2_r1(vec(i, : , j))
                end do
            end do
        else if ( dim .eq. 3 ) then
            do i = 1, m
                do j = 1, n
                    out(i, j) = fp8_e4m3_norm2_r1(vec(i , j, :))
                end do
            end do
        end if

    end function


end submodule
