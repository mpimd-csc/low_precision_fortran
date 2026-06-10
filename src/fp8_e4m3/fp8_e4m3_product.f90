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

submodule (lpf_fp8_e4m3) lpf_fp8_e4m3_product
    use iso_fortran_env
    implicit none
contains
        module pure function fp8_e4m3_product_r1(vec, mask) result(out)
            type(fp8_e4m3), intent(in) :: vec(:)
            logical, intent(in), optional :: mask(..)
            type(fp8_e4m3) :: out

            integer(int64) :: j

            out = 1.0

            if ( present ( mask )) then
                select rank (mask)
                    rank (0)
                        if (mask) then
                            out = fp8_e4m3_product_r1(vec)
                        else
                            out = 1.0
                        end if
                    rank (1)
                        if ( size(vec, 1) .ne. size(mask,1)) then
                            error stop "product: size(vec, 1) <> size(mask,1)"
                        end if
                        do j = 1, size(vec, 1)
                            if ( mask(j) ) then
                                out = out * vec(j)
                            end if
                        end do
                    rank default
                        error stop "product: rank(mask) <> rank(vec)"
                end select

            else
                do j = 1, size(vec)
                    out = out * vec(j)
                end do
            end if

        end function

        module pure function fp8_e4m3_product_r2(vec, mask) result(out)
            type(fp8_e4m3), intent(in) :: vec(:,:)
            logical, intent(in), optional :: mask(..)
            type(fp8_e4m3) :: out
            integer(int64) :: j, k

            out = 1.0

            if ( present ( mask )) then
                select rank (mask)
                    rank (0)
                        if (mask) then
                            out = fp8_e4m3_product_r2(vec)
                        else
                            out = 1.0
                        end if
                    rank (2)
                        if ( size(vec, 1) .ne. size(mask,1) .or. size(vec,2) .ne. size(mask, 2)) then
                            error stop "product: size(vec, k) <> size(mask,k), k = 1 or 2"
                        end if
                        do k = 1, size(vec, 2)
                            do j = 1, size(vec, 1)
                                if ( mask(j,k) ) then
                                  out = out * vec(j,k)
                                end if
                            end do
                        end do
                    rank default
                        error stop "product: rank(mask) <> rank(vec)"
                end select

            else
                do k = 1, size(vec,2)
                    do j = 1, size(vec,1)
                        out = out * vec(j,k)
                    end do
                end do
            end if
        end function

        module pure function fp8_e4m3_product_r3(vec, mask) result(out)
            type(fp8_e4m3), intent(in) :: vec(:,:,:)
            logical, intent(in), optional :: mask(..)
            type(fp8_e4m3) :: out

            integer(int64) :: j, k, l

            out = 1.0

            if ( present ( mask )) then
                select rank (mask)
                    rank (0)
                        if (mask) then
                            out = fp8_e4m3_product_r3(vec)
                        else
                            out = 1.0
                        end if
                    rank (3)
                        if ( size(vec, 1) .ne. size(mask,1) .or. size(vec,2) .ne. size(mask, 2) &
                            & .or. size(vec, 3) .ne. size(mask,3)) then
                            error stop "product: size(vec, k) <> size(mask,k), k = 1 or 2 or 3"
                        end if
                        do l = 1, size(vec, 3)
                            do k = 1, size(vec, 2)
                                do j = 1, size(vec, 1)
                                    if ( mask(j,k,l) ) then
                                        out = out * vec(j,k,l )
                                    end if
                                end do
                            end do
                        end do
                    rank default
                        error stop "product: rank(mask) <> rank(vec)"
                end select

            else
                do l = 1, size(vec,3)
                    do k = 1, size(vec,2)
                        do j = 1, size(vec,1)
                            out = out * vec(j,k,l)
                        end do
                    end do
                end do
            end if

        end function

        module pure function fp8_e4m3_product_r1_dim(vec, dim, mask) result(out)
            type(fp8_e4m3), intent(in) :: vec(:)
            integer, intent(in) :: dim
            logical, intent(in), optional :: mask(..)
            type(fp8_e4m3) :: out

            if ( dim .ne. 1 ) then
                error stop "product: invalid dim"
            end if

            if ( present (mask)) then
                if ( size(vec, 1) .ne. size(mask,1)) then
                    error stop "product: size(vec, 1) <> size(mask,1)"
                end if

                out = fp8_e4m3_product_r1(vec, mask)
            else
                out = fp8_e4m3_product_r1(vec)
            end if

        end function

        module pure function fp8_e4m3_product_r2_dim(vec, dim, mask) result(out)
            type(fp8_e4m3), intent(in) :: vec(:, :)
            integer, intent(in) :: dim
            logical, intent(in), optional :: mask(..)
            type(fp8_e4m3), dimension(size(vec,merge(2, 1, dim == 1))) :: out

            integer(int64) :: j
            out = 1.0

            if ( present (mask) ) then
                select rank(mask)
                    rank(0)
                        if ( mask ) then
                            out = fp8_e4m3_product_r2_dim(vec, dim)
                        else
                            out = 1.0
                        end if
                    rank(2)
                        if ( size(vec, 1) .ne. size(mask,1) .or. size(vec,2) .ne. size(mask, 2)) then
                            error stop "product: size(vec, k) <> size(mask,k), k = 1 or 2"
                        end if

                        if ( dim .eq. 1 ) then
                            do j = 1, size(vec,2)
                                out(j) = fp8_e4m3_product_r1(vec(:,j), mask(:, j))
                            end do
                        else if ( dim .eq. 2 ) then
                            do j = 1, size(vec,1)
                                out(j) = fp8_e4m3_product_r1(vec(j,:), mask(j, :))
                            end do
                        end if
                    rank default
                        error stop "product: rank(vec) <> rank(mask)"
                end select
            else
                if ( dim .eq. 1 ) then
                    do j = 1, size(vec,2)
                        out(j) = fp8_e4m3_product_r1(vec(:,j))
                    end do
                else if ( dim .eq. 2 ) then
                    do j = 1, size(vec,1)
                        out(j) = fp8_e4m3_product_r1(vec(j,:))
                    end do
                end if
            end if
        end function

        module pure function fp8_e4m3_product_r3_dim(vec, dim, mask) result(out)
            type(fp8_e4m3), intent(in) :: vec(:, :, :)
            integer, intent(in) :: dim
            logical, intent(in), optional :: mask(..)
            type(fp8_e4m3), dimension( size(vec, merge(2, 1, dim == 1)), &
                & size(vec, merge(2, 3, dim == 3))) :: out

            integer(int64) :: i, j
            out = 1.0

            if ( present (mask) ) then
                select rank(mask)
                    rank(0)
                        if ( mask ) then
                            out = fp8_e4m3_product_r3_dim(vec, dim)
                        else
                            out = 1.0
                        end if
                     rank(3)
                        if ( size(vec, 1) .ne. size(mask,1) .or. size(vec,2) .ne. size(mask, 2) &
                            & .or. size(vec, 3) .ne. size(mask,3)) then
                            error stop "product: size(vec, k) <> size(mask,k), k = 1 or 2 or 3"
                        end if

                         if ( dim .eq. 1 ) then
                            do i = 1, size(vec,2)
                                do j = 1, size(vec,3)
                                    out(i, j) = fp8_e4m3_product_r1(vec(:, i , j), mask(:,i,j))
                                end do
                            end do
                        else if ( dim .eq. 2 ) then
                            do i = 1, size(vec,1)
                                do j = 1, size(vec,3)
                                    out(i, j) = fp8_e4m3_product_r1(vec(i, : , j), mask(i,:,j))
                                end do
                            end do
                        else if ( dim .eq. 3 ) then
                            do i = 1, size(vec,1)
                                do j = 1, size(vec,2)
                                    out(i, j) = fp8_e4m3_product_r1(vec(i , j, :), mask(i,j,:))
                                end do
                            end do
                        end if
                    rank default
                        error stop "product: rank(vec) <> rank(mask)"
                end select
            else
                if ( dim .eq. 1 ) then
                    do i = 1, size(vec,2)
                        do j = 1, size(vec,3)
                            out(i, j) = fp8_e4m3_product_r1(vec(:, i , j))
                        end do
                    end do
                else if ( dim .eq. 2 ) then
                    do i = 1, size(vec,1)
                        do j = 1, size(vec,3)
                            out(i, j) = fp8_e4m3_product_r1(vec(i, : , j))
                        end do
                    end do
                else if ( dim .eq. 3 ) then
                    do i = 1, size(vec,1)
                        do j = 1, size(vec,2)
                            out(i, j) = fp8_e4m3_product_r1(vec(i , j, :))
                        end do
                    end do
                end if
            end if

        end function

end submodule

