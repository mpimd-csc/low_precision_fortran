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
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

submodule (lpf_fp8_e5m2) lpf_fp8_e5m2_maxval
    use iso_c_binding
    use iso_fortran_env, only: real32, real64
    use lpf_types
    implicit none

contains
    ! Overall max in 1D
    module pure function maxval_fp8_e5m2_1d(array) result(max_value)
        type(fp8_e5m2), dimension(:), intent(in) :: array
        type(fp8_e5m2) :: max_value
        integer(lpf_default_int_kind) :: i
        type(fp8_e5m2) :: hval;


        if (size(array) .eq. 0) then
            hval = fp8_e5m2(1.0)
            max_value = huge(hval)
            return
        end if

        max_value = array(1)

        do i = 2, size(array)
            if (array(i) .gt. max_value) then
                max_value = array(i)
            end if
        end do
    end function maxval_fp8_e5m2_1d

    ! Overall max in 2D
    module pure function maxval_fp8_e5m2_2d(array) result(max_value)
        type(fp8_e5m2), dimension(:,:), intent(in) :: array
        type(fp8_e5m2) :: max_value
        integer(lpf_default_int_kind) :: i, j

        max_value = array(1, 1)

        do i = 1, size(array, 1)
            do j = 1, size(array, 2)
                if (array(i, j) .gt. max_value) then
                    max_value = array(i, j)
                end if
            end do
        end do
    end function maxval_fp8_e5m2_2d

    ! Overall max in 3D
    module pure function maxval_fp8_e5m2_3d(array) result(max_value)
        type(fp8_e5m2), dimension(:,:,:), intent(in) :: array
        type(fp8_e5m2) :: max_value
        integer(lpf_default_int_kind) :: i, j, k

        max_value = array(1, 1, 1)

        do i = 1, size(array, 1)
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    if (array(i, j, k) .gt. max_value) then
                        max_value = array(i, j, k)
                    end if
                end do
            end do
        end do
    end function maxval_fp8_e5m2_3d

    ! Overall max in 4D
    module pure function maxval_fp8_e5m2_4d(array) result(max_value)
        type(fp8_e5m2), dimension(:,:,:,:), intent(in) :: array
        type(fp8_e5m2) :: max_value
        integer(lpf_default_int_kind) :: i, j, k, l

        max_value = array(1, 1, 1,1)

        do i = 1, size(array, 1)
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    do l = 1, size(array, 4)
                        if (array(i, j, k,l) .gt. max_value) then
                            max_value = array(i, j, k, l)
                        end if
                    end do
                end do
            end do
        end do
    end function maxval_fp8_e5m2_4d


    ! Overall Max in 1D but with dim argument.
    module pure function maxval_fp8_e5m2_1d_dim(array, dim) result(max_value)
        type(fp8_e5m2), dimension(:), intent(in) :: array
        integer, intent(in) :: dim
        type(fp8_e5m2) :: max_value

        if (dim .ne. 1) then
            error stop 'Invalid dimension for 1D array'
        end if

        max_value = maxval_fp8_e5m2_1d(array)
    end function maxval_fp8_e5m2_1d_dim

    ! Max along dim in 2D
    module pure function maxval_fp8_e5m2_2d_dim(array, dim) result(max_value)
        type(fp8_e5m2), dimension(:,:), intent(in) :: array
        integer, intent(in) :: dim
        type(fp8_e5m2), dimension(size(array, merge(2, 1, dim == 1))) :: max_value
        integer(lpf_default_int_kind) :: i, j

        if (dim < 1 .or. dim > 2) then
            error stop 'Invalid dimension for 2D array'
        end if

        ! Initialisierung mit den Elementen des Arrays
        if (dim == 1) then
            do j = 1, size(array, 2)
                max_value(j) = array(1, j)
                do i = 2, size(array, 1)
                    if (array(i, j) .gt. max_value(j)) then
                        max_value(j) = array(i, j)
                    end if
                end do
            end do
        else
            do i = 1, size(array, 1)
                max_value(i) = array(i,1)
                do j = 2, size(array, 2)
                    if (array(i, j) .gt. max_value(i)) then
                        max_value(i) = array(i, j)
                    end if
                end do
            end do
        end if
    end function maxval_fp8_e5m2_2d_dim

    ! Max along dim in 3D
    module pure function maxval_fp8_e5m2_3d_dim(array, dim) result(max_value)
        type(fp8_e5m2), dimension(:,:,:), intent(in) :: array
        integer, intent(in) :: dim
        !
        ! The input is a m x n x k array.
        ! if dim == 1, the output is n x k, we need dimension 2 and 3
        ! if dim == 2, the output is m x k, we need dimension 1 and 3
        ! if dim == 3, the output is m x n, we need dimension 1 and 2
        !
        type(fp8_e5m2), dimension( size(array, merge(2, 1, dim == 1)), &
            & size(array, merge(2, 3, dim == 3))) :: max_value
        integer(lpf_default_int_kind) :: i, j, k

        if (dim < 1 .or. dim > 3) then
            error stop 'Invalid dimension for 3D array'
        end if

        if (dim == 1) then
            ! dim == 1, for each (j,k) iterate over i
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    max_value(j,k) = array(1, j, k)
                    do i = 2, size(array, 1)
                        if (array(i, j, k) .gt. max_value(j,k)) then
                            max_value(j,k) = array(i, j, k)
                        end if
                    end do
                end do
            end do
        elseif (dim == 2) then
            ! dim == 2, for each (i,k) iterate over j
            do i = 1, size(array, 1)
                do k = 1, size(array, 3)
                    max_value(i,k) = array(i, 1, k)
                    do j = 2, size(array, 2)
                        if (array(i, j, k) .gt. max_value(i,k)) then
                            max_value(i,k) = array(i, j, k)
                        end if
                    end do
                end do
            end do
        else
            ! dim == 3, for each (i,j) iterate over k
            do i = 1, size(array, 1)
                do j = 1, size(array, 2)
                    max_value(i,j) = array(i,j,1)
                    do k = 2, size(array, 3)
                        if (array(i, j, k) .gt. max_value(i,j)) then
                            max_value(i,j) = array(i, j, k)
                        end if
                    end do
                end do
            end do
        end if

    end function maxval_fp8_e5m2_3d_dim

    ! Max along dim in 4D
    module pure function maxval_fp8_e5m2_4d_dim(array, dim) result(max_value)
        type(fp8_e5m2), dimension(:,:,:,:), intent(in) :: array
        integer, intent(in) :: dim
        !
        ! The input is a m x n x k x l array.
        ! if dim == 1, the output is n x k x l, we need dimension 2, 3, 4
        ! if dim == 2, the output is m x k x l, we need dimension 1, 3, 4
        ! if dim == 3, the output is m x n x l, we need dimension 1, 2, 4
        ! if dim == 4, the output is m x n x k, we need dimension 1, 2, 3
        !
        type(fp8_e5m2), dimension( size(array, merge(2, 1, dim == 1)), &
            & size(array, merge(3, 2, dim < 3)), &
            & size(array, merge(4, 3, dim == 4))) :: max_value
        integer(lpf_default_int_kind) :: i, j, k, l

        if (dim < 1 .or. dim > 4) then
            error stop 'Invalid dimension for 3D array'
        end if

        if (dim == 1) then
            ! dim == 1, for each (j,k,l) iterate over i
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    do l = 1, size(array, 4)
                        max_value(j,k,l) = array(1, j, k, l)
                        do i = 2, size(array, 1)
                            if (array(i, j, k, l) .gt. max_value(j,k,l)) then
                                max_value(j,k,l) = array(i, j, k,l)
                            end if
                        end do
                    end do
                end do
            end do
        elseif (dim == 2) then
            ! dim == 2, for each (i,k,l) iterate over j
            do i = 1, size(array, 1)
                do k = 1, size(array, 3)
                    do l =1, size(array, 4)
                        max_value(i,k,l) = array(i, 1, k,l)
                        do j = 2, size(array, 2)
                            if (array(i, j, k,l) .gt. max_value(i,k,l)) then
                                max_value(i,k, l) = array(i, j, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        elseif ( dim == 3) then
            ! dim == 3, for each (i,j, l) iterate over k
            do i = 1, size(array, 1)
                do j = 1, size(array, 2)
                    do l = 1, size(array, 4)
                        max_value(i,j, l) = array(i,j,1, l)
                        do k = 2, size(array, 3)
                            if (array(i, j, k, l) .gt. max_value(i,j, l)) then
                                max_value(i,j,l) = array(i, j, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        else
            ! dim == 4, for each (i,j, k) iterate over l
            do i = 1, size(array, 1)
                do j = 1, size(array, 2)
                    do k = 2, size(array, 3)
                        max_value(i,j, k) = array(i,j,k,1)
                        do l = 1, size(array, 4)
                            if (array(i, j, k, l) .gt. max_value(i,j, k)) then
                                max_value(i,j,k) = array(i, j, k, l)
                            end if
                        end do
                    end do
                end do
            end do

        end if

    end function maxval_fp8_e5m2_4d_dim



end submodule
