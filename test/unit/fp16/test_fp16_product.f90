!  SPDX-License-Identifier: LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program; if not, write to the Free Software Foundation,
!  Inc. a copy of the GNU Lesser General Public License.
!
program test_fp16_product
    use lpf_fp16
    use iso_fortran_env, only: real32, real64
    use lpf_fp16_test_utils
    implicit none

    print *, "Running tests for fp16_product..."

    call test_product_1d()
    call test_product_2d()
    call test_product_3d()
    call test_product_2d_dim()
    call test_product_3d_dim()
    call test_product_mask()
    call test_product_edge_cases()

    print *, "All fp16_product tests passed!"

contains

    subroutine test_product_1d()
        implicit none
        type(fp16), dimension(4) :: vec
        type(fp16) :: res
        real(real32), dimension(4) :: ref_vec
        real(real32) :: ref_res

        ref_vec = [1.5, 2.0, 0.5, 2.0]
        vec = ref_vec

        res = product(vec)
        ref_res = product(ref_vec)
        call check_fp16_real64('product_1d', res, real(ref_res, real64) , FP16_TOL)
    end subroutine

    subroutine test_product_2d()
        implicit none
        type(fp16), dimension(2,2) :: vec
        type(fp16) :: res
        real(real32), dimension(2,2) :: ref_vec
        real(real32) :: ref_res

        ref_vec  = reshape([1.5, 2.0, 0.5, 2.0], [2,2])
        vec = ref_vec

        res = product(vec)
        ref_res = product(ref_vec)

        call check_fp16_real64('product_2d', res, real(ref_res, real64) , FP16_TOL)
    end subroutine

    subroutine test_product_3d()
        implicit none
        type(fp16), dimension(2,2,2) :: vec
        type(fp16) :: res
        real(real32), dimension(2,2,2) :: ref_vec
        real(real32) :: ref_res

        ref_vec  = reshape([1.1, 2.0, 0.5, 2.0, 1.5, 1.0, 0.2, 3.0], [2,2,2])
        vec = ref_vec

        res = product(vec)
        ref_res = product(ref_vec)
        call check_fp16_real64('product_3d', res, real(ref_res, real64) , FP16_TOL)
    end subroutine

    subroutine test_product_2d_dim()
        implicit none
        type(fp16), dimension(2,2) :: vec
        type(fp16), dimension(2) :: res_dim1, res_dim2
        real(real32), dimension(2,2) :: ref_vec
        real(real32), dimension(2) :: ref_dim1, ref_dim2
        integer :: i

        ref_vec  = reshape([1.5, 2.0, 0.5, 3.0], [2,2])
        vec = ref_vec

        res_dim1 = product(vec, dim=1)
        ref_dim1 = product(ref_vec, dim=1)
        do i = 1, 2
            call check_fp16_real64('product_2d_dim1_' // integer_to_char(i), res_dim1(i), real(ref_dim1(i), real64) &
                & , FP16_TOL)
        end do

        res_dim2 = product(vec, dim=2)
        ref_dim2 = product(ref_vec, dim=2)
        do i = 1, 2
            call check_fp16_real64('product_2d_dim2_' // integer_to_char(i), res_dim2(i), real(ref_dim2(i), real64) &
                & , FP16_TOL)
        end do
    end subroutine

    subroutine test_product_3d_dim()
        implicit none
        type(fp16), dimension(2,2,2) :: vec
        type(fp16), dimension(2,2) :: res_dim1, res_dim2, res_dim3
        real(real32), dimension(2,2,2) :: ref_vec
        real(real32), dimension(2,2) :: ref_dim1, ref_dim2, ref_dim3
        integer :: i, j

        ref_vec  = reshape([1.1, 2.0, 0.5, 2.0, 1.5, 1.0, 0.2, 3.0], [2,2,2])
        vec = ref_vec

        res_dim1 = product(vec, dim=1)
        ref_dim1 = product(ref_vec, dim=1)
        do j = 1, 2
            do i = 1, 2
                call check_fp16_real64('product_3d_dim1_' // trim(integer_to_char(i)) // '_' &
                    & // integer_to_char(j), res_dim1(i,j), real(ref_dim1(i,j), real64), FP16_TOL)
            end do
        end do

        res_dim2 = product(vec, dim=2)
        ref_dim2 = product(ref_vec, dim=2)
        do j = 1, 2
            do i = 1, 2
                call check_fp16_real64('product_3d_dim2_' // trim(integer_to_char(i)) // '_' // integer_to_char(j), &
                    & res_dim2(i,j), real(ref_dim2(i,j), real64), FP16_TOL)
            end do
        end do

        res_dim3 = product(vec, dim=3)
        ref_dim3 = product(ref_vec, dim=3)
        do j = 1, 2
            do i = 1, 2
                call check_fp16_real64('product_3d_dim3_' // trim(integer_to_char(i)) // '_' &
                    & // integer_to_char(j), res_dim3(i,j), real(ref_dim3(i,j), real64), FP16_TOL)
            end do
        end do
    end subroutine

    subroutine test_product_mask()
        implicit none
        type(fp16), dimension(4) :: vec
        logical, dimension(4) :: mask
        type(fp16) :: res
        real(real32), dimension(4) :: ref_vec
        real(real32) :: ref_res
        integer :: i

        ref_vec  = [1.5, 2.0, 0.5, 2.0]
        mask  = [.TRUE., .FALSE., .TRUE., .FALSE.]
        vec = ref_vec

        res = product(vec, mask=mask)
        ref_res = product(ref_vec, mask=mask)

        call check_fp16_real64('product_mask_1d', res, real(ref_res, real64), FP16_TOL)
    end subroutine

    subroutine test_product_edge_cases()
        implicit none
        type(fp16) :: res
        real(real32) :: ref_res

        ! Single element
        res = product([FP16(1.5)])
        ref_res = product([1.5_real32])
        call check_fp16_real64('product_edge_single', res, real(ref_res, real64), FP16_TOL)

        ! Array with zero
        res = product([FP16(1.5), FP16(0.0), FP16(2.0)])
        ref_res = product([1.5_real32, 0.0_real32, 2.0_real32])
        call check_fp16_real64('product_edge_zero', res, real(ref_res, real64), FP16_TOL)

        ! Array with negative numbers
        res = product([FP16(-1.5), FP16(-2.0), FP16(0.5)])
        ref_res = product([-1.5_real32, -2.0_real32, 0.5_real32])
        call check_fp16_real64('product_edge_neg', res, real(ref_res, real64), FP16_TOL)

        ! Scalar mask - all true
        res = product([FP16(1.5), FP16(2.0)], mask=.TRUE.)
        ref_res = product([1.5_real32, 2.0_real32])
        call check_fp16_real64('product_edge_mask_T', res, real(ref_res, real64), FP16_TOL)

        ! Scalar mask - all false
        res = product([FP16(1.5), FP16(2.0)], mask=.FALSE.)
        ref_res = 1.0_real32
        call check_fp16_real64('product_edge_mask_F', res, real(ref_res, real64), FP16_TOL)
    end subroutine

end program test_fp16_product
