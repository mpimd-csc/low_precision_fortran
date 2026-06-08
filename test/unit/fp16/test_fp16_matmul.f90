!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 3 of the License, or any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with the program, and/or a separate License file
!  associated with the software.

PROGRAM test_fp16_matmul
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
    IMPLICIT NONE

    write(*, '(A)') 'Running fp16 matmul tests...'

    CALL test_matmul_vm()
    CALL test_matmul_mv()
    CALL test_matmul_mm()

    CALL test_summary()

CONTAINS

    subroutine test_matmul_vm()
        IMPLICIT NONE
        write(*, '(A)') 'Testing fp16_matmul_vm...'
        call run_vm_test(4, 4, 'identity', FP16_TOL_TIGHT, 'matmul_vm_identity')
        call run_vm_test(4, 4, 'zero', FP16_TOL_TIGHT, 'matmul_vm_zero')
        call run_vm_test(4, 4, 'det', FP16_TOL_LOOSE, 'matmul_vm_det_small')
        call run_vm_test(16, 16, 'det', FP16_TOL_LOOSE, 'matmul_vm_det_medium')
    end subroutine

    subroutine test_matmul_mv()
        IMPLICIT NONE
        write(*, '(A)') 'Testing fp16_matmul_mv...'
        call run_mv_test(4, 4, 'identity', FP16_TOL_TIGHT, 'matmul_mv_identity')
        call run_mv_test(4, 4, 'zero', FP16_TOL_TIGHT, 'matmul_mv_zero')
        call run_mv_test(4, 4, 'det', FP16_TOL_LOOSE, 'matmul_mv_det_small')
        call run_mv_test(16, 16, 'det', FP16_TOL_LOOSE, 'matmul_mv_det_medium')
    end subroutine

    subroutine test_matmul_mm()
        IMPLICIT NONE
        write(*, '(A)') 'Testing fp16_matmul_mm...'
        call run_mm_test(4, 4, 4, 'identity', FP16_TOL_TIGHT, 'matmul_mm_identity')
        call run_mm_test(4, 4, 4, 'zero', FP16_TOL_TIGHT, 'matmul_mm_zero')
        call run_mm_test(4, 4, 4, 'det', FP16_TOL_LOOSE, 'matmul_mm_det_small')
        call run_mm_test(16, 16, 16, 'det', FP16_TOL_LOOSE, 'matmul_mm_det_medium')
    end subroutine

    subroutine run_vm_test(n, m, case_type, tol, name)
        IMPLICIT NONE
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: case_type, name
        real(real64), intent(in) :: tol
        real(real64), allocatable :: a_r64(:), b_r64(:,:), c_r64(:)
        type(FP16), allocatable :: a_fp16(:), b_fp16(:,:), c_fp16(:)
        integer :: i

        allocate(a_r64(n), b_r64(n, m), c_r64(m), a_fp16(n), b_fp16(n, m), c_fp16(m))

        if (case_type == 'identity') then
            a_r64 = 1.0_real64
            b_r64 = 0.0_real64
            do i = 1, n
                b_r64(i, i) = 1.0_real64
            end do
        else if (case_type == 'zero') then
            a_r64 = 1.0_real64
            b_r64 = 0.0_real64
        else if (case_type == 'det') then
            call fill_deterministic_1d(a_r64)
            call fill_deterministic_2d(b_r64)
        end if

        a_fp16 = FP16(a_r64)
        b_fp16 = FP16(b_r64)
        c_fp16 = matmul(a_fp16, b_fp16)
        c_r64 = matmul(a_r64, b_r64)

        call check_array_fp16_real64_1d(name, c_fp16, c_r64, tol)

        deallocate(a_r64, b_r64, c_r64, a_fp16, b_fp16, c_fp16)
    end subroutine

    subroutine run_mv_test(n, m, case_type, tol, name)
        IMPLICIT NONE
        integer, intent(in) :: n, m
        character(len=*), intent(in) :: case_type, name
        real(real64), intent(in) :: tol
        real(real64), allocatable :: a_r64(:,:), b_r64(:), c_r64(:)
        type(FP16), allocatable :: a_fp16(:,:), b_fp16(:), c_fp16(:)
        integer :: i

        allocate(a_r64(n, m), b_r64(m), c_r64(n), a_fp16(n, m), b_fp16(m), c_fp16(n))

        if (case_type == 'identity') then
            b_r64 = 1.0_real64
            a_r64 = 0.0_real64
            do i = 1, n
                a_r64(i, i) = 1.0_real64
            end do
        else if (case_type == 'zero') then
            b_r64 = 1.0_real64
            a_r64 = 0.0_real64
        else if (case_type == 'det') then
            call fill_deterministic_2d(a_r64)
            call fill_deterministic_1d(b_r64)
        end if

        a_fp16 = FP16(a_r64)
        b_fp16 = FP16(b_r64)
        c_fp16 = matmul(a_fp16, b_fp16)
        c_r64 = matmul(a_r64, b_r64)

        call check_array_fp16_real64_1d(name, c_fp16, c_r64, tol)

        deallocate(a_r64, b_r64, c_r64, a_fp16, b_fp16, c_fp16)
    end subroutine

    subroutine run_mm_test(n, m, k, case_type, tol, name)
        IMPLICIT NONE
        integer, intent(in) :: n, m, k
        character(len=*), intent(in) :: case_type, name
        real(real64), intent(in) :: tol
        real(real64), allocatable :: a_r64(:,:), b_r64(:,:), c_r64(:,:)
        type(FP16), allocatable :: a_fp16(:,:), b_fp16(:,:), c_fp16(:,:)
        integer :: i

        allocate(a_r64(n, k), b_r64(k, m), c_r64(n, m), a_fp16(n, k), b_fp16(k, m), c_fp16(n, m))

        if (case_type == 'identity') then
            a_r64 = 1.0_real64
            b_r64 = 0.0_real64
            do i = 1, k
                b_r64(i, i) = 1.0_real64
            end do
        else if (case_type == 'zero') then
            a_r64 = 1.0_real64
            b_r64 = 0.0_real64
        else if (case_type == 'det') then
            call fill_deterministic_2d(a_r64)
            call fill_deterministic_2d(b_r64)
        end if

        a_fp16 = FP16(a_r64)
        b_fp16 = FP16(b_r64)
        c_fp16 = matmul(a_fp16, b_fp16)
        c_r64 = matmul(a_r64, b_r64)

        call check_array_fp16_real64_2d(name, c_fp16, c_r64, tol)

        deallocate(a_r64, b_r64, c_r64, a_fp16, b_fp16, c_fp16)
    end subroutine

    subroutine check_array_fp16_real64_1d(name, got, expected, tol)
        IMPLICIT NONE
        character(len=*), intent(in) :: name
        type(FP16), intent(in) :: got(:)
        real(real64), intent(in) :: expected(:)
        real(real64), intent(in) :: tol
        integer :: i, n
        character(len=20) :: idx_str

        n = size(got)
        do i = 1, n
            write(idx_str, '(I0)') i
            call check_fp16_real64(trim(name) // '_' // trim(idx_str), got(i), expected(i), tol)
        end do
    end subroutine

    subroutine check_array_fp16_real64_2d(name, got, expected, tol)
        IMPLICIT NONE
        character(len=*), intent(in) :: name
        type(FP16), intent(in) :: got(:,:)
        real(real64), intent(in) :: expected(:,:)
        real(real64), intent(in) :: tol
        integer :: i, j, n, m
        character(len=20) :: idx_i, idx_j

        n = size(got, 1)
        m = size(got, 2)
        do i = 1, n
            write(idx_i, '(I0)') i
            do j = 1, m
                write(idx_j, '(I0)') j
                call check_fp16_real64(trim(name) // '_' // trim(idx_i) // '_' // trim(idx_j), got(i,j), expected(i,j), tol)
            end do
        end do
    end subroutine

    subroutine fill_deterministic_1d(arr)
        IMPLICIT NONE
        real(real64), intent(out) :: arr(:)
        integer :: i, sz

        sz = size(arr)
        do i = 1, sz
            arr(i) = real(i, real64) * 0.1_real64
        end do
    end subroutine

    subroutine fill_deterministic_2d(arr)
        IMPLICIT NONE
        real(real64), intent(out) :: arr(:,:)
        integer :: i, j, n, m

        n = size(arr, 1)
        m = size(arr, 2)
        do i = 1, n
            do j = 1, m
                arr(i,j) = real(i + j * n, real64) * 0.01_real64
            end do
        end do
    end subroutine

END PROGRAM test_fp16_matmul
