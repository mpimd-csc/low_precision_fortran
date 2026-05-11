!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 array reduction functions.
!
!  Tests: maxval, minval, maxloc, minloc for 1D-4D arrays

PROGRAM test_bf16_reductions
    USE iso_fortran_env, only: real32, real64, int32
    USE LPF_BF16
    USE lpf_bf16_test_utils
    IMPLICIT NONE

    CALL test_maxval_1d()
    CALL test_minval_1d()
    CALL test_maxval_2d()
    CALL test_minval_2d()
    CALL test_maxval_2d_dim()
    CALL test_minval_2d_dim()
    CALL test_maxloc_1d()
    CALL test_minloc_1d()
    CALL test_maxval_1d_single()
    CALL test_minval_1d_allneg()

    CALL test_summary()

CONTAINS

    subroutine test_maxval_1d()
        type(BF16), dimension(5) :: arr
        arr = [BF16(1.0_real32), BF16(5.0_real32), BF16(3.0_real32), BF16(-2.0_real32), BF16(4.0_real32)]
        call check_bf16_real64('maxval_1d', maxval(arr), 5.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_minval_1d()
        type(BF16), dimension(5) :: arr
        arr = [BF16(1.0_real32), BF16(5.0_real32), BF16(3.0_real32), BF16(-2.0_real32), BF16(4.0_real32)]
        call check_bf16_real64('minval_1d', minval(arr), -2.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_maxval_2d()
        type(BF16), dimension(2,3) :: arr
        arr = reshape([BF16(1.0_real32), BF16(2.0_real32), &
                       BF16(3.0_real32), BF16(10.0_real32), &
                       BF16(5.0_real32), BF16(-1.0_real32)], [2,3])
        call check_bf16_real64('maxval_2d', maxval(arr), 10.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_minval_2d()
        type(BF16), dimension(2,3) :: arr
        arr = reshape([BF16(1.0_real32), BF16(2.0_real32), &
                       BF16(3.0_real32), BF16(10.0_real32), &
                       BF16(5.0_real32), BF16(-1.0_real32)], [2,3])
        call check_bf16_real64('minval_2d', minval(arr), -1.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_maxval_2d_dim()
        type(BF16), dimension(2,3) :: arr
        type(BF16), dimension(3) :: row_max
        type(BF16), dimension(2) :: col_max
        arr = reshape([BF16(1.0_real32), BF16(2.0_real32), &
                       BF16(3.0_real32), BF16(10.0_real32), &
                       BF16(5.0_real32), BF16(-1.0_real32)], [2,3])
        ! maxval along dim=1: max of each column (column-major, dim=1 gives row of results)
        row_max = maxval(arr, dim=1)
        call check_bf16_real64('maxval_2d_dim1_c1', row_max(1), 2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('maxval_2d_dim1_c2', row_max(2), 10.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('maxval_2d_dim1_c3', row_max(3), 5.0_real64, BF16_TOL_TIGHT)

        ! maxval along dim=2: max of each row
        col_max = maxval(arr, dim=2)
        call check_bf16_real64('maxval_2d_dim2_r1', col_max(1), 5.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('maxval_2d_dim2_r2', col_max(2), 10.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_minval_2d_dim()
        type(BF16), dimension(2,3) :: arr
        type(BF16), dimension(3) :: row_min
        arr = reshape([BF16(1.0_real32), BF16(-2.0_real32), &
                       BF16(3.0_real32), BF16(0.5_real32), &
                       BF16(5.0_real32), BF16(-1.0_real32)], [2,3])
        row_min = minval(arr, dim=1)
        call check_bf16_real64('minval_2d_dim1_c1', row_min(1), -2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('minval_2d_dim1_c2', row_min(2), 0.5_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('minval_2d_dim1_c3', row_min(3), -1.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_maxloc_1d()
        type(BF16), dimension(5) :: arr
        arr = [BF16(1.0_real32), BF16(5.0_real32), BF16(3.0_real32), BF16(-2.0_real32), BF16(4.0_real32)]
        call check_integer('maxloc_1d', maxloc(arr), 2)
    end subroutine

    subroutine test_minloc_1d()
        type(BF16), dimension(5) :: arr
        arr = [BF16(1.0_real32), BF16(5.0_real32), BF16(3.0_real32), BF16(-2.0_real32), BF16(4.0_real32)]
        call check_integer('minloc_1d', minloc(arr), 4)
    end subroutine

    subroutine test_maxval_1d_single()
        type(BF16), dimension(1) :: arr
        arr = [BF16(7.0_real32)]
        call check_bf16_real64('maxval_1d_single', maxval(arr), 7.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_minval_1d_allneg()
        type(BF16), dimension(3) :: arr
        arr = [BF16(-1.0_real32), BF16(-5.0_real32), BF16(-3.0_real32)]
        call check_bf16_real64('minval_1d_allneg', minval(arr), -5.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('maxval_1d_allneg', maxval(arr), -1.0_real64, BF16_TOL_TIGHT)
    end subroutine

END PROGRAM test_bf16_reductions