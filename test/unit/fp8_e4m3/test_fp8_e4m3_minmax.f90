!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 min/max functions.
!
!  Tests: min(FP8_E4M3,FP8_E4M3), max(FP8_E4M3,FP8_E4M3), min3, max3
!  Known bug: max3 has condition 'x .gt. y .and. x .gt. y' instead of
!             'x .gt. y .and. x .gt. z' (same as bf16/fp16)

PROGRAM test_fp8_e4m3_minmax
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
    IMPLICIT NONE

    CALL test_min2()
    CALL test_max2()
    CALL test_min3()
    CALL test_max3()

    CALL test_summary()

CONTAINS

    subroutine test_min2()
        call check_fp8_e4m3_real64('min_1_2',    min(FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32)), 1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min_neg1_1', &
            min(FP8_E4M3(-1.0_real32), FP8_E4M3(1.0_real32)), -1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min_eq',     min(FP8_E4M3(3.0_real32), FP8_E4M3(3.0_real32)), 3.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min_neg',    &
            min(FP8_E4M3(-5.0_real32), FP8_E4M3(-2.0_real32)), -5.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_max2()
        call check_fp8_e4m3_real64('max_1_2',    max(FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32)), 2.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('max_neg1_1', max(FP8_E4M3(-1.0_real32), FP8_E4M3(1.0_real32)), 1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('max_eq',     max(FP8_E4M3(3.0_real32), FP8_E4M3(3.0_real32)), 3.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('max_neg',    &
            max(FP8_E4M3(-5.0_real32), FP8_E4M3(-2.0_real32)), -2.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_min3()
        ! min3 should work correctly (condition: x.lt.y.and.x.lt.z)
        call check_fp8_e4m3_real64('min3_first',  &
            min(FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(3.0_real32)), 1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min3_middle', &
            min(FP8_E4M3(3.0_real32), FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32)), 1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min3_last',   &
            min(FP8_E4M3(3.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(1.0_real32)), 1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min3_neg',    &
            min(FP8_E4M3(-1.0_real32), FP8_E4M3(5.0_real32), FP8_E4M3(2.0_real32)), -1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('min3_eq',     &
            min(FP8_E4M3(2.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(2.0_real32)), 2.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_max3()
        ! max3 has a known bug: condition 'x .gt. y .and. x .gt. y' should be
        ! 'x .gt. y .and. x .gt. z'. When x > y but x < z, the condition still
        ! passes (because x.gt.y is checked twice) and x is returned instead of z.

        ! Cases where max3 works (first arg is the maximum):
        call check_fp8_e4m3_real64('max3_first',  &
            max(FP8_E4M3(7.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(2.0_real32)), 7.0_real64, FP8_E4M3_TOL_TIGHT)
        ! Cases where the else branch catches it (z is max, x <= y):
        call check_fp8_e4m3_real64('max3_last_else', &
            max(FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(3.0_real32)), 3.0_real64, FP8_E4M3_TOL_TIGHT)
        ! Equal values:
        call check_fp8_e4m3_real64('max3_eq',     &
            max(FP8_E4M3(2.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(2.0_real32)), 2.0_real64, FP8_E4M3_TOL_TIGHT)

        call check_fp8_e4m3_real64('max3_x_gt_y_lt_z', &
            max(FP8_E4M3(5.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(7.0_real32)), 7.0_real64, &
            FP8_E4M3_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp8_e4m3_minmax
