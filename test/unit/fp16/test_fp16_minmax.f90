!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 min/max functions.
!
!  Tests: min(FP16,FP16), max(FP16,FP16), min3, max3
!  Known bug: max3_fp16 has condition 'x .gt. y .and. x .gt. y' instead of
!             'x .gt. y .and. x .gt. z' (fp16.F90 line 1421)

PROGRAM test_fp16_minmax
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
    IMPLICIT NONE

    CALL test_min2()
    CALL test_max2()
    CALL test_min3()
    CALL test_max3()

    CALL test_summary()

CONTAINS

    subroutine test_min2()
        call check_fp16_real64('min_1_2',    min(FP16(1.0_real32), FP16(2.0_real32)), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min_neg1_1', min(FP16(-1.0_real32), FP16(1.0_real32)), -1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min_eq',     min(FP16(3.0_real32), FP16(3.0_real32)), 3.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min_neg',    min(FP16(-5.0_real32), FP16(-2.0_real32)), -5.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_max2()
        call check_fp16_real64('max_1_2',    max(FP16(1.0_real32), FP16(2.0_real32)), 2.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('max_neg1_1', max(FP16(-1.0_real32), FP16(1.0_real32)), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('max_eq',     max(FP16(3.0_real32), FP16(3.0_real32)), 3.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('max_neg',    max(FP16(-5.0_real32), FP16(-2.0_real32)), -2.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_min3()
        ! min3 should work correctly
        call check_fp16_real64('min3_first',  &
            min(FP16(1.0_real32), FP16(2.0_real32), FP16(3.0_real32)), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min3_middle', &
            min(FP16(3.0_real32), FP16(1.0_real32), FP16(2.0_real32)), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min3_last',   &
            min(FP16(3.0_real32), FP16(2.0_real32), FP16(1.0_real32)), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min3_neg',    &
            min(FP16(-1.0_real32), FP16(5.0_real32), FP16(2.0_real32)), -1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('min3_eq',     &
            min(FP16(2.0_real32), FP16(2.0_real32), FP16(2.0_real32)), 2.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_max3()
        ! Cases where max3 works (first arg is the maximum):
        call check_fp16_real64('max3_first',  &
            max(FP16(7.0_real32), FP16(3.0_real32), FP16(2.0_real32)), 7.0_real64, FP16_TOL_TIGHT)
        ! Cases where the else branch catches it (z is max, x <= y):
        call check_fp16_real64('max3_last_else', &
            max(FP16(1.0_real32), FP16(2.0_real32), FP16(3.0_real32)), 3.0_real64, FP16_TOL_TIGHT)
        ! Equal values:
        call check_fp16_real64('max3_eq',     &
            max(FP16(2.0_real32), FP16(2.0_real32), FP16(2.0_real32)), 2.0_real64, FP16_TOL_TIGHT)

        ! KNOWN BUG: max3_fp16 condition 'x.gt.y .and. x.gt.y' duplicates y check
        ! x=5 > y=3 but x=5 < z=7, condition passes, returns 5 (wrong)
        call check_fp16_real64('max3_x_gt_y_lt_z', &
            max(FP16(5.0_real32), FP16(3.0_real32), FP16(7.0_real32)), 7.0_real64, &
            FP16_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp16_minmax
