!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 min/max functions.
!
!  Tests: min(BF16,BF16), max(BF16,BF16), min3, max3
!  Known bug: max3_bf16 has condition 'x .gt. y .and. x .gt. y' instead of
!             'x .gt. y .and. x .gt. z' (bf16.F90 line 1421)

PROGRAM test_bf16_minmax
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
    IMPLICIT NONE

    CALL test_min2()
    CALL test_max2()
    CALL test_min3()
    CALL test_max3()

    CALL test_summary()

CONTAINS

    subroutine test_min2()
        call check_bf16_real64('min_1_2',    min(BF16(1.0_real32), BF16(2.0_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min_neg1_1', min(BF16(-1.0_real32), BF16(1.0_real32)), -1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min_eq',     min(BF16(3.0_real32), BF16(3.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min_neg',    min(BF16(-5.0_real32), BF16(-2.0_real32)), -5.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_max2()
        call check_bf16_real64('max_1_2',    max(BF16(1.0_real32), BF16(2.0_real32)), 2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('max_neg1_1', max(BF16(-1.0_real32), BF16(1.0_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('max_eq',     max(BF16(3.0_real32), BF16(3.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('max_neg',    max(BF16(-5.0_real32), BF16(-2.0_real32)), -2.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_min3()
        ! min3 should work correctly (verified: line 1397 has x.lt.y.and.x.lt.z)
        call check_bf16_real64('min3_first',  &
            min(BF16(1.0_real32), BF16(2.0_real32), BF16(3.0_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min3_middle', &
            min(BF16(3.0_real32), BF16(1.0_real32), BF16(2.0_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min3_last',   &
            min(BF16(3.0_real32), BF16(2.0_real32), BF16(1.0_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min3_neg',    &
            min(BF16(-1.0_real32), BF16(5.0_real32), BF16(2.0_real32)), -1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('min3_eq',     &
            min(BF16(2.0_real32), BF16(2.0_real32), BF16(2.0_real32)), 2.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_max3()
        ! max3_bf16 has a known bug: condition 'x .gt. y .and. x .gt. y' should be
        ! 'x .gt. y .and. x .gt. z'. When x > y but x < z, the condition still
        ! passes (because x.gt.y is checked twice) and x is returned instead of z.

        ! Cases where max3 works (first arg is the maximum):
        call check_bf16_real64('max3_first',  &
            max(BF16(7.0_real32), BF16(3.0_real32), BF16(2.0_real32)), 7.0_real64, BF16_TOL_TIGHT)
        ! Cases where the else branch catches it (z is max, x <= y):
        call check_bf16_real64('max3_last_else', &
            max(BF16(1.0_real32), BF16(2.0_real32), BF16(3.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        ! Equal values:
        call check_bf16_real64('max3_eq',     &
            max(BF16(2.0_real32), BF16(2.0_real32), BF16(2.0_real32)), 2.0_real64, BF16_TOL_TIGHT)

        ! : x=5 > y=3 but x=5 < z=7, condition passes, returns 5 (wrong)
        call check_bf16_real64('max3_x_gt_y_lt_z', &
            max(BF16(5.0_real32), BF16(3.0_real32), BF16(7.0_real32)), 7.0_real64, &
            BF16_TOL_TIGHT)
    end subroutine

END PROGRAM test_bf16_minmax
