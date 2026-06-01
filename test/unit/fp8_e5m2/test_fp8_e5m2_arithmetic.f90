!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E5M2 arithmetic operators.
!
!  Tests: +, -, *, /, ** for FP8_E5M2-FP8_E5M2, FP8_E5M2-real32, real32-FP8_E5M2 variants,
!         unary minus, edge cases (zero, identity, div-by-zero)
!  Note: FP8 E4M3 max value is 448, no infinity. 1/0 produces NaN.

PROGRAM test_fp8_e5m2_arithmetic
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E5M2
    USE lpf_fp8_e5m2_test_utils
    IMPLICIT NONE

    CALL test_add()
    CALL test_subtract()
    CALL test_multiply()
    CALL test_divide()
    CALL test_power()
    CALL test_unary_minus()
    CALL test_edge_cases()

    CALL test_summary()

CONTAINS

    ! ---- Addition ----
    subroutine test_add()
        call check_fp8_e5m2_real64('add_fp8_fp8',  FP8_E5M2(2.0_real32) + FP8_E5M2(3.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('add_fp8_real',  FP8_E5M2(2.0_real32) + 3.0_real32, 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('add_real_fp8',  2.0_real32 + FP8_E5M2(3.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('add_fp8_real64',  FP8_E5M2(2.0_real32) + 3.0_real64, 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('add_real64_fp8',  3.0_real64 + FP8_E5M2(2.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('add_neg',       FP8_E5M2(-2.0_real32) + FP8_E5M2(3.0_real32), 1.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('add_zero',      FP8_E5M2(5.0_real32) + FP8_E5M2(0.0_real32), 5.0_real64, FP8_E5M2_TOL)
    end subroutine

    ! ---- Subtraction ----
    subroutine test_subtract()
        call check_fp8_e5m2_real64('sub_fp8_fp8',  FP8_E5M2(5.0_real32) - FP8_E5M2(3.0_real32), 2.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_fp8_real',  FP8_E5M2(5.0_real32) - 3.0_real32, 2.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_real_fp8',  5.0_real32 - FP8_E5M2(3.0_real32), 2.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_fp8_real64', FP8_E5M2(5.0_real32) - 3.0_real64, 2.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_real64_fp8', 5.0_real64 - FP8_E5M2(3.0_real32), 2.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_neg',       FP8_E5M2(1.0_real32) - FP8_E5M2(5.0_real32), -4.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_zero',      FP8_E5M2(5.0_real32) - FP8_E5M2(0.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('sub_same',      FP8_E5M2(3.0_real32) - FP8_E5M2(3.0_real32), 0.0_real64, FP8_E5M2_TOL_TIGHT)
    end subroutine

    ! ---- Multiplication ----
    subroutine test_multiply()
        call check_fp8_e5m2_real64('mul_fp8_fp8',  FP8_E5M2(4.0_real32) * FP8_E5M2(3.0_real32), 12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_fp8_real',  FP8_E5M2(4.0_real32) * 3.0_real32, 12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_real_fp8',  4.0_real32 * FP8_E5M2(3.0_real32), 12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_fp8_real64', FP8_E5M2(4.0_real32) * 3.0_real64, 12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_real64_fp8', 4.0_real64 * FP8_E5M2(3.0_real32), 12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_neg_neg',   FP8_E5M2(-3.0_real32) * FP8_E5M2(-4.0_real32), 12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_neg_pos',   FP8_E5M2(-3.0_real32) * FP8_E5M2(4.0_real32), -12.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_zero',      FP8_E5M2(5.0_real32) * FP8_E5M2(0.0_real32), 0.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('mul_one',       FP8_E5M2(7.0_real32) * FP8_E5M2(1.0_real32), 7.0_real64, FP8_E5M2_TOL)
    end subroutine

    ! ---- Division ----
    subroutine test_divide()
        call check_fp8_e5m2_real64('div_fp8_fp8',  FP8_E5M2(10.0_real32) / FP8_E5M2(2.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('div_fp8_real',  FP8_E5M2(10.0_real32) / 2.0_real32, 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('div_real_fp8',  10.0_real32 / FP8_E5M2(2.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('div_fp8_real64', FP8_E5M2(10.0_real32) / 2.0_real64, 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('div_real64_fp8', 10.0_real64 / FP8_E5M2(2.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('div_neg',       FP8_E5M2(-6.0_real32) / FP8_E5M2(3.0_real32), -2.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('div_one',       FP8_E5M2(7.0_real32) / FP8_E5M2(1.0_real32), 7.0_real64, FP8_E5M2_TOL)
    end subroutine

    ! ---- Power ----
    subroutine test_power()
        call check_fp8_e5m2_real64('pow_fp8_fp8',  FP8_E5M2(2.0_real32) ** FP8_E5M2(3.0_real32), 8.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('pow_fp8_real',  FP8_E5M2(2.0_real32) ** 3.0_real32, 8.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('pow_fp8_real64', FP8_E5M2(2.0_real32) ** 3.0_real64, 8.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('pow_fp8_int',   FP8_E5M2(2.0_real32) ** 3, 8.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('pow_half',      FP8_E5M2(4.0_real32) ** FP8_E5M2(0.5_real32), 2.0_real64, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('pow_zero_exp',  FP8_E5M2(5.0_real32) ** FP8_E5M2(0.0_real32), 1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('pow_one_base',  FP8_E5M2(1.0_real32) ** FP8_E5M2(10.0_real32), 1.0_real64, FP8_E5M2_TOL_TIGHT)
    end subroutine

    ! ---- Unary minus ----
    subroutine test_unary_minus()
        call check_fp8_e5m2_real64('unary_minus_pos',  -FP8_E5M2(5.0_real32), -5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('unary_minus_neg',  -FP8_E5M2(-5.0_real32), 5.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('unary_minus_zero', -FP8_E5M2(0.0_real32), 0.0_real64, FP8_E5M2_TOL_TIGHT)
    end subroutine

    ! ---- Edge cases ----
    subroutine test_edge_cases()
        ! Zero arithmetic
        call check_fp8_e5m2_real64('add_zero_lhs',  FP8_E5M2(0.0_real32) + FP8_E5M2(3.0_real32), 3.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_zero_lhs',  FP8_E5M2(0.0_real32) * FP8_E5M2(3.0_real32), 0.0_real64, FP8_E5M2_TOL_TIGHT)

        ! Identity operations
        call check_fp8_e5m2_real64('add_zero_rhs',  FP8_E5M2(3.0_real32) + FP8_E5M2(0.0_real32), 3.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('mul_one_rhs',   FP8_E5M2(3.0_real32) * FP8_E5M2(1.0_real32), 3.0_real64, FP8_E5M2_TOL)

        ! Division by zero produces Infinity in FP8 E5M2
        ! Verify the result is Inf
        call check_logical('div_by_zero_is_inf', isinf(FP8_E5M2(1.0_real32) / FP8_E5M2(0.0_real32)), .TRUE.)

        ! Commutativity check
        call check_fp8_e5m2_real64('add_commutes', &
            FP8_E5M2(2.0_real32) + FP8_E5M2(3.0_real32), &
            dble(FP8_E5M2(3.0_real32) + FP8_E5M2(2.0_real32)), FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('mul_commutes', &
            FP8_E5M2(4.0_real32) * FP8_E5M2(5.0_real32), &
            dble(FP8_E5M2(5.0_real32) * FP8_E5M2(4.0_real32)), FP8_E5M2_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp8_e5m2_arithmetic
