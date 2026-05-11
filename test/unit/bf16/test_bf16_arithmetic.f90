!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 arithmetic operators.
!
!  Tests: +, -, *, /, ** for BF16-BF16, BF16-real32, real32-BF16 variants,
!         unary minus, edge cases (zero, negative, div-by-zero, pow)

PROGRAM test_bf16_arithmetic
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
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
        call check_bf16_real64('add_bf16_bf16',  BF16(2.0_real32) + BF16(3.0_real32), 5.0_real64, BF16_TOL)
        call check_bf16_real64('add_bf16_real',  BF16(2.0_real32) + 3.0_real32, 5.0_real64, BF16_TOL)
        call check_bf16_real64('add_real_bf16',  2.0_real32 + BF16(3.0_real32), 5.0_real64, BF16_TOL)
        call check_bf16_real64('add_neg',        BF16(-2.0_real32) + BF16(3.0_real32), 1.0_real64, BF16_TOL)
        call check_bf16_real64('add_zero',       BF16(5.0_real32) + BF16(0.0_real32), 5.0_real64, BF16_TOL)
    end subroutine

    ! ---- Subtraction ----
    subroutine test_subtract()
        call check_bf16_real64('sub_bf16_bf16',  BF16(5.0_real32) - BF16(3.0_real32), 2.0_real64, BF16_TOL)
        call check_bf16_real64('sub_bf16_real',  BF16(5.0_real32) - 3.0_real32, 2.0_real64, BF16_TOL)
        call check_bf16_real64('sub_real_bf16' , 5.0_real32 - BF16(3.0_real32), 2.0_real64, BF16_TOL)
        call check_bf16_real64('sub_neg',        BF16(1.0_real32) - BF16(5.0_real32), -4.0_real64, BF16_TOL)
        call check_bf16_real64('sub_zero',       BF16(5.0_real32) - BF16(0.0_real32), 5.0_real64, BF16_TOL)
        call check_bf16_real64('sub_same',       BF16(3.0_real32) - BF16(3.0_real32), 0.0_real64, BF16_TOL_TIGHT)
    end subroutine

    ! ---- Multiplication ----
    subroutine test_multiply()
        call check_bf16_real64('mul_bf16_bf16',  BF16(4.0_real32) * BF16(3.0_real32), 12.0_real64, BF16_TOL)
        call check_bf16_real64('mul_bf16_real',  BF16(4.0_real32) * 3.0_real32, 12.0_real64, BF16_TOL)
        call check_bf16_real64('mul_real_bf16',  4.0_real32 * BF16(3.0_real32), 12.0_real64, BF16_TOL)
        call check_bf16_real64('mul_neg_neg',    BF16(-3.0_real32) * BF16(-4.0_real32), 12.0_real64, BF16_TOL)
        call check_bf16_real64('mul_neg_pos',    BF16(-3.0_real32) * BF16(4.0_real32), -12.0_real64, BF16_TOL)
        call check_bf16_real64('mul_zero',       BF16(5.0_real32) * BF16(0.0_real32), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('mul_one',        BF16(7.0_real32) * BF16(1.0_real32), 7.0_real64, BF16_TOL)
    end subroutine

    ! ---- Division ----
    subroutine test_divide()
        call check_bf16_real64('div_bf16_bf16',  BF16(10.0_real32) / BF16(2.0_real32), 5.0_real64, BF16_TOL)
        call check_bf16_real64('div_bf16_real',  BF16(10.0_real32) / 2.0_real32, 5.0_real64, BF16_TOL)
        call check_bf16_real64('div_real_bf16',  10.0_real32 / BF16(2.0_real32), 5.0_real64, BF16_TOL)
        call check_bf16_real64('div_neg',        BF16(-6.0_real32) / BF16(3.0_real32), -2.0_real64, BF16_TOL)
        call check_bf16_real64('div_one',        BF16(7.0_real32) / BF16(1.0_real32), 7.0_real64, BF16_TOL)
    end subroutine

    ! ---- Power ----
    subroutine test_power()
        call check_bf16_real64('pow_bf16_bf16',  BF16(2.0_real32) ** BF16(3.0_real32), 8.0_real64, BF16_TOL)
        call check_bf16_real64('pow_bf16_real',  BF16(2.0_real32) ** 3.0_real32, 8.0_real64, BF16_TOL)
        call check_bf16_real64('pow_bf16_int',   BF16(2.0_real32) ** 3, 8.0_real64, BF16_TOL)
        call check_bf16_real64('pow_half',       BF16(4.0_real32) ** BF16(0.5_real32), 2.0_real64, BF16_TOL)
        call check_bf16_real64('pow_zero_exp',   BF16(5.0_real32) ** BF16(0.0_real32), 1.0_real64, BF16_TOL)
        call check_bf16_real64('pow_one_base',   BF16(1.0_real32) ** BF16(100.0_real32), 1.0_real64, BF16_TOL)
    end subroutine

    ! ---- Unary minus ----
    subroutine test_unary_minus()
        call check_bf16_real64('unary_minus_pos',  -BF16(5.0_real32), -5.0_real64, BF16_TOL)
        call check_bf16_real64('unary_minus_neg',  -BF16(-5.0_real32), 5.0_real64, BF16_TOL)
        call check_bf16_real64('unary_minus_zero', -BF16(0.0_real32), 0.0_real64, BF16_TOL_TIGHT)
    end subroutine

    ! ---- Edge cases ----
    subroutine test_edge_cases()
        ! Zero arithmetic
        call check_bf16_real64('add_zero_lhs',  BF16(0.0_real32) + BF16(3.0_real32), 3.0_real64, BF16_TOL)
        call check_bf16_real64('mul_zero_lhs',  BF16(0.0_real32) * BF16(3.0_real32), 0.0_real64, BF16_TOL_TIGHT)

        ! Identity operations
        call check_bf16_real64('add_zero_rhs',  BF16(3.0_real32) + BF16(0.0_real32), 3.0_real64, BF16_TOL)
        call check_bf16_real64('mul_one_rhs',   BF16(3.0_real32) * BF16(1.0_real32), 3.0_real64, BF16_TOL)

        ! Division by zero produces inf - verify the result is inf/large
        ! huge(BF16) is the max finite value, but 1/0 gives +inf in bf16
        ! Just verify the result is positive and very large
        call check_logical('div_by_zero_positive', &
            dble(BF16(1.0_real32) / BF16(0.0_real32)) .gt. 0.0_real64, .TRUE.)

        ! Commutativity check
        call check_bf16_real64('add_commutes', &
            BF16(2.0_real32) + BF16(3.0_real32), &
            dble(BF16(3.0_real32) + BF16(2.0_real32)), BF16_TOL_TIGHT)
        call check_bf16_real64('mul_commutes', &
            BF16(4.0_real32) * BF16(5.0_real32), &
            dble(BF16(5.0_real32) * BF16(4.0_real32)), BF16_TOL_TIGHT)
    end subroutine

END PROGRAM test_bf16_arithmetic
