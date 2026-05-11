!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 arithmetic operators.

PROGRAM test_fp16_arithmetic
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
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

    subroutine test_add()
        call check_fp16_real64('add_fp16_fp16',  FP16(2.0_real32) + FP16(3.0_real32), 5.0_real64, FP16_TOL)
        call check_fp16_real64('add_fp16_real',  FP16(2.0_real32) + 3.0_real32, 5.0_real64, FP16_TOL)
        call check_fp16_real64('add_real_fp16',  2.0_real32 + FP16(3.0_real32), 5.0_real64, FP16_TOL)
        call check_fp16_real64('add_neg',        FP16(-2.0_real32) + FP16(3.0_real32), 1.0_real64, FP16_TOL)
        call check_fp16_real64('add_zero',       FP16(5.0_real32) + FP16(0.0_real32), 5.0_real64, FP16_TOL)
    end subroutine

    subroutine test_subtract()
        call check_fp16_real64('sub_fp16_fp16',  FP16(5.0_real32) - FP16(3.0_real32), 2.0_real64, FP16_TOL)
        call check_fp16_real64('sub_fp16_real',  FP16(5.0_real32) - 3.0_real32, 2.0_real64, FP16_TOL)
        call check_fp16_real64('sub_real_fp16',  5.0_real32 - FP16(3.0_real32), 2.0_real64, FP16_TOL)
        call check_fp16_real64('sub_neg',        FP16(1.0_real32) - FP16(5.0_real32), -4.0_real64, FP16_TOL)
        call check_fp16_real64('sub_zero',       FP16(5.0_real32) - FP16(0.0_real32), 5.0_real64, FP16_TOL)
        call check_fp16_real64('sub_same',       FP16(3.0_real32) - FP16(3.0_real32), 0.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_multiply()
        call check_fp16_real64('mul_fp16_fp16',  FP16(4.0_real32) * FP16(3.0_real32), 12.0_real64, FP16_TOL)
        call check_fp16_real64('mul_fp16_real',  FP16(4.0_real32) * 3.0_real32, 12.0_real64, FP16_TOL)
        call check_fp16_real64('mul_real_fp16',  4.0_real32 * FP16(3.0_real32), 12.0_real64, FP16_TOL)
        call check_fp16_real64('mul_neg_neg',    FP16(-3.0_real32) * FP16(-4.0_real32), 12.0_real64, FP16_TOL)
        call check_fp16_real64('mul_neg_pos',    FP16(-3.0_real32) * FP16(4.0_real32), -12.0_real64, FP16_TOL)
        call check_fp16_real64('mul_zero',       FP16(5.0_real32) * FP16(0.0_real32), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('mul_one',        FP16(7.0_real32) * FP16(1.0_real32), 7.0_real64, FP16_TOL)
    end subroutine

    subroutine test_divide()
        call check_fp16_real64('div_fp16_fp16',  FP16(10.0_real32) / FP16(2.0_real32), 5.0_real64, FP16_TOL)
        call check_fp16_real64('div_fp16_real',  FP16(10.0_real32) / 2.0_real32, 5.0_real64, FP16_TOL)
        call check_fp16_real64('div_real_fp16',  10.0_real32 / FP16(2.0_real32), 5.0_real64, FP16_TOL)
        call check_fp16_real64('div_neg',        FP16(-6.0_real32) / FP16(3.0_real32), -2.0_real64, FP16_TOL)
        call check_fp16_real64('div_one',        FP16(7.0_real32) / FP16(1.0_real32), 7.0_real64, FP16_TOL)
    end subroutine

    subroutine test_power()
        call check_fp16_real64('pow_fp16_fp16',  FP16(2.0_real32) ** FP16(3.0_real32), 8.0_real64, FP16_TOL)
        call check_fp16_real64('pow_fp16_real',  FP16(2.0_real32) ** 3.0_real32, 8.0_real64, FP16_TOL)
        call check_fp16_real64('pow_fp16_int',   FP16(2.0_real32) ** 3, 8.0_real64, FP16_TOL)
        call check_fp16_real64('pow_half',       FP16(4.0_real32) ** FP16(0.5_real32), 2.0_real64, FP16_TOL)
        call check_fp16_real64('pow_zero_exp',   FP16(5.0_real32) ** FP16(0.0_real32), 1.0_real64, FP16_TOL)
        call check_fp16_real64('pow_one_base',   FP16(1.0_real32) ** FP16(100.0_real32), 1.0_real64, FP16_TOL)
    end subroutine

    subroutine test_unary_minus()
        call check_fp16_real64('unary_minus_pos',  -FP16(5.0_real32), -5.0_real64, FP16_TOL)
        call check_fp16_real64('unary_minus_neg',  -FP16(-5.0_real32), 5.0_real64, FP16_TOL)
        call check_fp16_real64('unary_minus_zero', -FP16(0.0_real32), 0.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_edge_cases()
        call check_fp16_real64('add_zero_lhs',  FP16(0.0_real32) + FP16(3.0_real32), 3.0_real64, FP16_TOL)
        call check_fp16_real64('mul_zero_lhs',  FP16(0.0_real32) * FP16(3.0_real32), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('add_zero_rhs',  FP16(3.0_real32) + FP16(0.0_real32), 3.0_real64, FP16_TOL)
        call check_fp16_real64('mul_one_rhs',   FP16(3.0_real32) * FP16(1.0_real32), 3.0_real64, FP16_TOL)
        call check_logical('div_by_zero_positive', &
            dble(FP16(1.0_real32) / FP16(0.0_real32)) .gt. 0.0_real64, .TRUE.)
        call check_fp16_real64('add_commutes', &
            FP16(2.0_real32) + FP16(3.0_real32), &
            dble(FP16(3.0_real32) + FP16(2.0_real32)), FP16_TOL_TIGHT)
        call check_fp16_real64('mul_commutes', &
            FP16(4.0_real32) * FP16(5.0_real32), &
            dble(FP16(5.0_real32) * FP16(4.0_real32)), FP16_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp16_arithmetic