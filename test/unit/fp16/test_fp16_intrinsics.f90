!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 numeric intrinsic functions.
!
!  FP16 (IEEE half) constants: digits=11, epsilon=9.765625e-4, huge=65504,
!  tiny=6.103515625e-5, minexponent=-13, maxexponent=16, precision=3, range=4

PROGRAM test_fp16_intrinsics
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
    IMPLICIT NONE

    CALL test_digits()
    CALL test_epsilon()
    CALL test_huge()
    CALL test_tiny()
    CALL test_radix()
    CALL test_exponent()
    CALL test_fraction()
    CALL test_minmax_exponent()
    CALL test_precision()
    CALL test_range()
    CALL test_erfc_scaled()
    CALL test_mod()
    CALL test_modulo()
    CALL test_nearest()
    CALL test_nint()
    CALL test_scale()
    CALL test_sign()
    CALL test_isnan_isinf()

    CALL test_summary()

CONTAINS

    subroutine test_digits()
        ! __FLT16_MANT_DIG__ = 11 for IEEE half (10 explicit + 1 implicit)
        call check_integer('digits_fp16', digits(FP16(1.0_real32)), 11)
    end subroutine

    subroutine test_epsilon()
        call check_fp16_real64('epsilon_fp16', epsilon(FP16(1.0_real32)), 9.765625e-4_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_huge()
        call check_fp16_real64('huge_fp16', huge(FP16(1.0_real32)), 65504.0_real64, FP16_TOL_LOOSE)
    end subroutine

    subroutine test_tiny()
        call check_fp16_real64('tiny_fp16', tiny(FP16(1.0_real32)), 6.103515625e-5_real64, FP16_TOL_LOOSE)
    end subroutine

    subroutine test_radix()
        call check_integer('radix_fp16', radix(FP16(1.0_real32)), 2)
    end subroutine

    subroutine test_exponent()
        call check_integer('exponent_1', exponent(FP16(1.0_real32)), 1)
        call check_integer('exponent_0', exponent(FP16(0.0_real32)), 0)
        call check_integer('exponent_4', exponent(FP16(4.0_real32)), 3)
        call check_integer('exponent_2', exponent(FP16(2.0_real32)), 2)
    end subroutine

    subroutine test_fraction()
        ! For FP16(4.0): 4.0 = 0.5 * 2^3, so fraction = 0.5
        call check_fp16_real64('fraction_4', fraction(FP16(4.0_real32)), 0.5_real64, FP16_TOL)
        ! For FP16(1.0): 1.0 = 0.5 * 2^1, so fraction = 0.5
        call check_fp16_real64('fraction_1', fraction(FP16(1.0_real32)), 0.5_real64, FP16_TOL)
    end subroutine

    subroutine test_minmax_exponent()
        call check_integer('minexponent_fp16', minexponent(FP16(1.0_real32)), -13)
        call check_integer('maxexponent_fp16', maxexponent(FP16(1.0_real32)), 16)
    end subroutine

    subroutine test_precision()
        call check_integer('precision_fp16', precision(FP16(1.0_real32)), 3)
    end subroutine

    subroutine test_range()
        ! __FLT16_MIN_10_EXP__ = -4, __FLT16_MAX_10_EXP__ = 4
        ! range = max(|-4|, 4) = 4
        call check_integer('range_fp16', range(FP16(1.0_real32)), 4)
    end subroutine

    subroutine test_erfc_scaled()
        call check_fp16_real64('erfc_scaled_0', erfc_scaled(FP16(0.0_real32)), 1.0_real64, FP16_TOL)
        call check_fp16_real64('erfc_scaled_1', erfc_scaled(FP16(1.0_real32)), &
            erfc(1.0_real64) * exp(1.0_real64 * 1.0_real64), FP16_TOL_LOOSE)
    end subroutine

    subroutine test_mod()
        call check_fp16_real64('mod_7_3',  mod(FP16(7.0_real32), FP16(3.0_real32)), 1.0_real64, FP16_TOL)
        call check_fp16_real64('mod_7_2',  mod(FP16(7.0_real32), FP16(2.0_real32)), 1.0_real64, FP16_TOL)
        call check_fp16_real64('mod_neg7_3', mod(FP16(-7.0_real32), FP16(3.0_real32)), -1.0_real64, FP16_TOL)
    end subroutine

    subroutine test_modulo()
        call check_fp16_real64('modulo_7_3',   modulo(FP16(7.0_real32), FP16(3.0_real32)), 1.0_real64, FP16_TOL)
        call check_fp16_real64('modulo_neg7_3', modulo(FP16(-7.0_real32), FP16(3.0_real32)), 2.0_real64, FP16_TOL)
    end subroutine

    subroutine test_nearest()
        type(FP16) :: x_up, x_down, x
        x = FP16(1.0_real32)
        x_up = nearest(x, FP16(1.0_real32))
        x_down = nearest(x, FP16(-1.0_real32))
        call check_logical('nearest_up_gt_1', dble(x_up) .gt. 1.0_real64, .TRUE.)
        call check_logical('nearest_down_lt_1', dble(x_down) .lt. 1.0_real64, .TRUE.)
    end subroutine

    subroutine test_nint()
        call check_fp16_real64('nint_1.5',  nint(FP16(1.5_real32)), 2.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('nint_0.4',  nint(FP16(0.4_real32)), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('nint_neg0.4', nint(FP16(-0.4_real32)), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('nint_neg1.5', nint(FP16(-1.5_real32)), -2.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_scale()
        call check_fp16_real64('scale_1_2',  scale(FP16(1.0_real32), 2), 4.0_real64, FP16_TOL)
        call check_fp16_real64('scale_1_neg1', scale(FP16(1.0_real32), -1), 0.5_real64, FP16_TOL)
        call check_fp16_real64('scale_3_2',  scale(FP16(3.0_real32), 2), 12.0_real64, FP16_TOL)
    end subroutine

    subroutine test_sign()
        call check_fp16_real64('sign_pos_pos',  sign(FP16(3.0_real32), FP16(1.0_real32)), 3.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('sign_neg_pos',  sign(FP16(-3.0_real32), FP16(1.0_real32)), 3.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('sign_pos_neg',  sign(FP16(3.0_real32), FP16(-1.0_real32)), -3.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('sign_neg_neg',  sign(FP16(-3.0_real32), FP16(-1.0_real32)), -3.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_isnan_isinf()
        call check_logical('isnan_finite', isnan(FP16(1.0_real32)), .FALSE.)
        call check_logical('isnan_zero',  isnan(FP16(0.0_real32)), .FALSE.)
        call check_logical('isinf_finite', isinf(FP16(1.0_real32)), .FALSE.)
        call check_logical('isinf_zero',  isinf(FP16(0.0_real32)), .FALSE.)
    end subroutine

END PROGRAM test_fp16_intrinsics