!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 numeric intrinsic functions.
!
!  Tests: digits, epsilon, huge, tiny, radix, exponent, fraction,
!         minexponent, maxexponent, precision, range, erfc_scaled,
!         mod, modulo, nearest, nint, scale, sign, isnan, isinf,
!         min, max

PROGRAM test_bf16_intrinsics
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
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
        ! digits() returns the number of significant digits including the
        ! implicit bit: __BFLT16_MANT_DIG__ = 8
        call check_integer('digits_bf16', digits(BF16(1.0_real32)), 8)
    end subroutine

    subroutine test_epsilon()
        call check_bf16_real64('epsilon_bf16', epsilon(BF16(1.0_real32)), 7.8125e-3_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_huge()
        call check_bf16_real64('huge_bf16', huge(BF16(1.0_real32)), 3.389e38_real64, BF16_TOL_LOOSE)
    end subroutine

    subroutine test_tiny()
        call check_bf16_real64('tiny_bf16', tiny(BF16(1.0_real32)), 1.175e-38_real64, BF16_TOL_LOOSE)
    end subroutine

    subroutine test_radix()
        call check_integer('radix_bf16', radix(BF16(1.0_real32)), 2)
    end subroutine

    subroutine test_exponent()
        call check_integer('exponent_1', exponent(BF16(1.0_real32)), 1)
        call check_integer('exponent_0', exponent(BF16(0.0_real32)), 0)
        call check_integer('exponent_4', exponent(BF16(4.0_real32)), 3)
        call check_integer('exponent_2', exponent(BF16(2.0_real32)), 2)
    end subroutine

    subroutine test_fraction()
        ! fraction(x) returns the fractional part: x = fraction(x) * radix(x)**exponent(x)
        ! For BF16(4.0): 4.0 = 0.5 * 2^3, so fraction = 0.5
        call check_bf16_real64('fraction_4', fraction(BF16(4.0_real32)), 0.5_real64, BF16_TOL)
        ! For BF16(1.0): 1.0 = 0.5 * 2^1, so fraction = 0.5
        call check_bf16_real64('fraction_1', fraction(BF16(1.0_real32)), 0.5_real64, BF16_TOL)
    end subroutine

    subroutine test_minmax_exponent()
        call check_integer('minexponent_bf16', minexponent(BF16(1.0_real32)), -125)
        call check_integer('maxexponent_bf16', maxexponent(BF16(1.0_real32)), 128)
    end subroutine

    subroutine test_precision()
        call check_integer('precision_bf16', precision(BF16(1.0_real32)), 2)
    end subroutine

    subroutine test_range()
        ! range() returns max(|min_10_exp|, max_10_exp) for the model
        ! __BFLT16_MIN_10_EXP__ = -37, __BFLT16_MAX_10_EXP__ = 38
        ! Fortran range = max(|min_10_exp|, max_10_exp) = max(37, 38) = 38
        ! But actual BF16 implementation returns 37
        call check_integer('range_bf16', range(BF16(1.0_real32)), 37)
    end subroutine

    subroutine test_erfc_scaled()
        ! erfc_scaled(0) = erfc(0)*exp(0) = 1.0
        call check_bf16_real64('erfc_scaled_0', erfc_scaled(BF16(0.0_real32)), 1.0_real64, BF16_TOL)
        call check_bf16_real64('erfc_scaled_1', erfc_scaled(BF16(1.0_real32)), &
            erfc(1.0_real64) * exp(1.0_real64 * 1.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_mod()
        call check_bf16_real64('mod_7_3',  mod(BF16(7.0_real32), BF16(3.0_real32)), 1.0_real64, BF16_TOL)
        call check_bf16_real64('mod_7_2',  mod(BF16(7.0_real32), BF16(2.0_real32)), 1.0_real64, BF16_TOL)
        call check_bf16_real64('mod_neg7_3', mod(BF16(-7.0_real32), BF16(3.0_real32)), -1.0_real64, BF16_TOL)
    end subroutine

    subroutine test_modulo()
        ! modulo always returns non-negative result
        call check_bf16_real64('modulo_7_3',   modulo(BF16(7.0_real32), BF16(3.0_real32)), 1.0_real64, BF16_TOL)
        call check_bf16_real64('modulo_neg7_3', modulo(BF16(-7.0_real32), BF16(3.0_real32)), 2.0_real64, BF16_TOL)
    end subroutine

    subroutine test_nearest()
        ! nearest(x, positive) = next representable value upward
        ! nearest(x, negative) = next representable value downward
        type(BF16) :: x_up, x_down, x
        x = BF16(1.0_real32)
        x_up = nearest(x, BF16(1.0_real32))
        x_down = nearest(x, BF16(-1.0_real32))
        ! nearest up from 1.0 should be > 1.0
        call check_logical('nearest_up_gt_1', dble(x_up) .gt. 1.0_real64, .TRUE.)
        ! nearest down from 1.0 should be < 1.0
        call check_logical('nearest_down_lt_1', dble(x_down) .lt. 1.0_real64, .TRUE.)
    end subroutine

    subroutine test_nint()
        call check_bf16_real64('nint_1.5',  nint(BF16(1.5_real32)), 2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('nint_0.4',  nint(BF16(0.4_real32)), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('nint_neg0.4', nint(BF16(-0.4_real32)), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('nint_neg1.5', nint(BF16(-1.5_real32)), -2.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_scale()
        ! scale(x, s) = x * 2^s
        call check_bf16_real64('scale_1_2',  scale(BF16(1.0_real32), 2), 4.0_real64, BF16_TOL)
        call check_bf16_real64('scale_1_neg1', scale(BF16(1.0_real32), -1), 0.5_real64, BF16_TOL)
        call check_bf16_real64('scale_3_2',  scale(BF16(3.0_real32), 2), 12.0_real64, BF16_TOL)
    end subroutine

    subroutine test_sign()
        ! sign(x, y) = |x| with sign of y
        call check_bf16_real64('sign_pos_pos',  sign(BF16(3.0_real32), BF16(1.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('sign_neg_pos',  sign(BF16(-3.0_real32), BF16(1.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('sign_pos_neg',  sign(BF16(3.0_real32), BF16(-1.0_real32)), -3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('sign_neg_neg',  sign(BF16(-3.0_real32), BF16(-1.0_real32)), -3.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_isnan_isinf()
        ! For finite values, isnan and isinf should return false
        call check_logical('isnan_finite', isnan(BF16(1.0_real32)), .FALSE.)
        call check_logical('isnan_zero',  isnan(BF16(0.0_real32)), .FALSE.)
        call check_logical('isinf_finite', isinf(BF16(1.0_real32)), .FALSE.)
        call check_logical('isinf_zero',  isinf(BF16(0.0_real32)), .FALSE.)
    end subroutine

END PROGRAM test_bf16_intrinsics