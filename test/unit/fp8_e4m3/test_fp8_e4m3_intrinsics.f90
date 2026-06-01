!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 numeric intrinsic functions.
!
!  Tests: digits, epsilon, huge, tiny, radix, exponent, fraction,
!         minexponent, maxexponent, precision, range, erfc_scaled,
!         mod, modulo, nearest, nint, scale, sign, isnan, isinf
!
!  FP8 E4M3 format: 4 exponent bits, 3 mantissa bits, bias=7
!  digits=4, epsilon=0.125, huge=448, tiny=0.015625 (min normal)
!  minexponent=-6, maxexponent=8, precision=0, range=2
!  No infinity encoding; 1/0 produces NaN (0x7F)

PROGRAM test_fp8_e4m3_intrinsics
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
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
        ! digits() = M_BITS + 1 = 3 + 1 = 4
        call check_integer('digits_fp8_e4m3', digits(FP8_E4M3(1.0_real32)), 4)
    end subroutine

    subroutine test_epsilon()
        ! epsilon = 2^(-3) = 0.125
        call check_fp8_e4m3_real64('epsilon_fp8_e4m3', epsilon(FP8_E4M3(1.0_real32)), 0.125_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_huge()
        ! huge = 448 (max positive normal value)
        call check_fp8_e4m3_real64('huge_fp8_e4m3', huge(FP8_E4M3(1.0_real32)), 448.0_real64, FP8_E4M3_TOL_LOOSE)
    end subroutine

    subroutine test_tiny()
        ! tiny = 2^(-6) = 0.015625 (min positive normal value)
        call check_fp8_e4m3_real64('tiny_fp8_e4m3', tiny(FP8_E4M3(1.0_real32)), 0.015625_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_radix()
        call check_integer('radix_fp8_e4m3', radix(FP8_E4M3(1.0_real32)), 2)
    end subroutine

    subroutine test_exponent()
        call check_integer('exponent_0.25', exponent(FP8_E4M3(0.25_real32)), -1)
        call check_integer('exponent_0.5', exponent(FP8_E4M3(0.5_real32)), 0)
        call check_integer('exponent_0', exponent(FP8_E4M3(0.0_real32)), 0)
        call check_integer('exponent_1', exponent(FP8_E4M3(1.0_real32)), 1)
        call check_integer('exponent_2', exponent(FP8_E4M3(2.0_real32)), 2)
        call check_integer('exponent_4', exponent(FP8_E4M3(4.0_real32)), 3)
        call check_integer('exponent_subnormal_small', exponent(FP8_E4M3(0.001953125_real32)), -6)
        call check_integer('exponent_subnormal_large', exponent(FP8_E4M3(0.013671875_real32)), -6)
    end subroutine

    subroutine test_fraction()
        call check_fp8_e4m3_real64('fraction_0.375', fraction(FP8_E4M3(0.375_real32)), 0.75_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('fraction_0', fraction(FP8_E4M3(0.0_real32)), 0.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('fraction_1', fraction(FP8_E4M3(1.0_real32)), 0.5_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('fraction_2', fraction(FP8_E4M3(2.0_real32)), 0.5_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('fraction_3', fraction(FP8_E4M3(3.0_real32)), 0.75_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('fraction_4', fraction(FP8_E4M3(4.0_real32)), 0.5_real64, FP8_E4M3_TOL)
    end subroutine

    subroutine test_minmax_exponent()
        ! minexponent = -6, maxexponent = 8 (E4M3 format)
        call check_integer('minexponent_fp8_e4m3', minexponent(FP8_E4M3(1.0_real32)), -6)
        call check_integer('maxexponent_fp8_e4m3', maxexponent(FP8_E4M3(1.0_real32)), 8)
    end subroutine

    subroutine test_precision()
        ! precision = INT((digits-1) * LOG10(2)) = INT(3 * 0.301) = 0
        call check_integer('precision_fp8_e4m3', precision(FP8_E4M3(1.0_real32)), 0)
    end subroutine

    subroutine test_range()
        ! range = INT(LOG10(huge)) = INT(LOG10(448)) = INT(2.651) = 2
        call check_integer('range_fp8_e4m3', range(FP8_E4M3(1.0_real32)), 2)
    end subroutine

    subroutine test_erfc_scaled()
        ! erfc_scaled(0) = erfc(0)*exp(0) = 1.0
        call check_fp8_e4m3_real64('erfc_scaled_0', erfc_scaled(FP8_E4M3(0.0_real32)), 1.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('erfc_scaled_1', erfc_scaled(FP8_E4M3(1.0_real32)), &
            erfc(1.0_real64) * exp(1.0_real64 * 1.0_real64), FP8_E4M3_TOL_LOOSE)
    end subroutine

    subroutine test_mod()
        call check_fp8_e4m3_real64('mod_7_3',  mod(FP8_E4M3(7.0_real32), FP8_E4M3(3.0_real32)), 1.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('mod_7_2',  mod(FP8_E4M3(7.0_real32), FP8_E4M3(2.0_real32)), 1.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('mod_neg7_3', mod(FP8_E4M3(-7.0_real32), FP8_E4M3(3.0_real32)), -1.0_real64, FP8_E4M3_TOL)
    end subroutine

    subroutine test_modulo()
        ! modulo always returns non-negative result
        call check_fp8_e4m3_real64('modulo_7_3',   modulo(FP8_E4M3(7.0_real32), FP8_E4M3(3.0_real32)), 1.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('modulo_neg7_3', modulo(FP8_E4M3(-7.0_real32), FP8_E4M3(3.0_real32)), 2.0_real64, FP8_E4M3_TOL)
    end subroutine

    subroutine test_nearest()
        ! nearest(x, positive) = next representable value upward
        ! nearest(x, negative) = next representable value downward
        type(FP8_E4M3) :: x_up, x_down, x
        x = FP8_E4M3(1.0_real32)
        x_up = nearest(x, FP8_E4M3(1.0_real32))
        x_down = nearest(x, FP8_E4M3(-1.0_real32))

        ! nearest up from 1.0 should be > 1.0
        call check_logical('nearest_up_gt_1', dble(x_up) .eq. 1.125_real64, .TRUE.)
        ! nearest down from 1.0 should be < 1.0
        call check_logical('nearest_down_lt_1', dble(x_down) .eq. 0.9375_real64, .TRUE.)

        ! Transitions between subnormal and normal
        x = FP8_E4M3(0.015625_real32) ! min normal (2^-6)
        x_down = nearest(x, FP8_E4M3(-1.0_real32))
        call check_logical('nearest_down_from_min_normal', dble(x_down) .eq. 0.013671875_real64, .TRUE.)

        x = FP8_E4M3(0.013671875_real32) ! max subnormal
        x_up = nearest(x, FP8_E4M3(1.0_real32))
        call check_logical('nearest_up_to_min_normal', dble(x_up) .eq. 0.015625_real64, .TRUE.)

        ! Transition from zero
        x = FP8_E4M3(0.0_real32)
        x_up = nearest(x, FP8_E4M3(1.0_real32))
        call check_logical('nearest_up_from_zero', dble(x_up) .eq. 0.001953125_real64, .TRUE.)
    end subroutine

    subroutine test_nint()
        call check_fp8_e4m3_real64('nint_1.5',  nint(FP8_E4M3(1.5_real32)), 2.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('nint_0.4',  nint(FP8_E4M3(0.4_real32)), 0.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('nint_neg0.4', nint(FP8_E4M3(-0.4_real32)), 0.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('nint_neg1.5', nint(FP8_E4M3(-1.5_real32)), -2.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_scale()
        call check_fp8_e4m3_real64('scale_1_2', &
            scale(FP8_E4M3(1.0_real32), 2), 4.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('scale_1_neg1', &
            scale(FP8_E4M3(1.0_real32), -1), 0.5_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('scale_3_2', &
            scale(FP8_E4M3(3.0_real32), 2), 12.0_real64, FP8_E4M3_TOL)
    end subroutine

    subroutine test_sign()
        ! Known bug: sign() has fp8_e4m3_from_float and fp8_e4m3_to_float swapped,
        ! causing garbage values (e.g., sign(3,1) returns 2.0 instead of 3.0)
        call check_fp8_e4m3_real64('sign_pos_pos', &
            sign(FP8_E4M3(3.0_real32), FP8_E4M3(1.0_real32)), 3.0_real64, &
            FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('sign_neg_pos', &
            sign(FP8_E4M3(-3.0_real32), FP8_E4M3(1.0_real32)), 3.0_real64, &
            FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('sign_pos_neg', &
            sign(FP8_E4M3(3.0_real32), FP8_E4M3(-1.0_real32)), -3.0_real64, &
            FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('sign_neg_neg', &
            sign(FP8_E4M3(-3.0_real32), FP8_E4M3(-1.0_real32)), -3.0_real64, &
            FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_isnan_isinf()
        ! For finite values, isnan and isinf should return false
        call check_logical('isnan_finite', isnan(FP8_E4M3(1.0_real32)), .FALSE.)
        call check_logical('isnan_zero',  isnan(FP8_E4M3(0.0_real32)), .FALSE.)
        call check_logical('isinf_finite', isinf(FP8_E4M3(1.0_real32)), .FALSE.)
        call check_logical('isinf_zero',  isinf(FP8_E4M3(0.0_real32)), .FALSE.)
        ! In FP8 E4M3, 1/0 produces NaN (no infinity encoding)
        ! isnan and isinf both return true for NaN bit patterns
        call check_logical('isnan_div_by_zero', isnan(FP8_E4M3(1.0_real32) / FP8_E4M3(0.0_real32)), .TRUE.)
    end subroutine

END PROGRAM test_fp8_e4m3_intrinsics
