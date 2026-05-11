!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 single-argument math functions.
!
!  Tests: abs, acos, acosh, asin, asinh, atan, atanh, bessel_j0, bessel_j1,
!         bessel_y0, bessel_y1, ceiling, cos, cosh, erf, erfc, exp, floor,
!         gamma, log, log10, log_gamma, sin, sinh, sqrt, tan, tanh

PROGRAM test_bf16_math_one_arg
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
    IMPLICIT NONE

    CALL test_abs()
    CALL test_acos()
    CALL test_acosh()
    CALL test_asin()
    CALL test_asinh()
    CALL test_atan()
    CALL test_atanh()
    CALL test_bessel_j0()
    CALL test_bessel_j1()
    CALL test_bessel_y0()
    CALL test_bessel_y1()
    CALL test_ceiling()
    CALL test_cos()
    CALL test_cosh()
    CALL test_erf()
    CALL test_erfc()
    CALL test_exp()
    CALL test_floor()
    CALL test_gamma()
    CALL test_log()
    CALL test_log10()
    CALL test_log_gamma()
    CALL test_sin()
    CALL test_sinh()
    CALL test_sqrt()
    CALL test_tan()
    CALL test_tanh()

    CALL test_summary()

CONTAINS

    subroutine test_abs()
        call check_bf16_real64('abs_pos',  abs(BF16(3.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('abs_neg',  abs(BF16(-3.0_real32)), 3.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('abs_zero', abs(BF16(0.0_real32)), 0.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_acos()
        call check_bf16_real64('acos_0',    acos(BF16(0.0_real32)), acos(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('acos_0.5',  acos(BF16(0.5_real32)), acos(0.5_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('acos_1',    acos(BF16(1.0_real32)), acos(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('acos_neg05', acos(BF16(-0.5_real32)), acos(-0.5_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_acosh()
        call check_bf16_real64('acosh_1',  acosh(BF16(1.0_real32)), acosh(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('acosh_2',  acosh(BF16(2.0_real32)), acosh(2.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('acosh_5',  acosh(BF16(5.0_real32)), acosh(5.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_asin()
        call check_bf16_real64('asin_0',    asin(BF16(0.0_real32)), asin(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('asin_0.5',  asin(BF16(0.5_real32)), asin(0.5_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('asin_1',    asin(BF16(1.0_real32)), asin(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('asin_neg05', asin(BF16(-0.5_real32)), asin(-0.5_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_asinh()
        call check_bf16_real64('asinh_0',  asinh(BF16(0.0_real32)), asinh(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('asinh_1',  asinh(BF16(1.0_real32)), asinh(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('asinh_neg1', asinh(BF16(-1.0_real32)), asinh(-1.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_atan()
        call check_bf16_real64('atan_0',   atan(BF16(0.0_real32)), atan(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('atan_1',   atan(BF16(1.0_real32)), atan(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('atan_neg1', atan(BF16(-1.0_real32)), atan(-1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('atan_0.5', atan(BF16(0.5_real32)), atan(0.5_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_atanh()
        call check_bf16_real64('atanh_0',   atanh(BF16(0.0_real32)), atanh(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('atanh_0.5', atanh(BF16(0.5_real32)), atanh(0.5_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('atanh_neg05', atanh(BF16(-0.5_real32)), atanh(-0.5_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_bessel_j0()
        call check_bf16_real64('bessel_j0_0',  bessel_j0(BF16(0.0_real32)), 1.0_real64, BF16_TOL_LOOSE)
        call check_bf16_real64('bessel_j0_1',  bessel_j0(BF16(1.0_real32)), &
            real(bessel_j0(1.0_real64), kind=real64), BF16_TOL_LOOSE)
        call check_bf16_real64('bessel_j0_5',  bessel_j0(BF16(5.0_real32)), &
            real(bessel_j0(5.0_real64), kind=real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_bessel_j1()
        call check_bf16_real64('bessel_j1_0',  bessel_j1(BF16(0.0_real32)), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('bessel_j1_1',  bessel_j1(BF16(1.0_real32)), &
            real(bessel_j1(1.0_real64), kind=real64), BF16_TOL_LOOSE)
        call check_bf16_real64('bessel_j1_5',  bessel_j1(BF16(5.0_real32)), &
            real(bessel_j1(5.0_real64), kind=real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_bessel_y0()
        ! bessel_y0 is undefined at x=0, start at x=1
        call check_bf16_real64('bessel_y0_1',  bessel_y0(BF16(1.0_real32)), &
            real(bessel_y0(1.0_real64), kind=real64), BF16_TOL_LOOSE)
        call check_bf16_real64('bessel_y0_5',  bessel_y0(BF16(5.0_real32)), &
            real(bessel_y0(5.0_real64), kind=real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_bessel_y1()
        call check_bf16_real64('bessel_y1_1',  bessel_y1(BF16(1.0_real32)), &
            real(bessel_y1(1.0_real64), kind=real64), BF16_TOL_LOOSE)
        call check_bf16_real64('bessel_y1_5',  bessel_y1(BF16(5.0_real32)), &
            real(bessel_y1(5.0_real64), kind=real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_ceiling()
        call check_bf16_real64('ceiling_1.5',  ceiling(BF16(1.5_real32)), 2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ceiling_neg1.5', ceiling(BF16(-1.5_real32)), -1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ceiling_0.1',   ceiling(BF16(0.1_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ceiling_100.7', ceiling(BF16(100.7_real32)), 101.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_cos()
        call check_bf16_real64('cos_0',  cos(BF16(0.0_real32)), cos(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('cos_1',  cos(BF16(1.0_real32)), cos(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('cos_neg1', cos(BF16(-1.0_real32)), cos(-1.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_cosh()
        call check_bf16_real64('cosh_0',  cosh(BF16(0.0_real32)), cosh(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('cosh_1',  cosh(BF16(1.0_real32)), cosh(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('cosh_2',  cosh(BF16(2.0_real32)), cosh(2.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_erf()
        call check_bf16_real64('erf_0',   erf(BF16(0.0_real32)), erf(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('erf_0.5', erf(BF16(0.5_real32)), erf(0.5_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('erf_1',   erf(BF16(1.0_real32)), erf(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('erf_3',   erf(BF16(3.0_real32)), erf(3.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_erfc()
        call check_bf16_real64('erfc_0',   erfc(BF16(0.0_real32)), erfc(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('erfc_0.5', erfc(BF16(0.5_real32)), erfc(0.5_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('erfc_1',   erfc(BF16(1.0_real32)), erfc(1.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_exp()
        call check_bf16_real64('exp_0',   exp(BF16(0.0_real32)), exp(0.0_real64), BF16_TOL)
        call check_bf16_real64('exp_1',   exp(BF16(1.0_real32)), exp(1.0_real64), BF16_TOL)
        call check_bf16_real64('exp_neg1', exp(BF16(-1.0_real32)), exp(-1.0_real64), BF16_TOL)
        call check_bf16_real64('exp_2',   exp(BF16(2.0_real32)), exp(2.0_real64), BF16_TOL)
    end subroutine

    subroutine test_floor()
        call check_bf16_real64('floor_1.5',  floor(BF16(1.5_real32)), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('floor_neg1.5', floor(BF16(-1.5_real32)), -2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('floor_0.9',   floor(BF16(0.9_real32)), 0.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_gamma()
        call check_bf16_real64('gamma_1',  gamma(BF16(1.0_real32)), gamma(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('gamma_2',  gamma(BF16(2.0_real32)), gamma(2.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('gamma_5',  gamma(BF16(5.0_real32)), gamma(5.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('gamma_0.5', gamma(BF16(0.5_real32)), gamma(0.5_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_log()
        call check_bf16_real64('log_1',   log(BF16(1.0_real32)), log(1.0_real64), BF16_TOL)
        call check_bf16_real64('log_e',   log(BF16(2.718_real32)), log(2.718_real64), BF16_TOL)
        call check_bf16_real64('log_10',  log(BF16(10.0_real32)), log(10.0_real64), BF16_TOL)
    end subroutine

    subroutine test_log10()
        call check_bf16_real64('log10_1',   log10(BF16(1.0_real32)), log10(1.0_real64), BF16_TOL)
        call check_bf16_real64('log10_10',  log10(BF16(10.0_real32)), log10(10.0_real64), BF16_TOL)
        call check_bf16_real64('log10_100', log10(BF16(100.0_real32)), log10(100.0_real64), BF16_TOL)
    end subroutine

    subroutine test_log_gamma()
        call check_bf16_real64('log_gamma_1',  log_gamma(BF16(1.0_real32)), log_gamma(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('log_gamma_2',  log_gamma(BF16(2.0_real32)), log_gamma(2.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('log_gamma_5',  log_gamma(BF16(5.0_real32)), log_gamma(5.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('log_gamma_10', log_gamma(BF16(10.0_real32)), log_gamma(10.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_sin()
        call check_bf16_real64('sin_0',  sin(BF16(0.0_real32)), sin(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('sin_1',  sin(BF16(1.0_real32)), sin(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('sin_0.5', sin(BF16(0.5_real32)), sin(0.5_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_sinh()
        call check_bf16_real64('sinh_0',  sinh(BF16(0.0_real32)), sinh(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('sinh_1',  sinh(BF16(1.0_real32)), sinh(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('sinh_neg1', sinh(BF16(-1.0_real32)), sinh(-1.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_sqrt()
        call check_bf16_real64('sqrt_0',  sqrt(BF16(0.0_real32)), sqrt(0.0_real64), BF16_TOL_TIGHT)
        call check_bf16_real64('sqrt_1',  sqrt(BF16(1.0_real32)), sqrt(1.0_real64), BF16_TOL)
        call check_bf16_real64('sqrt_2',  sqrt(BF16(2.0_real32)), sqrt(2.0_real64), BF16_TOL)
        call check_bf16_real64('sqrt_4',  sqrt(BF16(4.0_real32)), sqrt(4.0_real64), BF16_TOL)
    end subroutine

    subroutine test_tan()
        call check_bf16_real64('tan_0',  tan(BF16(0.0_real32)), tan(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('tan_1',  tan(BF16(1.0_real32)), tan(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('tan_neg1', tan(BF16(-1.0_real32)), tan(-1.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_tanh()
        call check_bf16_real64('tanh_0',  tanh(BF16(0.0_real32)), tanh(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('tanh_1',  tanh(BF16(1.0_real32)), tanh(1.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('tanh_10', tanh(BF16(10.0_real32)), tanh(10.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('tanh_neg1', tanh(BF16(-1.0_real32)), tanh(-1.0_real64), BF16_TOL_LOOSE)
    end subroutine

END PROGRAM test_bf16_math_one_arg