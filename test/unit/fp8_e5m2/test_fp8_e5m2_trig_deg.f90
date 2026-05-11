!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E5M2 degree-based trigonometric functions.
!
!  Tests: cosd, sind, tand, cotan, cotand, acosd, asind, atand, atan2d
!  Note: Many trig_deg functions have a systematic bug where fp8_e5m2_from_float
!  and fp8_e5m2_to_float are swapped in fp8_e5m2_math_trig.c. Failing tests
!  are marked as known bugs until the implementation is fixed.

PROGRAM test_fp8_e5m2_trig_deg
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E5M2
    USE lpf_fp8_e5m2_test_utils
    IMPLICIT NONE

    real(real64), parameter :: PI = 3.14159265358979323846_real64
    real(real64), parameter :: DEG2RAD = PI / 180.0_real64
    real(real64), parameter :: RAD2DEG = 180.0_real64 / PI

    CALL test_cosd()
    CALL test_sind()
    CALL test_tand()
    CALL test_cotan()
    CALL test_cotand()
    CALL test_acosd()
    CALL test_asind()
    CALL test_atand()
    CALL test_atan2d()

    CALL test_summary()

CONTAINS

    subroutine test_cosd()
        ! Known bug: cosd has from_float/to_float swapped, returns 0 for most inputs
        call check_fp8_e5m2_real64('cosd_0',   &
            cosd(FP8_E5M2(0.0_real32)), cos(0.0_real64), FP8_E5M2_TOL_LOOSE)
        ! Since 90_fp8_e5m2 is 96, we compare with cosd(96.0D0)
        call check_fp8_e5m2_real64('cosd_90',  &
            cosd(FP8_E5M2(90.0_real32)), cos(96.0_real64 * DEG2RAD), FP8_E5M2_TOL_COARSE*3.0D0)
        call check_fp8_e5m2_real64('cosd_60',  &
            cosd(FP8_E5M2(60.0_real32)), cos(60.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('cosd_180', &
            cosd(FP8_E5M2(180.0_real32)), cos(180.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_sind()
        call check_fp8_e5m2_real64('sind_0',   &
            sind(FP8_E5M2(0.0_real32)), sin(0.0_real64), FP8_E5M2_TOL_LOOSE)
        ! Known bug: sind has from_float/to_float swapped
        call check_fp8_e5m2_real64('sind_90',  &
            sind(FP8_E5M2(90.0_real32)), sin(90.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('sind_30',  &
            sind(FP8_E5M2(30.0_real32)), sin(30.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('sind_45',  &
            sind(FP8_E5M2(45.0_real32)), sin(45.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)

    end subroutine

    subroutine test_tand()
        call check_fp8_e5m2_real64('tand_0',   &
            tand(FP8_E5M2(0.0_real32)), tan(0.0_real64), FP8_E5M2_TOL_LOOSE)
        ! Known bug: tand has from_float/to_float swapped, returns wrong values
        call check_fp8_e5m2_real64('tand_45',  &
            tand(FP8_E5M2(45.0_real32)), tan(45.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('tand_30',  &
            tand(FP8_E5M2(30.0_real32)), tan(30.0_real64 * DEG2RAD), FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_cotan()
        ! Known bug: cotan has from_float/to_float swapped
        call check_fp8_e5m2_real64('cotan_1',  &
            cotan(FP8_E5M2(1.0_real32)), -tan(1.0_real64 + PI/2.0_real64), FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('cotan_0.5', &
            cotan(FP8_E5M2(0.5_real32)), -tan(0.5_real64 + PI/2.0_real64), FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_cotand()
        ! Known bug: cotand has from_float/to_float swapped
        call check_fp8_e5m2_real64('cotand_45', &
            cotand(FP8_E5M2(45.0_real32)), &
            -tan((45.0_real64 + 90.0_real64) * DEG2RAD), FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('cotand_90', &
            cotand(FP8_E5M2(90.0_real32)), &
            -tan((96.0_real64 + 90.0_real64) * DEG2RAD), FP8_E5M2_TOL_COARSE*3.0D0)
    end subroutine

    subroutine test_acosd()
        ! Known bug: acosd has from_float/to_float swapped
        call check_fp8_e5m2_real64('acosd_0',   &
            acosd(FP8_E5M2(0.0_real32)), acos(0.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('acosd_1',   &
            acosd(FP8_E5M2(1.0_real32)), acos(1.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('acosd_0.5', &
            acosd(FP8_E5M2(0.5_real32)), acos(0.5_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_asind()
        call check_fp8_e5m2_real64('asind_0',   &
            asind(FP8_E5M2(0.0_real32)), asin(0.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        ! Known bug: asind has from_float/to_float swapped
        call check_fp8_e5m2_real64('asind_1',   &
            asind(FP8_E5M2(1.0_real32)), asin(1.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('asind_0.5', &
            asind(FP8_E5M2(0.5_real32)), asin(0.5_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_atand()
        call check_fp8_e5m2_real64('atand_0',   &
            atand(FP8_E5M2(0.0_real32)), atan(0.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        ! Known bug: atand has from_float/to_float swapped
        call check_fp8_e5m2_real64('atand_1',   &
            atand(FP8_E5M2(1.0_real32)), atan(1.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('atand_neg1', &
            atand(FP8_E5M2(-1.0_real32)), atan(-1.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_atan2d()
        call check_fp8_e5m2_real64('atan2d_1_1', &
            atan2d(FP8_E5M2(1.0_real32), FP8_E5M2(1.0_real32)), &
            atan2(1.0_real64, 1.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('atan2d_0_1', &
            atan2d(FP8_E5M2(0.0_real32), FP8_E5M2(1.0_real32)), &
            atan2(0.0_real64, 1.0_real64) * RAD2DEG, FP8_E5M2_TOL_LOOSE)
    end subroutine

END PROGRAM test_fp8_e5m2_trig_deg
