!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 degree-based trigonometric functions.
!
!  Tests: cosd, sind, tand, cotan, cotand, acosd, asind, atand, atan2d

PROGRAM test_bf16_trig_deg
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
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
        call check_bf16_real64('cosd_0',   &
            cosd(BF16(0.0_real32)), cos(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('cosd_90',  &
            cosd(BF16(90.0_real32)), cos(90.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
        call check_bf16_real64('cosd_60',  &
            cosd(BF16(60.0_real32)), cos(60.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
        call check_bf16_real64('cosd_180', &
            cosd(BF16(180.0_real32)), cos(180.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_sind()
        call check_bf16_real64('sind_0',   &
            sind(BF16(0.0_real32)), sin(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('sind_90',  &
            sind(BF16(90.0_real32)), sin(90.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
        call check_bf16_real64('sind_30',  &
            sind(BF16(30.0_real32)), sin(30.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_tand()
        call check_bf16_real64('tand_0',   &
            tand(BF16(0.0_real32)), tan(0.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('tand_45',  &
            tand(BF16(45.0_real32)), tan(45.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
        call check_bf16_real64('tand_30',  &
            tand(BF16(30.0_real32)), tan(30.0_real64 * DEG2RAD), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_cotan()
        ! cotan uses radians: cotan(x) = -tan(x + pi/2)
        call check_bf16_real64('cotan_1',  &
            cotan(BF16(1.0_real32)), -tan(1.0_real64 + PI/2.0_real64), BF16_TOL_LOOSE)
        call check_bf16_real64('cotan_0.5', &
            cotan(BF16(0.5_real32)), -tan(0.5_real64 + PI/2.0_real64), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_cotand()
        ! cotand uses degrees: cotand(x) = -tand(x + 90)
        call check_bf16_real64('cotand_45', &
            cotand(BF16(45.0_real32)), &
            -tan((45.0_real64 + 90.0_real64) * DEG2RAD), BF16_TOL_LOOSE)
        call check_bf16_real64('cotand_90', &
            cotand(BF16(90.0_real32)), &
            -tan((90.0_real64 + 90.0_real64) * DEG2RAD), BF16_TOL_LOOSE)
    end subroutine

    subroutine test_acosd()
        call check_bf16_real64('acosd_0',   &
            acosd(BF16(0.0_real32)), acos(0.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('acosd_1',   &
            acosd(BF16(1.0_real32)), acos(1.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('acosd_0.5', &
            acosd(BF16(0.5_real32)), acos(0.5_real64) * RAD2DEG, BF16_TOL_LOOSE)
    end subroutine

    subroutine test_asind()
        call check_bf16_real64('asind_0',   &
            asind(BF16(0.0_real32)), asin(0.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('asind_1',   &
            asind(BF16(1.0_real32)), asin(1.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('asind_0.5', &
            asind(BF16(0.5_real32)), asin(0.5_real64) * RAD2DEG, BF16_TOL_LOOSE)
    end subroutine

    subroutine test_atand()
        call check_bf16_real64('atand_0',   &
            atand(BF16(0.0_real32)), atan(0.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('atand_1',   &
            atand(BF16(1.0_real32)), atan(1.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('atand_neg1', &
            atand(BF16(-1.0_real32)), atan(-1.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
    end subroutine

    subroutine test_atan2d()
        call check_bf16_real64('atan2d_1_1', &
            atan2d(BF16(1.0_real32), BF16(1.0_real32)), &
            atan2(1.0_real64, 1.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
        call check_bf16_real64('atan2d_0_1', &
            atan2d(BF16(0.0_real32), BF16(1.0_real32)), &
            atan2(0.0_real64, 1.0_real64) * RAD2DEG, BF16_TOL_LOOSE)
    end subroutine

END PROGRAM test_bf16_trig_deg