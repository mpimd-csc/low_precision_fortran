!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 degree-based trigonometric functions.

PROGRAM test_fp16_trig_deg
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
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
        call check_fp16_real64('cosd_0',   &
            cosd(FP16(0.0_real32)), cos(0.0_real64), FP16_TOL_LOOSE)
        call check_fp16_real64('cosd_90',  &
            cosd(FP16(90.0_real32)), cos(90.0_real64 * DEG2RAD), FP16_TOL_LOOSE)
        call check_fp16_real64('cosd_60',  &
            cosd(FP16(60.0_real32)), cos(60.0_real64 * DEG2RAD), FP16_TOL_LOOSE)
    end subroutine

    subroutine test_sind()
        call check_fp16_real64('sind_0',   &
            sind(FP16(0.0_real32)), sin(0.0_real64), FP16_TOL_LOOSE)
        call check_fp16_real64('sind_90',  &
            sind(FP16(90.0_real32)), sin(90.0_real64 * DEG2RAD), FP16_TOL_LOOSE)
        call check_fp16_real64('sind_30',  &
            sind(FP16(30.0_real32)), sin(30.0_real64 * DEG2RAD), FP16_TOL_LOOSE)
    end subroutine

    subroutine test_tand()
        call check_fp16_real64('tand_0',   &
            tand(FP16(0.0_real32)), tan(0.0_real64), FP16_TOL_LOOSE)
        call check_fp16_real64('tand_45',  &
            tand(FP16(45.0_real32)), tan(45.0_real64 * DEG2RAD), FP16_TOL_LOOSE)
    end subroutine

    subroutine test_cotan()
        call check_fp16_real64('cotan_1',  &
            cotan(FP16(1.0_real32)), -tan(1.0_real64 + PI/2.0_real64), FP16_TOL_LOOSE)
        call check_fp16_real64('cotan_0.5', &
            cotan(FP16(0.5_real32)), -tan(0.5_real64 + PI/2.0_real64), FP16_TOL_LOOSE)
    end subroutine

    subroutine test_cotand()
        call check_fp16_real64('cotand_45', &
            cotand(FP16(45.0_real32)), &
            -tan((45.0_real64 + 90.0_real64) * DEG2RAD), FP16_TOL_LOOSE)
        call check_fp16_real64('cotand_90', &
            cotand(FP16(90.0_real32)), &
            -tan((90.0_real64 + 90.0_real64) * DEG2RAD), FP16_TOL_LOOSE)
    end subroutine

    subroutine test_acosd()
        call check_fp16_real64('acosd_0',   &
            acosd(FP16(0.0_real32)), acos(0.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('acosd_1',   &
            acosd(FP16(1.0_real32)), acos(1.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('acosd_0.5', &
            acosd(FP16(0.5_real32)), acos(0.5_real64) * RAD2DEG, FP16_TOL_LOOSE)
    end subroutine

    subroutine test_asind()
        call check_fp16_real64('asind_0',   &
            asind(FP16(0.0_real32)), asin(0.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('asind_1',   &
            asind(FP16(1.0_real32)), asin(1.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('asind_0.5', &
            asind(FP16(0.5_real32)), asin(0.5_real64) * RAD2DEG, FP16_TOL_LOOSE)
    end subroutine

    subroutine test_atand()
        call check_fp16_real64('atand_0',   &
            atand(FP16(0.0_real32)), atan(0.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('atand_1',   &
            atand(FP16(1.0_real32)), atan(1.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('atand_neg1', &
            atand(FP16(-1.0_real32)), atan(-1.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
    end subroutine

    subroutine test_atan2d()
        call check_fp16_real64('atan2d_1_1', &
            atan2d(FP16(1.0_real32), FP16(1.0_real32)), &
            atan2(1.0_real64, 1.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
        call check_fp16_real64('atan2d_0_1', &
            atan2d(FP16(0.0_real32), FP16(1.0_real32)), &
            atan2(0.0_real64, 1.0_real64) * RAD2DEG, FP16_TOL_LOOSE)
    end subroutine

END PROGRAM test_fp16_trig_deg