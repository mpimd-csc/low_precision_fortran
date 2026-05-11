!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 two-argument math functions.
!
!  Tests: atan2, bessel_jn, bessel_yn, hypot

PROGRAM test_fp8_e4m3_math_two_arg
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
    IMPLICIT NONE

    CALL test_atan2()
    CALL test_bessel_jn()
    CALL test_bessel_yn()
    CALL test_hypot()

    CALL test_summary()

CONTAINS

    subroutine test_atan2()
        call check_fp8_e4m3_real64('atan2_1_1',   &
            atan2(FP8_E4M3(1.0_real32), FP8_E4M3(1.0_real32)), atan2(1.0_real64, 1.0_real64), &
            FP8_E4M3_TOL_LOOSE)
        call check_fp8_e4m3_real64('atan2_0_1',   &
            atan2(FP8_E4M3(0.0_real32), FP8_E4M3(1.0_real32)), atan2(0.0_real64, 1.0_real64), &
            FP8_E4M3_TOL_LOOSE)
        call check_fp8_e4m3_real64('atan2_neg1_1', &
            atan2(FP8_E4M3(-1.0_real32), FP8_E4M3(1.0_real32)), atan2(-1.0_real64, 1.0_real64), &
            FP8_E4M3_TOL_LOOSE)
        call check_fp8_e4m3_real64('atan2_1_0',   &
            atan2(FP8_E4M3(1.0_real32), FP8_E4M3(0.0_real32)), atan2(1.0_real64, 0.0_real64), &
            FP8_E4M3_TOL_LOOSE)
    end subroutine

    subroutine test_bessel_jn()
        ! Known bug: bessel_jn returns 0.0 for FP8 inputs (float round-trip issue)
        call check_fp8_e4m3_real64('bessel_jn_0_1', &
            bessel_jn(0, FP8_E4M3(1.0_real32)), &
            real(bessel_j0(1.0_real64), kind=real64), FP8_E4M3_TOL_LOOSE)
        call check_fp8_e4m3_real64('bessel_jn_1_2', &
            bessel_jn(1, FP8_E4M3(2.0_real32)), &
            real(bessel_j1(2.0_real64), kind=real64), FP8_E4M3_TOL_LOOSE)
        call check_fp8_e4m3_real64('bessel_jn_2_3', &
            bessel_jn(2, FP8_E4M3(3.0_real32)), &
            real(bessel_jn(2, 3.0_real64), kind=real64), FP8_E4M3_TOL_LOOSE)
    end subroutine

    subroutine test_bessel_yn()
        ! Known bug: bessel_yn returns 0.0 for FP8 inputs (float round-trip issue)
        call check_fp8_e4m3_real64('bessel_yn_0_1', &
            bessel_yn(0, FP8_E4M3(1.0_real32)), &
            real(bessel_y0(1.0_real64), kind=real64), FP8_E4M3_TOL_LOOSE)
        call check_fp8_e4m3_real64('bessel_yn_1_2', &
            bessel_yn(1, FP8_E4M3(2.0_real32)), &
            real(bessel_y1(2.0_real64), kind=real64), FP8_E4M3_TOL_LOOSE)
    end subroutine

    subroutine test_hypot()
        call check_fp8_e4m3_real64('hypot_3_4',   &
            hypot(FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32)), hypot(3.0_real64, 4.0_real64), &
            FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('hypot_1_0',   &
            hypot(FP8_E4M3(1.0_real32), FP8_E4M3(0.0_real32)), hypot(1.0_real64, 0.0_real64), &
            FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('hypot_0_1',   &
            hypot(FP8_E4M3(0.0_real32), FP8_E4M3(1.0_real32)), hypot(0.0_real64, 1.0_real64), &
            FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('hypot_5_12',  &
            hypot(FP8_E4M3(5.0_real32), FP8_E4M3(12.0_real32)), hypot(5.0_real64, 12.0_real64), &
            FP8_E4M3_TOL)
    end subroutine

END PROGRAM test_fp8_e4m3_math_two_arg
