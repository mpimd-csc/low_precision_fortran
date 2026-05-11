!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E5M2 constructors and conversions.
!
!  Tests: FP8_E5M2(real32), FP8_E5M2(real64), FP8_E5M2(int32), FP8_E5M2(int64),
!         real(FP8_E5M2), dble(FP8_E5M2)

PROGRAM test_fp8_e5m2_constructors
    USE iso_fortran_env, only: real32, real64, int32, int64
    USE LPF_FP8_E5M2
    USE lpf_fp8_e5m2_test_utils
    IMPLICIT NONE

    CALL test_constructors_real32()
    CALL test_constructors_real64()
    CALL test_constructors_int32()
    CALL test_constructors_int64()
    CALL test_conversions()

    CALL test_summary()

CONTAINS

    subroutine test_constructors_real32()
        call check_fp8_e5m2_real64('ctor_real32_zero', FP8_E5M2(0.0_real32), 0.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real32_one',  FP8_E5M2(1.0_real32), 1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real32_neg',  FP8_E5M2(-1.0_real32), -1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real32_half',  FP8_E5M2(0.5_real32), 0.5_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real32_two',  FP8_E5M2(2.0_real32), 2.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real32_small', FP8_E5M2(0.002_real32), 0.002_real64, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('ctor_real32_large', FP8_E5M2(448.0_real32), 448.0_real64, FP8_E5M2_TOL_LOOSE)
        call check_fp8_e5m2_real64('ctor_real32_neg_large', FP8_E5M2(-448.0_real32), -448.0_real64, FP8_E5M2_TOL_LOOSE)
    end subroutine

    subroutine test_constructors_real64()
        call check_fp8_e5m2_real64('ctor_real64_zero',  FP8_E5M2(0.0_real64), 0.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real64_one',   FP8_E5M2(1.0_real64), 1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real64_neg',   FP8_E5M2(-1.0_real64), -1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real64_half',  FP8_E5M2(0.5_real64), 0.5_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_real64_third', FP8_E5M2(1.0_real64/3.0_real64), 1.0_real64/3.0_real64, FP8_E5M2_TOL)
    end subroutine

    subroutine test_constructors_int32()
        call check_fp8_e5m2_real64('ctor_int32_zero',  FP8_E5M2(0_int32), 0.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int32_one',   FP8_E5M2(1_int32), 1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int32_neg',   FP8_E5M2(-5_int32), -5.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int32_42',    FP8_E5M2(42_int32), 42.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int32_100',   FP8_E5M2(100_int32), 100.0_real64, FP8_E5M2_TOL)
        call check_fp8_e5m2_real64('ctor_int32_7',     FP8_E5M2(7_int32), 7.0_real64, FP8_E5M2_TOL_TIGHT)
    end subroutine

    subroutine test_constructors_int64()
        call check_fp8_e5m2_real64('ctor_int64_zero',  FP8_E5M2(0_int64), 0.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int64_one',   FP8_E5M2(1_int64), 1.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int64_neg',   FP8_E5M2(-5_int64), -5.0_real64, FP8_E5M2_TOL_TIGHT)
        call check_fp8_e5m2_real64('ctor_int64_42',    FP8_E5M2(42_int64), 42.0_real64, FP8_E5M2_TOL_TIGHT)
    end subroutine

    subroutine test_conversions()
        type(FP8_E5M2) :: x
        real(real32) :: r32
        real(real64) :: r64

        ! real(FP8_E5M2) conversion
        x = FP8_E5M2(2.0_real32)
        r32 = real(x)
        call check_fp8_e5m2_real64('real_conv_2.0', FP8_E5M2(r32), 2.0_real64, FP8_E5M2_TOL_TIGHT)

        ! dble(FP8_E5M2) conversion
        x = FP8_E5M2(2.0_real32)
        r64 = dble(x)
        call check_fp8_e5m2_real64('dble_conv_2.0', FP8_E5M2(real(r64, kind=real32)), 2.0_real64, FP8_E5M2_TOL_TIGHT)

        ! Round-trip: FP8_E5M2 -> real32 -> FP8_E5M2
        x = FP8_E5M2(1.0_real32)
        r32 = real(x)
        call check_fp8_e5m2_real64('roundtrip_real32', FP8_E5M2(r32), 1.0_real64, FP8_E5M2_TOL_TIGHT)

        ! Round-trip: FP8_E5M2 -> real64 -> FP8_E5M2
        x = FP8_E5M2(-3.0_real32)
        r64 = dble(x)
        call check_fp8_e5m2_real64('roundtrip_real64', FP8_E5M2(real(r64, kind=real32)), -3.0_real64, FP8_E5M2_TOL_TIGHT)

        ! Consistency: same value via real32 and real64
        x = FP8_E5M2(1.0_real32)
        call check_fp8_e5m2_real64('ctor_consistency', FP8_E5M2(1.0_real64), dble(x), FP8_E5M2_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp8_e5m2_constructors