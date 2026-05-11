!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 constructors and conversions.

PROGRAM test_fp16_constructors
    USE iso_fortran_env, only: real32, real64, int32, int64
    USE LPF_FP16
    USE lpf_fp16_test_utils
    IMPLICIT NONE

    CALL test_constructors_real32()
    CALL test_constructors_real64()
    CALL test_constructors_int32()
    CALL test_constructors_int64()
    CALL test_conversions()

    CALL test_summary()

CONTAINS

    subroutine test_constructors_real32()
        call check_fp16_real64('ctor_real32_zero', FP16(0.0_real32), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real32_one',  FP16(1.0_real32), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real32_neg',  FP16(-1.0_real32), -1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real32_half',  FP16(0.5_real32), 0.5_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real32_pi',   FP16(3.14159_real32), 3.14159_real64, FP16_TOL)
        call check_fp16_real64('ctor_real32_two',  FP16(2.0_real32), 2.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real32_small', FP16(6.104e-5_real32), 6.104e-5_real64, FP16_TOL_LOOSE)
        call check_fp16_real64('ctor_real32_large', FP16(65504.0_real32), 65504.0_real64, FP16_TOL_LOOSE)
    end subroutine

    subroutine test_constructors_real64()
        call check_fp16_real64('ctor_real64_zero',  FP16(0.0_real64), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real64_one',   FP16(1.0_real64), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real64_neg',   FP16(-1.0_real64), -1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real64_half',  FP16(0.5_real64), 0.5_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_real64_third', FP16(1.0_real64/3.0_real64), 1.0_real64/3.0_real64, FP16_TOL)
    end subroutine

    subroutine test_constructors_int32()
        call check_fp16_real64('ctor_int32_zero',  FP16(0_int32), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int32_one',   FP16(1_int32), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int32_neg',   FP16(-5_int32), -5.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int32_42',    FP16(42_int32), 42.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int32_100',   FP16(100_int32), 100.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_constructors_int64()
        call check_fp16_real64('ctor_int64_zero',  FP16(0_int64), 0.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int64_one',   FP16(1_int64), 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int64_neg',   FP16(-5_int64), -5.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('ctor_int64_42',    FP16(42_int64), 42.0_real64, FP16_TOL_TIGHT)
    end subroutine

    subroutine test_conversions()
        type(FP16) :: x
        real(real32) :: r32
        real(real64) :: r64

        x = FP16(2.5_real32)
        r32 = real(x)
        call check_fp16_real64('real_conv_2.5', FP16(r32), 2.5_real64, FP16_TOL_TIGHT)

        x = FP16(2.5_real32)
        r64 = dble(x)
        call check_fp16_real64('dble_conv_2.5', FP16(real(r64, kind=real32)), 2.5_real64, FP16_TOL_TIGHT)

        x = FP16(1.0_real32)
        r32 = real(x)
        call check_fp16_real64('roundtrip_real32', FP16(r32), 1.0_real64, FP16_TOL_TIGHT)

        x = FP16(-3.0_real32)
        r64 = dble(x)
        call check_fp16_real64('roundtrip_real64', FP16(real(r64, kind=real32)), -3.0_real64, FP16_TOL_TIGHT)

        x = FP16(1.0_real32)
        call check_fp16_real64('ctor_consistency', FP16(1.0_real64), dble(x), FP16_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp16_constructors