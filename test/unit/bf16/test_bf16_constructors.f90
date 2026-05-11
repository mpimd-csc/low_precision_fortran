!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 constructors and conversions.
!
!  Tests: BF16(real32), BF16(real64), BF16(int32), BF16(int64),
!         real(BF16), dble(BF16)

PROGRAM test_bf16_constructors
    USE iso_fortran_env, only: real32, real64, int32, int64
    USE LPF_BF16
    USE lpf_bf16_test_utils
    IMPLICIT NONE

    CALL test_constructors_real32()
    CALL test_constructors_real64()
    CALL test_constructors_int32()
    CALL test_constructors_int64()
    CALL test_conversions()

    CALL test_summary()

CONTAINS

    subroutine test_constructors_real32()
        call check_bf16_real64('ctor_real32_zero', BF16(0.0_real32), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real32_one',  BF16(1.0_real32), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real32_neg',  BF16(-1.0_real32), -1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real32_half',  BF16(0.5_real32), 0.5_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real32_pi',   BF16(3.14159_real32), 3.14159_real64, BF16_TOL)
        call check_bf16_real64('ctor_real32_two',  BF16(2.0_real32), 2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real32_small', BF16(1.175e-38_real32), 1.175e-38_real64, BF16_TOL_LOOSE)
        call check_bf16_real64('ctor_real32_large', BF16(3.389e38_real32), 3.389e38_real64, BF16_TOL_LOOSE)
        call check_bf16_real64('ctor_real32_neg_large', BF16(-3.389e38_real32), -3.389e38_real64, BF16_TOL_LOOSE)
    end subroutine

    subroutine test_constructors_real64()
        call check_bf16_real64('ctor_real64_zero',  BF16(0.0_real64), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real64_one',   BF16(1.0_real64), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real64_neg',   BF16(-1.0_real64), -1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real64_half',  BF16(0.5_real64), 0.5_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_real64_third', BF16(1.0_real64/3.0_real64), 1.0_real64/3.0_real64, BF16_TOL)
        call check_bf16_real64('ctor_real64_pi',    BF16(3.141592653589793_real64), 3.141592653589793_real64, BF16_TOL)
    end subroutine

    subroutine test_constructors_int32()
        call check_bf16_real64('ctor_int32_zero',  BF16(0_int32), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int32_one',   BF16(1_int32), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int32_neg',   BF16(-5_int32), -5.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int32_42',    BF16(42_int32), 42.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int32_1000',  BF16(1000_int32), 1000.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int32_127',   BF16(127_int32), 127.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_constructors_int64()
        call check_bf16_real64('ctor_int64_zero',  BF16(0_int64), 0.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int64_one',   BF16(1_int64), 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int64_neg',   BF16(-5_int64), -5.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('ctor_int64_42',    BF16(42_int64), 42.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_conversions()
        type(BF16) :: x
        real(real32) :: r32
        real(real64) :: r64

        ! real(BF16) conversion
        x = BF16(2.5_real32)
        r32 = real(x)
        call check_bf16_real64('real_conv_2.5', BF16(r32), 2.5_real64, BF16_TOL_TIGHT)

        ! dble(BF16) conversion
        x = BF16(2.5_real32)
        r64 = dble(x)
        call check_bf16_real64('dble_conv_2.5', BF16(real(r64, kind=real32)), 2.5_real64, BF16_TOL_TIGHT)

        ! Round-trip: BF16 -> real32 -> BF16
        x = BF16(1.0_real32)
        r32 = real(x)
        call check_bf16_real64('roundtrip_real32', BF16(r32), 1.0_real64, BF16_TOL_TIGHT)

        ! Round-trip: BF16 -> real64 -> BF16
        x = BF16(-3.0_real32)
        r64 = dble(x)
        call check_bf16_real64('roundtrip_real64', BF16(real(r64, kind=real32)), -3.0_real64, BF16_TOL_TIGHT)

        ! Consistency: same value via real32 and real64
        x = BF16(1.0_real32)
        call check_bf16_real64('ctor_consistency', BF16(1.0_real64), dble(x), BF16_TOL_TIGHT)
    end subroutine

END PROGRAM test_bf16_constructors