!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 assignment operators.
!
!  Tests: assign_int32, assign_int64, assign_real, assign_double,
!         assign_bf16, assign_to_real

PROGRAM test_bf16_assignment
    USE iso_fortran_env, only: real32, real64, int32, int64
    USE LPF_BF16
    USE lpf_bf16_test_utils
    IMPLICIT NONE

    CALL test_assign_int32()
    CALL test_assign_int64()
    CALL test_assign_real()
    CALL test_assign_double()
    CALL test_assign_bf16()
    CALL test_assign_to_real()
    CALL test_assign_chain()

    CALL test_summary()

CONTAINS

    subroutine test_assign_int32()
        type(BF16) :: x
        x = 42
        call check_bf16_real64('assign_int32_42', x, 42.0_real64, BF16_TOL_TIGHT)

        x = 0
        call check_bf16_real64('assign_int32_0', x, 0.0_real64, BF16_TOL_TIGHT)

        x = -5
        call check_bf16_real64('assign_int32_neg5', x, -5.0_real64, BF16_TOL_TIGHT)

        x = 1000
        call check_bf16_real64('assign_int32_1000', x, 1000.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_assign_int64()
        type(BF16) :: x
        x = 42_int64
        call check_bf16_real64('assign_int64_42', x, 42.0_real64, BF16_TOL_TIGHT)

        x = 0_int64
        call check_bf16_real64('assign_int64_0', x, 0.0_real64, BF16_TOL_TIGHT)

        x = -5_int64
        call check_bf16_real64('assign_int64_neg5', x, -5.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_assign_real()
        type(BF16) :: x
        x = 3.14
        call check_bf16_real64('assign_real32_pi', x, 3.14_real64, BF16_TOL)

        x = 0.0
        call check_bf16_real64('assign_real32_zero', x, 0.0_real64, BF16_TOL_TIGHT)

        x = -1.0
        call check_bf16_real64('assign_real32_neg1', x, -1.0_real64, BF16_TOL_TIGHT)

        x = 100.0
        call check_bf16_real64('assign_real32_100', x, 100.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_assign_double()
        type(BF16) :: x
        x = 3.14d0
        call check_bf16_real64('assign_real64_pi', x, 3.14_real64, BF16_TOL)

        x = 0.0d0
        call check_bf16_real64('assign_real64_zero', x, 0.0_real64, BF16_TOL_TIGHT)

        x = -1.0d0
        call check_bf16_real64('assign_real64_neg1', x, -1.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_assign_bf16()
        type(BF16) :: x, y
        y = BF16(7.5_real32)
        x = y
        call check_bf16_real64('assign_bf16_7.5', x, 7.5_real64, BF16_TOL_TIGHT)

        y = BF16(-2.0_real32)
        x = y
        call check_bf16_real64('assign_bf16_neg2', x, -2.0_real64, BF16_TOL_TIGHT)
    end subroutine

    subroutine test_assign_to_real()
        type(BF16) :: x
        real(real32) :: r32
        x = BF16(4.2_real32)
        r32 = x
        call check_bf16_real64('assign_to_real32', BF16(r32), 4.2_real64, BF16_TOL)
    end subroutine

    subroutine test_assign_chain()
        type(BF16) :: x
        x = 42
        call check_bf16_real64('chain_int_to_bf16', x, 42.0_real64, BF16_TOL_TIGHT)

        x = 3.14
        call check_bf16_real64('chain_real_to_bf16', x, 3.14_real64, BF16_TOL)

        x = BF16(2.0_real32)
        call check_bf16_real64('chain_bf16_overwrite', x, 2.0_real64, BF16_TOL_TIGHT)
    end subroutine

END PROGRAM test_bf16_assignment