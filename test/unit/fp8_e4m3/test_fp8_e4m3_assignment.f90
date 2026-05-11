!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 assignment.
!
!  Tests: assignment from int32, int64, real32, real64, same-type;
!         assignment to real32; chained assignments

PROGRAM test_fp8_e4m3_assignment
    USE iso_fortran_env, only: real32, real64, int32, int64
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
    IMPLICIT NONE

    CALL test_assign_int32()
    CALL test_assign_int64()
    CALL test_assign_real32()
    CALL test_assign_real64()
    CALL test_assign_same_type()
    CALL test_assign_to_real32()
    CALL test_assign_chained()

    CALL test_summary()

CONTAINS

    subroutine test_assign_int32()
        type(FP8_E4M3) :: x
        x = 5_int32
        call check_fp8_e4m3_real64('assign_int32_5', x, 5.0_real64, FP8_E4M3_TOL_TIGHT)
        x = -3_int32
        call check_fp8_e4m3_real64('assign_int32_neg3', x, -3.0_real64, FP8_E4M3_TOL_TIGHT)
        x = 0_int32
        call check_fp8_e4m3_real64('assign_int32_0', x, 0.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_assign_int64()
        type(FP8_E4M3) :: x
        x = 10_int64
        call check_fp8_e4m3_real64('assign_int64_10', x, 10.0_real64, FP8_E4M3_TOL_TIGHT)
        x = -7_int64
        call check_fp8_e4m3_real64('assign_int64_neg7', x, -7.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_assign_real32()
        type(FP8_E4M3) :: x
        x = 2.0_real32
        call check_fp8_e4m3_real64('assign_real32_2', x, 2.0_real64, FP8_E4M3_TOL_TIGHT)
        x = -1.5_real32
        call check_fp8_e4m3_real64('assign_real32_neg1.5', x, -1.5_real64, FP8_E4M3_TOL)
        x = 0.0_real32
        call check_fp8_e4m3_real64('assign_real32_0', x, 0.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_assign_real64()
        type(FP8_E4M3) :: x
        x = 3.0_real64
        call check_fp8_e4m3_real64('assign_real64_3', x, 3.0_real64, FP8_E4M3_TOL_TIGHT)
        x = -0.5_real64
        call check_fp8_e4m3_real64('assign_real64_neg0.5', x, -0.5_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_assign_same_type()
        type(FP8_E4M3) :: x, y
        x = FP8_E4M3(7.0_real32)
        y = x
        call check_fp8_e4m3_real64('assign_same_type', y, 7.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_assign_to_real32()
        type(FP8_E4M3) :: x
        real(real32) :: r
        x = FP8_E4M3(4.0_real32)
        r = x
        call check_fp8_e4m3_real64('assign_to_real32', FP8_E4M3(r), 4.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_assign_chained()
        type(FP8_E4M3) :: x, y, z
        x = FP8_E4M3(6.0_real32)
        y = x
        z = y
        call check_fp8_e4m3_real64('assign_chained', z, 6.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp8_e4m3_assignment