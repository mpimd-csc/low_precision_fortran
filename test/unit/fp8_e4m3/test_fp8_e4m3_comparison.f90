!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 comparison operators.
!
!  Tests: .lt., .le., .gt., .ge., .eq., .ne. for FP8_E4M3-FP8_E4M3,
!         FP8_E4M3-real32, real32-FP8_E4M3 variants

PROGRAM test_fp8_e4m3_comparison
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
    IMPLICIT NONE

    CALL test_lt()
    CALL test_le()
    CALL test_gt()
    CALL test_ge()
    CALL test_eq()
    CALL test_ne()

    CALL test_summary()

CONTAINS

    subroutine test_lt()
        ! FP8_E4M3-FP8_E4M3
        call check_logical('lt_fp8_fp8_2_3',    FP8_E4M3(2.0_real32) .lt. FP8_E4M3(3.0_real32), .TRUE.)
        call check_logical('lt_fp8_fp8_3_2',    FP8_E4M3(3.0_real32) .lt. FP8_E4M3(2.0_real32), .FALSE.)
        call check_logical('lt_fp8_fp8_eq',     FP8_E4M3(2.0_real32) .lt. FP8_E4M3(2.0_real32), .FALSE.)
        call check_logical('lt_fp8_fp8_neg',    FP8_E4M3(-1.0_real32) .lt. FP8_E4M3(1.0_real32), .TRUE.)

        ! FP8_E4M3-real32
        call check_logical('lt_fp8_real_2_3',   FP8_E4M3(2.0_real32) .lt. 3.0_real32, .TRUE.)
        call check_logical('lt_fp8_real_3_2',   FP8_E4M3(3.0_real32) .lt. 2.0_real32, .FALSE.)

        ! real32-FP8_E4M3
        call check_logical('lt_real_fp8_2_3',   2.0_real32 .lt. FP8_E4M3(3.0_real32), .TRUE.)
        call check_logical('lt_real_fp8_3_2',   3.0_real32 .lt. FP8_E4M3(2.0_real32), .FALSE.)
    end subroutine

    subroutine test_le()
        call check_logical('le_fp8_fp8_2_3',    FP8_E4M3(2.0_real32) .le. FP8_E4M3(3.0_real32), .TRUE.)
        call check_logical('le_fp8_fp8_3_2',    FP8_E4M3(3.0_real32) .le. FP8_E4M3(2.0_real32), .FALSE.)
        call check_logical('le_fp8_fp8_eq',     FP8_E4M3(2.0_real32) .le. FP8_E4M3(2.0_real32), .TRUE.)
        call check_logical('le_fp8_real_2_2',   FP8_E4M3(2.0_real32) .le. 2.0_real32, .TRUE.)
        call check_logical('le_real_fp8_2_2',   2.0_real32 .le. FP8_E4M3(2.0_real32), .TRUE.)
    end subroutine

    subroutine test_gt()
        call check_logical('gt_fp8_fp8_3_2',    FP8_E4M3(3.0_real32) .gt. FP8_E4M3(2.0_real32), .TRUE.)
        call check_logical('gt_fp8_fp8_2_3',    FP8_E4M3(2.0_real32) .gt. FP8_E4M3(3.0_real32), .FALSE.)
        call check_logical('gt_fp8_fp8_eq',     FP8_E4M3(2.0_real32) .gt. FP8_E4M3(2.0_real32), .FALSE.)
        call check_logical('gt_fp8_real_3_2',   FP8_E4M3(3.0_real32) .gt. 2.0_real32, .TRUE.)
        call check_logical('gt_real_fp8_3_2',   3.0_real32 .gt. FP8_E4M3(2.0_real32), .TRUE.)
    end subroutine

    subroutine test_ge()
        call check_logical('ge_fp8_fp8_3_2',    FP8_E4M3(3.0_real32) .ge. FP8_E4M3(2.0_real32), .TRUE.)
        call check_logical('ge_fp8_fp8_2_3',    FP8_E4M3(2.0_real32) .ge. FP8_E4M3(3.0_real32), .FALSE.)
        call check_logical('ge_fp8_fp8_eq',     FP8_E4M3(2.0_real32) .ge. FP8_E4M3(2.0_real32), .TRUE.)
        call check_logical('ge_fp8_real_2_2',   FP8_E4M3(2.0_real32) .ge. 2.0_real32, .TRUE.)
        call check_logical('ge_real_fp8_2_2',   2.0_real32 .ge. FP8_E4M3(2.0_real32), .TRUE.)
    end subroutine

    subroutine test_eq()
        call check_logical('eq_fp8_fp8_eq',     FP8_E4M3(2.0_real32) .eq. FP8_E4M3(2.0_real32), .TRUE.)
        call check_logical('eq_fp8_fp8_neq',    FP8_E4M3(2.0_real32) .eq. FP8_E4M3(3.0_real32), .FALSE.)
        call check_logical('eq_fp8_real_eq',    FP8_E4M3(2.0_real32) .eq. 2.0_real32, .TRUE.)
        call check_logical('eq_real_fp8_eq',    2.0_real32 .eq. FP8_E4M3(2.0_real32), .TRUE.)
        call check_logical('eq_zero',           FP8_E4M3(0.0_real32) .eq. FP8_E4M3(0.0_real32), .TRUE.)
    end subroutine

    subroutine test_ne()
        call check_logical('ne_fp8_fp8_neq',    FP8_E4M3(2.0_real32) .ne. FP8_E4M3(3.0_real32), .TRUE.)
        call check_logical('ne_fp8_fp8_eq',     FP8_E4M3(2.0_real32) .ne. FP8_E4M3(2.0_real32), .FALSE.)
        call check_logical('ne_fp8_real_neq',   FP8_E4M3(2.0_real32) .ne. 3.0_real32, .TRUE.)
        call check_logical('ne_real_fp8_neq',   2.0_real32 .ne. FP8_E4M3(3.0_real32), .TRUE.)
    end subroutine

END PROGRAM test_fp8_e4m3_comparison