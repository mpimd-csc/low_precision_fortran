!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 comparison operators.

PROGRAM test_fp16_comparison
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
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
        call check_logical('lt_fp16_fp16_true',   FP16(1.0_real32) .lt. FP16(2.0_real32), .TRUE.)
        call check_logical('lt_fp16_fp16_false',  FP16(2.0_real32) .lt. FP16(1.0_real32), .FALSE.)
        call check_logical('lt_fp16_fp32',        FP16(1.0_real32) .lt. 2.0_real32, .TRUE.)
        call check_logical('lt_fp32_fp16',        1.0_real32 .lt. FP16(2.0_real32), .TRUE.)
        call check_logical('lt_negative',         FP16(-1.0_real32) .lt. FP16(0.0_real32), .TRUE.)
        call check_logical('lt_equal',            FP16(1.0_real32) .lt. FP16(1.0_real32), .FALSE.)
    end subroutine

    subroutine test_le()
        call check_logical('le_fp16_fp16_true',   FP16(1.0_real32) .le. FP16(2.0_real32), .TRUE.)
        call check_logical('le_equal',            FP16(1.0_real32) .le. FP16(1.0_real32), .TRUE.)
        call check_logical('le_fp16_fp32',        FP16(1.0_real32) .le. 2.0_real32, .TRUE.)
        call check_logical('le_fp32_fp16',        1.0_real32 .le. FP16(2.0_real32), .TRUE.)
        call check_logical('le_false',            FP16(3.0_real32) .le. FP16(1.0_real32), .FALSE.)
    end subroutine

    subroutine test_gt()
        call check_logical('gt_fp16_fp16_true',   FP16(3.0_real32) .gt. FP16(2.0_real32), .TRUE.)
        call check_logical('gt_fp16_fp16_false',  FP16(1.0_real32) .gt. FP16(2.0_real32), .FALSE.)
        call check_logical('gt_fp16_fp32',        FP16(3.0_real32) .gt. 2.0_real32, .TRUE.)
        call check_logical('gt_fp32_fp16',        3.0_real32 .gt. FP16(2.0_real32), .TRUE.)
        call check_logical('gt_negative',         FP16(0.0_real32) .gt. FP16(-1.0_real32), .TRUE.)
    end subroutine

    subroutine test_ge()
        call check_logical('ge_fp16_fp16_true',   FP16(3.0_real32) .ge. FP16(2.0_real32), .TRUE.)
        call check_logical('ge_equal',            FP16(1.0_real32) .ge. FP16(1.0_real32), .TRUE.)
        call check_logical('ge_fp16_fp32',        FP16(3.0_real32) .ge. 2.0_real32, .TRUE.)
        call check_logical('ge_fp32_fp16',        3.0_real32 .ge. FP16(2.0_real32), .TRUE.)
        call check_logical('ge_false',            FP16(1.0_real32) .ge. FP16(2.0_real32), .FALSE.)
    end subroutine

    subroutine test_eq()
        call check_logical('eq_fp16_fp16',  FP16(5.0_real32) .eq. FP16(5.0_real32), .TRUE.)
        call check_logical('eq_fp16_fp32',  FP16(5.0_real32) .eq. 5.0_real32, .TRUE.)
        call check_logical('eq_fp32_fp16',  5.0_real32 .eq. FP16(5.0_real32), .TRUE.)
        call check_logical('eq_zero',       FP16(0.0_real32) .eq. FP16(0.0_real32), .TRUE.)
        call check_logical('eq_not_equal',  FP16(1.0_real32) .eq. FP16(2.0_real32), .FALSE.)
    end subroutine

    subroutine test_ne()
        call check_logical('ne_fp16_fp16',  FP16(1.0_real32) .ne. FP16(2.0_real32), .TRUE.)
        call check_logical('ne_fp16_fp32',  FP16(1.0_real32) .ne. 2.0_real32, .TRUE.)
        call check_logical('ne_fp32_fp16',  1.0_real32 .ne. FP16(2.0_real32), .TRUE.)
        call check_logical('ne_equal',      FP16(5.0_real32) .ne. FP16(5.0_real32), .FALSE.)
        call check_logical('ne_zero',       FP16(0.0_real32) .ne. FP16(0.0_real32), .FALSE.)
    end subroutine

END PROGRAM test_fp16_comparison