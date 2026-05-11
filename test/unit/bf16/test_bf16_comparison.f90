!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 comparison operators.
!
!  Tests: .lt., .le., .gt., .ge., .eq., .ne. with BF16-BF16, BF16-fp32, fp32-BF16

PROGRAM test_bf16_comparison
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
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
        call check_logical('lt_bf16_bf16_true',   BF16(1.0_real32) .lt. BF16(2.0_real32), .TRUE.)
        call check_logical('lt_bf16_bf16_false',  BF16(2.0_real32) .lt. BF16(1.0_real32), .FALSE.)
        call check_logical('lt_bf16_fp32',        BF16(1.0_real32) .lt. 2.0_real32, .TRUE.)
        call check_logical('lt_fp32_bf16',        1.0_real32 .lt. BF16(2.0_real32), .TRUE.)
        call check_logical('lt_negative',         BF16(-1.0_real32) .lt. BF16(0.0_real32), .TRUE.)
        call check_logical('lt_equal',            BF16(1.0_real32) .lt. BF16(1.0_real32), .FALSE.)
    end subroutine

    subroutine test_le()
        call check_logical('le_bf16_bf16_true',   BF16(1.0_real32) .le. BF16(2.0_real32), .TRUE.)
        call check_logical('le_equal',            BF16(1.0_real32) .le. BF16(1.0_real32), .TRUE.)
        call check_logical('le_bf16_fp32',        BF16(1.0_real32) .le. 2.0_real32, .TRUE.)
        call check_logical('le_fp32_bf16',        1.0_real32 .le. BF16(2.0_real32), .TRUE.)
        call check_logical('le_false',            BF16(3.0_real32) .le. BF16(1.0_real32), .FALSE.)
    end subroutine

    subroutine test_gt()
        call check_logical('gt_bf16_bf16_true',   BF16(3.0_real32) .gt. BF16(2.0_real32), .TRUE.)
        call check_logical('gt_bf16_bf16_false',  BF16(1.0_real32) .gt. BF16(2.0_real32), .FALSE.)
        call check_logical('gt_bf16_fp32',        BF16(3.0_real32) .gt. 2.0_real32, .TRUE.)
        call check_logical('gt_fp32_bf16',        3.0_real32 .gt. BF16(2.0_real32), .TRUE.)
        call check_logical('gt_negative',         BF16(0.0_real32) .gt. BF16(-1.0_real32), .TRUE.)
    end subroutine

    subroutine test_ge()
        call check_logical('ge_bf16_bf16_true',   BF16(3.0_real32) .ge. BF16(2.0_real32), .TRUE.)
        call check_logical('ge_equal',            BF16(1.0_real32) .ge. BF16(1.0_real32), .TRUE.)
        call check_logical('ge_bf16_fp32',        BF16(3.0_real32) .ge. 2.0_real32, .TRUE.)
        call check_logical('ge_fp32_bf16',        3.0_real32 .ge. BF16(2.0_real32), .TRUE.)
        call check_logical('ge_false',            BF16(1.0_real32) .ge. BF16(2.0_real32), .FALSE.)
    end subroutine

    subroutine test_eq()
        call check_logical('eq_bf16_bf16',  BF16(5.0_real32) .eq. BF16(5.0_real32), .TRUE.)
        call check_logical('eq_bf16_fp32',  BF16(5.0_real32) .eq. 5.0_real32, .TRUE.)
        call check_logical('eq_fp32_bf16',  5.0_real32 .eq. BF16(5.0_real32), .TRUE.)
        call check_logical('eq_zero',       BF16(0.0_real32) .eq. BF16(0.0_real32), .TRUE.)
        call check_logical('eq_not_equal',  BF16(1.0_real32) .eq. BF16(2.0_real32), .FALSE.)
        call check_logical('eq_negative',   BF16(-1.0_real32) .eq. BF16(-1.0_real32), .TRUE.)
    end subroutine

    subroutine test_ne()
        call check_logical('ne_bf16_bf16',  BF16(1.0_real32) .ne. BF16(2.0_real32), .TRUE.)
        call check_logical('ne_bf16_fp32',  BF16(1.0_real32) .ne. 2.0_real32, .TRUE.)
        call check_logical('ne_fp32_bf16',  1.0_real32 .ne. BF16(2.0_real32), .TRUE.)
        call check_logical('ne_equal',      BF16(5.0_real32) .ne. BF16(5.0_real32), .FALSE.)
        call check_logical('ne_zero',       BF16(0.0_real32) .ne. BF16(0.0_real32), .FALSE.)
    end subroutine

END PROGRAM test_bf16_comparison