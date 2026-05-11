!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP16 formatted I/O.

PROGRAM test_fp16_io
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_fp16_test_utils
    IMPLICIT NONE

    CALL test_write_listdirected()
    CALL test_read_roundtrip()
    CALL test_write_read_multiple()

    CALL test_summary()

CONTAINS

    subroutine test_write_listdirected()
        type(FP16) :: x
        character(len=40) :: buffer
        real(real64) :: val

        x = FP16(3.14_real32)
        write(buffer, '(DT)') x
        read(buffer, '(DT)') x
        val = dble(x)
        call check_fp16_real64('write_listdirected_3.14', FP16(real(val, kind=real32)), 3.14_real64, FP16_TOL)
    end subroutine

    subroutine test_read_roundtrip()
        type(FP16) :: x, y
        character(len=40) :: buffer

        x = FP16(1.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp16_real64('roundtrip_1.0', y, 1.0_real64, FP16_TOL_TIGHT)

        x = FP16(-2.5_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp16_real64('roundtrip_neg2.5', y, -2.5_real64, FP16_TOL)

        x = FP16(0.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp16_real64('roundtrip_0.0', y, 0.0_real64, FP16_TOL_TIGHT)

        x = FP16(100.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp16_real64('roundtrip_100', y, 100.0_real64, FP16_TOL)
    end subroutine

    subroutine test_write_read_multiple()
        type(FP16) :: a, b, c
        character(len=40) :: buf1, buf2, buf3
        type(FP16) :: ra, rb, rc

        a = FP16(1.0_real32)
        b = FP16(2.0_real32)
        c = FP16(3.0_real32)

        write(buf1, '(DT)') a
        write(buf2, '(DT)') b
        write(buf3, '(DT)') c

        read(buf1, *) ra
        read(buf2, *) rb
        read(buf3, *) rc

        call check_fp16_real64('multi_read_a', ra, 1.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('multi_read_b', rb, 2.0_real64, FP16_TOL_TIGHT)
        call check_fp16_real64('multi_read_c', rc, 3.0_real64, FP16_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp16_io