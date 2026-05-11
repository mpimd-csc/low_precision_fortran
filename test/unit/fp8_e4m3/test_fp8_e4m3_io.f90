!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 formatted I/O.
!
!  Tests: write (list-directed, DT format), read round-trip

PROGRAM test_fp8_e4m3_io
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
    IMPLICIT NONE

    CALL test_write_listdirected()
    CALL test_read_roundtrip()
    CALL test_write_read_multiple()

    CALL test_summary()

CONTAINS

    subroutine test_write_listdirected()
        type(FP8_E4M3) :: x
        character(len=40) :: buffer
        real(real64) :: val

        x = FP8_E4M3(3.14_real32)
        write(buffer, '(DT)') x
        ! Verify we can read the buffer back
        read(buffer, '(DT)') x
        val = dble(x)
        call check_fp8_e4m3_real64('write_listdirected_3.14', FP8_E4M3(real(val, kind=real32)), 3.14_real64, FP8_E4M3_TOL)
    end subroutine

    subroutine test_read_roundtrip()
        type(FP8_E4M3) :: x, y
        character(len=40) :: buffer

        ! Test round-trip for various values
        x = FP8_E4M3(1.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp8_e4m3_real64('roundtrip_1.0', y, 1.0_real64, FP8_E4M3_TOL_TIGHT)

        x = FP8_E4M3(-2.5_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp8_e4m3_real64('roundtrip_neg2.5', y, -2.5_real64, FP8_E4M3_TOL)

        x = FP8_E4M3(0.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp8_e4m3_real64('roundtrip_0.0', y, 0.0_real64, FP8_E4M3_TOL_TIGHT)

        x = FP8_E4M3(100.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_fp8_e4m3_real64('roundtrip_100', y, 100.0_real64, FP8_E4M3_TOL)
    end subroutine

    subroutine test_write_read_multiple()
        type(FP8_E4M3) :: a, b, c
        character(len=40) :: buf1, buf2, buf3
        type(FP8_E4M3) :: ra, rb, rc

        a = FP8_E4M3(1.0_real32)
        b = FP8_E4M3(2.0_real32)
        c = FP8_E4M3(3.0_real32)

        ! Write/read each value separately to avoid list-directed parsing issues
        write(buf1, '(DT)') a
        write(buf2, '(DT)') b
        write(buf3, '(DT)') c

        read(buf1, *) ra
        read(buf2, *) rb
        read(buf3, *) rc

        call check_fp8_e4m3_real64('multi_read_a', ra, 1.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('multi_read_b', rb, 2.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('multi_read_c', rc, 3.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

END PROGRAM test_fp8_e4m3_io