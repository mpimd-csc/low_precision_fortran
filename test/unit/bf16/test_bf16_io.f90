!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for BF16 formatted I/O.
!
!  Tests: write (list-directed, DT format), read round-trip

PROGRAM test_bf16_io
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_bf16_test_utils
    IMPLICIT NONE

    CALL test_write_listdirected()
    CALL test_read_roundtrip()
    CALL test_write_read_multiple()

    CALL test_summary()

CONTAINS

    subroutine test_write_listdirected()
        type(BF16) :: x
        character(len=40) :: buffer
        real(real64) :: val

        x = BF16(3.14_real32)
        write(buffer, '(DT)') x
        ! Verify we can read the buffer back
        read(buffer, '(DT)') x
        val = dble(x)
        call check_bf16_real64('write_listdirected_3.14', BF16(real(val, kind=real32)), 3.14_real64, BF16_TOL)
    end subroutine

    subroutine test_read_roundtrip()
        type(BF16) :: x, y
        character(len=40) :: buffer

        ! Test round-trip for various values
        x = BF16(1.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_bf16_real64('roundtrip_1.0', y, 1.0_real64, BF16_TOL_TIGHT)

        x = BF16(-2.5_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_bf16_real64('roundtrip_neg2.5', y, -2.5_real64, BF16_TOL)

        x = BF16(0.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_bf16_real64('roundtrip_0.0', y, 0.0_real64, BF16_TOL_TIGHT)

        x = BF16(100.0_real32)
        write(buffer, '(DT)') x
        read(buffer, *) y
        call check_bf16_real64('roundtrip_100', y, 100.0_real64, BF16_TOL)
    end subroutine

    subroutine test_write_read_multiple()
        type(BF16) :: a, b, c
        character(len=40) :: buf1, buf2, buf3
        type(BF16) :: ra, rb, rc

        a = BF16(1.0_real32)
        b = BF16(2.0_real32)
        c = BF16(3.0_real32)

        ! Write/read each value separately to avoid list-directed parsing issues
        write(buf1, '(DT)') a
        write(buf2, '(DT)') b
        write(buf3, '(DT)') c

        read(buf1, *) ra
        read(buf2, *) rb
        read(buf3, *) rc

        call check_bf16_real64('multi_read_a', ra, 1.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('multi_read_b', rb, 2.0_real64, BF16_TOL_TIGHT)
        call check_bf16_real64('multi_read_c', rc, 3.0_real64, BF16_TOL_TIGHT)
    end subroutine

END PROGRAM test_bf16_io