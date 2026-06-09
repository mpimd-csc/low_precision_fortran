!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM test_fp16_dot_product
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_FP16
    USE lpf_fp16_test_utils
    IMPLICIT NONE

    CALL run_tests()
    CALL test_summary()

CONTAINS

    SUBROUTINE run_tests()
        INTEGER, PARAMETER :: NTEST = 5
        INTEGER(int64), dimension(NTEST) :: n = [3, 4, 5, 10, 0]
        TYPE(FP16), dimension(10, NTEST) :: x
        TYPE(FP16), dimension(10, NTEST) :: y
        TYPE(FP16) :: result
        REAL(real32) :: ref_result
        REAL(real32), dimension(10) :: xr, yr
        INTEGER :: i

        ! Test Case 1: Simple positive values
        x(1, 1) = FP16(1.0_real32)
        x(2, 1) = FP16(2.0_real32)
        x(3, 1) = FP16(3.0_real32)
        y(1, 1) = FP16(4.0_real32)
        y(2, 1) = FP16(5.0_real32)
        y(3, 1) = FP16(6.0_real32)

        ! Test Case 2: Positive and negative values
        x(1, 2) = FP16(1.0_real32)
        x(2, 2) = FP16(-2.0_real32)
        x(3, 2) = FP16(3.0_real32)
        x(4, 2) = FP16(-4.0_real32)
        y(1, 2) = FP16(1.0_real32)
        y(2, 2) = FP16(1.0_real32)
        y(3, 2) = FP16(1.0_real32)
        y(4, 2) = FP16(1.0_real32)

        ! Test Case 3: Including zeros
        x(1, 3) = FP16(0.0_real32)
        x(2, 3) = FP16(1.0_real32)
        x(3, 3) = FP16(0.0_real32)
        x(4, 3) = FP16(-1.0_real32)
        x(5, 3) = FP16(2.0_real32)
        y(1, 3) = FP16(1.0_real32)
        y(2, 3) = FP16(0.0_real32)
        y(3, 3) = FP16(1.0_real32)
        y(4, 3) = FP16(1.0_real32)
        y(5, 3) = FP16(1.0_real32)

        ! Test Case 4: Precision test (more elements)
        do i = 1, 10
            x(i, 4) = FP16(real(1, real32))
            y(i, 4) = FP16(real(i, real32)/10.0)
        end do

        ! Test Case 5: Empty vectors
        ! n(5) = 0, arrays are size 10 but we only use n(i) elements.

        do i = 1, NTEST
            ! Use slices to match the intended size n(i)
            result = dot_product(x(1:n(i), i), y(1:n(i), i))

            ! Reference result using real32
            xr(1:n(i)) = real(x(1:n(i), i))
            yr(1:n(i)) = real(y(1:n(i), i))
            ref_result = dot_product(xr(1:n(i)), yr(1:n(i)))

            CALL check_fp16_real64('dot_product case ' // trim(integer_to_char(i)), &
                                  result, real(ref_result, real64), 1.0d-2)
        end do
    end SUBROUTINE
END PROGRAM test_fp16_dot_product
