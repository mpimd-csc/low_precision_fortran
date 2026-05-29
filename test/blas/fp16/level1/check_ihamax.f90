!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM check_ihamax
    USE iso_fortran_env, only: real32, int64
    USE LPF_FP16
    USE lpf_blas_fp16
    IMPLICIT NONE

    LOGICAL :: ok
    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "IHAMAX tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 6
        INTEGER(int64), dimension(NTEST) :: n = [5, 3, 4, 4, 1, 0]
        INTEGER(int64), dimension(NTEST) :: incx = [1, 1, 1, 2, 1, 1]
        TYPE(FP16), dimension(10, NTEST) :: data
        INTEGER(int64), dimension(NTEST) :: expected = [3, 2, 1, 1, 1, 0]
        INTEGER(int64) :: result
        INTEGER :: i
        LOGICAL :: lok

        data(1, 1) = FP16(1.0_real32)
        data(2, 1) = FP16(-2.0_real32)
        data(3, 1) = FP16(33.0_real32)
        data(4, 1) = FP16(-4.0_real32)
        data(5, 1) = FP16(5.0_real32)
        data(6, 1) = FP16(0.0_real32)
        data(7, 1) = FP16(0.0_real32)
        data(8, 1) = FP16(0.0_real32)
        data(9, 1) = FP16(0.0_real32)
        data(10, 1) = FP16(0.0_real32)
        data(1, 2) = FP16(1.5_real32)
        data(2, 2) = FP16(-4.5_real32)
        data(3, 2) = FP16(3.5_real32)
        data(4, 2) = FP16(0.0_real32)
        data(5, 2) = FP16(0.0_real32)
        data(6, 2) = FP16(0.0_real32)
        data(7, 2) = FP16(0.0_real32)
        data(8, 2) = FP16(0.0_real32)
        data(9, 2) = FP16(0.0_real32)
        data(10, 2) = FP16(0.0_real32)
        data(:, 3) = 0.0_real32
        data(1, 4) = FP16(1.0_real32)
        data(2, 4) = FP16(2.0_real32)
        data(3, 4) = FP16(-1.0_real32)
        data(4, 4) = FP16(2.0_real32)
        data(5, 4) = FP16(1.0_real32)
        data(6, 4) = FP16(2.0_real32)
        data(7, 4) = FP16(-1.0_real32)
        data(8, 4) = FP16(2.0_real32)
        data(9, 4) = FP16(1.0_real32)
        data(10, 4) = FP16(2.0_real32)
        data(1, 5) = FP16(10.0_real32)
        data(2, 5) = FP16(0.0_real32)
        data(3, 5) = FP16(0.0_real32)
        data(4, 5) = FP16(0.0_real32)
        data(5, 5) = FP16(0.0_real32)
        data(6, 5) = FP16(0.0_real32)
        data(7, 5) = FP16(0.0_real32)
        data(8, 5) = FP16(0.0_real32)
        data(9, 5) = FP16(0.0_real32)
        data(10, 5) = FP16(0.0_real32)
        data(:, 6) = 9.0_real32

        do i = 1, NTEST
            result = iamax(n(i), data(:, i), incx(i))

            lok = .TRUE.
            IF (result /= expected(i)) THEN
                print "(A, I2, A, I3, A, I3, A, I2, A, I2)", &
                      "IHAMAX  -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", Result = ", result, ", Expected = ", expected(i)
                ok = .FALSE.
                lok = .FALSE.
            END IF
            IF (lok) THEN
                print "(A, I2, A, I3, A, I3, A, I2, A, I2)", &
                      "IHAMAX  -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", Result = ", result, ", Expected = ", expected(i)
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_ihamax
