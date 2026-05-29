!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM check_bcopy
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_BF16
    USE lpf_blas_bf16
    IMPLICIT NONE

    LOGICAL :: ok
    INTEGER :: i

    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HCOPY tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 8
        INTEGER(int64), dimension(NTEST) :: n = [5, 3, 4, 4, 4, 3, 3, 0]
        INTEGER(int64), dimension(NTEST) :: incx = [1, 1, 1, 2, 2, 1, 1, 1]
        INTEGER(int64), dimension(NTEST) :: incy = [1, 1, 1, 2, 1, 2, 1, 1]
        TYPE(BF16), dimension(10, NTEST) :: x
        TYPE(BF16), dimension(10, NTEST) :: y
        TYPE(BF16), dimension(10, NTEST) :: expected
        TYPE(BF16), dimension(10) :: ytemp
        INTEGER :: k
        LOGICAL :: lok
        ! X data
        x(1, 1) = BF16(1.0_real32)
        x(2, 1) = BF16(-2.0_real32)
        x(3, 1) = BF16(3.0_real32)
        x(4, 1) = BF16(-4.0_real32)
        x(5, 1) = BF16(5.0_real32)
        x(6, 1) = BF16(0.0_real32)
        x(7, 1) = BF16(0.0_real32)
        x(8, 1) = BF16(0.0_real32)
        x(9, 1) = BF16(0.0_real32)
        x(10, 1) = BF16(0.0_real32)
        x(1, 2) = BF16(1.5_real32)
        x(2, 2) = BF16(-2.5_real32)
        x(3, 2) = BF16(3.5_real32)
        x(4, 2) = BF16(0.0_real32)
        x(5, 2) = BF16(0.0_real32)
        x(6, 2) = BF16(0.0_real32)
        x(7, 2) = BF16(0.0_real32)
        x(8, 2) = BF16(0.0_real32)
        x(9, 2) = BF16(0.0_real32)
        x(10, 2) = BF16(0.0_real32)
        x(1, 3) = BF16(3.0_real32)
        x(2, 3) = BF16(-4.0_real32)
        x(3, 3) = BF16(0.0_real32)
        x(4, 3) = BF16(0.0_real32)
        x(5, 3) = BF16(0.0_real32)
        x(6, 3) = BF16(0.0_real32)
        x(7, 3) = BF16(0.0_real32)
        x(8, 3) = BF16(0.0_real32)
        x(9, 3) = BF16(0.0_real32)
        x(10, 3) = BF16(0.0_real32)
        x(1, 4) = BF16(1.0_real32)
        x(2, 4) = BF16(2.0_real32)
        x(3, 4) = BF16(-1.0_real32)
        x(4, 4) = BF16(2.0_real32)
        x(5, 4) = BF16(1.0_real32)
        x(6, 4) = BF16(2.0_real32)
        x(7, 4) = BF16(-1.0_real32)
        x(8, 4) = BF16(2.0_real32)
        x(9, 4) = BF16(1.0_real32)
        x(10, 4) = BF16(2.0_real32)
        x(1, 5) = BF16(10.0_real32)
        x(2, 5) = BF16(4.0_real32)
        x(3, 5) = BF16(-4.0_real32)
        x(4, 5) = BF16(-5.0_real32)
        x(5, 5) = BF16(0.4_real32)
        x(6, 5) = BF16(12.0_real32)
        x(7, 5) = BF16(0.0_real32)
        x(8, 5) = BF16(0.0_real32)
        x(9, 5) = BF16(0.0_real32)
        x(10, 5) = BF16(0.0_real32)
        x(:, 6) = 9.0_real32
        x(1, 7) = BF16(1.0_real32)
        x(2, 7) = BF16(2.0_real32)
        x(3, 7) = BF16(-1.0_real32)
        x(4, 7) = BF16(2.0_real32)
        x(5, 7) = BF16(1.0_real32)
        x(6, 7) = BF16(2.0_real32)
        x(7, 7) = BF16(-1.0_real32)
        x(8, 7) = BF16(2.0_real32)
        x(9, 7) = BF16(1.0_real32)
        x(10, 7) = BF16(2.0_real32)
        x(1, 8) = BF16(1.0_real32)
        x(2, 8) = BF16(2.0_real32)
        x(3, 8) = BF16(-1.0_real32)
        x(4, 8) = BF16(2.0_real32)
        x(5, 8) = BF16(1.0_real32)
        x(6, 8) = BF16(2.0_real32)
        x(7, 8) = BF16(-1.0_real32)
        x(8, 8) = BF16(2.0_real32)
        x(9, 8) = BF16(1.0_real32)
        x(10, 8) = BF16(2.0_real32)

        ! Y data
        y(1, 1) = BF16(1.0_real32)
        y(2, 1) = BF16(-2.0_real32)
        y(3, 1) = BF16(3.0_real32)
        y(4, 1) = BF16(-4.0_real32)
        y(5, 1) = BF16(5.0_real32)
        y(6, 1) = BF16(0.0_real32)
        y(7, 1) = BF16(0.0_real32)
        y(8, 1) = BF16(0.0_real32)
        y(9, 1) = BF16(0.0_real32)
        y(10, 1) = BF16(0.0_real32)
        y(1, 2) = BF16(1.5_real32)
        y(2, 2) = BF16(-2.5_real32)
        y(3, 2) = BF16(3.5_real32)
        y(4, 2) = BF16(0.0_real32)
        y(5, 2) = BF16(0.0_real32)
        y(6, 2) = BF16(0.0_real32)
        y(7, 2) = BF16(0.0_real32)
        y(8, 2) = BF16(0.0_real32)
        y(9, 2) = BF16(0.0_real32)
        y(10, 2) = BF16(0.0_real32)
        y(1, 3) = BF16(1.0_real32)
        y(2, 3) = BF16(2.0_real32)
        y(3, 3) = BF16(-1.0_real32)
        y(4, 3) = BF16(2.0_real32)
        y(5, 3) = BF16(1.0_real32)
        y(6, 3) = BF16(2.0_real32)
        y(7, 3) = BF16(-1.0_real32)
        y(8, 3) = BF16(2.0_real32)
        y(9, 3) = BF16(1.0_real32)
        y(10, 3) = BF16(2.0_real32)
        y(1, 4) = BF16(0.0_real32)
        y(2, 4) = BF16(3.0_real32)
        y(3, 4) = BF16(3.0_real32)
        y(4, 4) = BF16(0.0_real32)
        y(5, 4) = BF16(0.0_real32)
        y(6, 4) = BF16(0.0_real32)
        y(7, 4) = BF16(0.0_real32)
        y(8, 4) = BF16(0.0_real32)
        y(9, 4) = BF16(0.0_real32)
        y(10, 4) = BF16(0.0_real32)
        y(1, 5) = BF16(1.0_real32)
        y(2, 5) = BF16(2.0_real32)
        y(3, 5) = BF16(-1.0_real32)
        y(4, 5) = BF16(2.0_real32)
        y(5, 5) = BF16(1.0_real32)
        y(6, 5) = BF16(2.0_real32)
        y(7, 5) = BF16(-1.0_real32)
        y(8, 5) = BF16(2.0_real32)
        y(9, 5) = BF16(1.0_real32)
        y(10, 5) = BF16(2.0_real32)
        y(:, 6) = [BF16(0.0_real32), BF16(0.0_real32), BF16(0.0_real32), &
                   BF16(0.0_real32), BF16(0.0_real32), BF16(-9.0_real32), &
                   BF16(-9.0_real32), BF16(-9.0_real32), BF16(-9.0_real32), BF16(-9.0_real32)]
        y(1, 7) = BF16(1.0_real32)
        y(2, 7) = BF16(2.0_real32)
        y(3, 7) = BF16(-1.0_real32)
        y(4, 7) = BF16(2.0_real32)
        y(5, 7) = BF16(1.0_real32)
        y(6, 7) = BF16(2.0_real32)
        y(7, 7) = BF16(-1.0_real32)
        y(8, 7) = BF16(2.0_real32)
        y(9, 7) = BF16(1.0_real32)
        y(10, 7) = BF16(2.0_real32)
        y(1, 8) = BF16(10.0_real32)
        y(2, 8) = BF16(4.0_real32)
        y(3, 8) = BF16(-4.0_real32)
        y(4, 8) = BF16(-5.0_real32)
        y(5, 8) = BF16(0.4_real32)
        y(6, 8) = BF16(12.0_real32)
        y(7, 8) = BF16(0.0_real32)
        y(8, 8) = BF16(0.0_real32)
        y(9, 8) = BF16(0.0_real32)
        y(10, 8) = BF16(0.0_real32)

        ! Expected results
        expected(1, 1) = BF16(1.0_real32)
        expected(2, 1) = BF16(-2.0_real32)
        expected(3, 1) = BF16(3.0_real32)
        expected(4, 1) = BF16(-4.0_real32)
        expected(5, 1) = BF16(5.0_real32)
        expected(6, 1) = BF16(0.0_real32)
        expected(7, 1) = BF16(0.0_real32)
        expected(8, 1) = BF16(0.0_real32)
        expected(9, 1) = BF16(0.0_real32)
        expected(10, 1) = BF16(0.0_real32)
        expected(1, 2) = BF16(1.5_real32)
        expected(2, 2) = BF16(-2.5_real32)
        expected(3, 2) = BF16(3.5_real32)
        expected(4, 2) = BF16(0.0_real32)
        expected(5, 2) = BF16(0.0_real32)
        expected(6, 2) = BF16(0.0_real32)
        expected(7, 2) = BF16(0.0_real32)
        expected(8, 2) = BF16(0.0_real32)
        expected(9, 2) = BF16(0.0_real32)
        expected(10, 2) = BF16(0.0_real32)
        expected(1, 3) = BF16(3.0_real32)
        expected(2, 3) = BF16(-4.0_real32)
        expected(3, 3) = BF16(0.0_real32)
        expected(4, 3) = BF16(0.0_real32)
        expected(5, 3) = BF16(1.0_real32)
        expected(6, 3) = BF16(2.0_real32)
        expected(7, 3) = BF16(-1.0_real32)
        expected(8, 3) = BF16(2.0_real32)
        expected(9, 3) = BF16(1.0_real32)
        expected(10, 3) = BF16(2.0_real32)
        expected(1, 4) = BF16(1.0_real32)
        expected(2, 4) = BF16(3.0_real32)
        expected(3, 4) = BF16(-1.0_real32)
        expected(4, 4) = BF16(0.0_real32)
        expected(5, 4) = BF16(1.0_real32)
        expected(6, 4) = BF16(0.0_real32)
        expected(7, 4) = BF16(-1.0_real32)
        expected(8, 4) = BF16(0.0_real32)
        expected(9, 4) = BF16(0.0_real32)
        expected(10, 4) = BF16(0.0_real32)
        expected(1, 5) = BF16(10.0_real32)
        expected(2, 5) = BF16(-4.0_real32)
        expected(3, 5) = BF16(0.399902_real32)
        expected(4, 5) = BF16(0.0_real32)
        expected(5, 5) = BF16(1.0_real32)
        expected(6, 5) = BF16(2.0_real32)
        expected(7, 5) = BF16(-1.0_real32)
        expected(8, 5) = BF16(2.0_real32)
        expected(9, 5) = BF16(1.0_real32)
        expected(10, 5) = BF16(2.0_real32)
        expected(1, 6) = BF16(9.0_real32)
        expected(2, 6) = BF16(0.0_real32)
        expected(3, 6) = BF16(9.0_real32)
        expected(4, 6) = BF16(0.0_real32)
        expected(5, 6) = BF16(9.0_real32)
        expected(6, 6) = BF16(-9.0_real32)
        expected(7, 6) = BF16(-9.0_real32)
        expected(8, 6) = BF16(-9.0_real32)
        expected(9, 6) = BF16(-9.0_real32)
        expected(10, 6) = BF16(-9.0_real32)
        expected(1, 7) = BF16(1.0_real32)
        expected(2, 7) = BF16(2.0_real32)
        expected(3, 7) = BF16(-1.0_real32)
        expected(4, 7) = BF16(2.0_real32)
        expected(5, 7) = BF16(1.0_real32)
        expected(6, 7) = BF16(2.0_real32)
        expected(7, 7) = BF16(-1.0_real32)
        expected(8, 7) = BF16(2.0_real32)
        expected(9, 7) = BF16(1.0_real32)
        expected(10, 7) = BF16(2.0_real32)
        expected(1, 8) = BF16(10.0_real32)
        expected(2, 8) = BF16(4.0_real32)
        expected(3, 8) = BF16(-4.0_real32)
        expected(4, 8) = BF16(-5.0_real32)
        expected(5, 8) = BF16(0.399902_real32)
        expected(6, 8) = BF16(12.0_real32)
        expected(7, 8) = BF16(0.0_real32)
        expected(8, 8) = BF16(0.0_real32)
        expected(9, 8) = BF16(0.0_real32)
        expected(10, 8) = BF16(0.0_real32)

        do i = 1, NTEST
            ytemp = y(:, i)
            CALL copy(n(i), x(:, i), incx(i), ytemp, incy(i))

            lok = .TRUE.
            do k = 1, 10
                IF (abs(dble(ytemp(k)) - dble(expected(k, i))) > 1.0d-2 * max(abs(dble(expected(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I3, A, I3, A, I3, A, I2, A, F10.6, A, F10.6)", &
                          "HCOPY   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", INCY = ", incy(i), ", K = ", k, ", Result = ", dble(ytemp(k)), &
                          ", Expected = ", dble(expected(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
            end do
            IF (lok) THEN
                print "(A, I2, A, I3, A, I3, A, I3)", &
                      "HCOPY   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), ", INCY = ", incy(i)
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_bcopy
