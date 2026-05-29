!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM check_hswap
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_FP16
    USE lpf_blas_fp16
    IMPLICIT NONE

    LOGICAL :: ok
    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HSWAP tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 8
        INTEGER(int64), dimension(NTEST) :: n = [5, 3, 4, 4, 4, 3, 3, 0]
        INTEGER(int64), dimension(NTEST) :: incx = [1, 1, 1, 2, 2, 1, 1, 1]
        INTEGER(int64), dimension(NTEST) :: incy = [1, 1, 1, 2, 1, 2, 1, 1]
        TYPE(FP16), dimension(10, NTEST) :: x
        TYPE(FP16), dimension(10, NTEST) :: y
        TYPE(FP16), dimension(10) :: xtemp, ytemp
        INTEGER :: k, i, j
        LOGICAL :: lok

        x(1, 1) = FP16(1.0_real32)
        x(2, 1) = FP16(-2.0_real32)
        x(3, 1) = FP16(3.0_real32)
        x(4, 1) = FP16(-4.0_real32)
        x(5, 1) = FP16(5.0_real32)
        x(6, 1) = FP16(0.0_real32)
        x(7, 1) = FP16(0.0_real32)
        x(8, 1) = FP16(0.0_real32)
        x(9, 1) = FP16(0.0_real32)
        x(10, 1) = FP16(0.0_real32)
        x(1, 2) = FP16(1.5_real32)
        x(2, 2) = FP16(-2.5_real32)
        x(3, 2) = FP16(3.5_real32)
        x(4, 2) = FP16(0.0_real32)
        x(5, 2) = FP16(0.0_real32)
        x(6, 2) = FP16(0.0_real32)
        x(7, 2) = FP16(0.0_real32)
        x(8, 2) = FP16(0.0_real32)
        x(9, 2) = FP16(0.0_real32)
        x(10, 2) = FP16(0.0_real32)
        x(1, 3) = FP16(3.0_real32)
        x(2, 3) = FP16(-4.0_real32)
        x(3, 3) = FP16(0.0_real32)
        x(4, 3) = FP16(0.0_real32)
        x(5, 3) = FP16(0.0_real32)
        x(6, 3) = FP16(0.0_real32)
        x(7, 3) = FP16(0.0_real32)
        x(8, 3) = FP16(0.0_real32)
        x(9, 3) = FP16(0.0_real32)
        x(10, 3) = FP16(0.0_real32)
        x(1, 4) = FP16(1.0_real32)
        x(2, 4) = FP16(2.0_real32)
        x(3, 4) = FP16(-1.0_real32)
        x(4, 4) = FP16(2.0_real32)
        x(5, 4) = FP16(1.0_real32)
        x(6, 4) = FP16(2.0_real32)
        x(7, 4) = FP16(-1.0_real32)
        x(8, 4) = FP16(2.0_real32)
        x(9, 4) = FP16(1.0_real32)
        x(10, 4) = FP16(2.0_real32)
        x(1, 5) = FP16(10.0_real32)
        x(2, 5) = FP16(4.0_real32)
        x(3, 5) = FP16(-4.0_real32)
        x(4, 5) = FP16(-5.0_real32)
        x(5, 5) = FP16(0.4_real32)
        x(6, 5) = FP16(12.0_real32)
        x(7, 5) = FP16(0.0_real32)
        x(8, 5) = FP16(0.0_real32)
        x(9, 5) = FP16(0.0_real32)
        x(10, 5) = FP16(0.0_real32)
        x(:, 6) = 9.0_real32
        x(1, 7) = FP16(1.0_real32)
        x(2, 7) = FP16(2.0_real32)
        x(3, 7) = FP16(-1.0_real32)
        x(4, 7) = FP16(2.0_real32)
        x(5, 7) = FP16(1.0_real32)
        x(6, 7) = FP16(2.0_real32)
        x(7, 7) = FP16(-1.0_real32)
        x(8, 7) = FP16(2.0_real32)
        x(9, 7) = FP16(1.0_real32)
        x(10, 7) = FP16(2.0_real32)
        x(1, 8) = FP16(1.0_real32)
        x(2, 8) = FP16(2.0_real32)
        x(3, 8) = FP16(-1.0_real32)
        x(4, 8) = FP16(2.0_real32)
        x(5, 8) = FP16(1.0_real32)
        x(6, 8) = FP16(2.0_real32)
        x(7, 8) = FP16(-1.0_real32)
        x(8, 8) = FP16(2.0_real32)
        x(9, 8) = FP16(1.0_real32)
        x(10, 8) = FP16(2.0_real32)

        y(1, 1) = FP16(1.0_real32)
        y(2, 1) = FP16(-2.0_real32)
        y(3, 1) = FP16(3.0_real32)
        y(4, 1) = FP16(-4.0_real32)
        y(5, 1) = FP16(5.0_real32)
        y(6, 1) = FP16(0.0_real32)
        y(7, 1) = FP16(0.0_real32)
        y(8, 1) = FP16(0.0_real32)
        y(9, 1) = FP16(0.0_real32)
        y(10, 1) = FP16(0.0_real32)
        y(1, 2) = FP16(1.5_real32)
        y(2, 2) = FP16(-2.5_real32)
        y(3, 2) = FP16(3.5_real32)
        y(4, 2) = FP16(0.0_real32)
        y(5, 2) = FP16(0.0_real32)
        y(6, 2) = FP16(0.0_real32)
        y(7, 2) = FP16(0.0_real32)
        y(8, 2) = FP16(0.0_real32)
        y(9, 2) = FP16(0.0_real32)
        y(10, 2) = FP16(0.0_real32)
        y(1, 3) = FP16(1.0_real32)
        y(2, 3) = FP16(2.0_real32)
        y(3, 3) = FP16(-1.0_real32)
        y(4, 3) = FP16(2.0_real32)
        y(5, 3) = FP16(1.0_real32)
        y(6, 3) = FP16(2.0_real32)
        y(7, 3) = FP16(-1.0_real32)
        y(8, 3) = FP16(2.0_real32)
        y(9, 3) = FP16(1.0_real32)
        y(10, 3) = FP16(2.0_real32)
        y(1, 4) = FP16(0.0_real32)
        y(2, 4) = FP16(3.0_real32)
        y(3, 4) = FP16(3.0_real32)
        y(4, 4) = FP16(0.0_real32)
        y(5, 4) = FP16(0.0_real32)
        y(6, 4) = FP16(0.0_real32)
        y(7, 4) = FP16(0.0_real32)
        y(8, 4) = FP16(0.0_real32)
        y(9, 4) = FP16(0.0_real32)
        y(10, 4) = FP16(0.0_real32)
        y(1, 5) = FP16(1.0_real32)
        y(2, 5) = FP16(2.0_real32)
        y(3, 5) = FP16(-1.0_real32)
        y(4, 5) = FP16(2.0_real32)
        y(5, 5) = FP16(1.0_real32)
        y(6, 5) = FP16(2.0_real32)
        y(7, 5) = FP16(-1.0_real32)
        y(8, 5) = FP16(2.0_real32)
        y(9, 5) = FP16(1.0_real32)
        y(10, 5) = FP16(2.0_real32)
        y(1, 6) = FP16(0.0_real32)
        y(2, 6) = FP16(0.0_real32)
        y(3, 6) = FP16(0.0_real32)
        y(4, 6) = FP16(0.0_real32)
        y(5, 6) = FP16(0.0_real32)
        y(6, 6) = FP16(-9.0_real32)
        y(7, 6) = FP16(-9.0_real32)
        y(8, 6) = FP16(-9.0_real32)
        y(9, 6) = FP16(-9.0_real32)
        y(10, 6) = FP16(-9.0_real32)
        y(1, 7) = FP16(1.0_real32)
        y(2, 7) = FP16(2.0_real32)
        y(3, 7) = FP16(-1.0_real32)
        y(4, 7) = FP16(2.0_real32)
        y(5, 7) = FP16(1.0_real32)
        y(6, 7) = FP16(2.0_real32)
        y(7, 7) = FP16(-1.0_real32)
        y(8, 7) = FP16(2.0_real32)
        y(9, 7) = FP16(1.0_real32)
        y(10, 7) = FP16(2.0_real32)
        y(1, 8) = FP16(10.0_real32)
        y(2, 8) = FP16(4.0_real32)
        y(3, 8) = FP16(-4.0_real32)
        y(4, 8) = FP16(-5.0_real32)
        y(5, 8) = FP16(0.4_real32)
        y(6, 8) = FP16(12.0_real32)
        y(7, 8) = FP16(0.0_real32)
        y(8, 8) = FP16(0.0_real32)
        y(9, 8) = FP16(0.0_real32)
        y(10, 8) = FP16(0.0_real32)

        do i = 1, NTEST
            xtemp = x(:, i)
            ytemp = y(:, i)
            CALL swap(n(i), xtemp, incx(i), ytemp, incy(i))

            lok = .TRUE.
            do k = 1, n(i)
                IF (xtemp(1+(k-1)*incx(i)) .ne. y((1+(k-1)*incy(i)), i) .OR. &
                    & ytemp(1+(k-1)*incy(i)) .ne. x(1+(k-1)*incx(i), i)) THEN
                    print "(A, I2, A, I3, A, I3, A, I3, A, I2)", &
                          "HSWAP   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", INCY = ", incy(i), ", K = ", k
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
            end do
            IF (lok) THEN
                print "(A, I2, A, I3, A, I3, A, I3)", &
                      "HSWAP   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), ", INCY = ", incy(i)
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_hswap
