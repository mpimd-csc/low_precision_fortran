!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM check_hasum
    USE iso_fortran_env, only: real32, real64
    USE iso_c_binding, only: c_int64_t
    USE LPF_FP16
    USE lpf_blas_fp16
    IMPLICIT NONE

    LOGICAL :: mok
    INTEGER :: i

    mok = .TRUE.

    CALL run_tests(mok)

    IF (.NOT. mok) THEN
        PRINT *, "HASUM tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 6
        INTEGER(c_int64_t), dimension(NTEST) :: n = [5, 3, 4, 4, 1, 0]
        INTEGER(c_int64_t), dimension(NTEST) :: incx = [1, 1, 1, 2, 1, 1]
        TYPE(FP16), dimension(10, NTEST) :: data
        TYPE(FP16), dimension(NTEST) :: hasum_expected
        TYPE(FP16) :: result

        ! Initialize data
        ! Testcase 1: {1.0, -2.0, 3.0, -4.0, 5.0}
        data(1, 1) = FP16(1.0_real32)
        data(2, 1) = FP16(-2.0_real32)
        data(3, 1) = FP16(3.0_real32)
        data(4, 1) = FP16(-4.0_real32)
        data(5, 1) = FP16(5.0_real32)
        data(6, 1) = FP16(0.0_real32)
        data(7, 1) = FP16(0.0_real32)
        data(8, 1) = FP16(0.0_real32)
        data(9, 1) = FP16(0.0_real32)
        data(10, 1) = FP16(0.0_real32)
        ! Testcase 2: {1.5, -2.5, 3.5}
        data(1, 2) = FP16(1.5_real32)
        data(2, 2) = FP16(-2.5_real32)
        data(3, 2) = FP16(3.5_real32)
        data(4, 2) = FP16(0.0_real32)
        data(5, 2) = FP16(0.0_real32)
        data(6, 2) = FP16(0.0_real32)
        data(7, 2) = FP16(0.0_real32)
        data(8, 2) = FP16(0.0_real32)
        data(9, 2) = FP16(0.0_real32)
        data(10, 2) = FP16(0.0_real32)
        ! Testcase 3: {0, 0, 0, 0}
        data(:, 3) = 0.0_real32
        ! Testcase 4: {1, 2, -1, 2} (with INCX=2)
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
        ! Testcase 5: {10}
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
        ! Testcase 6: {9...} but N=0
        data(:, 6) = 9.0_real32


        hasum_expected(1) = FP16(15.0_real32)
        hasum_expected(2) = FP16(7.5_real32)
        hasum_expected(3) = FP16(0.0_real32)
        hasum_expected(4) = FP16(4.0_real32)
        hasum_expected(5) = FP16(10.0_real32)
        hasum_expected(6) = FP16(0.0_real32)

        do i = 1, NTEST
            result = asum(n(i), data(:, i), incx(i))
            if (abs(dble(result) - dble(hasum_expected(i))) < 1.0d-2) then
                print "(A, I2, A, I3, A, I3, A, F10.6, A, F10.6)", &
                      "HASUM   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", Result = ", dble(result), ", Expected = ", dble(hasum_expected(i))
            else
                print "(A, I2, A, I3, A, I3, A, F10.6, A, F10.6)", &
                      "HASUM   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", Result = ", dble(result), ", Expected = ", dble(hasum_expected(i))
                ok = .FALSE.
            end if
        end do
    end SUBROUTINE
END PROGRAM check_hasum
