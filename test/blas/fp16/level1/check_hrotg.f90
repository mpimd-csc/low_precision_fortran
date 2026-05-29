!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM check_hrotg
    USE iso_fortran_env, only: real32, real64
    USE LPF_FP16
    USE lpf_blas_fp16
    IMPLICIT NONE

    LOGICAL :: ok

    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HROTG tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 8
        TYPE(FP16) :: sa, sb, sc, ss
        TYPE(FP16), dimension(NTEST) :: da1, db1, dc1, ds1, datrue, dbtrue
        INTEGER :: i
        LOGICAL :: lok

        da1(1) = FP16(0.3_real32)
        da1(2) = FP16(0.4_real32)
        da1(3) = FP16(-0.3_real32)
        da1(4) = FP16(-0.4_real32)
        da1(5) = FP16(-0.3_real32)
        da1(6) = FP16(0.0_real32)
        da1(7) = FP16(0.0_real32)
        da1(8) = FP16(1.0_real32)
        db1(1) = FP16(0.4_real32)
        db1(2) = FP16(0.3_real32)
        db1(3) = FP16(0.4_real32)
        db1(4) = FP16(0.3_real32)
        db1(5) = FP16(-0.4_real32)
        db1(6) = FP16(0.0_real32)
        db1(7) = FP16(1.0_real32)
        db1(8) = FP16(0.0_real32)
        dc1(1) = FP16(0.6_real32)
        dc1(2) = FP16(0.8_real32)
        dc1(3) = FP16(-0.6_real32)
        dc1(4) = FP16(0.8_real32)
        dc1(5) = FP16(0.6_real32)
        dc1(6) = FP16(1.0_real32)
        dc1(7) = FP16(0.0_real32)
        dc1(8) = FP16(1.0_real32)
        ds1(1) = FP16(0.8_real32)
        ds1(2) = FP16(0.6_real32)
        ds1(3) = FP16(0.8_real32)
        ds1(4) = FP16(-0.6_real32)
        ds1(5) = FP16(0.8_real32)
        ds1(6) = FP16(0.0_real32)
        ds1(7) = FP16(1.0_real32)
        ds1(8) = FP16(0.0_real32)
        datrue(1) = FP16(0.5_real32)
        datrue(2) = FP16(0.5_real32)
        datrue(3) = FP16(0.5_real32)
        datrue(4) = FP16(-0.5_real32)
        datrue(5) = FP16(-0.5_real32)
        datrue(6) = FP16(0.0_real32)
        datrue(7) = FP16(1.0_real32)
        datrue(8) = FP16(1.0_real32)
        dbtrue(1) = FP16(1.0_real32/0.6_real32)
        dbtrue(2) = FP16(0.6_real32)
        dbtrue(3) = FP16(-1.0_real32/0.6_real32)
        dbtrue(4) = FP16(-0.6_real32)
        dbtrue(5) = FP16(1.0_real32/0.6_real32)
        dbtrue(6) = FP16(0.0_real32)
        dbtrue(7) = FP16(1.0_real32)
        dbtrue(8) = FP16(0.0_real32)

        do i = 1, NTEST
            sa = da1(i)
            sb = db1(i)
            CALL rotg(sa, sb, sc, ss)

            lok = .TRUE.
            IF (abs(dble(sa) - dble(datrue(i))) > 1.0d-2 * max(abs(dble(datrue(i))), 1.0d-5)) THEN
                print "(A, I2, A, F10.6, A, F10.6)", &
                      "HROTG   -- FAIL -- Testcase ", i, ": SA : Result = ", dble(sa), &
                      ", Expected = ", dble(datrue(i))
                ok = .FALSE.
                lok = .FALSE.
            END IF
            IF (abs(dble(sb) - dble(dbtrue(i))) > 1.0d-2 * max(abs(dble(dbtrue(i))), 1.0d-5)) THEN
                print "(A, I2, A, F10.6, A, F10.6)", &
                      "HROTG   -- FAIL -- Testcase ", i, ": SB : Result = ", dble(sb), &
                      ", Expected = ", dble(dbtrue(i))
                ok = .FALSE.
                lok = .FALSE.
            END IF
            IF (abs(dble(sc) - dble(dc1(i))) > 1.0d-2 * max(abs(dble(dc1(i))), 1.0d-5)) THEN
                print "(A, I2, A, F10.6, A, F10.6)", &
                      "HROTG   -- FAIL -- Testcase ", i, ": SC : Result = ", dble(sc), &
                      ", Expected = ", dble(dc1(i))
                ok = .FALSE.
                lok = .FALSE.
            END IF
            IF (abs(dble(ss) - dble(ds1(i))) > 1.0d-2 * max(abs(dble(ds1(i))), 1.0d-5)) THEN
                print "(A, I2, A, F10.6, A, F10.6)", &
                      "HROTG   -- FAIL -- Testcase ", i, ": SS : Result = ", dble(ss), &
                      ", Expected = ", dble(ds1(i))
                ok = .FALSE.
                lok = .FALSE.
            END IF
            IF (lok) THEN
                print "(A, I2, A)", "HROTG   -- PASS -- Testcase ", i, "."
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_hrotg
