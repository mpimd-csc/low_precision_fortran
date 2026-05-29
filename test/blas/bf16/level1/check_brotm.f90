!  SPDX-License-Identifier: LGPL-3.0-or-later
PROGRAM check_brotm
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_BF16
    USE lpf_blas_bf16
    IMPLICIT NONE

    LOGICAL :: ok

    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HROTM tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 8
        INTEGER(int64), dimension(NTEST) :: n = [1, 1, 2, 4, 1, 1, 3, 7]
        INTEGER(int64), dimension(NTEST) :: incx = [1, 1, 1, 2, 2, 1, 1, 1]
        INTEGER(int64), dimension(NTEST) :: incy = [1, 1, 1, 2, 1, 2, 1, 1]
        TYPE(BF16), dimension(10) :: dx, dy
        TYPE(BF16), dimension(5, 4) :: dpar
        TYPE(BF16), dimension(10, NTEST) :: xresult
        TYPE(BF16), dimension(10, NTEST) :: yresult
        TYPE(BF16), dimension(10) :: xtemp, ytemp
        TYPE(BF16), dimension(5) :: partemp
        INTEGER :: k, i
        LOGICAL :: lok

        dx(1) = BF16(0.6_real32)
        dx(2) = BF16(0.1_real32)
        dx(3) = BF16(-0.5_real32)
        dx(4) = BF16(0.8_real32)
        dx(5) = BF16(0.9_real32)
        dx(6) = BF16(-0.3_real32)
        dx(7) = BF16(-0.4_real32)
        dx(8) = BF16(0.0_real32)
        dx(9) = BF16(0.0_real32)
        dx(10) = BF16(0.0_real32)
        dy(1) = BF16(0.5_real32)
        dy(2) = BF16(-0.9_real32)
        dy(3) = BF16(0.3_real32)
        dy(4) = BF16(0.7_real32)
        dy(5) = BF16(-0.6_real32)
        dy(6) = BF16(0.2_real32)
        dy(7) = BF16(0.8_real32)
        dy(8) = BF16(0.0_real32)
        dy(9) = BF16(0.0_real32)
        dy(10) = BF16(0.0_real32)

        dpar(1, 1) = BF16(-2.0_real32)
        dpar(2, 1) = BF16(0.0_real32)
        dpar(3, 1) = BF16(0.0_real32)
        dpar(4, 1) = BF16(0.0_real32)
        dpar(5, 1) = BF16(0.0_real32)
        dpar(1, 2) = BF16(-1.0_real32)
        dpar(2, 2) = BF16(2.0_real32)
        dpar(3, 2) = BF16(-3.0_real32)
        dpar(4, 2) = BF16(-4.0_real32)
        dpar(5, 2) = BF16(5.0_real32)
        dpar(1, 3) = BF16(0.0_real32)
        dpar(2, 3) = BF16(0.0_real32)
        dpar(3, 3) = BF16(2.0_real32)
        dpar(4, 3) = BF16(-3.0_real32)
        dpar(5, 3) = BF16(0.0_real32)
        dpar(1, 4) = BF16(1.0_real32)
        dpar(2, 4) = BF16(5.0_real32)
        dpar(3, 4) = BF16(2.0_real32)
        dpar(4, 4) = BF16(0.0_real32)
        dpar(5, 4) = BF16(-4.0_real32)

        xresult(1, 1) = BF16(0.600098_real32)
        xresult(2, 1) = BF16(0.099976_real32)
        xresult(3, 1) = BF16(-0.500000_real32)
        xresult(4, 1) = BF16(0.799805_real32)
        xresult(5, 1) = BF16(0.899902_real32)
        xresult(6, 1) = BF16(-0.300049_real32)
        xresult(7, 1) = BF16(-0.399902_real32)
        xresult(8, 1) = BF16(0.0_real32)
        xresult(9, 1) = BF16(0.0_real32)
        xresult(10, 1) = BF16(0.0_real32)
        xresult(1, 2) = BF16(-0.799805_real32)
        xresult(2, 2) = BF16(0.099976_real32)
        xresult(3, 2) = BF16(-0.500000_real32)
        xresult(4, 2) = BF16(0.799805_real32)
        xresult(5, 2) = BF16(0.899902_real32)
        xresult(6, 2) = BF16(-0.300049_real32)
        xresult(7, 2) = BF16(-0.399902_real32)
        xresult(8, 2) = BF16(0.0_real32)
        xresult(9, 2) = BF16(0.0_real32)
        xresult(10, 2) = BF16(0.0_real32)
        xresult(1, 3) = BF16(-0.899902_real32)
        xresult(2, 3) = BF16(2.798828_real32)
        xresult(3, 3) = BF16(-0.500000_real32)
        xresult(4, 3) = BF16(0.799805_real32)
        xresult(5, 3) = BF16(0.899902_real32)
        xresult(6, 3) = BF16(-0.300049_real32)
        xresult(7, 3) = BF16(-0.399902_real32)
        xresult(8, 3) = BF16(0.0_real32)
        xresult(9, 3) = BF16(0.0_real32)
        xresult(10, 3) = BF16(0.0_real32)
        xresult(1, 4) = BF16(3.500000_real32)
        xresult(2, 4) = BF16(0.099976_real32)
        xresult(3, 4) = BF16(-2.199219_real32)
        xresult(4, 4) = BF16(0.799805_real32)
        xresult(5, 4) = BF16(3.898438_real32)
        xresult(6, 4) = BF16(-0.300049_real32)
        xresult(7, 4) = BF16(-1.199219_real32)
        xresult(8, 4) = BF16(0.0_real32)
        xresult(9, 4) = BF16(0.0_real32)
        xresult(10, 4) = BF16(0.0_real32)
        xresult(:, 5) = xresult(:, 1)
        xresult(:, 6) = xresult(:, 2)
        xresult(1, 7) = BF16(-0.899902_real32)
        xresult(2, 7) = BF16(2.798828_real32)
        xresult(3, 7) = BF16(-1.400391_real32)
        xresult(4, 7) = BF16(0.799805_real32)
        xresult(5, 7) = BF16(0.899902_real32)
        xresult(6, 7) = BF16(-0.300049_real32)
        xresult(7, 7) = BF16(-0.399902_real32)
        xresult(8, 7) = BF16(0.0_real32)
        xresult(9, 7) = BF16(0.0_real32)
        xresult(10, 7) = BF16(0.0_real32)
        xresult(1, 8) = BF16(3.500000_real32)
        xresult(2, 8) = BF16(-0.399902_real32)
        xresult(3, 8) = BF16(-2.199219_real32)
        xresult(4, 8) = BF16(4.699219_real32)
        xresult(5, 8) = BF16(3.898438_real32)
        xresult(6, 8) = BF16(-1.300781_real32)
        xresult(7, 8) = BF16(-1.199219_real32)
        xresult(8, 8) = BF16(0.0_real32)
        xresult(9, 8) = BF16(0.0_real32)
        xresult(10, 8) = BF16(0.0_real32)

        yresult(1, 1) = BF16(0.500000_real32)
        yresult(2, 1) = BF16(-0.899902_real32)
        yresult(3, 1) = BF16(0.300049_real32)
        yresult(4, 1) = BF16(0.700195_real32)
        yresult(5, 1) = BF16(-0.600098_real32)
        yresult(6, 1) = BF16(0.199951_real32)
        yresult(7, 1) = BF16(0.799805_real32)
        yresult(8, 1) = BF16(0.0_real32)
        yresult(9, 1) = BF16(0.0_real32)
        yresult(10, 1) = BF16(0.0_real32)
        yresult(1, 2) = BF16(0.699707_real32)
        yresult(2, 2) = BF16(-0.899902_real32)
        yresult(3, 2) = BF16(0.300049_real32)
        yresult(4, 2) = BF16(0.700195_real32)
        yresult(5, 2) = BF16(-0.600098_real32)
        yresult(6, 2) = BF16(0.199951_real32)
        yresult(7, 2) = BF16(0.799805_real32)
        yresult(8, 2) = BF16(0.0_real32)
        yresult(9, 2) = BF16(0.0_real32)
        yresult(10, 2) = BF16(0.0_real32)
        yresult(1, 3) = BF16(1.700195_real32)
        yresult(2, 3) = BF16(-0.700195_real32)
        yresult(3, 3) = BF16(0.300049_real32)
        yresult(4, 3) = BF16(0.700195_real32)
        yresult(5, 3) = BF16(-0.600098_real32)
        yresult(6, 3) = BF16(0.199951_real32)
        yresult(7, 3) = BF16(0.799805_real32)
        yresult(8, 3) = BF16(0.0_real32)
        yresult(9, 3) = BF16(0.0_real32)
        yresult(10, 3) = BF16(0.0_real32)
        yresult(1, 4) = BF16(-2.599609_real32)
        yresult(2, 4) = BF16(-0.899902_real32)
        yresult(3, 4) = BF16(-0.700195_real32)
        yresult(4, 4) = BF16(0.700195_real32)
        yresult(5, 4) = BF16(1.500000_real32)
        yresult(6, 4) = BF16(0.199951_real32)
        yresult(7, 4) = BF16(-2.798828_real32)
        yresult(8, 4) = BF16(0.0_real32)
        yresult(9, 4) = BF16(0.0_real32)
        yresult(10, 4) = BF16(0.0_real32)
        yresult(:, 5) = yresult(:, 1)
        yresult(:, 6) = yresult(:, 2)
        yresult(1, 7) = BF16(1.700195_real32)
        yresult(2, 7) = BF16(-0.700195_real32)
        yresult(3, 7) = BF16(-0.700195_real32)
        yresult(4, 7) = BF16(0.700195_real32)
        yresult(5, 7) = BF16(-0.600098_real32)
        yresult(6, 7) = BF16(0.199951_real32)
        yresult(7, 7) = BF16(0.799805_real32)
        yresult(8, 7) = BF16(0.0_real32)
        yresult(9, 7) = BF16(0.0_real32)
        yresult(10, 7) = BF16(0.0_real32)
        yresult(1, 8) = BF16(-2.599609_real32)
        yresult(2, 8) = BF16(3.500000_real32)
        yresult(3, 8) = BF16(-0.700195_real32)
        yresult(4, 8) = BF16(-3.601562_real32)
        yresult(5, 8) = BF16(1.500000_real32)
        yresult(6, 8) = BF16(-0.499756_real32)
        yresult(7, 8) = BF16(-2.798828_real32)
        yresult(8, 8) = BF16(0.0_real32)
        yresult(9, 8) = BF16(0.0_real32)
        yresult(10, 8) = BF16(0.0_real32)

        do i = 1, NTEST
            xtemp = dx
            ytemp = dy
            partemp = dpar(:,mod(i-1, 4)+1)
            CALL rotm(n(i), xtemp, incx(i), ytemp, incy(i), partemp)

            lok = .TRUE.
            do k = 1, n(i)
                IF (abs(dble(xtemp(k)) - dble(xresult(k, i))) > 1.0d-2 * max(abs(dble(xresult(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I3, A, I3, A, I3, A, I2, A, F10.6, A, F10.6)", &
                          "HROTM   -- FAIL -- Testcase X ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", INCY = ", incy(i), ", K = ", k, ", X, computed = ", dble(xtemp(k)), &
                          ", expected = ", dble(xresult(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
                IF (abs(dble(ytemp(k)) - dble(yresult(k, i))) > 1.0d-2 * max(abs(dble(yresult(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I3, A, I3, A, I3, A, I2, A, F10.6, A, F10.6)", &
                          "HROTM   -- FAIL -- Testcase Y ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", INCY = ", incy(i), ", K = ", k, ", Y, computed = ", dble(ytemp(k)), &
                          ", expected = ", dble(yresult(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
            end do
            IF (lok) THEN
                print "(A, I2, A, I3, A, I3, A, I3)", &
                      "HROTM   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), ", INCY = ", incy(i)
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_brotm
