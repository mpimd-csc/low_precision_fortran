!  SPDX-License-Identifier: LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program; if not, write to the Free Software Foundation,
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

PROGRAM check_brot
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_BF16
    USE lpf_blas_bf16
    IMPLICIT NONE

    LOGICAL :: ok
    INTEGER :: i

    ok = .TRUE.
    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HROT tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 8
        TYPE(BF16) :: sc, ss
        INTEGER(c_int64_t), dimension(NTEST) :: n = [1, 1, 2, 4, 1, 1, 3, 7]
        INTEGER(c_int64_t), dimension(NTEST) :: incx = [1, 1, 1, 2, 2, 1, 1, 1]
        INTEGER(c_int64_t), dimension(NTEST) :: incy = [1, 1, 1, 2, 1, 2, 1, 1]
        TYPE(BF16), dimension(10) :: dx, dy
        TYPE(BF16), dimension(10, NTEST) :: xresult
        TYPE(BF16), dimension(10, NTEST) :: yresult
        TYPE(BF16), dimension(10) :: xtemp, ytemp
        INTEGER :: k
        LOGICAL :: lok

        sc = BF16(0.8_real32)
        ss = BF16(0.6_real32)

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

        xresult(1, 1) = BF16(0.779785_real32)
        xresult(2, 1) = BF16(0.099976_real32)
        xresult(3, 1) = BF16(-0.500000_real32)
        xresult(4, 1) = BF16(0.799805_real32)
        xresult(5, 1) = BF16(0.899902_real32)
        xresult(6, 1) = BF16(-0.300049_real32)
        xresult(7, 1) = BF16(-0.399902_real32)
        xresult(8, 1) = BF16(0.0_real32)
        xresult(9, 1) = BF16(0.0_real32)
        xresult(10, 1) = BF16(0.0_real32)
        xresult(:, 2) = xresult(:, 1)
        xresult(1, 3) = BF16(0.779785_real32)
        xresult(2, 3) = BF16(-0.459961_real32)
        xresult(3, 3) = BF16(-0.500000_real32)
        xresult(4, 3) = BF16(0.799805_real32)
        xresult(5, 3) = BF16(0.899902_real32)
        xresult(6, 3) = BF16(-0.300049_real32)
        xresult(7, 3) = BF16(-0.399902_real32)
        xresult(8, 3) = BF16(0.0_real32)
        xresult(9, 3) = BF16(0.0_real32)
        xresult(10, 3) = BF16(0.0_real32)
        xresult(1, 4) = BF16(0.779785_real32)
        xresult(2, 4) = BF16(0.099976_real32)
        xresult(3, 4) = BF16(-0.219849_real32)
        xresult(4, 4) = BF16(0.799805_real32)
        xresult(5, 4) = BF16(0.359619_real32)
        xresult(6, 4) = BF16(-0.300049_real32)
        xresult(7, 4) = BF16(0.160156_real32)
        xresult(8, 4) = BF16(0.0_real32)
        xresult(9, 4) = BF16(0.0_real32)
        xresult(10, 4) = BF16(0.0_real32)
        xresult(:, 5) = xresult(:, 1)
        xresult(:, 6) = xresult(:, 2)
        xresult(1, 7) = BF16(0.779785_real32)
        xresult(2, 7) = BF16(-0.459961_real32)
        xresult(3, 7) = BF16(-0.219849_real32)
        xresult(4, 7) = BF16(0.799805_real32)
        xresult(5, 7) = BF16(0.899902_real32)
        xresult(6, 7) = BF16(-0.300049_real32)
        xresult(7, 7) = BF16(-0.399902_real32)
        xresult(8, 7) = BF16(0.0_real32)
        xresult(9, 7) = BF16(0.0_real32)
        xresult(10, 7) = BF16(0.0_real32)
        xresult(1, 8) = BF16(0.779785_real32)
        xresult(2, 8) = BF16(-0.459961_real32)
        xresult(3, 8) = BF16(-0.219849_real32)
        xresult(4, 8) = BF16(1.059570_real32)
        xresult(5, 8) = BF16(0.359619_real32)
        xresult(6, 8) = BF16(-0.119995_real32)
        xresult(7, 8) = BF16(0.162109_real32)
        xresult(8, 8) = BF16(0.0_real32)
        xresult(9, 8) = BF16(0.0_real32)
        xresult(10, 8) = BF16(0.0_real32)

        yresult(1, 1) = BF16(0.039062_real32)
        yresult(2, 1) = BF16(-0.899902_real32)
        yresult(3, 1) = BF16(0.300049_real32)
        yresult(4, 1) = BF16(0.700195_real32)
        yresult(5, 1) = BF16(-0.600098_real32)
        yresult(6, 1) = BF16(0.199951_real32)
        yresult(7, 1) = BF16(0.799805_real32)
        yresult(8, 1) = BF16(0.0_real32)
        yresult(9, 1) = BF16(0.0_real32)
        yresult(10, 1) = BF16(0.0_real32)
        yresult(:, 2) = yresult(:, 1)
        yresult(1, 3) = BF16(0.039062_real32)
        yresult(2, 3) = BF16(-0.779785_real32)
        yresult(3, 3) = BF16(0.300049_real32)
        yresult(4, 3) = BF16(0.700195_real32)
        yresult(5, 3) = BF16(-0.600098_real32)
        yresult(6, 3) = BF16(0.199951_real32)
        yresult(7, 3) = BF16(0.799805_real32)
        yresult(8, 3) = BF16(0.0_real32)
        yresult(9, 3) = BF16(0.0_real32)
        yresult(10, 3) = BF16(0.0_real32)
        yresult(1, 4) = BF16(0.039062_real32)
        yresult(2, 4) = BF16(-0.899902_real32)
        yresult(3, 4) = BF16(0.540039_real32)
        yresult(4, 4) = BF16(0.700195_real32)
        yresult(5, 4) = BF16(-1.019531_real32)
        yresult(6, 4) = BF16(0.199951_real32)
        yresult(7, 4) = BF16(0.879883_real32)
        yresult(8, 4) = BF16(0.0_real32)
        yresult(9, 4) = BF16(0.0_real32)
        yresult(10, 4) = BF16(0.0_real32)
        yresult(:, 5) = yresult(:, 1)
        yresult(:, 6) = yresult(:, 2)
        yresult(1, 7) = BF16(0.039062_real32)
        yresult(2, 7) = BF16(-0.779785_real32)
        yresult(3, 7) = BF16(0.540039_real32)
        yresult(4, 7) = BF16(0.700195_real32)
        yresult(5, 7) = BF16(-0.600098_real32)
        yresult(6, 7) = BF16(0.199951_real32)
        yresult(7, 7) = BF16(0.799805_real32)
        yresult(8, 7) = BF16(0.0_real32)
        yresult(9, 7) = BF16(0.0_real32)
        yresult(10, 7) = BF16(0.0_real32)
        yresult(1, 8) = BF16(0.039062_real32)
        yresult(2, 8) = BF16(-0.779785_real32)
        yresult(3, 8) = BF16(0.540039_real32)
        yresult(4, 8) = BF16(0.076172_real32)
        yresult(5, 8) = BF16(-1.019531_real32)
        yresult(6, 8) = BF16(0.340088_real32)
        yresult(7, 8) = BF16(0.879883_real32)
        yresult(8, 8) = BF16(0.0_real32)
        yresult(9, 8) = BF16(0.0_real32)
        yresult(10, 8) = BF16(0.0_real32)

        do i = 1, NTEST
            xtemp = dx
            ytemp = dy
            CALL rot(n(i), xtemp, incx(i), ytemp, incy(i), sc, ss)

            lok = .TRUE.
            do k = 1, n(i)
                IF (abs(dble(xtemp(k)) - dble(xresult(k, i))) > 1.0d-1 * max(abs(dble(xresult(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I3, A, I3, A, I3, A, I2, A, F10.6, A, F10.6)", &
                          "HROT   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", INCY = ", incy(i), ", K = ", k, ", X, computed = ", dble(xtemp(k)), &
                          ", expected = ", dble(xresult(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
                IF (abs(dble(ytemp(k)) - dble(yresult(k, i))) > 1.0d-1 * max(abs(dble(yresult(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I3, A, I3, A, I3, A, I2, A, F10.6, A, F10.6)", &
                          "HROT   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", INCY = ", incy(i), ", K = ", k, ", Y, computed = ", dble(ytemp(k)), &
                          ", expected = ", dble(yresult(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
            end do
            IF (lok) THEN
                print "(A, I2, A, I3, A, I3, A, I3)", &
                      "HROT   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), ", INCY = ", incy(i)
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_brot
