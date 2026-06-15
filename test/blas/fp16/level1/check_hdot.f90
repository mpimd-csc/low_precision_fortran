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

PROGRAM check_hdot
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_FP16
    USE lpf_blas_fp16
    IMPLICIT NONE

    LOGICAL :: ok
    INTEGER :: i

    ok = .TRUE.
    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HDOT tests failed!"
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
        TYPE(FP16), dimension(NTEST) :: expected
        TYPE(FP16) :: result
        ! X data
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

        ! Y data
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
        y(:, 6) = 9.0_real32
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

        expected(1) = FP16(55.0_real32)
        expected(2) = FP16(20.75_real32)
        expected(3) = FP16(-5.0_real32)
        expected(4) = FP16(-3.0_real32)
        expected(5) = FP16(1.599609_real32)
        expected(6) = FP16(243.0_real32)
        expected(7) = FP16(6.0_real32)
        expected(8) = FP16(0.0_real32)

        do i = 1, NTEST
            result = dot(n(i), x(:, i), incx(i), y(:, i), incy(i))
            if (abs(dble(result) - dble(expected(i))) < 1.0d-2) then
                print "(A, I2, A, I3, A, I3, A, I3, A, F10.6, A, F10.6)", &
                      "HDOT   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", INCY = ", incy(i), ", Result = ", dble(result), ", Expected = ", dble(expected(i))
            else
                print "(A, I2, A, I3, A, I3, A, I3, A, F10.6, A, F10.6)", &
                      "HDOT   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", INCY = ", incy(i), ", Result = ", dble(result), ", Expected = ", dble(expected(i))
                ok = .FALSE.
            end if
        end do
    end SUBROUTINE
END PROGRAM check_hdot
