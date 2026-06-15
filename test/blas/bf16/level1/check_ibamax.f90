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

PROGRAM check_ibamax
    USE iso_fortran_env, only: real32, int64
    USE LPF_BF16
    USE lpf_blas_bf16
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
        INTEGER(int64) :: n_val
        INTEGER(int64), dimension(NTEST) :: n
        INTEGER(int64), dimension(NTEST) :: incx
        TYPE(BF16), dimension(10, NTEST) :: data
        INTEGER(int64), dimension(NTEST) :: expected
        INTEGER(int64) :: result
        INTEGER :: i
        LOGICAL :: lok

        n(1) = 5
        n(2) = 3
        n(3) = 4
        n(4) = 4
        n(5) = 1
        n(6) = 0
        incx(1) = 1
        incx(2) = 1
        incx(3) = 1
        incx(4) = 2
        incx(5) = 1
        incx(6) = 1
        expected(1) = 3
        expected(2) = 2
        expected(3) = 1
        expected(4) = 1
        expected(5) = 1
        expected(6) = 0

        data(1, 1) = BF16(1.0_real32)
        data(2, 1) = BF16(-2.0_real32)
        data(3, 1) = BF16(33.0_real32)
        data(4, 1) = BF16(-4.0_real32)
        data(5, 1) = BF16(5.0_real32)
        data(6, 1) = BF16(0.0_real32)
        data(7, 1) = BF16(0.0_real32)
        data(8, 1) = BF16(0.0_real32)
        data(9, 1) = BF16(0.0_real32)
        data(10, 1) = BF16(0.0_real32)
        data(1, 2) = BF16(1.5_real32)
        data(2, 2) = BF16(-4.5_real32)
        data(3, 2) = BF16(3.5_real32)
        data(4, 2) = BF16(0.0_real32)
        data(5, 2) = BF16(0.0_real32)
        data(6, 2) = BF16(0.0_real32)
        data(7, 2) = BF16(0.0_real32)
        data(8, 2) = BF16(0.0_real32)
        data(9, 2) = BF16(0.0_real32)
        data(10, 2) = BF16(0.0_real32)
        data(:, 3) = 0.0_real32
        data(1, 4) = BF16(1.0_real32)
        data(2, 4) = BF16(2.0_real32)
        data(3, 4) = BF16(-1.0_real32)
        data(4, 4) = BF16(2.0_real32)
        data(5, 4) = BF16(1.0_real32)
        data(6, 4) = BF16(2.0_real32)
        data(7, 4) = BF16(-1.0_real32)
        data(8, 4) = BF16(2.0_real32)
        data(9, 4) = BF16(1.0_real32)
        data(10, 4) = BF16(2.0_real32)
        data(1, 5) = BF16(10.0_real32)
        data(2, 5) = BF16(0.0_real32)
        data(3, 5) = BF16(0.0_real32)
        data(4, 5) = BF16(0.0_real32)
        data(5, 5) = BF16(0.0_real32)
        data(6, 5) = BF16(0.0_real32)
        data(7, 5) = BF16(0.0_real32)
        data(8, 5) = BF16(0.0_real32)
        data(9, 5) = BF16(0.0_real32)
        data(10, 5) = BF16(0.0_real32)
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
END PROGRAM check_ibamax
