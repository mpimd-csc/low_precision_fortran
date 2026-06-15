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

PROGRAM check_hscal
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_FP16
    USE lpf_blas_fp16
    IMPLICIT NONE

    LOGICAL :: ok
    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HSCAL tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 6
        INTEGER(int64), dimension(NTEST) :: n = [5, 3, 4, 4, 1, 0]
        INTEGER(int64), dimension(NTEST) :: incx = [1, 1, 1, 2, 1, 1]
        TYPE(FP16), dimension(10, NTEST) :: data
        TYPE(FP16), dimension(NTEST) :: sa
        TYPE(FP16), dimension(10, NTEST) :: expected
        INTEGER :: k, i
        LOGICAL :: lok

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
        ! Testcase 4: {1, 2, -1, 2, 1, 2, -1, 2, 1, 2}
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

        sa(1) = FP16(0.1_real32)
        sa(2) = FP16(0.0625_real32)
        sa(3) = FP16(1.0_real32)
        sa(4) = FP16(-4.0_real32)
        sa(5) = FP16(10.0_real32)
        sa(6) = FP16(0.0_real32)

        expected(1, 1) = FP16(0.1_real32)
        expected(2, 1) = FP16(-0.2_real32)
        expected(3, 1) = FP16(0.3_real32)
        expected(4, 1) = FP16(-0.4_real32)
        expected(5, 1) = FP16(0.5_real32)
        expected(6, 1) = FP16(0.0_real32)
        expected(7, 1) = FP16(0.0_real32)
        expected(8, 1) = FP16(0.0_real32)
        expected(9, 1) = FP16(0.0_real32)
        expected(10, 1) = FP16(0.0_real32)
        expected(1, 2) = FP16(0.09375_real32)
        expected(2, 2) = FP16(-0.15625_real32)
        expected(3, 2) = FP16(0.21875_real32)
        expected(4, 2) = FP16(0.0_real32)
        expected(5, 2) = FP16(0.0_real32)
        expected(6, 2) = FP16(0.0_real32)
        expected(7, 2) = FP16(0.0_real32)
        expected(8, 2) = FP16(0.0_real32)
        expected(9, 2) = FP16(0.0_real32)
        expected(10, 2) = FP16(0.0_real32)
        expected(:, 3) = 0.0_real32
        expected(1, 4) = FP16(-4.0_real32)
        expected(2, 4) = FP16(2.0_real32)
        expected(3, 4) = FP16(4.0_real32)
        expected(4, 4) = FP16(2.0_real32)
        expected(5, 4) = FP16(-4.0_real32)
        expected(6, 4) = FP16(2.0_real32)
        expected(7, 4) = FP16(4.0_real32)
        expected(8, 4) = FP16(2.0_real32)
        expected(9, 4) = FP16(1.0_real32)
        expected(10, 4) = FP16(2.0_real32)
        expected(1, 5) = FP16(100.0_real32)
        expected(2, 5) = FP16(0.0_real32)
        expected(3, 5) = FP16(0.0_real32)
        expected(4, 5) = FP16(0.0_real32)
        expected(5, 5) = FP16(0.0_real32)
        expected(6, 5) = FP16(0.0_real32)
        expected(7, 5) = FP16(0.0_real32)
        expected(8, 5) = FP16(0.0_real32)
        expected(9, 5) = FP16(0.0_real32)
        expected(10, 5) = FP16(0.0_real32)
        expected(1, 6) = FP16(9.0_real32)
        expected(2, 6) = FP16(9.0_real32)
        expected(3, 6) = FP16(9.0_real32)
        expected(4, 6) = FP16(9.0_real32)
        expected(5, 6) = FP16(9.0_real32)
        expected(6, 6) = FP16(9.0_real32)
        expected(7, 6) = FP16(9.0_real32)
        expected(8, 6) = FP16(9.0_real32)
        expected(9, 6) = FP16(9.0_real32)
        expected(10, 6) = FP16(9.0_real32)

        do i = 1, NTEST
            CALL scal(n(i), sa(i), data(:, i), incx(i))

            lok = .TRUE.
            do k = 1, 10
                IF (abs(dble(data(k, i)) - dble(expected(k, i))) > 1.0d-2 * max(abs(dble(expected(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I3, A, I3, A, I2, A, F10.6, A, F10.6)", &
                          "HSCAL   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                          ", K = ", k, ", Result = ", dble(data(k, i)), &
                          ", Expected = ", dble(expected(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
            end do
            IF (lok) THEN
                print "(A, I2, A, I3, A, I3, A, F10.6)", &
                      "HSCAL   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", SA = ", dble(sa(i))
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_hscal
