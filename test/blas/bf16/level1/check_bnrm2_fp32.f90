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

PROGRAM check_bnrm2_fp32
    USE iso_fortran_env, only: real32, real64, int64
    USE LPF_BF16
    USE lpf_blas_bf16
    IMPLICIT NONE

    LOGICAL :: ok
    INTEGER :: i

    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HNRM2_FP32 tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 6
        INTEGER(c_int64_t), dimension(NTEST) :: n = [5, 3, 4, 4, 1, 0]
        INTEGER(c_int64_t), dimension(NTEST) :: incx = [1, 1, 1, 2, 1, 1]
        TYPE(BF16), dimension(10, NTEST) :: data
        TYPE(BF16), dimension(NTEST) :: bnrm2_expected
        TYPE(BF16) :: result

        ! Initialize data
        ! Testcase 1: {1.0, -2.0, 3.0, -4.0, 5.0}
        data(1, 1) = BF16(1.0_real32)
        data(2, 1) = BF16(-2.0_real32)
        data(3, 1) = BF16(3.0_real32)
        data(4, 1) = BF16(-4.0_real32)
        data(5, 1) = BF16(5.0_real32)
        data(6, 1) = BF16(0.0_real32)
        data(7, 1) = BF16(0.0_real32)
        data(8, 1) = BF16(0.0_real32)
        data(9, 1) = BF16(0.0_real32)
        data(10, 1) = BF16(0.0_real32)
        ! Testcase 2: {1.5, -2.5, 3.5}
        data(1, 2) = BF16(1.5_real32)
        data(2, 2) = BF16(-2.5_real32)
        data(3, 2) = BF16(3.5_real32)
        data(4, 2) = BF16(0.0_real32)
        data(5, 2) = BF16(0.0_real32)
        data(6, 2) = BF16(0.0_real32)
        data(7, 2) = BF16(0.0_real32)
        data(8, 2) = BF16(0.0_real32)
        data(9, 2) = BF16(0.0_real32)
        data(10, 2) = BF16(0.0_real32)
        ! Testcase 3: {0, 0, 0, 0}
        data(:, 3) = 0.0_real32
        ! Testcase 4: {1, 2, 1, 2, 1, 2, 1, 2, 1, 2}
        data(1, 4) = BF16(1.0_real32)
        data(2, 4) = BF16(2.0_real32)
        data(3, 4) = BF16(1.0_real32)
        data(4, 4) = BF16(2.0_real32)
        data(5, 4) = BF16(1.0_real32)
        data(6, 4) = BF16(2.0_real32)
        data(7, 4) = BF16(1.0_real32)
        data(8, 4) = BF16(2.0_real32)
        data(9, 4) = BF16(1.0_real32)
        data(10, 4) = BF16(2.0_real32)
        ! Testcase 5: {10}
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
        ! Testcase 6: {9...} but N=0
        data(:, 6) = 9.0_real32

        bnrm2_expected(1) = BF16(7.416198_real32)
        bnrm2_expected(2) = BF16(4.5552_real32)
        bnrm2_expected(3) = BF16(0.0_real32)
        bnrm2_expected(4) = BF16(2.0_real32)
        bnrm2_expected(5) = BF16(10.0_real32)
        bnrm2_expected(6) = BF16(0.0_real32)

        do i = 1, NTEST
            result = nrm2_fp32(n(i), data(:, i), incx(i))
            if (abs(dble(result) - dble(bnrm2_expected(i))) < 1.0d-2) then
                print "(A, I2, A, I3, A, I3, A, F10.6, A, F10.6)", &
                      "HNRM2_FP32   -- PASS -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", Result = ", dble(result), ", Expected = ", dble(bnrm2_expected(i))
            else
                print "(A, I2, A, I3, A, I3, A, F10.6, A, F10.6)", &
                      "HNRM2_FP32   -- FAIL -- Testcase ", i, ": N = ", n(i), ", INCX = ", incx(i), &
                      ", Result = ", dble(result), ", Expected = ", dble(bnrm2_expected(i))
                ok = .FALSE.
            end if
        end do
    end SUBROUTINE
END PROGRAM check_bnrm2_fp32
