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

PROGRAM check_brotmg
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    USE lpf_blas_bf16
    IMPLICIT NONE

    LOGICAL :: ok
    ok = .TRUE.

    CALL run_tests(ok)

    IF (.NOT. ok) THEN
        PRINT *, "HROTMG tests failed!"
        STOP 1
    END IF

CONTAINS

    SUBROUTINE run_tests(ok)
        LOGICAL, INTENT(INOUT) :: ok
        INTEGER, PARAMETER :: NTEST = 4
        TYPE(BF16), dimension(4, NTEST) :: dab
        TYPE(BF16), dimension(9, NTEST) :: dtrue
        TYPE(BF16), dimension(9) :: dtemp
        INTEGER :: i, k
        LOGICAL :: lok

        dab(1, 1) = BF16(0.1_real32)
        dab(2, 1) = BF16(0.3_real32)
        dab(3, 1) = BF16(1.2_real32)
        dab(4, 1) = BF16(0.2_real32)
        dab(1, 2) = BF16(0.7_real32)
        dab(2, 2) = BF16(0.2_real32)
        dab(3, 2) = BF16(0.6_real32)
        dab(4, 2) = BF16(4.2_real32)
        dab(:, 3) = 0.0_real32
        dab(1, 4) = BF16(4.0_real32)
        dab(2, 4) = BF16(-1.0_real32)
        dab(3, 4) = BF16(2.0_real32)
        dab(4, 4) = BF16(4.0_real32)

        dtrue(1, 1) = BF16(0.0_real32)
        dtrue(2, 1) = BF16(0.0_real32)
        dtrue(3, 1) = BF16(1.3_real32)
        dtrue(4, 1) = BF16(0.2_real32)
        dtrue(5, 1) = BF16(0.0_real32)
        dtrue(6, 1) = BF16(0.0_real32)
        dtrue(7, 1) = BF16(0.0_real32)
        dtrue(8, 1) = BF16(0.5_real32)
        dtrue(9, 1) = BF16(0.0_real32)
        dtrue(1, 2) = BF16(0.0_real32)
        dtrue(2, 2) = BF16(0.0_real32)
        dtrue(3, 2) = BF16(4.5_real32)
        dtrue(4, 2) = BF16(4.2_real32)
        dtrue(5, 2) = BF16(1.0_real32)
        dtrue(6, 2) = BF16(0.5_real32)
        dtrue(7, 2) = BF16(0.0_real32)
        dtrue(8, 2) = BF16(0.0_real32)
        dtrue(9, 2) = BF16(0.0_real32)
        dtrue(1, 3) = BF16(0.0_real32)
        dtrue(2, 3) = BF16(0.0_real32)
        dtrue(3, 3) = BF16(0.0_real32)
        dtrue(4, 3) = BF16(0.0_real32)
        dtrue(5, 3) = BF16(-2.0_real32)
        dtrue(6, 3) = BF16(0.0_real32)
        dtrue(7, 3) = BF16(0.0_real32)
        dtrue(8, 3) = BF16(0.0_real32)
        dtrue(9, 3) = BF16(0.0_real32)
        dtrue(1, 4) = BF16(0.0_real32)
        dtrue(2, 4) = BF16(0.0_real32)
        dtrue(3, 4) = BF16(0.0_real32)
        dtrue(4, 4) = BF16(4.0_real32)
        dtrue(5, 4) = BF16(-1.0_real32)
        dtrue(6, 4) = BF16(0.0_real32)
        dtrue(7, 4) = BF16(0.0_real32)
        dtrue(8, 4) = BF16(0.0_real32)
        dtrue(9, 4) = BF16(0.0_real32)

        ! Initialize computed values for DTRUE
        dtrue(1, 1) = BF16(12.0_real32 / 130.0_real32)
        dtrue(2, 1) = BF16(36.0_real32 / 130.0_real32)
        dtrue(7, 1) = BF16(-1.0_real32 / 6.0_real32)
        dtrue(1, 2) = BF16(14.0_real32 / 75.0_real32)
        dtrue(2, 2) = BF16(49.0_real32 / 75.0_real32)
        dtrue(9, 2) = BF16(1.0_real32 / 7.0_real32)

        do i = 1, NTEST
            dtemp = 0.0_real32
            do k = 1, 4
                dtemp(k) = dab(k, i)
            end do

            CALL rotmg(dtemp(1), dtemp(2), dtemp(3), dtemp(4), dtemp(5))

            lok = .TRUE.
            do k = 1, 9
                IF (abs(dble(dtemp(k)) - dble(dtrue(k, i))) > 1.0d-2 * max(abs(dble(dtrue(k, i))), 1.0d-5)) THEN
                    print "(A, I2, A, I2, A, F10.6, A, F10.6)", &
                          "HROTMG  -- FAIL -- Testcase ", i, ": Value ", k, ": Result = ", dble(dtemp(k)), &
                          ", Expected = ", dble(dtrue(k, i))
                    ok = .FALSE.
                    lok = .FALSE.
                END IF
            end do
            IF (lok) THEN
                print "(A, I2, A)", "HROTMG  -- PASS -- Testcase ", i, "."
            END IF
        end do
    end SUBROUTINE
END PROGRAM check_brotmg
