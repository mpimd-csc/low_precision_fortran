!  SPDX-License-Identifier LGPL-3.0-or-later
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

MODULE lpf_bf16_test_utils
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: check_bf16_real64
    PUBLIC :: check_logical
    PUBLIC :: check_integer
    PUBLIC :: check_bf16_known_bug
    PUBLIC :: test_summary
    PUBLIC :: integer_to_char
    PUBLIC :: BF16_TOL
    PUBLIC :: BF16_TOL_TIGHT
    PUBLIC :: BF16_TOL_LOOSE

    real(real64), parameter :: BF16_TOL      = 1.0d-2
    real(real64), parameter :: BF16_TOL_TIGHT = 1.0d-3
    real(real64), parameter :: BF16_TOL_LOOSE = 5.0d-2

    integer :: n_pass = 0
    integer :: n_fail = 0
    integer :: n_known_bug = 0

CONTAINS

    subroutine check_bf16_real64(name, got_bf16, expected_r64, tol)
        character(len=*), intent(in) :: name
        type(BF16), intent(in) :: got_bf16
        real(real64), intent(in) :: expected_r64
        real(real64), intent(in) :: tol

        real(real64) :: got_r64, diff, threshold

        got_r64 = dble(got_bf16)

        if (abs(expected_r64) < 1.0d-10) then
            threshold = 1.0d-5
            diff = abs(got_r64 - expected_r64)
        else
            threshold = max(abs(expected_r64) * tol, 1.0d-5)
            diff = abs(got_r64 - expected_r64)
        end if

        if (diff <= threshold) then
            n_pass = n_pass + 1
            write(*, '(A, T40, A)') '  TEST: ' // trim(name), 'PASS'
        else
            n_fail = n_fail + 1
            write(*, '(A, T40, A, ES10.3, A, ES10.3, A)') &
                '  TEST: ' // trim(name), 'FAIL (got=', got_r64, &
                ', expected=', expected_r64, ')'
        end if
    end subroutine

    subroutine check_logical(name, got, expected)
        character(len=*), intent(in) :: name
        logical, intent(in) :: got, expected

        if (got .eqv. expected) then
            n_pass = n_pass + 1
            write(*, '(A, T40, A)') '  TEST: ' // trim(name), 'PASS'
        else
            n_fail = n_fail + 1
            write(*, '(A, T40, A, L3, A, L3, A)') &
                '  TEST: ' // trim(name), 'FAIL (got=', got, ', expected=', expected, ')'
        end if
    end subroutine

    subroutine check_integer(name, got, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: got, expected

        if (got == expected) then
            n_pass = n_pass + 1
            write(*, '(A, T40, A)') '  TEST: ' // trim(name), 'PASS'
        else
            n_fail = n_fail + 1
            write(*, '(A, T40, A, I0, A, I0, A)') &
                '  TEST: ' // trim(name), 'FAIL (got=', got, ', expected=', expected, ')'
        end if
    end subroutine

    subroutine check_bf16_known_bug(name, got_bf16, expected_r64, tol, bug_id)
        character(len=*), intent(in) :: name
        type(BF16), intent(in) :: got_bf16
        real(real64), intent(in) :: expected_r64
        real(real64), intent(in) :: tol
        character(len=*), intent(in) :: bug_id

        real(real64) :: got_r64, diff, threshold

        got_r64 = dble(got_bf16)

        if (abs(expected_r64) < 1.0d-10) then
            threshold = 1.0d-5
            diff = abs(got_r64 - expected_r64)
        else
            threshold = max(abs(expected_r64) * tol, 1.0d-5)
            diff = abs(got_r64 - expected_r64)
        end if

        n_known_bug = n_known_bug + 1
        if (diff <= threshold) then
            n_pass = n_pass + 1
            write(*, '(A, T40, A)') '  TEST: ' // trim(name) // ' [BUG-' // trim(bug_id) // ']', &
                'PASS (bug appears fixed)'
        else
            write(*, '(A, T40, A, ES10.3, A, ES10.3, A)') &
                '  TEST: ' // trim(name) // ' [BUG-' // trim(bug_id) // ']', &
                'KNOWN BUG (got=', got_r64, ', expected=', expected_r64, ')'
        end if
    end subroutine

    subroutine test_summary()
        integer :: total
        total = n_pass + n_fail
        write(*, '(A)', advance='no') repeat('-', 60)
        write(*, *)
        write(*, '(A, I5, A)') '  Passed:     ', n_pass, &
            ' (including ' // trim(integer_to_char(n_known_bug)) // ' known bugs)'
        write(*, '(A, I5, A)') '  Failed:     ', n_fail
        write(*, '(A, I5)')    '  Total:      ', total
        write(*, '(A)', advance='no') repeat('-', 60)
        write(*, *)

        if (n_fail > 0) then
            write(*, '(A)') '  SOME TESTS FAILED'
            stop 1
        else
            write(*, '(A)') '  ALL TESTS PASSED'
            stop 0
        end if
    end subroutine

    function integer_to_char(n) result(str)
        integer, intent(in) :: n
        character(len=20) :: str
        write(str, *) n
        str = trim(ADJUSTL(str))
    end function

END MODULE lpf_bf16_test_utils
