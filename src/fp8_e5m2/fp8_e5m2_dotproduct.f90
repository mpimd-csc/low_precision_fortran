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
!

submodule (lpf_fp8_e5m2) lpf_fp8_e5m2_dot_product
    use iso_fortran_env
    implicit none

    contains
        module pure function fp8_e5m2_dot_product(a, b) result(c)
            type(fp8_e5m2), intent(in) :: a(:)
            type(fp8_e5m2), intent(in) :: b(:)
            type(fp8_e5m2) :: c

            integer(int64) ::  i

            if ( size(a,1) .ne. size(b,1) ) error stop 'dot_product: dimension of A and B not matching'

            c = 0.0
            do i = 1, size(a,1)
                c = c +  a(i) * b(i)
            end do

        end function

end submodule
