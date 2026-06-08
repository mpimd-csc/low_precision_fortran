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

module lpf_xerbla
    use iso_c_binding
    use lpf_types
    implicit none

    interface
        subroutine xerbla_set_function_internal_old(func) bind(C, name = "lpf_blas_xerbla_set_function_fortran_old")
            use, intrinsic :: iso_c_binding
            type(c_funptr), value, intent(in) :: func
        end subroutine
        subroutine xerbla_set_function_internal_new(func) bind(C, name = "lpf_blas_xerbla_set_function_fortran_new")
            use, intrinsic :: iso_c_binding
            type(c_funptr), value, intent(in) :: func
        end subroutine

    end interface

    abstract interface
        subroutine xerbla_function(msg, info)
            use, intrinsic :: iso_c_binding
            import :: lpf_default_int_kind
            character(len=:), pointer, intent(in) :: msg
            integer(lpf_default_int_kind), intent(in) :: info
        end subroutine
    end interface
contains
    subroutine xerbla_set_function_f77(func)
        procedure(xerbla_function) :: func
        type(c_funptr) :: cfun

        cfun = c_funloc(func)
        call xerbla_set_function_internal_old(cfun)

    end subroutine

    subroutine xerbla_set_function_f90(func)
        procedure(xerbla_function) :: func

        type(c_funptr) :: cfun

        cfun = c_funloc(func)
        call xerbla_set_function_internal_new(cfun)

    end subroutine



end module

