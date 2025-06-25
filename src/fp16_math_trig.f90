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

submodule (fp16_support) fp16_math_trig
    use iso_c_binding
    use iso_fortran_env
    implicit none

    ! C Interfaces
    INTERFACE
        pure subroutine helper_cosd_fp16(out, in) bind(c, name = "__fp16_helper_cosd")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_sind_fp16(out, in) bind(c, name = "__fp16_helper_sind")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_tand_fp16(out, in) bind(c, name = "__fp16_helper_tand")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_cotan_fp16(out, in) bind(c, name = "__fp16_helper_cotan")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_cotand_fp16(out, in) bind(c, name = "__fp16_helper_cotand")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_acosd_fp16(out, in) bind(c, name = "__fp16_helper_acosd")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_asind_fp16(out, in) bind(c, name = "__fp16_helper_asind")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_atand_fp16(out, in) bind(c, name = "__fp16_helper_atand")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_atan2d_fp16(out, in1, in2) bind(c, name = "__fp16_helper_atan2d")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine




    END INTERFACE

contains
    module elemental function cosd_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_cosd_fp16(ret_val%value, x%value)
    end function cosd_fp16

    module elemental function sind_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_sind_fp16(ret_val%value, x%value)
    end function sind_fp16

    module elemental function tand_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_tand_fp16(ret_val%value, x%value)
    end function tand_fp16

    module elemental function cotan_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_cotan_fp16(ret_val%value, x%value)
    end function cotan_fp16

    module elemental function cotand_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_cotand_fp16(ret_val%value, x%value)
    end function cotand_fp16

    module elemental function acosd_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_acosd_fp16(ret_val%value, x%value)
    end function acosd_fp16

    module elemental function asind_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_asind_fp16(ret_val%value, x%value)
    end function asind_fp16

    module elemental function atand_fp16(x) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16) :: ret_val
        call helper_atand_fp16(ret_val%value, x%value)
    end function atand_fp16

    module elemental function atan2d_fp16(x,y) result(ret_val)
        type(FP16), intent(in) :: x
        type(FP16), intent(in) :: y
        type(FP16) :: ret_val
        call helper_atan2d_fp16(ret_val%value, x%value, y%value)
    end function atan2d_fp16

end submodule fp16_math_trig
