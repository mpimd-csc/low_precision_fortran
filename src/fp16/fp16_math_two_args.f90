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

submodule (lpf_fp16) lpf_fp16_math_two_args
    use iso_c_binding
    use iso_fortran_env, only: real32, real64
    implicit none

    interface
        pure subroutine helper_fp16_atan2(out, in1, in2) bind(c, name="__fp16_helper_atan2")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine
        pure subroutine helper_fp16_bessel_jn(out, in1, in2) bind(c, name="__fp16_helper_bessel_jn")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine
        pure subroutine helper_fp16_bessel_yn(out, in1, in2) bind(c, name="__fp16_helper_bessel_yn")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine
        pure subroutine helper_fp16_hypot(out, in1, in2) bind(c, name="__fp16_helper_hypot")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine
    end interface
contains
    module elemental function atan2_fp16(in1, in2) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2

            call helper_fp16_atan2(out%value, in1%value, in2%value)
    end function

    module elemental function bessel_jn_fp16(in1, in2) result(out)
            type(FP16) :: out
            integer, intent(in) :: in1
            type(FP16), intent(in) :: in2
            integer(c_int) :: cin1
            cin1 = int(in1, c_int)
            call helper_fp16_bessel_jn(out%value, cin1, in2%value)
    end function

    module elemental function bessel_yn_fp16(in1, in2) result(out)
            type(FP16) :: out
            integer, intent(in) :: in1
            type(FP16), intent(in) :: in2
            integer(c_int) :: cin1

            cin1 = int(in1, c_int )
            call helper_fp16_bessel_yn(out%value, cin1, in2%value)
    end function

    module elemental function hypot_fp16(in1, in2) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2

            call helper_fp16_hypot(out%value, in1%value, in2%value)
    end function
end submodule
