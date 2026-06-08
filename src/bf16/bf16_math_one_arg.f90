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

submodule (lpf_bf16) lpf_bf16_math_one_arg
    use iso_c_binding
    use iso_fortran_env, only: real32, real64
    implicit none

    interface
        pure subroutine helper_bf16_acos(out, in) bind(c, name="__bf16_helper_acos")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_acosh(out, in) bind(c, name="__bf16_helper_acosh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_asin(out, in) bind(c, name="__bf16_helper_asin")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_asinh(out, in) bind(c, name="__bf16_helper_asinh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_atan(out, in) bind(c, name="__bf16_helper_atan")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_atanh(out, in) bind(c, name="__bf16_helper_atanh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_bessel_j0(out, in) bind(c, name="__bf16_helper_bessel_j0")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_bessel_j1(out, in) bind(c, name="__bf16_helper_bessel_j1")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_bessel_y0(out, in) bind(c, name="__bf16_helper_bessel_y0")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_bessel_y1(out, in) bind(c, name="__bf16_helper_bessel_y1")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_ceiling(out, in) bind(c, name="__bf16_helper_ceiling")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_cos(out, in) bind(c, name="__bf16_helper_cos")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_cosh(out, in) bind(c, name="__bf16_helper_cosh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_erf(out, in) bind(c, name="__bf16_helper_erf")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_erfc(out, in) bind(c, name="__bf16_helper_erfc")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_exp(out, in) bind(c, name="__bf16_helper_exp")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_floor(out, in) bind(c, name="__bf16_helper_floor")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_gamma(out, in) bind(c, name="__bf16_helper_gamma")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_log(out, in) bind(c, name="__bf16_helper_log")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_log10(out, in) bind(c, name="__bf16_helper_log10")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_log_gamma(out, in) bind(c, name="__bf16_helper_log_gamma")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_sin(out, in) bind(c, name="__bf16_helper_sin")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_sinh(out, in) bind(c, name="__bf16_helper_sinh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_sqrt(out, in) bind(c, name="__bf16_helper_sqrt")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_tan(out, in) bind(c, name="__bf16_helper_tan")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_bf16_tanh(out, in) bind(c, name="__bf16_helper_tanh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
    end interface
contains
    module elemental function acos_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_acos(out%value, in%value)
    end function
    module elemental function acosh_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_acosh(out%value, in%value)
    end function
    module elemental function asin_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_asin(out%value, in%value)
    end function
    module elemental function asinh_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_asinh(out%value, in%value)
    end function
    module elemental function atan_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_atan(out%value, in%value)
    end function
    module elemental function atanh_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_atanh(out%value, in%value)
    end function
    module elemental function bessel_j0_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_bessel_j0(out%value, in%value)
    end function
    module elemental function bessel_j1_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_bessel_j1(out%value, in%value)
    end function
    module elemental function bessel_y0_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_bessel_y0(out%value, in%value)
    end function
    module elemental function bessel_y1_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_bessel_y1(out%value, in%value)
    end function
    module elemental function ceiling_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_ceiling(out%value, in%value)
    end function
    module elemental function cos_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_cos(out%value, in%value)
    end function
    module elemental function cosh_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_cosh(out%value, in%value)
    end function
    module elemental function erf_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_erf(out%value, in%value)
    end function
    module elemental function erfc_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_erfc(out%value, in%value)
    end function
    module elemental function exp_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_exp(out%value, in%value)
    end function
    module elemental function floor_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_floor(out%value, in%value)
    end function
    module elemental function gamma_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_gamma(out%value, in%value)
    end function
    module elemental function log_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_log(out%value, in%value)
    end function
    module elemental function log10_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_log10(out%value, in%value)
    end function
    module elemental function log_gamma_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_log_gamma(out%value, in%value)
    end function
    module elemental function sin_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_sin(out%value, in%value)
    end function
    module elemental function sinh_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_sinh(out%value, in%value)
    end function
    module elemental function sqrt_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_sqrt(out%value, in%value)
    end function
    module elemental function tan_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_tan(out%value, in%value)
    end function
    module elemental function tanh_bf16(in) result(out)
            type(BF16) :: out
            type(BF16), intent(in) :: in

            call helper_bf16_tanh(out%value, in%value)
    end function
end submodule
