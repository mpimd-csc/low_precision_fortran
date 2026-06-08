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

submodule (lpf_fp8_e5m2) lpf_fp8_e5m2_math_one_arg
    use iso_c_binding
    use iso_fortran_env, only: real32, real64
    implicit none

    interface
        pure subroutine helper_fp8_e5m2_acos(out, in) bind(c, name="__fp8_e5m2_helper_acos")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_acosh(out, in) bind(c, name="__fp8_e5m2_helper_acosh")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_asin(out, in) bind(c, name="__fp8_e5m2_helper_asin")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_asinh(out, in) bind(c, name="__fp8_e5m2_helper_asinh")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_atan(out, in) bind(c, name="__fp8_e5m2_helper_atan")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_atanh(out, in) bind(c, name="__fp8_e5m2_helper_atanh")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_bessel_j0(out, in) bind(c, name="__fp8_e5m2_helper_bessel_j0")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_bessel_j1(out, in) bind(c, name="__fp8_e5m2_helper_bessel_j1")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_bessel_y0(out, in) bind(c, name="__fp8_e5m2_helper_bessel_y0")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_bessel_y1(out, in) bind(c, name="__fp8_e5m2_helper_bessel_y1")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_ceiling(out, in) bind(c, name="__fp8_e5m2_helper_ceiling")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_cos(out, in) bind(c, name="__fp8_e5m2_helper_cos")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_cosh(out, in) bind(c, name="__fp8_e5m2_helper_cosh")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_erf(out, in) bind(c, name="__fp8_e5m2_helper_erf")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_erfc(out, in) bind(c, name="__fp8_e5m2_helper_erfc")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_exp(out, in) bind(c, name="__fp8_e5m2_helper_exp")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_floor(out, in) bind(c, name="__fp8_e5m2_helper_floor")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_gamma(out, in) bind(c, name="__fp8_e5m2_helper_gamma")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_log(out, in) bind(c, name="__fp8_e5m2_helper_log")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_log10(out, in) bind(c, name="__fp8_e5m2_helper_log10")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_log_gamma(out, in) bind(c, name="__fp8_e5m2_helper_log_gamma")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_sin(out, in) bind(c, name="__fp8_e5m2_helper_sin")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_sinh(out, in) bind(c, name="__fp8_e5m2_helper_sinh")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_sqrt(out, in) bind(c, name="__fp8_e5m2_helper_sqrt")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_tan(out, in) bind(c, name="__fp8_e5m2_helper_tan")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp8_e5m2_tanh(out, in) bind(c, name="__fp8_e5m2_helper_tanh")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine
    end interface
contains
    module elemental function acos_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_acos(out%value, in%value)
    end function
    module elemental function acosh_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_acosh(out%value, in%value)
    end function
    module elemental function asin_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_asin(out%value, in%value)
    end function
    module elemental function asinh_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_asinh(out%value, in%value)
    end function
    module elemental function atan_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_atan(out%value, in%value)
    end function
    module elemental function atanh_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_atanh(out%value, in%value)
    end function
    module elemental function bessel_j0_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_bessel_j0(out%value, in%value)
    end function
    module elemental function bessel_j1_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_bessel_j1(out%value, in%value)
    end function
    module elemental function bessel_y0_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_bessel_y0(out%value, in%value)
    end function
    module elemental function bessel_y1_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_bessel_y1(out%value, in%value)
    end function
    module elemental function ceiling_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_ceiling(out%value, in%value)
    end function
    module elemental function cos_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_cos(out%value, in%value)
    end function
    module elemental function cosh_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_cosh(out%value, in%value)
    end function
    module elemental function erf_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_erf(out%value, in%value)
    end function
    module elemental function erfc_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_erfc(out%value, in%value)
    end function
    module elemental function exp_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_exp(out%value, in%value)
    end function
    module elemental function floor_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_floor(out%value, in%value)
    end function
    module elemental function gamma_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_gamma(out%value, in%value)
    end function
    module elemental function log_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_log(out%value, in%value)
    end function
    module elemental function log10_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_log10(out%value, in%value)
    end function
    module elemental function log_gamma_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_log_gamma(out%value, in%value)
    end function
    module elemental function sin_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_sin(out%value, in%value)
    end function
    module elemental function sinh_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_sinh(out%value, in%value)
    end function
    module elemental function sqrt_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_sqrt(out%value, in%value)
    end function
    module elemental function tan_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_tan(out%value, in%value)
    end function
    module elemental function tanh_fp8_e5m2(in) result(out)
            type(FP8_E5M2) :: out
            type(FP8_E5M2), intent(in) :: in

            call helper_fp8_e5m2_tanh(out%value, in%value)
    end function
end submodule
