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

submodule (lpf_fp8_e5m2) lpf_fp8_e5m2_math_generic
    use iso_c_binding
    use iso_fortran_env, only: real32, real64
    implicit none

    ! C Interfaces
    INTERFACE
        pure subroutine helper_abs_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_abs")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
            integer(c_int8_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_epsilon_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_epsilon")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_digits_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_digits")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine


        pure subroutine helper_exponent_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_exponent")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_fraction_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_fraction")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_erfc_scaled_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_erfc_scaled")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_huge_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_huge")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_tiny_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_tiny")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(out) :: out
        end subroutine


        pure subroutine helper_minexponent_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_minexponent")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_maxexponent_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_maxexponent")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_mod_fp8_e5m2(out, in1, in2) bind(c, name = "__fp8_e5m2_helper_mod")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in1, in2
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_modulo_fp8_e5m2(out, in1, in2) bind(c, name = "__fp8_e5m2_helper_modulo")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in1, in2
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_nearest_fp8_e5m2(out, in1, in2) bind(c, name = "__fp8_e5m2_helper_nearest")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in1, in2
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_nint_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_nint")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_precision_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_precision")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_range_fp8_e5m2(out) bind(c, name = "__fp8_e5m2_helper_range")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_scale_fp8_e5m2(out, in, s) bind(c, name = "__fp8_e5m2_helper_scale")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int), intent(in), value :: s
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_sign_fp8_e5m2(out, in1, in2) bind(c, name = "__fp8_e5m2_helper_sign")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in1
            integer(c_int8_t), intent(in), value :: in2
            integer(c_int8_t), intent(out) :: out
        end subroutine

        pure subroutine helper_isinf_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_isinf")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_isnan_fp8_e5m2(out, in) bind(c, name = "__fp8_e5m2_helper_isnan")
            use, intrinsic :: iso_c_binding
            integer(c_int8_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

    END INTERFACE

contains
    module elemental function abs_fp8_e5m2(x) result(abs_val)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: abs_val
        call helper_abs_fp8_e5m2(abs_val%value, x%value)
    end function abs_fp8_e5m2

    module elemental function digits_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out
        integer(c_int) :: outi
        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        call helper_digits_fp8_e5m2(outi)
        out = outi
    end function digits_fp8_e5m2


    module elemental function epsilon_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out%value = x%value
        end if

        call helper_epsilon_fp8_e5m2(out%value)
    end function epsilon_fp8_e5m2

    module elemental function exponent_fp8_e5m2 (x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out
        integer(c_int) :: outi
        call helper_exponent_fp8_e5m2(outi, x%value)
        out = outi
    end function exponent_fp8_e5m2

    module elemental function fraction_fp8_e5m2 (x) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: out

        call helper_fraction_fp8_e5m2(out%value, x%value)
    end function fraction_fp8_e5m2


    module elemental function radix_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        out = 2
    end function radix_fp8_e5m2

    module elemental function huge_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out%value = INT(x%value, 2)
        end if

        call helper_huge_fp8_e5m2(out%value)
    end function huge_fp8_e5m2

    module elemental function tiny_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out%value = INT(x%value, 2)
        end if

        call helper_tiny_fp8_e5m2(out%value)
    end function tiny_fp8_e5m2


    module elemental function minexponent_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out
        integer(c_int) :: outt

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        call helper_minexponent_fp8_e5m2(outt)
        out = INT(outt)
    end function minexponent_fp8_e5m2

    module elemental function maxexponent_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out

        integer(c_int) :: outt
        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out= INT(x%value)
        end if

        call helper_maxexponent_fp8_e5m2(outt)
        out = INT(outt)
    end function maxexponent_fp8_e5m2


    module elemental function erfc_scaled_fp8_e5m2 (x) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: out

        call helper_erfc_scaled_fp8_e5m2(out%value, x%value)
    end function erfc_scaled_fp8_e5m2

    module elemental function mod_fp8_e5m2 (x, y) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2), intent(in) :: y
        type(FP8_E5M2) :: out

        call helper_mod_fp8_e5m2(out%value, x%value, y%value)
    end function mod_fp8_e5m2

    module elemental function modulo_fp8_e5m2 (x, y) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2), intent(in) :: y
        type(FP8_E5M2) :: out

        call helper_modulo_fp8_e5m2(out%value, x%value, y%value)
    end function modulo_fp8_e5m2

    module elemental function nearest_fp8_e5m2 (x, y) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2), intent(in) :: y
        type(FP8_E5M2) :: out

        call helper_nearest_fp8_e5m2(out%value, x%value, y%value)
    end function nearest_fp8_e5m2

    module elemental function nint_fp8_e5m2 (x) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2) :: out

        call helper_nint_fp8_e5m2(out%value, x%value)
    end function nint_fp8_e5m2

    module elemental function precision_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out

        integer(c_int) :: outt
        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        call helper_precision_fp8_e5m2(outt)
        out = INT(outt, kind = lpf_default_int_kind)
    end function precision_fp8_e5m2

    module elemental function range_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer :: out
        integer(c_int) :: lout

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if
        call helper_range_fp8_e5m2(lout)
        out = lout
    end function range_fp8_e5m2

    module elemental function scale_fp8_e5m2_32(x, s) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer(lpf_int32_kind), intent(in) :: s
        type(FP8_E5M2) :: out
        integer(c_int) :: tmp
        tmp  = int(s, kind = c_int )
        call helper_scale_fp8_e5m2(out%value, x%value, tmp)
    end function scale_fp8_e5m2_32

    module elemental function scale_fp8_e5m2_64(x, s) result(out)
        type(FP8_E5M2), intent(in) :: x
        integer(lpf_int64_kind), intent(in) :: s
        type(FP8_E5M2) :: out
        integer(c_int) :: tmp
        tmp  = int(s, kind = c_int )
        call helper_scale_fp8_e5m2(out%value, x%value, tmp)
    end function scale_fp8_e5m2_64


    module elemental function sign_fp8_e5m2(x, y) result(out)
        type(FP8_E5M2), intent(in) :: x
        type(FP8_E5M2), intent(in) :: y
        type(FP8_E5M2) :: out

        call helper_sign_fp8_e5m2(out%value, x%value, y%value)
    end function sign_fp8_e5m2

    module elemental function isnan_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        logical :: out

        integer(c_int) :: outi
        call helper_isnan_fp8_e5m2(outi, x%value)
        if ( outi.ne.0) then
            out = .TRUE.
        else
            out = .FALSE.
        end if
    end function

    module elemental function isinf_fp8_e5m2(x) result(out)
        type(FP8_E5M2), intent(in) :: x
        logical :: out

        integer(c_int) :: outi
        call helper_isinf_fp8_e5m2(outi, x%value)
        if ( outi.ne.0) then
            out = .TRUE.
        else
            out = .FALSE.
        end if
    end function


end submodule
