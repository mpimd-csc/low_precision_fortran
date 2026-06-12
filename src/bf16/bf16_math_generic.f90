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

submodule (lpf_bf16) lpf_bf16_math_generic
    use iso_c_binding
    use iso_fortran_env, only: real32, real64
    implicit none

    ! C Interfaces
    INTERFACE
        pure subroutine helper_abs_bf16(out, in) bind(c, name = "__bf16_helper_abs")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_epsilon_bf16(out) bind(c, name = "__bf16_helper_epsilon")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_digits_bf16(out) bind(c, name = "__bf16_helper_digits")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine


        pure subroutine helper_exponent_bf16(out, in) bind(c, name = "__bf16_helper_exponent")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_fraction_bf16(out, in) bind(c, name = "__bf16_helper_fraction")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_erfc_scaled_bf16(out, in) bind(c, name = "__bf16_helper_erfc_scaled")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_huge_bf16(out) bind(c, name = "__bf16_helper_huge")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_tiny_bf16(out) bind(c, name = "__bf16_helper_tiny")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
        end subroutine


        pure subroutine helper_minexponent_bf16(out) bind(c, name = "__bf16_helper_minexponent")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_maxexponent_bf16(out) bind(c, name = "__bf16_helper_maxexponent")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_mod_bf16(out, in1, in2) bind(c, name = "__bf16_helper_mod")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in1, in2
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_modulo_bf16(out, in1, in2) bind(c, name = "__bf16_helper_modulo")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in1, in2
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_nearest_bf16(out, in1, in2) bind(c, name = "__bf16_helper_nearest")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in1, in2
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_nint_bf16(out, in) bind(c, name = "__bf16_helper_nint")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_precision_bf16(out) bind(c, name = "__bf16_helper_precision")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_range_bf16(out) bind(c, name = "__bf16_helper_range")
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_scale_bf16(out, in, s) bind(c, name = "__bf16_helper_scale")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int), intent(in), value :: s
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_sign_bf16(out, in1, in2) bind(c, name = "__bf16_helper_sign")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_isinf_bf16(out, in) bind(c, name = "__bf16_helper_isinf")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_isnan_bf16(out, in) bind(c, name = "__bf16_helper_isnan")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

    END INTERFACE

contains
    module elemental function abs_bf16(x) result(abs_val)
        type(BF16), intent(in) :: x
        type(BF16) :: abs_val
        call helper_abs_bf16(abs_val%value, x%value)
    end function abs_bf16

    module elemental function digits_bf16(x) result(out)
        type(BF16), intent(in) :: x
        integer :: out
        integer(c_int) :: outi
        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        call helper_digits_bf16(outi)
        out = outi
    end function digits_bf16


    module elemental function epsilon_bf16(x) result(out)
        type(BF16), intent(in) :: x
        type(BF16) :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out%value = x%value
        end if

        call helper_epsilon_bf16(out%value)
    end function epsilon_bf16

    module elemental function exponent_bf16 (x) result(out)
        type(BF16), intent(in) :: x
        integer :: out
        integer(c_int) :: outi
        call helper_exponent_bf16(outi, x%value)
        out = outi
    end function exponent_bf16

    module elemental function fraction_bf16 (x) result(out)
        type(BF16), intent(in) :: x
        type(BF16) :: out

        call helper_fraction_bf16(out%value, x%value)
    end function fraction_bf16


    module elemental function radix_bf16(x) result(out)
        type(BF16), intent(in) :: x
        integer :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        out = 2
    end function radix_bf16

    module elemental function huge_bf16(x) result(out)
        type(BF16), intent(in) :: x
        type(BF16) :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out%value = INT(x%value, 2)
        end if

        call helper_huge_bf16(out%value)
    end function huge_bf16

    module elemental function tiny_bf16(x) result(out)
        type(BF16), intent(in) :: x
        type(BF16) :: out

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out%value = INT(x%value, 2)
        end if

        call helper_tiny_bf16(out%value)
    end function tiny_bf16


    module elemental function minexponent_bf16(x) result(out)
        type(BF16), intent(in) :: x
        integer :: out
        integer(c_int) :: outt

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        call helper_minexponent_bf16(outt)
        out = INT(outt)
    end function minexponent_bf16

    module elemental function maxexponent_bf16(x) result(out)
        type(BF16), intent(in) :: x
        integer :: out

        integer(c_int) :: outt
        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out= INT(x%value)
        end if

        call helper_maxexponent_bf16(outt)
        out = INT(outt)
    end function maxexponent_bf16


    module elemental function erfc_scaled_bf16 (x) result(out)
        type(BF16), intent(in) :: x
        type(BF16) :: out

        call helper_erfc_scaled_bf16(out%value, x%value)
    end function erfc_scaled_bf16

    module elemental function mod_bf16 (x, y) result(out)
        type(BF16), intent(in) :: x
        type(BF16), intent(in) :: y
        type(BF16) :: out

        call helper_mod_bf16(out%value, x%value, y%value)
    end function mod_bf16

    module elemental function modulo_bf16 (x, y) result(out)
        type(BF16), intent(in) :: x
        type(BF16), intent(in) :: y
        type(BF16) :: out

        call helper_modulo_bf16(out%value, x%value, y%value)
    end function modulo_bf16

    module elemental function nearest_bf16 (x, y) result(out)
        type(BF16), intent(in) :: x
        type(BF16), intent(in) :: y
        type(BF16) :: out

        call helper_nearest_bf16(out%value, x%value, y%value)
    end function nearest_bf16

    module elemental function nint_bf16 (x) result(out)
        type(BF16), intent(in) :: x
        type(BF16) :: out

        call helper_nint_bf16(out%value, x%value)
    end function nint_bf16

    module elemental function precision_bf16(x) result(out)
        type(BF16), intent(in) :: x
        integer :: out

        integer(c_int) :: outt
        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if

        call helper_precision_bf16(outt)
        out = INT(outt)
    end function precision_bf16

    module elemental function range_bf16(x) result(out)
        type(BF16), intent(in) :: x
        integer :: out
        integer(c_int) :: lout

        ! dummy to avoid compiler warnings
        if ( 0 .eq. 1) then
            out = INT(x%value)
        end if
        call helper_range_bf16(lout)
        out = lout
    end function range_bf16

    module elemental function scale_bf16_32(x, s) result(out)
        type(BF16), intent(in) :: x
        integer(int32), intent(in) :: s
        type(BF16) :: out
        integer(c_int) :: tmp
        tmp  = int(s, kind = c_int )
        call helper_scale_bf16(out%value, x%value, tmp)
    end function scale_bf16_32

    module elemental function scale_bf16_64(x, s) result(out)
        type(BF16), intent(in) :: x
        integer(int64), intent(in) :: s
        type(BF16) :: out
        integer(c_int) :: tmp
        tmp  = int(s, kind = c_int )
        call helper_scale_bf16(out%value, x%value, tmp)
    end function scale_bf16_64


    module elemental function sign_bf16(x, y) result(out)
        type(BF16), intent(in) :: x
        type(BF16), intent(in) :: y
        type(BF16) :: out

        call helper_sign_bf16(out%value, x%value, y%value)
    end function sign_bf16

    module elemental function isnan_bf16(x) result(out)
        type(BF16), intent(in) :: x
        logical :: out

        integer(c_int) :: outi
        call helper_isnan_bf16(outi, x%value)
        if ( outi.ne.0) then
            out = .TRUE.
        else
            out = .FALSE.
        end if
    end function

    module elemental function isinf_bf16(x) result(out)
        type(BF16), intent(in) :: x
        logical :: out

        integer(c_int) :: outi
        call helper_isinf_bf16(outi, x%value)
        if ( outi.ne.0) then
            out = .TRUE.
        else
            out = .FALSE.
        end if
    end function


end submodule
