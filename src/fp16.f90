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
MODULE FP16_SUPPORT
    USE iso_c_binding
    USE iso_fortran_env
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: FP16, ASSIGNMENT(=)
    PUBLIC :: write(formatted)
    PUBLIC :: operator(+), operator(-), operator(*), operator(/)
    PUBLIC :: operator(.lt.), operator(.le.), operator(.gt.), operator(.ge.)
    PUBLIC :: operator(.eq.), operator(.ne.)
    PUBLIC :: operator(**)
    PUBLIC :: abs
    PUBLIC :: acos
    PUBLIC :: acosh
    PUBLIC :: asin
    PUBLIC :: asinh
    PUBLIC :: atan
    PUBLIC :: atanh
    PUBLIC :: bessel_j0
    PUBLIC :: bessel_j1
    PUBLIC :: bessel_y0
    PUBLIC :: bessel_y1
    PUBLIC :: ceiling
    PUBLIC :: cos
    PUBLIC :: cosh
    PUBLIC :: erf
    PUBLIC :: erfc
    PUBLIC :: exp
    PUBLIC :: floor
    PUBLIC :: gamma
    PUBLIC :: log
    PUBLIC :: log10
    PUBLIC :: log_gamma
    PUBLIC :: sin
    PUBLIC :: sinh
    PUBLIC :: sqrt
    PUBLIC :: tan
    PUBLIC :: tanh
    PUBLIC :: atan2
    PUBLIC :: bessel_jn
    PUBLIC :: bessel_yn
    PUBLIC :: hypot
    PUBLIC :: radix
    PUBLIC :: epsilon
    PUBLIC :: exponent
    PUBLIC :: fraction
    PUBLIC :: cosd
    PUBLIC :: sind
    PUBLIC :: tand
    PUBLIC :: cotan
    PUBLIC :: cotand
    PUBLIC :: acosd
    PUBLIC :: asind
    PUBLIC :: atand
    PUBLIC :: atan2d
    PUBLIC :: erfc_scaled
    PUBLIC :: huge
    PUBLIC :: tiny
    PUBLIC :: minexponent
    PUBLIC :: maxexponent
    PUBLIC :: mod
    PUBLIC :: modulo
    PUBLIC :: nearest
    PUBLIC :: nint
    PUBLIC :: precision
    PUBLIC :: range
    PUBLIC :: scale
    PUBLIC :: sign
    PUBLIC :: maxval
    PUBLIC :: maxloc
    PUBLIC :: minval
    PUBLIC :: minloc
    PUBLIC :: isnan
    PUBLIC :: isinf
    PUBLIC :: min
    PUBLIC :: max

    TYPE, BIND(C) :: FP16
        INTEGER(c_int16_t) :: value
    END TYPE

    interface FP16
        module procedure :: construct_real
        module procedure :: construct_double
        module procedure :: construct_int
    end interface

    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE assign_int, assign_real, assign_double, assign_fp16, assign_to_real
    END INTERFACE

    INTERFACE write(formatted)
        MODULE PROCEDURE  write_formatted
    END INTERFACE

    ! Interface für den + Operator
    interface operator(+)
        module procedure add_fp16_fp16, add_fp16_real, add_real_fp16
    end interface operator(+)

    ! Interface für den - Operator
    interface operator(-)
        module procedure subtract_fp16_fp16, subtract_fp16_real, subtract_real_fp16
    end interface operator(-)

    ! Interface für den * Operator
    interface operator(*)
        module procedure multiply_fp16_fp16, multiply_fp16_real, multiply_real_fp16
    end interface operator(*)

    ! Interface für den / Operator
    interface operator(/)
        module procedure divide_fp16_fp16, divide_fp16_real, divide_real_fp16
    end interface operator(/)

    ! Interface für den ** Operator
    interface operator(**)
        module procedure power_fp16_fp16, power_fp16_real, power_fp16_int
    end interface operator(**)

    ! .lt. operator
    interface operator(.lt.)
        module procedure lt_fp16_fp16, lt_fp16_fp32, lt_fp32_fp16
    end interface operator(.lt.)

    ! .le. operator
    interface operator(.le.)
        module procedure le_fp16_fp16, le_fp16_fp32, le_fp32_fp16
    end interface operator(.le.)

    ! .gt. operator
    interface operator(.gt.)
        module procedure gt_fp16_fp16, gt_fp16_fp32, gt_fp32_fp16
    end interface operator(.gt.)

    ! .ge. operator
    interface operator(.ge.)
        module procedure ge_fp16_fp16, ge_fp16_fp32, ge_fp32_fp16
    end interface operator(.ge.)

    ! .eq. operator
    interface operator(.eq.)
        module procedure eq_fp16_fp16, eq_fp16_fp32, eq_fp32_fp16
    end interface operator(.eq.)

    ! .ne. operator
    interface operator(.ne.)
        module procedure ne_fp16_fp16, ne_fp16_fp32, ne_fp32_fp16
    end interface operator(.ne.)



    ! Math functions
    interface abs
        module elemental function abs_fp16(x) result(abs_val)
            type(FP16), intent(in) :: x
            type(FP16) :: abs_val
        end function abs_fp16
    end interface

    interface epsilon
        module elemental function epsilon_fp16(x) result(abs_val)
            type(FP16), intent(in) :: x
            type(FP16) :: abs_val
        end function epsilon_fp16
    end interface

    interface radix
        module elemental function radix_fp16(x) result(out)
            type(FP16), intent(in) :: x
            integer :: out
        end function radix_fp16
    end interface

    interface exponent
        module elemental function exponent_fp16(x) result(out)
            type(FP16), intent(in) :: x
            integer :: out
        end function exponent_fp16
    end interface

    interface fraction
        module elemental function fraction_fp16(x) result(out)
            type(FP16), intent(in) :: x
            type(FP16) :: out
        end function fraction_fp16
    end interface

    interface erfc_scaled
        module elemental function erfc_scaled_fp16(x) result(out)
            type(FP16), intent(in) :: x
            type(FP16) :: out
        end function erfc_scaled_fp16
    end interface

    interface huge
        module elemental function huge_fp16(x) result(out)
            type(FP16), intent(in) :: x
            type(FP16) :: out
        end function huge_fp16
    end interface

    interface tiny
        module elemental function tiny_fp16(x) result(out)
            type(FP16), intent(in) :: x
            type(FP16) :: out
        end function tiny_fp16
    end interface


    interface minexponent
        module elemental function minexponent_fp16(x) result(out)
            type(FP16), intent(in) :: x
            integer :: out
        end function minexponent_fp16
    end interface

    interface maxexponent
        module elemental function maxexponent_fp16(x) result(out)
            type(FP16), intent(in) :: x
            integer :: out
        end function maxexponent_fp16
    end interface

    interface mod
        module elemental function mod_fp16(x, y) result(out)
            type(FP16), intent(in) :: x
            type(FP16), intent(in) :: y
            type(FP16) :: out
        end function mod_fp16
    end interface

    interface modulo
        module elemental function modulo_fp16(x, y) result(out)
            type(FP16), intent(in) :: x
            type(FP16), intent(in) :: y
            type(FP16) :: out
        end function modulo_fp16
    end interface

    interface nearest
        module elemental function nearest_fp16(x, y) result(out)
            type(FP16), intent(in) :: x
            type(FP16), intent(in) :: y
            type(FP16) :: out
        end function nearest_fp16
    end interface

    interface nint
        module elemental function nint_fp16(x) result(out)
            type(FP16), intent(in) :: x
            type(FP16) :: out
        end function nint_fp16
    end interface

    interface precision
        module elemental function precision_fp16(x) result(out)
            type(FP16), intent(in) :: x
            integer :: out
        end function precision_fp16
    end interface

    interface range
        module elemental function range_fp16(x) result(out)
            type(FP16), intent(in) :: x
            integer :: out
        end function range_fp16
    end interface

    interface scale
        module elemental function scale_fp16(x, s) result(out)
            type(FP16), intent(in) :: x
            integer, intent(in) :: s
            type(FP16) :: out
        end function scale_fp16
    end interface

    interface sign
        module elemental function sign_fp16(x, y) result(out)
            type(FP16), intent(in) :: x
            type(FP16), intent(in) :: y
            type(FP16) :: out
        end function sign_fp16
    end interface

    interface isnan
        module elemental function isnan_fp16(x) result(out)
            type(FP16), intent(in) :: x
            logical :: out
        end function isnan_fp16
    end interface

    interface isinf
        module elemental function isinf_fp16(x) result(out)
            type(FP16), intent(in) :: x
            logical :: out
        end function isinf_fp16
    end interface

    interface min
        module procedure min_fp16
        module procedure min3_fp16
    end interface min

    interface max
        module procedure max_fp16
        module procedure max3_fp16
    end interface max

    interface acos
        module elemental function acos_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function acos_fp16
    end interface
    interface acosh
        module elemental function acosh_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function acosh_fp16
    end interface
    interface asin
        module elemental function asin_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function asin_fp16
    end interface
    interface asinh
        module elemental function asinh_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function asinh_fp16
    end interface
    interface atan
        module elemental function atan_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function atan_fp16
    end interface
    interface atanh
        module elemental function atanh_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function atanh_fp16
    end interface
    interface bessel_j0
        module elemental function bessel_j0_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function bessel_j0_fp16
    end interface
    interface bessel_j1
        module elemental function bessel_j1_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function bessel_j1_fp16
    end interface
    interface bessel_y0
        module elemental function bessel_y0_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function bessel_y0_fp16
    end interface
    interface bessel_y1
        module elemental function bessel_y1_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function bessel_y1_fp16
    end interface
    interface ceiling
        module elemental function ceiling_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function ceiling_fp16
    end interface
    interface cos
        module elemental function cos_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function cos_fp16
    end interface
    interface cosh
        module elemental function cosh_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function cosh_fp16
    end interface
    interface erf
        module elemental function erf_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function erf_fp16
    end interface
    interface erfc
        module elemental function erfc_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function erfc_fp16
    end interface
    interface exp
        module elemental function exp_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function exp_fp16
    end interface
    interface floor
        module elemental function floor_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function floor_fp16
    end interface
    interface gamma
        module elemental function gamma_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function gamma_fp16
    end interface
    interface log
        module elemental function log_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function log_fp16
    end interface
    interface log10
        module elemental function log10_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function log10_fp16
    end interface
    interface log_gamma
        module elemental function log_gamma_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function log_gamma_fp16
    end interface
    interface sin
        module elemental function sin_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function sin_fp16
    end interface
    interface sinh
        module elemental function sinh_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function sinh_fp16
    end interface
    interface sqrt
        module elemental function sqrt_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function sqrt_fp16
    end interface
    interface tan
        module elemental function tan_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function tan_fp16
    end interface
    interface tanh
        module elemental function tanh_fp16(in) result(out)
            type(FP16), intent(in) :: in
            type(FP16) :: out
        end function tanh_fp16
    end interface

    interface atan2
        module elemental function atan2_fp16(in1, in2) result(out)
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2
            type(FP16) :: out
        end function atan2_fp16
    end interface
    interface bessel_jn
        module elemental function bessel_jn_fp16(in1, in2) result(out)
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2
            type(FP16) :: out
        end function bessel_jn_fp16
    end interface
    interface bessel_yn
        module elemental function bessel_yn_fp16(in1, in2) result(out)
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2
            type(FP16) :: out
        end function bessel_yn_fp16
    end interface
    interface hypot
        module elemental function hypot_fp16(in1, in2) result(out)
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2
            type(FP16) :: out
        end function hypot_fp16
    end interface

    interface cosd
        module elemental function cosd_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function cosd_fp16
    end interface

    interface sind
        module elemental function sind_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function sind_fp16
    end interface

    interface tand
        module elemental function tand_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function tand_fp16
    end interface

    interface cotan
        module elemental function cotan_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function cotan_fp16
    end interface

    interface cotand
        module elemental function cotand_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function cotand_fp16
    end interface

    interface acosd
        module elemental function acosd_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function acosd_fp16
    end interface

    interface asind
        module elemental function asind_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function asind_fp16
    end interface

    interface atand
        module elemental function atand_fp16(x) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16) :: ret_val
        end function atand_fp16
    end interface

    interface atan2d
        module elemental function atan2d_fp16(x,y) result(ret_val)
            type(FP16), intent(in) :: x
            type(FP16), intent(in) :: y
            type(FP16) :: ret_val
        end function atan2d_fp16
    end interface

    interface maxval
        module pure function maxval_fp16_1d(array) result(max_value)
            type(fp16), dimension(:), intent(in) :: array
            type(fp16) :: max_value
        end function
        module pure function maxval_fp16_2d(array) result(max_value)
            type(fp16), dimension(:,:), intent(in) :: array
            type(fp16) :: max_value
        end function
        module pure function maxval_fp16_3d(array) result(max_value)
            type(fp16), dimension(:,:,:), intent(in) :: array
            type(fp16) :: max_value
        end function
        module pure function maxval_fp16_4d(array) result(max_value)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            type(fp16) :: max_value
        end function
        module pure function maxval_fp16_1d_dim(array, dim) result(max_value)
            type(fp16), dimension(:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16) :: max_value
        end function
        module pure function maxval_fp16_2d_dim(array, dim) result(max_value)
            type(fp16), dimension(:,:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16), dimension(size(array, merge(2, 1, dim == 1))) :: max_value
        end function
        module pure function maxval_fp16_3d_dim(array, dim) result(max_value)
            type(fp16), dimension(:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16), dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(2, 3, dim == 3))) :: max_value
        end function
        module pure function maxval_fp16_4d_dim(array, dim) result(max_value)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16), dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(3, 2, dim < 3)), &
                & size(array, merge(4, 3, dim == 4))) :: max_value
        end function
    end interface

    interface minval
        module pure function minval_fp16_1d(array) result(min_value)
            type(fp16), dimension(:), intent(in) :: array
            type(fp16) :: min_value
        end function
        module pure function minval_fp16_2d(array) result(min_value)
            type(fp16), dimension(:,:), intent(in) :: array
            type(fp16) :: min_value
        end function
        module pure function minval_fp16_3d(array) result(min_value)
            type(fp16), dimension(:,:,:), intent(in) :: array
            type(fp16) :: min_value
        end function
        module pure function minval_fp16_4d(array) result(min_value)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            type(fp16) :: min_value
        end function
        module pure function minval_fp16_1d_dim(array, dim) result(min_value)
            type(fp16), dimension(:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16) :: min_value
        end function
        module pure function minval_fp16_2d_dim(array, dim) result(min_value)
            type(fp16), dimension(:,:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16), dimension(size(array, merge(2, 1, dim == 1))) :: min_value
        end function
        module pure function minval_fp16_3d_dim(array, dim) result(min_value)
            type(fp16), dimension(:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16), dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(2, 3, dim == 3))) :: min_value
        end function
        module pure function minval_fp16_4d_dim(array, dim) result(min_value)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            type(fp16), dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(3, 2, dim < 3)), &
                & size(array, merge(4, 3, dim == 4))) :: min_value
        end function
    end interface


    interface maxloc
        module pure function maxloc_fp16_1d(array) result(max_loc)
            type(fp16), dimension(:), intent(in) :: array
            integer :: max_loc
        end function
        module pure function maxloc_fp16_2d(array) result(max_loc)
            type(fp16), dimension(:,:), intent(in) :: array
            integer, dimension(2) :: max_loc
        end function
        module pure function maxloc_fp16_3d(array) result(max_loc)
            type(fp16), dimension(:,:,:), intent(in) :: array
            integer, dimension(3) :: max_loc
        end function
        module pure function maxloc_fp16_4d(array) result(max_loc)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            integer, dimension(4) :: max_loc
        end function
        module pure function maxloc_fp16_1d_dim(array, dim) result(max_loc)
            type(fp16), dimension(:), intent(in) :: array
            integer, intent(in) :: dim
            integer :: max_loc
        end function
        module pure function maxloc_fp16_2d_dim(array, dim) result(max_loc)
            type(fp16), dimension(:,:), intent(in) :: array
            integer, intent(in) :: dim
            integer, dimension(size(array, merge(2, 1, dim == 1))) :: max_loc
        end function
        module pure function maxloc_fp16_3d_dim(array, dim) result(max_loc)
            type(fp16), dimension(:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            integer, dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(2, 3, dim == 3))) :: max_loc
        end function
        module pure function maxloc_fp16_4d_dim(array, dim) result(max_loc)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            integer, dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(3, 2, dim < 3)), &
                & size(array, merge(4, 3, dim == 4))) :: max_loc
        end function
    end interface


    interface minloc
        module pure function minloc_fp16_1d(array) result(min_loc)
            type(fp16), dimension(:), intent(in) :: array
            integer :: min_loc
        end function
        module pure function minloc_fp16_2d(array) result(min_loc)
            type(fp16), dimension(:,:), intent(in) :: array
            integer, dimension(2) :: min_loc
        end function
        module pure function minloc_fp16_3d(array) result(min_loc)
            type(fp16), dimension(:,:,:), intent(in) :: array
            integer, dimension(3) :: min_loc
        end function
        module pure function minloc_fp16_4d(array) result(min_loc)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            integer, dimension(4) :: min_loc
        end function
        module pure function minloc_fp16_1d_dim(array, dim) result(min_loc)
            type(fp16), dimension(:), intent(in) :: array
            integer, intent(in) :: dim
            integer :: min_loc
        end function
        module pure function minloc_fp16_2d_dim(array, dim) result(min_loc)
            type(fp16), dimension(:,:), intent(in) :: array
            integer, intent(in) :: dim
            integer, dimension(size(array, merge(2, 1, dim == 1))) :: min_loc
        end function
        module pure function minloc_fp16_3d_dim(array, dim) result(min_loc)
            type(fp16), dimension(:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            integer, dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(2, 3, dim == 3))) :: min_loc
        end function
        module pure function minloc_fp16_4d_dim(array, dim) result(min_loc)
            type(fp16), dimension(:,:,:,:), intent(in) :: array
            integer, intent(in) :: dim
            integer, dimension( size(array, merge(2, 1, dim == 1)), &
                & size(array, merge(3, 2, dim < 3)), &
                & size(array, merge(4, 3, dim == 4))) :: min_loc
        end function
    end interface



    ! C Interfaces
    INTERFACE
        PURE SUBROUTINE SET_FP16_FROM_INT(OUT, IN) BIND(C, name = "__fp16_helper_set_from_int")
            USE, INTRINSIC:: iso_c_binding
            INTEGER(c_int16_t), INTENT(OUT) :: OUT
            INTEGER(c_int), VALUE, INTENT(IN) :: IN
        END SUBROUTINE
        PURE SUBROUTINE SET_FP16_FROM_REAL(OUT, IN) BIND(C, name = "__fp16_helper_set_from_float")
            USE, INTRINSIC:: iso_c_binding
            INTEGER(c_int16_t), INTENT(OUT) :: OUT
            REAL(c_float), VALUE, INTENT(IN) :: IN
        END SUBROUTINE
        PURE SUBROUTINE SET_FP16_FROM_DOUBLE(OUT, IN) BIND(C, name = "__fp16_helper_set_from_double")
            USE, INTRINSIC:: iso_c_binding
            INTEGER(c_int16_t), INTENT(OUT) :: OUT
            REAL(c_double), VALUE, INTENT(IN) :: IN
        END SUBROUTINE
        PURE REAL(c_float) FUNCTION GET_FP16(IN) BIND(C, name = "__fp16_helper_get_float")
            USE, INTRINSIC :: iso_c_binding
            INTEGER(c_int16_t), INTENT(IN), VALUE :: IN
        END FUNCTION

        pure subroutine helper_add_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_add_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine
        pure subroutine helper_add_fp16_real(out, a, b) bind(c, name = "__fp16_helper_add_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine

        pure subroutine helper_mul_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_mul_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine

        pure subroutine helper_mul_fp16_real(out, a, b) bind(c, name = "__fp16_helper_mul_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine

        pure subroutine helper_sub_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_sub_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine

        pure subroutine helper_sub_fp16_real(out, a, b) bind(c, name = "__fp16_helper_sub_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine

        pure subroutine helper_sub_real_fp16(out, a, b) bind(c, name = "__fp16_helper_sub_real_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            real(c_float), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine

        pure subroutine helper_div_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_div_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine

        pure subroutine helper_div_fp16_real(out, a, b) bind(c, name = "__fp16_helper_div_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine

        pure subroutine helper_div_real_fp16(out, a, b) bind(c, name = "__fp16_helper_div_real_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            real(c_float), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine

        pure subroutine helper_power_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_pow_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine

        pure subroutine helper_power_fp16_real(out, a, b) bind(c, name = "__fp16_helper_pow_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine

        pure subroutine helper_power_fp16_int(out, a, b) bind(c, name = "__fp16_helper_pow_fp16_int")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: a
            integer(c_int), intent(in), value :: b
        end subroutine
    END INTERFACE



CONTAINS
    ! Contructors
    elemental function construct_real(in) result(r)
        real(real32), intent(in) :: in
        type(fp16) :: r
        call SET_FP16_FROM_REAL(r%value, in)
    end function

    elemental function construct_double(in) result(r)
        real(real64), intent(in) :: in
        type(fp16) :: r
        call SET_FP16_FROM_DOUBLE(r%value, in)
    end function

    elemental function construct_int(in) result(r)
        integer, intent(in) :: in
        type(fp16) :: r
        call SET_FP16_FROM_INT(r%value, in)
    end function


    !
    ! Assignment Operator
    !
    elemental SUBROUTINE assign_int(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        INTEGER, INTENT(IN) :: that
        CALL SET_FP16_FROM_INT(this%value, that)
    END SUBROUTINE

    elemental SUBROUTINE assign_real(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        REAL(real32), INTENT(IN) :: that
        CALL SET_FP16_FROM_REAL(this%value, that)
    END SUBROUTINE

    elemental SUBROUTINE assign_double(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        REAL(real64), INTENT(IN) :: that
        CALL SET_FP16_FROM_DOUBLE(this%value, that)
    END SUBROUTINE

    elemental SUBROUTINE assign_fp16(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        TYPE(FP16), INTENT(IN) :: that
        this%value = that%value
    END SUBROUTINE

    elemental SUBROUTINE assign_to_real(this, that)
        REAL(real32), INTENT(OUT) :: this
        TYPE(FP16), INTENT(IN) :: that
        this = GET_FP16(that%value)
    END SUBROUTINE

    !
    ! Operator(+)
    !
    elemental function add_fp16_fp16(this, that) result(sum)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: sum

        call helper_add_fp16_fp16(sum%value, this%value, that%value)
    end function

    elemental function add_fp16_real(this, that) result(sum)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that
        type(fp16) :: sum
        call helper_add_fp16_real(sum%value, this%value, that)
    end function

    elemental function add_real_fp16(this, that) result(sum)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: sum
        call helper_add_fp16_real(sum%value, that%value, this)
    end function

    !
    ! Operator (-)
    !
    elemental function subtract_fp16_fp16(this, that) result(diff)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: diff

        call helper_sub_fp16_fp16(diff%value, this%value, that%value)
    end function

    elemental function subtract_fp16_real(this, that) result(diff)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that
        type(fp16) :: diff

        call helper_sub_fp16_real(diff%value, this%value, that)
    end function

    elemental function subtract_real_fp16(this, that) result(diff)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: diff

        call helper_sub_real_fp16(diff%value, this, that%value)
    end function

    !
    ! Operator(*)
    !
    elemental function multiply_fp16_fp16(this, that) result(prod)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: prod

        call helper_mul_fp16_fp16(prod%value, this%value, that%value)
    end function

    elemental function multiply_fp16_real(this, that) result(prod)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that
        type(fp16) :: prod

        call helper_mul_fp16_real(prod%value, this%value, that)
    end function

    elemental function multiply_real_fp16(this, that) result(prod)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: prod

        call helper_mul_fp16_real(prod%value, that%value, this)
    end function

    !
    ! Operator(/)
    !
    elemental function divide_fp16_fp16(this, that) result(quot)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: quot

        call helper_div_fp16_fp16(quot%value, this%value, that%value)
    end function

    elemental function divide_fp16_real(this, that) result(quot)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that
        type(fp16) :: quot

        call helper_div_fp16_real(quot%value, this%value, that)
    end function

    elemental function divide_real_fp16(this, that) result(quot)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: quot

        call helper_div_real_fp16(quot%value, this, that%value)
    end function

    !
    ! Operator (.lt.)
    !
    elemental function lt_fp16_fp16(x, y) result(out)
        type(fp16), intent(in)  :: x
        type(fp16), intent(in) :: y
        logical :: out

        real(real32) :: xr, yr

        xr = GET_FP16(x%value)
        yr = GET_FP16(y%value)

        out = xr .lt. yr
    end function

    elemental function lt_fp16_fp32(x, y) result(out)
        type(fp16), intent(in)   :: x
        real(real32), intent(in) :: y
        logical :: out

        real(real32) :: xr

        xr = GET_FP16(x%value)

        out = xr .lt. y
    end function

    elemental function lt_fp32_fp16(x, y) result(out)
        type(fp16), intent(in)  :: y
        real(real32), intent(in) :: x
        logical :: out

        real(real32) :: yr

        yr = GET_FP16(y%value)

        out = x .lt. yr
    end function

    !
    ! Operator (.le.)
    !
    elemental function le_fp16_fp16(x, y) result(out)
        type(fp16), intent(in)  :: x
        type(fp16), intent(in) :: y
        logical :: out

        real(real32) :: xr, yr

        xr = GET_FP16(x%value)
        yr = GET_FP16(y%value)

        out = xr .le. yr
    end function

    elemental function le_fp16_fp32(x, y) result(out)
        type(fp16), intent(in)   :: x
        real(real32), intent(in) :: y
        logical :: out

        real(real32) :: xr

        xr = GET_FP16(x%value)

        out = xr .le. y
    end function

    elemental function le_fp32_fp16(x, y) result(out)
        type(fp16), intent(in)  :: y
        real(real32), intent(in) :: x
        logical :: out

        real(real32) :: yr

        yr = GET_FP16(y%value)

        out = x .le. yr
    end function

    !
    ! Operator (.gt.)
    !
    elemental function gt_fp16_fp16(x, y) result(out)
        type(fp16), intent(in)  :: x
        type(fp16), intent(in) :: y
        logical :: out

        real(real32) :: xr, yr

        xr = GET_FP16(x%value)
        yr = GET_FP16(y%value)

        out = xr .gt. yr
    end function

    elemental function gt_fp16_fp32(x, y) result(out)
        type(fp16), intent(in)   :: x
        real(real32), intent(in) :: y
        logical :: out

        real(real32) :: xr

        xr = GET_FP16(x%value)

        out = xr .gt. y
    end function

    elemental function gt_fp32_fp16(x, y) result(out)
        type(fp16), intent(in)  :: y
        real(real32), intent(in) :: x
        logical :: out

        real(real32) :: yr

        yr = GET_FP16(y%value)

        out = x .gt. yr
    end function

    !
    ! Operator (.ge.)
    !
    elemental function ge_fp16_fp16(x, y) result(out)
        type(fp16), intent(in)  :: x
        type(fp16), intent(in) :: y
        logical :: out

        real(real32) :: xr, yr

        xr = GET_FP16(x%value)
        yr = GET_FP16(y%value)

        out = xr .ge. yr
    end function

    elemental function ge_fp16_fp32(x, y) result(out)
        type(fp16), intent(in)   :: x
        real(real32), intent(in) :: y
        logical :: out

        real(real32) :: xr

        xr = GET_FP16(x%value)

        out = xr .ge. y
    end function

    elemental function ge_fp32_fp16(x, y) result(out)
        type(fp16), intent(in)  :: y
        real(real32), intent(in) :: x
        logical :: out

        real(real32) :: yr

        yr = GET_FP16(y%value)

        out = x .ge. yr
    end function

    !
    ! Operator (.eq.)
    !
    elemental function eq_fp16_fp16(x, y) result(out)
        type(fp16), intent(in)  :: x
        type(fp16), intent(in) :: y
        logical :: out

        real(real32) :: xr, yr

        xr = GET_FP16(x%value)
        yr = GET_FP16(y%value)

        out = xr .eq. yr
    end function

    elemental function eq_fp16_fp32(x, y) result(out)
        type(fp16), intent(in)   :: x
        real(real32), intent(in) :: y
        logical :: out

        real(real32) :: xr

        xr = GET_FP16(x%value)

        out = xr .eq. y
    end function

    elemental function eq_fp32_fp16(x, y) result(out)
        type(fp16), intent(in)  :: y
        real(real32), intent(in) :: x
        logical :: out

        real(real32) :: yr

        yr = GET_FP16(y%value)

        out = x .eq. yr
    end function

    !
    ! Operator (.ne.)
    !
    elemental function ne_fp16_fp16(x, y) result(out)
        type(fp16), intent(in)  :: x
        type(fp16), intent(in) :: y
        logical :: out

        real(real32) :: xr, yr

        xr = GET_FP16(x%value)
        yr = GET_FP16(y%value)

        out = xr .ne. yr
    end function

    elemental function ne_fp16_fp32(x, y) result(out)
        type(fp16), intent(in)   :: x
        real(real32), intent(in) :: y
        logical :: out

        real(real32) :: xr

        xr = GET_FP16(x%value)

        out = xr .ne. y
    end function

    elemental function ne_fp32_fp16(x, y) result(out)
        type(fp16), intent(in)  :: y
        real(real32), intent(in) :: x
        logical :: out

        real(real32) :: yr

        yr = GET_FP16(y%value)

        out = x .ne. yr
    end function


    !
    ! Operator(**)
    !
    elemental function power_fp16_fp16(this, that) result(power)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that
        type(fp16) :: power

        call helper_power_fp16_fp16(power%value, this%value, that%value)
    end function

    elemental function power_fp16_real(this, that) result(power)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that
        type(fp16) :: power

        call helper_power_fp16_real(power%value, this%value, that)
    end function

    elemental function power_fp16_int(this, that) result(power)
        type(fp16), intent(in) :: this
        integer, intent(in) :: that
        type(fp16) :: power

        call helper_power_fp16_int(power%value, this%value, that)
    end function

    elemental function min_fp16(x, y) result(out)
        type(fp16), intent(in) :: x, y
        type(fp16) :: out

        if ( x .lt. y) then
            out = x
        else
            out = y
        end if
    end function

    elemental function min3_fp16(x, y, z) result(out)
        type(fp16), intent(in) :: x, y, z
        type(fp16) :: out

        if ( x .lt. y .and. x.lt.z) then
           out = x
        else if ( y.lt.x .and. y .lt. z) then
            out = y
        else
            out = z
        end if
    end function

    elemental function max_fp16(x, y) result(out)
        type(fp16), intent(in) :: x, y
        type(fp16) :: out

        if ( x .gt. y) then
            out = x
        else
            out = y
        end if
    end function

    elemental function max3_fp16(x, y, z) result(out)
        type(fp16), intent(in) :: x, y, z
        type(fp16) :: out

        if ( x .gt. y .and. x .gt. y) then
            out = x
        else if ( y.gt. x .and. y.gt.z) then
            out = y
        else
            out = z
        end if
    end function




    ! Formated output
    subroutine write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        type(FP16), intent(in) :: dtv
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
        character(10) :: pfmt

        REAL(real32) :: x
        x = GET_FP16(dtv%value)

        if (iotype == 'LISTDIRECTED') then
            write(unit, '(F0.4)', iostat=iostat, iomsg=iomsg) x
        else if (iotype == 'DT') then
            if (size(v_list) == 0 ) then
                write(unit, '(F0.4)', iostat=iostat, iomsg=iomsg) x
            else if (size(v_list) == 1 ) then
                write(pfmt, '(a,i2,a)') '(F0.', v_list(1),')'
                write(unit, pfmt, iostat = iostat, iomsg = iomsg) x
            else if (size(v_list) == 2 ) then
                write(pfmt, '(a,i0,a,i0,a)') '(F', v_list(1),'.',v_list(2),')'
                write(unit, pfmt, iostat = iostat, iomsg = iomsg) x
            else
                iostat = -1
                iomsg = 'Too many options in DT setting'
            endif
        else
            iostat = -1
            iomsg = 'Unsupported iotype'
        end if
    end subroutine write_formatted

END MODULE
