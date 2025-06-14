submodule (fp16_support) fp16_math_one_arg
    use iso_c_binding
    use iso_fortran_env
    implicit none

    interface
        pure subroutine helper_fp16_acos(out, in) bind(c, name="__fp16_helper_acos")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_acosh(out, in) bind(c, name="__fp16_helper_acosh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_asin(out, in) bind(c, name="__fp16_helper_asin")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_asinh(out, in) bind(c, name="__fp16_helper_asinh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_atan(out, in) bind(c, name="__fp16_helper_atan")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_atanh(out, in) bind(c, name="__fp16_helper_atanh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_bessel_j0(out, in) bind(c, name="__fp16_helper_bessel_j0")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_bessel_j1(out, in) bind(c, name="__fp16_helper_bessel_j1")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_bessel_y0(out, in) bind(c, name="__fp16_helper_bessel_y0")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_bessel_y1(out, in) bind(c, name="__fp16_helper_bessel_y1")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_ceiling(out, in) bind(c, name="__fp16_helper_ceiling")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_cos(out, in) bind(c, name="__fp16_helper_cos")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_cosh(out, in) bind(c, name="__fp16_helper_cosh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_erf(out, in) bind(c, name="__fp16_helper_erf")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_erfc(out, in) bind(c, name="__fp16_helper_erfc")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_exp(out, in) bind(c, name="__fp16_helper_exp")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_floor(out, in) bind(c, name="__fp16_helper_floor")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_gamma(out, in) bind(c, name="__fp16_helper_gamma")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_log(out, in) bind(c, name="__fp16_helper_log")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_log10(out, in) bind(c, name="__fp16_helper_log10")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_log_gamma(out, in) bind(c, name="__fp16_helper_log_gamma")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_sin(out, in) bind(c, name="__fp16_helper_sin")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_sinh(out, in) bind(c, name="__fp16_helper_sinh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_sqrt(out, in) bind(c, name="__fp16_helper_sqrt")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_tan(out, in) bind(c, name="__fp16_helper_tan")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
        pure subroutine helper_fp16_tanh(out, in) bind(c, name="__fp16_helper_tanh")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
    end interface
contains
    module elemental function acos_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_acos(out%value, in%value)
    end function
    module elemental function acosh_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_acosh(out%value, in%value)
    end function
    module elemental function asin_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_asin(out%value, in%value)
    end function
    module elemental function asinh_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_asinh(out%value, in%value)
    end function
    module elemental function atan_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_atan(out%value, in%value)
    end function
    module elemental function atanh_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_atanh(out%value, in%value)
    end function
    module elemental function bessel_j0_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_bessel_j0(out%value, in%value)
    end function
    module elemental function bessel_j1_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_bessel_j1(out%value, in%value)
    end function
    module elemental function bessel_y0_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_bessel_y0(out%value, in%value)
    end function
    module elemental function bessel_y1_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_bessel_y1(out%value, in%value)
    end function
    module elemental function ceiling_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_ceiling(out%value, in%value)
    end function
    module elemental function cos_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_cos(out%value, in%value)
    end function
    module elemental function cosh_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_cosh(out%value, in%value)
    end function
    module elemental function erf_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_erf(out%value, in%value)
    end function
    module elemental function erfc_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_erfc(out%value, in%value)
    end function
    module elemental function exp_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_exp(out%value, in%value)
    end function
    module elemental function floor_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_floor(out%value, in%value)
    end function
    module elemental function gamma_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_gamma(out%value, in%value)
    end function
    module elemental function log_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_log(out%value, in%value)
    end function
    module elemental function log10_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_log10(out%value, in%value)
    end function
    module elemental function log_gamma_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_log_gamma(out%value, in%value)
    end function
    module elemental function sin_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_sin(out%value, in%value)
    end function
    module elemental function sinh_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_sinh(out%value, in%value)
    end function
    module elemental function sqrt_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_sqrt(out%value, in%value)
    end function
    module elemental function tan_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_tan(out%value, in%value)
    end function
    module elemental function tanh_fp16(in) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in

            call helper_fp16_tanh(out%value, in%value)
    end function
end submodule fp16_math_one_arg
