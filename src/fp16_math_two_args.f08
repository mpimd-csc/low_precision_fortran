submodule (fp16_support) fp16_math_two_args
    use iso_c_binding
    use iso_fortran_env
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
            integer(c_int16_t), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine
        pure subroutine helper_fp16_bessel_yn(out, in1, in2) bind(c, name="__fp16_helper_bessel_yn")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in1
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
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2

            call helper_fp16_bessel_jn(out%value, in1%value, in2%value)
    end function
    module elemental function bessel_yn_fp16(in1, in2) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2

            call helper_fp16_bessel_yn(out%value, in1%value, in2%value)
    end function
    module elemental function hypot_fp16(in1, in2) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2

            call helper_fp16_hypot(out%value, in1%value, in2%value)
    end function
end submodule fp16_math_two_args
