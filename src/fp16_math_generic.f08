submodule (fp16_support) fp16_math_generic
    use iso_c_binding
    use iso_fortran_env
    implicit none

    ! C Interfaces
    INTERFACE
        pure subroutine helper_abs_fp16(out, in) bind(c, name = "__fp16_helper_abs")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine

        pure subroutine helper_epsilon_fp16(out) bind(c, name = "__fp16_helper_epsilon")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
        end subroutine

        pure subroutine helper_exponent_fp16(out, in) bind(c, name = "__fp16_helper_exponent")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int), intent(out) :: out
        end subroutine

        pure subroutine helper_fraction_fp16(out, in) bind(c, name = "__fp16_helper_fraction")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(in), value :: in
            integer(c_int16_t), intent(out) :: out
        end subroutine



    END INTERFACE

contains
    module elemental function abs_fp16(x) result(abs_val)
        type(FP16), intent(in) :: x
        type(FP16) :: abs_val
        call helper_abs_fp16(abs_val%value, x%value)
    end function abs_fp16

    module elemental function epsilon_fp16(x) result(out)
        type(FP16), intent(in) :: x
        type(FP16) :: out
        call helper_epsilon_fp16(out%value)
    end function epsilon_fp16

    module elemental function exponent_fp16 (x) result(out)
        type(FP16), intent(in) :: x
        integer :: out

        call helper_exponent_fp16(out, x%value)
    end function exponent_fp16

    module elemental function fraction_fp16 (x) result(out)
        type(FP16), intent(in) :: x
        type(FP16) :: out

        call helper_fraction_fp16(out%value, x%value)
    end function fraction_fp16


    module elemental function radix_fp16(x) result(out)
        type(FP16), intent(in) :: x
        integer :: out
        out = 2
    end function radix_fp16



end submodule fp16_math_generic
