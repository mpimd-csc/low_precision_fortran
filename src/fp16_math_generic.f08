submodule (fp16_support) fp16_math_generic
    use iso_c_binding
    use iso_fortran_env
    implicit none

    ! C Interfaces
    INTERFACE
        subroutine helper_abs_fp16(out, in) bind(c, name = "__fp16_helper_abs")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in
        end subroutine
    END INTERFACE

contains
    module function abs_fp16(x) result(abs_val)
        type(FP16), intent(in) :: x
        type(FP16) :: abs_val
        call helper_abs_fp16(abs_val%value, x%value)
    end function abs_fp16


end submodule fp16_math_generic
