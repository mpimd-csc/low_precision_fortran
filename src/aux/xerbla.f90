module lpf_xerbla
    use iso_c_binding
    use lpf_types
    implicit none

    interface
        subroutine xerbla_set_function_internal_old(func) bind(C, name = "lpf_blas_xerbla_set_function_fortran_old")
            use, intrinsic :: iso_c_binding
            type(c_funptr), value, intent(in) :: func
        end subroutine
        subroutine xerbla_set_function_internal_new(func) bind(C, name = "lpf_blas_xerbla_set_function_fortran_new")
            use, intrinsic :: iso_c_binding
            type(c_funptr), value, intent(in) :: func
        end subroutine

    end interface

    abstract interface
        subroutine xerbla_function(msg, info)
            use, intrinsic :: iso_c_binding
            import :: lpf_default_int_kind
            character(len=:), pointer, intent(in) :: msg
            integer(lpf_default_int_kind), intent(in) :: info
        end subroutine
    end interface
contains
    subroutine xerbla_set_function_f77(func)
        procedure(xerbla_function) :: func
        type(c_funptr) :: cfun

        cfun = c_funloc(func)
        call xerbla_set_function_internal_old(cfun)

    end subroutine

    subroutine xerbla_set_function_f90(func)
        procedure(xerbla_function) :: func

        type(c_funptr) :: cfun

        cfun = c_funloc(func)
        call xerbla_set_function_internal_new(cfun)

    end subroutine



end module

