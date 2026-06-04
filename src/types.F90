module lpf_types
    use iso_fortran_env, only: int32, int64
    use iso_c_binding, only: c_int32_t, c_int64_t, c_double
    implicit none

    integer, parameter :: lpf_default_int_kind = int32
    integer, parameter :: lpf_default_c_int_kind = c_int32_t

    integer, parameter :: lpf_int32_kind = int32
    integer, parameter :: lpf_int64_kind = int64

    interface
        real(c_double) function lpf_get_wtime() bind(C, name = "lpf_get_wtime")
            import :: c_double
        end function
    end interface

end module

