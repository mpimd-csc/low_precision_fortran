PROGRAM test_bf16_simple
    USE iso_fortran_env, only: real32, real64
    USE LPF_BF16
    IMPLICIT NONE
    type(BF16), allocatable :: a(:)
    integer :: n
    n = 4
    allocate(a(n))
    a = BF16(1.0_real64)
    write(*,*) "Success"
    deallocate(a)
END PROGRAM test_bf16_simple
