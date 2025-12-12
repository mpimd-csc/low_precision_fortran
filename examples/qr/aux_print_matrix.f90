module printmatrix
    use lpf_bf16
    use lpf_fp16
    implicit none

    interface print_matrix
        module procedure :: print_fp16, print_bf16, print_fp32
    end interface
contains
    subroutine print_fp32(A)
        use iso_fortran_env
        implicit none
        real(real32), DIMENSION(:,:) :: A
        integer :: k,l

        do k = 1, size(A, 1)
            do l = 1, size(A, 2)
                write(*, '(F16.9, " ")', advance = "NO" ) A(k,l)
            end do
            write(*,*) ""
        end do
    end subroutine

    subroutine print_fp16(A)
        use iso_fortran_env
        implicit none
        type(fp16), DIMENSION(:,:) :: A
        integer :: k,l

        do k = 1, size(A, 1)
            do l = 1, size(A, 2)
                write(*, '(DT(16,9), " ")', advance = "NO" ) A(k,l)
            end do
            write(*,*) ""
        end do
    end subroutine

    subroutine print_bf16(A)
        use iso_fortran_env
        implicit none
        type(bf16), DIMENSION(:,:) :: A
        integer :: k,l

        do k = 1, size(A, 1)
            do l = 1, size(A, 2)
                write(*, '(DT(16,9), " ")', advance = "NO" ) A(k,l)
            end do
            write(*,*) ""
        end do
    end subroutine
end module
