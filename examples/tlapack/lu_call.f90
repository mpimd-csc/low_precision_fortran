module solve_lu
    use lpf_fp16
    use lpf_bf16
    use LPF_FP8_E4M3
    use LPF_FP8_E5M2
    use lpf_blas_fp16
    use lpf_blas_bf16
    use lpf_blas_fp8_e4m3
    use lpf_blas_fp8_e5m2
    use iso_c_binding
    use iso_fortran_env

    implicit none
    interface
        subroutine solve_fp16(n, A, lda, x, incx) bind(c, name = "solve_fp16")
            import c_int64_t, fp16
            integer(c_int64_t), intent(in), value :: n, lda, incx
            type(fp16), intent(in) :: A(lda, *)
            type(fp16), intent(inout) :: x(*)
        end subroutine
        subroutine solve_bf16(n, A, lda, x, incx) bind(c, name = "solve_bf16")
            import c_int64_t, bf16
            integer(c_int64_t), intent(in), value :: n, lda, incx
            type(bf16), intent(in) :: A(lda, *)
            type(bf16), intent(inout) :: x(*)
        end subroutine
        subroutine solve_fp8_e4m3(n, A, lda, x, incx) bind(c, name = "solve_fp8_e4m3")
            import c_int64_t, FP8_E4M3
            integer(c_int64_t), intent(in), value :: n, lda, incx
            type(fp8_e4m3), intent(in) :: A(lda, *)
            type(fp8_e4m3), intent(inout) :: x(*)
        end subroutine
        subroutine solve_fp8_e5m2(n, A, lda, x, incx) bind(c, name = "solve_fp8_e5m2")
            import c_int64_t, FP8_E5M2
            integer(c_int64_t), intent(in), value :: n, lda, incx
            type(fp8_e5m2), intent(in) :: A(lda, *)
            type(fp8_e5m2), intent(inout) :: x(*)
        end subroutine



    end interface

contains

    subroutine lu_fp16(n, A)
        use lpf_fp16
        use lpf_blas_fp16
        use iso_fortran_env
        implicit none

        integer(int64), intent(in) :: n
        real(real32), dimension(:,:), intent(in) :: A
        type(fp16) :: A_lp(n,n)
        type(fp16) :: x(n), y(n), y_old(n)
        type(fp16) :: nrmA, nrmx, nrmr

        x = 1.0
        A_lp = A

        call gemv("N", n, n, FP16(1.0), A_lp, n, x, 1_int64, FP16(0.0), y, 1_int64)
        nrmx = nrm2_fp32(n, y, 1_int64)

        y_old = y
        call solve_fp16(n, A_lp, n, y, 1_int64)

        write(*,*) " Y = ", y

        call gemv("N", n, n, FP16(1.0), A_lp, n, y, 1_int64, FP16(-1.0), y_old, 1_int64);

        nrmA = nrm2_fp32(n*n, A_lp, 1_int64)
        nrmr = nrm2_fp32(n, y_old, 1_int64)

        write(*,*) "FP16: ||A * x - b || / ( ||A||_F + || b ||) = ", nrmr / (nrmA + nrmx)

    end subroutine

    subroutine lu_bf16(n, A)
        use lpf_bf16
        use lpf_blas_bf16
        use iso_fortran_env
        implicit none

        integer(int64), intent(in) :: n
        real(real32), dimension(:,:), intent(in) :: A
        type(bf16) :: A_lp(n,n)
        type(bf16) :: x(n), y(n), y_old(n)
        type(bf16) :: nrmA, nrmx, nrmr

        x = 1.0
        A_lp = A

        call gemv("N", n, n, BF16(1.0), A_lp, n, x, 1_int64, BF16(0.0), y, 1_int64)
        nrmx = nrm2_fp32(n, y, 1_int64)

        y_old = y
        call solve_bf16(n, A_lp, n, y, 1_int64)

        write(*,*) " Y = ", y
        call gemv("N", n, n, BF16(1.0), A_lp, n, y, 1_int64, BF16(-1.0), y_old, 1_int64);

        nrmA = nrm2_fp32(n*n, A_lp, 1_int64)
        nrmr = nrm2_fp32(n, y_old, 1_int64)

        write(*,*) "BF16: ||A * x - b || / ( ||A||_F + || b ||) = ", nrmr / (nrmA + nrmx)

    end subroutine

    subroutine lu_fp8_e4m3(n, A)
        use lpf_fp8_e4m3
        use lpf_blas_fp8_e4m3
        use iso_fortran_env
        implicit none

        integer(int64), intent(in) :: n
        real(real32), dimension(:,:), intent(in) :: A
        type(fp8_e4m3) :: A_lp(n,n)
        type(fp8_e4m3) :: x(n), y(n), y_old(n)
        type(fp8_e4m3) :: nrmA, nrmx, nrmr

        x = 1.0
        A_lp = A

        call gemv("N", n, n, FP8_E4M3(1.0), A_lp, n, x, 1_int64, FP8_E4M3(0.0), y, 1_int64)
        nrmx = nrm2_fp32(n, y, 1_int64)

        y_old = y
        call solve_fp8_e4m3(n, A_lp, n, y, 1_int64)

        write(*,*) " Y = ", y
        call gemv("N", n, n, FP8_E4M3(1.0), A_lp, n, y, 1_int64, FP8_E4M3(-1.0), y_old, 1_int64);

        nrmA = nrm2_fp32(n*n, A_lp, 1_int64)
        nrmr = nrm2_fp32(n, y_old, 1_int64)

        write(*,*) "FP8_E4M3: ||A * x - b || / ( ||A||_F + || b ||) = ", nrmr / (nrmA + nrmx)

    end subroutine

    subroutine lu_fp8_e5m2(n, A)
        use lpf_fp8_e5m2
        use lpf_blas_fp8_e5m2
        use iso_fortran_env
        implicit none

        integer(int64), intent(in) :: n
        real(real32), dimension(:,:), intent(in) :: A
        type(fp8_e5m2) :: A_lp(n,n)
        type(fp8_e5m2) :: x(n), y(n), y_old(n)
        type(fp8_e5m2) :: nrmA, nrmx, nrmr

        x = 1.0
        A_lp = A

        call gemv("N", n, n, FP8_E5M2(1.0), A_lp, n, x, 1_int64, FP8_E5M2(0.0), y, 1_int64)
        nrmx = nrm2_fp32(n, y, 1_int64)

        y_old = y
        call solve_fp8_e5m2(n, A_lp, n, y, 1_int64)

        write(*,*) " Y = ", y
        call gemv("N", n, n, FP8_E5M2(1.0), A_lp, n, y, 1_int64, FP8_E5M2(-1.0), y_old, 1_int64);

        nrmA = nrm2_fp32(n*n, A_lp, 1_int64)
        nrmr = nrm2_fp32(n, y_old, 1_int64)

        write(*,*) "FP8_E5M2: ||A * x - b || / ( ||A||_F + || b ||) = ", nrmr / (nrmA + nrmx)

    end subroutine



end module
program lu_call
    use solve_lu
    use lpf_fp16
    use lpf_blas_fp16
    use iso_c_binding
    use iso_fortran_env
    implicit none

    integer(int64), parameter :: n = 3
    real :: A(n, n)
    integer(int64) :: i

    A = 0.0
    do i = 1, n
        A(i,i) = 2.0
        if (i .lt. n ) then
            A(i,i+1) = -1.0
            A(i+1, i) = -1.0
        end if
    end do


    ! call RANDOM_NUMBER(A);

    call lu_fp16(n, A)
    call lu_bf16(n, A)
    call lu_fp8_e4m3(n, A)
    call lu_fp8_e5m2(n, A)


end program

