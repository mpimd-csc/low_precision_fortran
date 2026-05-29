module cg
    use lpf_bf16
    implicit none

contains
    subroutine cg_bf16(n, A, lda, x, maxit, tol)
        use lpf_bf16
        use lpf_blas_bf16
        implicit none

        integer, intent(in) :: n
        type(bf16), intent(in) :: A(lda, n)
        integer, intent(in) :: lda
        type(bf16), intent(inout) :: x(n)
        integer, intent(in) :: maxit
        type(bf16), intent(in) :: tol

        type(bf16), allocatable :: r(:), p(:), Ap(:), b(:)
        type(bf16) :: rho, rho_new, alpha, beta
        integer :: k
        type(bf16) :: one, zero

        one = bf16(1.0)
        zero = bf16(0.0)

        allocate(b(n), r(n), p(n), Ap(n))

        b = x
        x = zero

        r = b
        call gemv('N', n, n, -one, A, lda, x, 1, one, r, 1)
        p = r

        rho = dot(n, r, 1, r, 1)

        do k = 1, maxit
            ! Ap = A * p
            call gemv('N', n, n, one, A, lda, p, 1, zero, Ap, 1)

            ! alpha = rho / (p' * Ap)
            alpha = rho / dot(n, p, 1, Ap, 1)

            ! x = x + alpha * p
            x = x + alpha * p

            ! r = r - alpha * Ap

            call axpy(n, -alpha, Ap, 1, r, 1)

            ! rho_new = r' * r
            rho_new = dot(n, r, 1, r, 1)

            write(*,'(A, I5, A, DT"E"(10,5),A, DT"E"(10,5))') "IT = ", k, " RHO_NEW = ", rho_new, " TOL = ", tol
            if (sqrt(rho_new) .lt. tol) exit

            ! beta = rho_new / rho
            beta = rho_new / rho

            p = r + beta * p

            rho = rho_new
        end do

        deallocate(b, r, p, Ap)
    end subroutine

    subroutine cg_fp8_e4m3(n, A, lda, x, maxit, tol)
        use lpf_fp8_e4m3
        use lpf_blas_fp8_e4m3
        implicit none

        integer, intent(in) :: n
        type(fp8_e4m3), intent(in) :: A(lda, n)
        integer, intent(in) :: lda
        type(fp8_e4m3), intent(inout) :: x(n)
        integer, intent(in) :: maxit
        type(fp8_e4m3), intent(in) :: tol

        type(fp8_e4m3), allocatable :: r(:), p(:), Ap(:), b(:)
        type(fp8_e4m3) :: rho, rho_new, alpha, beta
        integer :: k
        type(fp8_e4m3) :: one, zero

        one = fp8_e4m3(1.0)
        zero = fp8_e4m3(0.0)

        allocate(b(n), r(n), p(n), Ap(n))

        b = x
        x = zero

        r = b
        call gemv('N', n, n, -one, A, lda, x, 1, one, r, 1)
        p = r

        rho = dot(n, r, 1, r, 1)

        do k = 1, maxit
            ! Ap = A * p
            call gemv('N', n, n, one, A, lda, p, 1, zero, Ap, 1)

            ! alpha = rho / (p' * Ap)
            alpha = rho / dot(n, p, 1, Ap, 1)

            ! x = x + alpha * p
            x = x + alpha * p

            ! r = r - alpha * Ap

            call axpy(n, -alpha, Ap, 1, r, 1)

            ! rho_new = r' * r
            rho_new = dot(n, r, 1, r, 1)

            write(*,'(A, I5, A, DT"E"(10,5),A, DT"E"(10,5))') "IT = ", k, " RHO_NEW = ", rho_new, " TOL = ", tol
            if (sqrt(rho_new) .lt. tol) exit

            ! beta = rho_new / rho
            beta = rho_new / rho

            p = r + beta * p

            rho = rho_new
        end do

        deallocate(b, r, p, Ap)
    end subroutine

    subroutine cg_fp8_e5m2(n, A, lda, x, maxit, tol)
        use lpf_fp8_e5m2
        use lpf_blas_fp8_e5m2
        implicit none

        integer, intent(in) :: n
        type(fp8_e5m2), intent(in) :: A(lda, n)
        integer, intent(in) :: lda
        type(fp8_e5m2), intent(inout) :: x(n)
        integer, intent(in) :: maxit
        type(fp8_e5m2), intent(in) :: tol

        type(fp8_e5m2), allocatable :: r(:), p(:), Ap(:), b(:)
        type(fp8_e5m2) :: rho, rho_new, alpha, beta
        integer :: k
        type(fp8_e5m2) :: one, zero

        one = fp8_e5m2(1.0)
        zero = fp8_e5m2(0.0)

        allocate(b(n), r(n), p(n), Ap(n))

        b = x
        x = zero

        r = b
        call gemv('N', n, n, -one, A, lda, x, 1, one, r, 1)
        p = r

        rho = dot(n, r, 1, r, 1)

        do k = 1, maxit
            ! Ap = A * p
            call gemv('N', n, n, one, A, lda, p, 1, zero, Ap, 1)

            ! alpha = rho / (p' * Ap)
            alpha = rho / dot(n, p, 1, Ap, 1)

            ! x = x + alpha * p
            x = x + alpha * p

            ! r = r - alpha * Ap

            call axpy(n, -alpha, Ap, 1, r, 1)

            ! rho_new = r' * r
            rho_new = dot(n, r, 1, r, 1)

            write(*,'(A, I5, A, DT"E"(10,5),A, DT"E"(10,5))') "IT = ", k, " RHO_NEW = ", rho_new, " TOL = ", tol
            if (sqrt(rho_new) .lt. tol) exit

            ! beta = rho_new / rho
            beta = rho_new / rho

            p = r + beta * p

            rho = rho_new
        end do

        deallocate(b, r, p, Ap)
    end subroutine

    subroutine cg_fp16(n, A, lda, x, maxit, tol)
        use lpf_fp16
        use lpf_blas_fp16
        implicit none

        integer, intent(in) :: n
        type(fp16), intent(in) :: A(lda, n)
        integer, intent(in) :: lda
        type(fp16), intent(inout) :: x(n)
        integer, intent(in) :: maxit
        type(fp16), intent(in) :: tol

        type(fp16), allocatable :: r(:), p(:), Ap(:), b(:)
        type(fp16) :: rho, rho_new, alpha, beta
        integer :: k
        type(fp16) :: one, zero

        one = fp16(1.0)
        zero = fp16(0.0)

        allocate(b(n), r(n), p(n), Ap(n))

        b = x
        x = zero

        r = b
        call gemv('N', n, n, -one, A, lda, x, 1, one, r, 1)
        p = r

        rho = dot(n, r, 1, r, 1)

        do k = 1, maxit
            ! Ap = A * p
            call gemv('N', n, n, one, A, lda, p, 1, zero, Ap, 1)

            ! alpha = rho / (p' * Ap)
            alpha = rho / dot(n, p, 1, Ap, 1)

            ! x = x + alpha * p
            x = x + alpha * p

            ! r = r - alpha * Ap

            call axpy(n, -alpha, Ap, 1, r, 1)

            ! rho_new = r' * r
            rho_new = dot(n, r, 1, r, 1)

            write(*,'(A, I5, A, DT"E"(10,5),A, DT"E"(10,5))') "IT = ", k, " RHO_NEW = ", rho_new, " TOL = ", tol
            if (sqrt(rho_new) .lt. tol) exit

            ! beta = rho_new / rho
            beta = rho_new / rho

            p = r + beta * p

            rho = rho_new
        end do

        deallocate(b, r, p, Ap)
    end subroutine




end module


program cg_test
    use lpf_bf16
    use lpf_fp16
    use lpf_fp8_e5m2
    use lpf_fp8_e4m3
    use lpf_blas_bf16
    use lpf_blas_fp16
    use lpf_blas_fp8_e5m2
    use lpf_blas_fp8_e4m3
    use iso_fortran_env
    use cg
    implicit none

    integer, parameter :: n = 20
    real(real32) :: A(n, n), x (n), nrmA, nrmb, nrmr
    type(fp16) :: A_fp16(n,n), x_fp16(n), b_fp16(n)
    type(bf16) :: A_bf16(n,n), x_bf16(n), b_bf16(n)
    type(fp8_e4m3) :: A_fp8_e4m3(n,n), x_fp8_e4m3(n), b_fp8_e4m3(n)
    type(fp8_e5m2) :: A_fp8_e5m2(n,n), x_fp8_e5m2(n), b_fp8_e5m2(n)

    integer :: i

    nrmA = 0.0
    A = 0.0
    do i = 1, n
        A(i,i) = 2.0
        nrmA = nrmA + 4.0
        if (i .lt. n ) then
            A(i,i+1) = -1.0
            A(i+1, i) = -1.0
            nrmA = nrmA + 2.0
        end if
    end do
    x = 0.0
    x(1) = 1;
    x(n) = 1;
    nrmb = NORM2(x)
    nrmA = sqrt(nrmA)

    A_fp16 = A
    A_bf16 = A
    A_fp8_e4m3 = A
    A_fp8_e5m2 = A

    x_fp16 = x
    x_bf16 = x
    x_fp8_e4m3 = x
    x_fp8_e5m2 = x
    b_fp16 = x
    b_bf16 = x
    b_fp8_e4m3 = x
    b_fp8_e5m2 = x







    call cg_bf16(n, A_bf16, n, x_bf16, 100, bf16(0.001))
    write(*,*) x_bf16
    call gemv("N", n, n, bf16(-1.0), A_bf16, n, x_bf16, 1, bf16(1.0), b_bf16, 1)
    nrmr = real(nrm2(n, b_bf16, 1))
    write(*,*) "BF16 Res = ", nrmr/(nrmA + nrmb)

    call cg_fp16(n, A_fp16, n, x_fp16, 100, fp16(0.001))
    write(*,*) x_fp16
    call gemv("N", n, n, fp16(-1.0), A_fp16, n, x_fp16, 1, fp16(1.0), b_fp16, 1)
    nrmr = real(nrm2(n, b_fp16, 1))
    write(*,*) "FP16 Res = ", nrmr/(nrmA + nrmb)


    call cg_fp8_e4m3(n, A_fp8_e4m3, n, x_fp8_e4m3, 100, fp8_e4m3(0.001))
    write(*,*) x_fp8_e4m3
    call gemv("N", n, n, fp8_e4m3(-1.0), A_fp8_e4m3, n, x_fp8_e4m3, 1, fp8_e4m3(1.0), b_fp8_e4m3, 1)
    nrmr = real(nrm2(n, b_fp8_e4m3, 1))
    write(*,*) "FP8_E4M3 = ", nrmr/(nrmA + nrmb)


    call cg_fp8_e5m2(n, A_fp8_e5m2, n, x_fp8_e5m2, 100, fp8_e5m2(0.001))
    write(*,*) x_fp8_e5m2
    call gemv("N", n, n, fp8_e4m3(-1.0), A_fp8_e4m3, n, x_fp8_e4m3, 1, fp8_e4m3(1.0), b_fp8_e4m3, 1)
    nrmr = real(nrm2(n, b_fp8_e4m3, 1))
    write(*,*) "FP8_E5M2 Res = ", nrmr/(nrmA + nrmb)





end program
