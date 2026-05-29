#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_gemm_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_gemm_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine gemm_impl(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        implicit none
        character, intent(in) :: transa, transb
        integer(int64), intent(in) :: k, lda, ldb, ldc, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(in) :: b(ldb,*)
        type(DT), intent(inout) :: c(ldc,*)

        type(DT) :: temp, zero
        integer(int64) :: i, j, l, nrowa, nrowb
        logical :: nota, notb

        zero = 0.0

        nota = (transa == 'N' .or. transa == 'n')
        notb = (transb == 'N' .or. transb == 'n')

        if (nota) then
            nrowa = m
        else
            nrowa = k
        end if
        if (notb) then
            nrowb = k
        else
            nrowb = n
        end if

        if ((m == 0) .or. (n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == 1.0))) return

        if (alpha == zero) then
            if (beta == zero) then
                do j = 1, n
                    do i = 1, m
                        c(i,j) = zero
                    end do
                end do
            else
                do j = 1, n
                    do i = 1, m
                        c(i,j) = beta*c(i,j)
                    end do
                end do
            end if
            return
        end if

        if (notb) then
            if (nota) then
                do j = 1, n
                    if (beta == zero) then
                        do i = 1, m
                            c(i,j) = zero
                        end do
                    else if (beta /= 1.0) then
                        do i = 1, m
                            c(i,j) = beta*c(i,j)
                        end do
                    end if
                    do l = 1, k
                        temp = alpha*b(l,j)
                        do i = 1, m
                            c(i,j) = c(i,j) + temp*a(i,l)
                        end do
                    end do
                end do
            else
                do j = 1, n
                    do i = 1, m
                        temp = zero
                        do l = 1, k
                            temp = temp + a(l,i)*b(l,j)
                        end do
                        if (beta == zero) then
                            c(i,j) = alpha*temp
                        else
                            c(i,j) = alpha*temp + beta*c(i,j)
                        end if
                    end do
                end do
            end if
        else
            if (nota) then
                do j = 1, n
                    if (beta == zero) then
                        do i = 1, m
                            c(i,j) = zero
                        end do
                    else if (beta /= 1.0) then
                        do i = 1, m
                            c(i,j) = beta*c(i,j)
                        end do
                    end if
                    do l = 1, k
                        temp = alpha*b(j,l)
                        do i = 1, m
                            c(i,j) = c(i,j) + temp*a(i,l)
                        end do
                    end do
                end do
            else
                do j = 1, n
                    do i = 1, m
                        temp = zero
                        do l = 1, k
                            temp = temp + a(l,i)*b(j,l)
                        end do
                        if (beta == zero) then
                            c(i,j) = alpha*temp
                        else
                            c(i,j) = alpha*temp + beta*c(i,j)
                        end if
                    end do
                end do
            end if
        end if
    end subroutine

    module subroutine gemm_64(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        implicit none
        character, intent(in) :: transa, transb
        integer(int64), intent(in) :: k, lda, ldb, ldc, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: b(..)
        type(DT), target, intent(inout) :: c(..)

        type(DT), CONTIGUOUS, pointer :: pa(:)
        type(DT), CONTIGUOUS, pointer :: pb(:)
        type(DT), CONTIGUOUS, pointer :: pc(:)
        integer(int64) :: la(1), lb(1), lc(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        total_size = size(b)
        lb(1) = total_size
        ptr = c_loc(b)
        call c_f_pointer(ptr, pb, lb)

        total_size = size(c)
        lc(1) = total_size
        ptr = c_loc(c)
        call c_f_pointer(ptr, pc, lc)

        call gemm_impl(transa, transb, m, n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc)
    end subroutine

    module subroutine gemm_32(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        implicit none
        character, intent(in) :: transa, transb
        integer(int32), intent(in) :: k, lda, ldb, ldc, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: b(..)
        type(DT), target, intent(inout) :: c(..)

        type(DT), CONTIGUOUS, pointer :: pa(:)
        type(DT), CONTIGUOUS, pointer :: pb(:)
        type(DT), CONTIGUOUS, pointer :: pc(:)
        integer(int64) :: la(1), lb(1), lc(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        total_size = size(b)
        lb(1) = total_size
        ptr = c_loc(b)
        call c_f_pointer(ptr, pb, lb)

        total_size = size(c)
        lc(1) = total_size
        ptr = c_loc(c)
        call c_f_pointer(ptr, pc, lc)

        call gemm_impl(transa, transb, int(m, int64), int(n, int64), &
            & int(k, int64), alpha, pa, int(lda, int64), pb, &
            & int(ldb, int64), beta, pc, int(ldc, int64))
    end subroutine

end submodule
