! SPDX-License-Identifier: LGPL-3.0-or-later
!
! \brief Symmetric Matrix-Matrix Multiplication (SYMM)
!
! This routine performs the operation:
! C := alpha * A * B + beta * C  (if side == 'L')
! C := alpha * B * A + beta * C  (if side == 'R')
! where A is a symmetric matrix.
!
! \param[in] side Character specifying the side on which the symmetric matrix A is multiplied: 'L' for left, 'R' for right.
! \param[in] uplo Character specifying the triangular part of the symmetric matrix A to be used: 'U' for upper, 'L' for lower.
! \param[in] m Number of rows of matrix B.
! \param[in] n Number of columns of matrix B.
! \param[in] alpha Scalar multiplier for the matrix product.
! \param[in] a The symmetric matrix A, stored in column-major order.
! \param[in] lda Leading dimension of matrix A.
! \param[in] b The matrix B, stored in column-major order.
! \param[in] ldb Leading dimension of matrix B.
! \param[in] beta Scalar multiplier for matrix C.
! \param[in,out] c The matrix C, stored in column-major order.
! \param[in] ldc Leading dimension of matrix C.

#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_symm_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_symm_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

    contains

        pure subroutine symm_impl(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            implicit none
            character, intent(in) :: side, uplo
            integer(int64), intent(in) :: lda, ldb, ldc, m, n
            type(DT), intent(in) :: alpha, beta
            type(DT), intent(in) :: a(lda,*)
            type(DT), intent(in) :: b(ldb,*)
            type(DT), intent(inout) :: c(ldc,*)

            type(DT) :: temp1, temp2, zero
            integer(int64) :: i, j, k, nrowa
            logical :: upper, left

            zero = 0.0
            upper = (uplo == 'U' .or. uplo == 'u')
            left = ( side == 'L' .or. side == 'l')

            if (side == 'L' .or. side == 'l') then
                nrowa = m
            else
                nrowa = n
            end if

            if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == 1.0))) return
            if (alpha.eq.zero) then
                if (beta.eq.zero) then
                    do j = 1,n
                        do  i = 1,m
                            c(i,j) = zero
                        end do
                    end do
                else
                    do j = 1,n
                        do i = 1,m
                            c(i,j) = beta*c(i,j)
                        end do
                    end do
                end if
                return
            end if
            !
            !     start the operations.
            !
            if (left) then
                !
                !        form  c := alpha*a*b + beta*c.
                !
                if (upper) then
                    do j = 1,n
                        do i = 1,m
                            temp1 = alpha*b(i,j)
                            temp2 = zero
                            do k = 1,i - 1
                                c(k,j) = c(k,j) + temp1*a(k,i)
                                temp2 = temp2 + b(k,j)*a(k,i)
                            end do
                            if (beta.eq.zero) then
                                c(i,j) = temp1*a(i,i) + alpha*temp2
                            else
                                c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                            end if
                        end do
                    end do
                else
                    do j = 1,n
                        do i = m,1,-1
                            temp1 = alpha*b(i,j)
                            temp2 = zero
                            do k = i + 1,m
                                c(k,j) = c(k,j) + temp1*a(k,i)
                                temp2 = temp2 + b(k,j)*a(k,i)
                            end do
                            if (beta.eq.zero) then
                                c(i,j) = temp1*a(i,i) + alpha*temp2
                            else
                                c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                            end if
                        end do
                    end do
                end if
            else
                !
                !        form  c := alpha*b*a + beta*c.
                !
                do  j = 1,n
                    temp1 = alpha*a(j,j)
                    if (beta.eq.zero) then
                        do i = 1,m
                            c(i,j) = temp1*b(i,j)
                        end do
                    else
                        do i = 1,m
                            c(i,j) = beta*c(i,j) + temp1*b(i,j)
                        end do
                    end if
                    do  k = 1,j - 1
                        if (upper) then
                            temp1 = alpha*a(k,j)
                        else
                            temp1 = alpha*a(j,k)
                        end if
                        do i = 1,m
                            c(i,j) = c(i,j) + temp1*b(i,k)
                        end do
                    end do
                    do k = j + 1,n
                        if (upper) then
                            temp1 = alpha*a(j,k)
                        else
                            temp1 = alpha*a(k,j)
                        end if
                        do i = 1,m
                            c(i,j) = c(i,j) + temp1*b(i,k)
                        end do
                    end do
                end do
            end if
            !
            return
            !

        end subroutine

        module subroutine symm_64(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            implicit none
            character, intent(in) :: side, uplo
            integer(int64), intent(in) :: lda, ldb, ldc, m, n
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

            call symm_impl(side, uplo, m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc)
        end subroutine

        module subroutine symm_32(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            implicit none
            character, intent(in) :: side, uplo
            integer(int32), intent(in) :: lda, ldb, ldc, m, n
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

            call symm_impl(side, uplo, int(m, int64), int(n, int64), &
                & alpha, pa, int(lda, int64), pb, int(ldb, int64), beta, pc, int(ldc, int64))
        end subroutine

    end submodule
