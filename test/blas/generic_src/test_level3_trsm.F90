#ifdef LPF_FP8_E4M3
#define DT fp8_e4m3
#define BLASMOD lpf_blas_fp8_e4m3
#define TYPEMOD lpf_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
#define DT fp8_e5m2
#define BLASMOD lpf_blas_fp8_e5m2
#define TYPEMOD lpf_fp8_e5m2
#endif
#ifdef LPF_FP16
#define DT fp16
#define BLASMOD lpf_blas_fp16
#define TYPEMOD lpf_fp16
#endif
#ifdef LPF_BF16
#define DT bf16
#define BLASMOD lpf_blas_bf16
#define TYPEMOD lpf_bf16
#endif

program test_level3_trsm
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_combinations()
    call test_summary()

contains

    subroutine test_combinations()
        integer(int64) :: m, n
        integer(int64) :: lda, ldb
        type(DT) :: alpha
        type(DT), allocatable :: a(:,:), b(:,:)
        real(real32), allocatable :: a32(:,:), b32(:,:)
        real(real32) :: alpha32 

        character(len=1) :: side, uplo, transa, diag
        character(len=1) :: s_opts(2), u_opts(2), t_opts(2), d_opts(2)
        integer(int64) :: i, l_idx, t_idx, j_idx, k_idx, m_idx

        s_opts = ['L', 'R']
        u_opts = ['U', 'L']
        t_opts = ['N', 'T']
        d_opts = ['U', 'N']
        alpha = 1.0

        do l_idx = 1, 2
            side = s_opts(l_idx)
            do j_idx = 1, 2
                uplo = u_opts(j_idx)
                do t_idx = 1, 2
                    transa = t_opts(t_idx)
                    do m_idx = 1, 2
                        diag = d_opts(m_idx)

                        if (side == 'L') then
                            m = 4
                            n = 3
                            lda = m
                            ldb = m
                        else ! 'R'
                            m = 3
                            n = 4
                            lda = n
                            ldb = n
                        end if

                        allocate(a(lda, lda))
                        allocate(b(ldb, n))
                        allocate(a32(lda, lda))
                        allocate(b32(ldb, n))

                        do i = 1, lda
                            do k_idx = 1, lda
                                a(i,k_idx) = 1 
                            end do
                        end do
                        do i = 1, ldb
                            do k_idx = 1, n
                                b(i,k_idx) = 1
                            end do
                        end do

                        a32 = a 
                        b32 = b
                        alpha32 = alpha 


                        call trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
                        call strsm(side, uplo, transa, diag, m, n, alpha32, a32, lda, b32, ldb)

                        do i = 1, ldb
                            do k_idx = 1, n 
                                call check_dt_real64("trsm_combinations_" // side // uplo // transa // diag, & 
                                    int(i*100 + k_idx, int64), b(i,k_idx), real(b32(i,k_idx), real64) , GENERIC_TOL)
                            end do 
                        end do 
 
                        deallocate(a, b)
                        deallocate(a32, b32)
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        implicit none
        character, intent(in) :: side, uplo, transa, diag
        integer(int64), intent(in) :: m, n, lda, ldb
        real(real32), intent(in) :: alpha
        real(real32), intent(in) :: a(lda,*)
        real(real32), intent(inout) :: b(ldb,*)

        real(real32) :: temp, zero, one
        integer(int64) :: i, j, k, nrowa
        logical :: lside, nounit, upper

        zero = 0.0
        one = 1.0
        lside = (side == 'L' .or. side == 'l')
        nounit = (diag == 'N' .or. diag == 'n')
        upper = (uplo == 'U' .or. uplo == 'u')

        if (lside) then
            nrowa = m
        else
            nrowa = n
        end if

        if ((m == 0) .or. (n == 0)) return

        if (alpha == zero) then
            do j = 1, n
                do i = 1, m
                    b(i,j) = zero
                end do
            end do
            return
        end if

        if (lside) then
            if (transa == 'N' .or. transa == 'n') then
                if (upper) then
                    do j = 1, n
                        if (alpha /= one) then
                            do i = 1, m
                                b(i,j) = alpha*b(i,j)
                            end do
                        end if
                        do k = m, 1, -1
                            if (b(k,j) /= zero) then
                                if (nounit) b(k,j) = b(k,j)/a(k,k)
                                do i = 1, k - 1
                                    b(i,j) = b(i,j) - b(k,j)*a(i,k)
                                end do
                            end if
                        end do
                    end do
                else
                    do j = 1, n
                        if (alpha /= one) then
                            do i = 1, m
                                b(i,j) = alpha*b(i,j)
                            end do
                        end if
                        do k = 1, m
                            if (b(k,j) /= zero) then
                                if (nounit) b(k,j) = b(k,j)/a(k,k)
                                do i = k + 1, m
                                    b(i,j) = b(i,j) - b(k,j)*a(i,k)
                                end do
                            end if
                        end do
                    end do
                end if
            else
                if (upper) then
                    do j = 1, n
                        do i = 1, m
                            temp = alpha*b(i,j)
                            do k = 1, i - 1
                                temp = temp - a(k,i)*b(k,j)
                            end do
                            if (nounit) temp = temp/a(i,i)
                            b(i,j) = temp
                        end do
                    end do
                else
                    do j = 1, n
                        do i = m, 1, -1
                            temp = alpha*b(i,j)
                            do k = i + 1, m
                                temp = temp - a(k,i)*b(k,j)
                            end do
                            if (nounit) temp = temp/a(i,i)
                            b(i,j) = temp
                        end do
                    end do
                end if
            end if
        else
            if (transa == 'N' .or. transa == 'n') then
                if (upper) then
                    do j = 1, n
                        if (alpha /= one) then
                            do i = 1, m
                                b(i,j) = alpha*b(i,j)
                            end do
                        end if
                        do k = 1, j - 1
                            if (a(k,j) /= zero) then
                                do i = 1, m
                                    b(i,j) = b(i,j) - a(k,j)*b(i,k)
                                end do
                            end if
                        end do
                        if (nounit) then
                            temp = one/a(j,j)
                            do i = 1, m
                                b(i,j) = temp*b(i,j)
                            end do
                        end if
                    end do
                else
                    do j = n, 1, -1
                        if (alpha /= one) then
                            do i = 1, m
                                b(i,j) = alpha*b(i,j)
                            end do
                        end if
                        do k = j + 1, n
                            if (a(k,j) /= zero) then
                                do i = 1, m
                                    b(i,j) = b(i,j) - a(k,j)*b(i,k)
                                end do
                            end if
                        end do
                        if (nounit) then
                            temp = one/a(j,j)
                            do i = 1, m
                                b(i,j) = temp*b(i,j)
                            end do
                        end if
                    end do
                end if
            else
                if (upper) then
                    do k = n, 1, -1
                        if (nounit) then
                            temp = one/a(k,k)
                            do i = 1, m
                                b(i,k) = temp*b(i,k)
                            end do
                        end if
                        do j = 1, k - 1
                            if (a(j,k) /= zero) then
                                temp = a(j,k)
                                do i = 1, m
                                    b(i,j) = b(i,j) - temp*b(i,k)
                                end do
                            end if
                        end do
                        if (alpha /= one) then
                            do i = 1, m
                                b(i,k) = alpha*b(i,k)
                            end do
                        end if
                    end do
                else
                    do k = 1, n
                        if (nounit) then
                            temp = one/a(k,k)
                            do i = 1, m
                                b(i,k) = temp*b(i,k)
                            end do
                        end if
                        do j = k + 1, n
                            if (a(j,k) /= zero) then
                                temp = a(j,k)
                                do i = 1, m
                                    b(i,j) = b(i,j) - temp*b(i,k)
                                end do
                            end if
                        end do
                        if (alpha /= one) then
                            do i = 1, m
                                b(i,k) = alpha*b(i,k)
                            end do
                        end if
                    end do
                end if
            end if
        end if
    end subroutine

end program test_level3_trsm
