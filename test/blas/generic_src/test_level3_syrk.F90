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

program test_level3_syrk
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_combinations()
    call test_summary()

contains

    subroutine test_combinations()
        integer(int64) :: n, k
        integer(int64) :: lda, ldc
        type(DT) :: alpha, beta
        type(DT), allocatable :: a(:,:), c(:,:)
        character(len=1) :: uplo, trans
        character(len=1) :: u_opts(2), t_opts(2)
        integer(int64) :: i, j, cols, l_idx, j_idx, k_idx, js, je

        u_opts = ['U', 'L']
        t_opts = ['N', 'T']
        alpha = 1.0
        beta = 0.0

        do l_idx = 1, 2
            uplo = u_opts(l_idx)
            do j_idx = 1, 2
                trans = t_opts(j_idx)

                n = 4
                k = 3
                if (trans == 'N') then
                    lda = n
                    cols = k
                else ! 'T'
                    lda = k
                    cols = n
                end if
                ldc = n

                allocate(a(lda, cols))
                allocate(c(ldc, n))

                do i = 1, size(a, 1)
                    do k_idx = 1, size(a, 2)
                        a(i,k_idx) = dble(i + k_idx)
                    end do
                end do
                c = 0.0

                call syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)

                ! Reference implementation: C = alpha*A*A' + beta*C
                do i = 1, n
                    if (uplo.eq.'U') then
                        js = i
                        je = n
                    else
                        js = 1
                        je = i
                    end if

                    do j = js, je
                        block
                            real(real64) :: ref
                            ref = 0.0_real64
                            do k_idx = 1, k
                                block
                                    real(real64) :: vala, valb
                                    if (trans == 'N') then
                                        vala = dble(a(i, k_idx))
                                        valb = dble(a(j, k_idx))
                                    else ! 'T'
                                        vala = dble(a(k_idx, i))
                                        valb = dble(a(k_idx, j))
                                    end if
                                    ref = ref + dble(alpha) * vala * valb
                                end block
                            end do
                            call check_dt_real64("syrk_combinations" // uplo //trans &
                                & , int(i*100 + j, int64), c(i,j), ref, GENERIC_TOL)
                        end block
                    end do
                end do

                deallocate(a, c)
                end do
            end do
        end subroutine

end program test_level3_syrk
