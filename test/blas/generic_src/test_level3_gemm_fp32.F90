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

program test_level3_gemm_fp32
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_combinations()
    call test_summary()

contains

    subroutine test_combinations()
        integer(int64) :: m, n, k
        integer(int64) :: lda, ldb, ldc
        real(real32) :: alpha
        real(real32) :: beta
        type(DT), allocatable :: a(:,:), b(:,:)
        real(real32), allocatable :: c(:,:)
        character(len=1) :: transa, transb
        character(len=1) :: t_opts(3)
        integer(int64) :: i, j, l_idx, j_idx, k_idx

        t_opts = ['N', 'T', 'C']
        alpha = 1.0
        beta = 0.0

        do l_idx = 1, 3
            transa = t_opts(l_idx)
            do j_idx = 1, 3
                transb = t_opts(j_idx)

                ! Set dimensions based on trans
                if (transa == 'N') then
                    m = 4
                    k = 3
                    lda = m
                else if (transa == 'T') then
                    m = 3
                    k = 4
                    lda = k
                else ! 'C'
                    m = 3
                    k = 4
                    lda = k
                end if

                if (transb == 'N') then
                    n = 4
                    ldb = k
                else if (transb == 'T') then
                    n = 3
                    ldb = k
                else ! 'C'
                    n = 3
                    ldb = k
                end if

                ldc = m

                allocate(a(lda, k))
                allocate(b(ldb, max(n, k)))
                allocate(c(ldc, n))

                do i = 1, size(a, 1)
                    do k_idx = 1, size(a, 2)
                        a(i,k_idx) = dble(i + k_idx)
                    end do
                end do
                do i = 1, size(b, 1)
                    do k_idx = 1, size(b, 2)
                        b(i,k_idx) = dble(i * k_idx)
                    end do
                end do
                c = 0.0

                call gemm_fp32(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

                ! Reference implementation
                do i = 1, m
                    do j = 1, n
                        block
                            real(real64) :: ref
                            ref = 0.0_real64
                            do k_idx = 1, k
                                block
                                    real(real64) :: vala, valb
                                    if (transa == 'N') then
                                        vala = dble(a(i, k_idx))
                                    else if (transa == 'T') then
                                        vala = dble(a(k_idx, i))
                                    else ! 'C'
                                        vala = dble(a(k_idx, i))
                                    end if
                                    if (transb == 'N') then
                                        valb = dble(b(k_idx, j))
                                    else if (transb == 'T') then
                                        valb = dble(b(j, k_idx))
                                    else ! 'C'
                                        valb = dble(b(j, k_idx))
                                    end if
                                    ref = ref + dble(alpha) * vala * valb
                                end block
                            end do
                            call check_real32_real64("gemm_fp32_combinations", int(i*100 + j, int64), c(i,j), ref, GENERIC_TOL)
                        end block
                    end do
                end do

                deallocate(a, b, c)
                end do
            end do
        end subroutine

end program test_level3_gemm_fp32
