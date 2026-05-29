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

program test_level2_symv
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical('U')
    call test_typical('L')
    call test_edge()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical(uplo)
        character(len=1) :: uplo
        integer(int64) :: n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(3, 3) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 2.0
        beta = 1.0
        a = reshape([1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0], [3, 3])
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]

        call symv(uplo, n, alpha, a, n, x, incx, beta, y, incy)

        do i = 1, 3
            block
                real(real64) :: expected
                expected = 1.0_real64
                do j = 1, 3
                    expected = expected + 2.0_real64 * dble(a(j,i)) * dble(x(j))
                end do
                call check_dt_real64("symv_typical_" // uplo, int(i,int64), y(i), expected, GENERIC_TOL)
            end block
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(0, 0) :: a
        type(DT), dimension(0) :: x, y

        alpha = 1.0
        beta = 1.0
        call symv('U', n, alpha, a, int(1, int64), x, incx, beta, y, incy)
    end subroutine

    subroutine test_stride()
        integer(int64) :: n = 2
        integer(int64) :: incx = 2, incy = 2
        type(DT) :: alpha, beta
        type(DT), dimension(2, 2) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 1.0
        beta = 1.0
        a = reshape([1.0, 2.0, 2.0, 3.0], [2, 2])
        x = [1.0, 0.0, 2.0]
        y = [1.0, 0.0, 1.0]

        call symv('U', n, alpha, a, n, x, incx, beta, y, incy)

        do i = 1, 2
            block
                real(real64) :: expected
                expected = 1.0_real64
                do j = 1, 2
                    expected = expected + 1.0_real64 * dble(a(j,i)) * dble(x(2*j-1))
                end do
                call check_dt_real64("symv_stride", int(i,int64), y(2*i-1), expected, GENERIC_TOL)
            end block
        end do
    end subroutine

end program test_level2_symv
