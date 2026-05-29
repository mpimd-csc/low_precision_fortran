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

program test_level2_gemv
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_edge()
    call test_options()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical()
        integer(int64) :: m = 3, n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(3, 3) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 2.0
        beta = 1.0
        a = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]

        call gemv('N', m, n, alpha, a, m, x, incx, beta, y, incy)

        do i = 1, 3
            ! Expected: y(i) = 2.0 * sum_j(a(i,j) * x(j)) + 1.0
            block
                real(real64) :: expected
                expected = 1.0_real64
                do j = 1, 3
                    expected = expected + 2.0_real64 * dble(a(i,j)) * dble(x(j))
                end do
                call check_dt_real64("gemv_typical", int(i,int64), y(i), expected, GENERIC_TOL)
            end block
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: m = 0, n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(0, 3) :: a
        type(DT), dimension(3) :: x, y

        alpha = 1.0
        beta = 1.0
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]

        call gemv('N', m, n, alpha, a, int(1, int64), x, incx, beta, y, incy)
        ! No-op, check for crash.
    end subroutine

    subroutine test_options()
        integer(int64) :: m = 3, n = 3
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(3, 3) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 1.0
        beta = 0.0
        a = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        x = [1.0, 2.0, 3.0]
        y = [0.0, 0.0, 0.0]

        call gemv('T', m, n, alpha, a, m, x, incx, beta, y, incy)

        do i = 1, 3
            block
                real(real64) :: expected
                expected = 0.0_real64
                do j = 1, 3
                    expected = expected + 1.0_real64 * dble(a(j,i)) * dble(x(j))
                end do
                call check_dt_real64("gemv_options", int(i,int64), y(i), expected, GENERIC_TOL)
            end block
        end do
        ! I'll rewrite the loop for correctness.
    end subroutine

    subroutine test_stride()
        integer(int64) :: m = 2, n = 2
        integer(int64) :: incx = 2, incy = 2
        type(DT) :: alpha, beta
        type(DT), dimension(2, 2) :: a
        type(DT), dimension(3) :: x, y
        integer :: i, j

        alpha = 1.0
        beta = 1.0
        a = reshape([1.0, 2.0, 3.0, 4.0], [2, 2])
        x = [1.0, 0.0, 2.0] ! x(1)=1, x(2)=2
        y = [1.0, 0.0, 1.0] ! y(1)=1, y(2)=1

        call gemv('N', m, n, alpha, a, m, x, incx, beta, y, incy)

        do i = 1, 2
            block
                real(real64) :: expected
                expected = 1.0_real64
                do j = 1, 2
                    expected = expected + 1.0_real64 * dble(a(i,j)) * dble(x(2*j-1))
                end do
                call check_dt_real64("gemv_stride", int(i,int64), y(2*i-1), expected, GENERIC_TOL)
            end block
        end do
    end subroutine

end program test_level2_gemv
