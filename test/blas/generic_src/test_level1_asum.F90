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

program test_level1_asum
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_n0()
    call test_n1()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical()
        integer(int64) :: n = 10
        integer(int64) :: incx = 1
        type(DT), dimension(10) :: x
        type(DT) :: result
        real(real64) :: expected = 0.0_real64
        integer :: i

        x = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0, -10.0]
        do i = 1, 10
            expected = expected + abs(dble(x(i)))
        end do

        result = asum(n, x, incx)
        call check_dt_real64("asum_typical", 0_int64, result, expected, GENERIC_TOL)
    end subroutine

    subroutine test_n0()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1
        type(DT), dimension(0) :: x
        type(DT) :: result

        result = asum(n, x, incx)
        call check_dt_real64("asum_n0", 0_int64, result, 0.0_real64, GENERIC_TOL)
    end subroutine

    subroutine test_n1()
        integer(int64) :: n = 1
        integer(int64) :: incx = 1
        type(DT), dimension(1) :: x
        type(DT) :: result

        x(1) = -5.0
        result = asum(n, x, incx)
        call check_dt_real64("asum_n1", 0_int64 , result, 5.0_real64, GENERIC_TOL)
    end subroutine

    subroutine test_stride()
        integer(int64) :: n = 5
        integer(int64) :: incx = 2
        type(DT), dimension(9) :: x
        type(DT) :: result
        real(real64) :: expected = 0.0_real64
        integer :: i

        x = [1.0, 0.0, -2.0, 0.0, 3.0, 0.0, -4.0, 0.0, 5.0]
        do i = 1, 5
            expected = expected + abs(dble(x(2*i-1)))
        end do

        result = asum(n, x, incx)
        call check_dt_real64("asum_stride", 0_int64,  result, expected, GENERIC_TOL)
    end subroutine

end program test_level1_asum
