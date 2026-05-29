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

program test_level1_copy
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_n0()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical()
        integer(int64) :: n = 10
        integer(int64) :: incx = 1, incy = 1
        type(DT), dimension(10) :: sx, sy
        integer :: i

        sx = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        sy = 0.0

        call copy(n, sx, incx, sy, incy)

        do i = 1, 10
            call check_dt_real64("copy_typical", int(i, int64), sy(i), dble(sx(i)), GENERIC_TOL)
        end do
    end subroutine

    subroutine test_n0()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1, incy = 1
        type(DT), dimension(0) :: sx, sy

        call copy(n, sx, incx, sy, incy)
    end subroutine

    subroutine test_stride()
        integer(int64) :: n = 5
        integer(int64) :: incx = 2, incy = 2
        type(DT), dimension(9) :: sx, sy
        integer :: i

        sx = [1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0]
        sy = 0.0

        call copy(n, sx, incx, sy, incy)

        do i = 1, 5
            call check_dt_real64("copy_stride", int(i, int64), sy(2*i-1), dble(sx(2*i-1)), GENERIC_TOL)
        end do
    end subroutine

end program test_level1_copy
