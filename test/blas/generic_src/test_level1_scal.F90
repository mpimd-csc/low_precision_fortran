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

program test_level1_scal
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
        integer(int64) :: incx = 1
        type(DT) :: sa
        type(DT), dimension(10) :: sx
        real(real64), dimension(10) :: orig
        integer :: i

        sa = 2.0
        sx = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        do i = 1, 10
            orig(i) = dble(sx(i))
        end do

        call scal(n, sa, sx, incx)

        do i = 1, 10
            call check_dt_real64("scal_typical", int(i, int64), sx(i), 2.0_real64 * orig(i), GENERIC_TOL)
        end do
    end subroutine

    subroutine test_n0()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1
        type(DT) :: sa
        type(DT), dimension(0) :: sx

        sa = 1.0
        call scal(n, sa, sx, incx)
    end subroutine

    subroutine test_stride()
        integer(int64) :: n = 5
        integer(int64) :: incx = 2
        type(DT) :: sa
        type(DT), dimension(9) :: sx
        real(real64), dimension(5) :: orig
        integer :: i

        sa = 3.0
        sx = [1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0]
        do i = 1, 5
            orig(i) = dble(sx(2*i-1))
        end do

        call scal(n, sa, sx, incx)

        do i = 1, 5
            call check_dt_real64("scal_stride", int(i, int64), sx(2*i-1), 3.0_real64 * orig(i), GENERIC_TOL)
        end do
    end subroutine

end program test_level1_scal
