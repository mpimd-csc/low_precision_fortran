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

program test_level1_rot
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_n0()

    call test_summary()

contains

    subroutine test_typical()
        integer(int64) :: n = 2
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: c, s
        type(DT), dimension(2) :: sx, sy

        c = 0.866025_real64 ! cos(30 deg)
        s = 0.5_real64      ! sin(30 deg)
        sx = [1.0, 2.0]
        sy = [3.0, 4.0]

        call rot(n, sx, incx, sy, incy, c, s)

       ! G = [ cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]
       ! x = [ 1, 2; 3, 4]
       ! G * x
       !  ans =
       !   2.3660   3.7321
       !   2.0981   2.4641

        call check_dt_real64("rot_typical_sx1", 0_int64, sx(1), 2.3660_real64, GENERIC_TOL)
        call check_dt_real64("rot_typical_sx2", 0_int64, sx(2), 3.7321_real64, GENERIC_TOL)
        call check_dt_real64("rot_typical_sy1", 0_int64, sy(1), 2.0981_real64, GENERIC_TOL)
        call check_dt_real64("rot_typical_sy2", 0_int64, sy(2), 2.4641_real64, GENERIC_TOL)
    end subroutine

    subroutine test_n0()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: c, s
        type(DT), dimension(0) :: sx, sy

        c = 1.0
        s = 0.0
        call rot(n, sx, incx, sy, incy, c, s)
    end subroutine

end program test_level1_rot
