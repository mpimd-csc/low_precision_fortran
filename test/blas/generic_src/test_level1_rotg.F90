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

program test_level1_rotg
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()
    call test_zero()

    call test_summary()

contains

    subroutine test_typical()
        type(DT) :: sa, sb, sc, ss

        sa = 3.0
        sb = 4.0

        call rotg(sa, sb, sc, ss)

        ! Expected: sa is modified to sqrt(sa^2+sb^2) = 5.0, and sb becomes z = 1/0.6
        ! sc = 3/5 = 0.6, ss = 4/5 = 0.8

        call check_dt_real64("rotg_typical_sa", 0_int64, sa, 5.0_real64, GENERIC_TOL)
        call check_dt_real64("rotg_typical_sb", 0_int64, sb, 1.6666_real64, GENERIC_TOL)
        call check_dt_real64("rotg_typical_sc", 0_int64, sc, 0.6_real64, GENERIC_TOL)
        call check_dt_real64("rotg_typical_ss", 0_int64, ss, 0.8_real64, GENERIC_TOL)
    end subroutine

    subroutine test_zero()
        type(DT) :: sa, sb, sc, ss

        sa = 0.0
        sb = 0.0

        call rotg(sa, sb, sc, ss)

        call check_dt_real64("rotg_zero_sa", 0_int64, sa, 0.0_real64, GENERIC_TOL)
        call check_dt_real64("rotg_zero_sb", 0_int64, sb, 0.0_real64, GENERIC_TOL)
        call check_dt_real64("rotg_zero_sc", 0_int64, sc, 1.0_real64, GENERIC_TOL)
        call check_dt_real64("rotg_zero_ss", 0_int64, ss, 0.0_real64, GENERIC_TOL)
    end subroutine

end program test_level1_rotg
