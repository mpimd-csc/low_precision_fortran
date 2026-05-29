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

program test_level1_rotmg
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical()

    call test_summary()

contains

    subroutine test_typical()
        type(DT) :: sa, sb, sc, sd
        type(DT), dimension(5) :: sparam
        integer :: i
        real(real64), dimension(5) :: sparam_ref

        sa = 1.0
        sb = 1.0
        sc = 1.0
        sd = 3.0

        sparam_ref(1) =  1.0000
        sparam_ref(2) = .334
        sparam_ref(3) = .0000
        sparam_ref(4) = .0000
        sparam_ref(5) = .334

        call rotmg(sa, sb, sc, sd, sparam)

        do i = 1, 5
            call check_dt_real64("rotmg_typical_sa", int(i, int64), sparam(i), sparam_ref(i), GENERIC_TOL)
        end do
    end subroutine

end program test_level1_rotmg
