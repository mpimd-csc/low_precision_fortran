#ifdef LPF_FP8_E4M3
#define DT fp8_e4m3
#define TYPEMOD lpf_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
#define DT fp8_e5m2
#define TYPEMOD lpf_fp8_e5m2
#endif
#ifdef LPF_FP16
#define DT fp16
#define TYPEMOD lpf_fp16
#endif
#ifdef LPF_BF16
#define DT bf16
#define TYPEMOD lpf_bf16
#endif

module generic_test_utils
    use TYPEMOD
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    integer :: n_pass = 0
    integer :: n_fail = 0

    real(real64), parameter :: GENERIC_TOL_FLOOR = 1.0d-4

#ifdef LPF_FP8_E4M3
    real(real64), parameter :: GENERIC_TOL = 1.0d-1
#endif
#ifdef LPF_FP8_E5M2
    real(real64), parameter :: GENERIC_TOL = 2.0d-1
#endif
#ifdef LPF_BF16
    real(real64), parameter :: GENERIC_TOL = 1.0d-2
#endif
#ifdef LPF_FP16
    real(real64), parameter :: GENERIC_TOL = 5.0d-3
#endif


contains

    subroutine check_dt_real64(name, idx, got_dt, expected_r64, tol)
        character(len=*), intent(in) :: name
        type(DT), intent(in) :: got_dt
        integer(int64), intent(in) :: idx
        real(real64), intent(in) :: expected_r64, tol
        real(real64) :: got_r64, diff, rel_err

        got_r64 = dble(got_dt)
        diff = abs(got_r64 - expected_r64)

        ! Relative error check with floor
        if (abs(expected_r64) > GENERIC_TOL_FLOOR) then
            rel_err = diff / abs(expected_r64)
        else
            rel_err = diff
        end if

        if (rel_err <= tol) then
            n_pass = n_pass + 1
        else
            print *, "FAIL: ", name, " (IDX = ", idx, ") Expected: ", expected_r64, &
                & " Got: ", got_r64, " RelErr: ", rel_err, " Tol: ", tol
            n_fail = n_fail + 1
        end if
    end subroutine

    subroutine check_real32_real64(name, idx, got_r32, expected_r64, tol)
        character(len=*), intent(in) :: name
        real(real32), intent(in) :: got_r32
        integer(int64), intent(in) :: idx
        real(real64), intent(in) :: expected_r64, tol
        real(real64) :: got_r64, diff, rel_err

        got_r64 = dble(got_r32)
        diff = abs(got_r64 - expected_r64)

        ! Relative error check with floor
        if (abs(expected_r64) > GENERIC_TOL_FLOOR) then
            rel_err = diff / abs(expected_r64)
        else
            rel_err = diff
        end if

        if (rel_err <= tol) then
            n_pass = n_pass + 1
        else
            print *, "FAIL: ", name, " (IDX = ", idx, ") Expected: ", expected_r64, &
                & " Got: ", got_r64, " RelErr: ", rel_err, " Tol: ", tol
            n_fail = n_fail + 1
        end if
    end subroutine

    subroutine check_integer(name, got, expected)
        character(len=*), intent(in) :: name
        integer(int64), intent(in) :: got, expected

        if (got .eq. expected) then
            n_pass = n_pass + 1
        else
            print *, "FAIL: ", name, " Expected: ", expected, " Got: ", got
            n_fail = n_fail + 1
        end if
    end subroutine

    subroutine check_logical(name, got, expected)
        character(len=*), intent(in) :: name
        logical, intent(in) :: got, expected

        if (got .eqv. expected) then
            n_pass = n_pass + 1
        else
            print *, "FAIL: ", name, " Expected: ", expected, " Got: ", got
            n_fail = n_fail + 1
        end if
    end subroutine

    subroutine test_summary()
        if (n_fail .eq. 0) then
            print *, "ALL TESTS PASSED (", n_pass, " tests)"
            stop 0
        else
            print *, "TESTS FAILED: ", n_fail, " failures out of ", n_pass + n_fail, " tests"
            stop 1
        end if
    end subroutine

end module generic_test_utils
