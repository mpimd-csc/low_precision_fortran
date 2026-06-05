!  SPDX-License-Identifier LGPL-3.0-or-later
!
!  Unit tests for FP8_E4M3 array reduction functions.
!
PROGRAM test_fp8_e4m3_reductions
    USE iso_fortran_env, only: real32, real64, int32
    USE LPF_FP8_E4M3
    USE lpf_fp8_e4m3_test_utils
    IMPLICIT NONE

    CALL test_maxval_1d()
    CALL test_minval_1d()
    CALL test_maxval_2d()
    CALL test_minval_2d()
    CALL test_maxval_2d_dim()
    CALL test_minval_2d_dim()
    CALL test_maxloc_1d()
    CALL test_minloc_1d()
    CALL test_maxval_1d_single()
    CALL test_minval_1d_allneg()

    CALL test_norm2_1d()
    CALL test_norm2_2d()
    CALL test_norm2_3d()
    CALL test_norm2_single()
    CALL test_norm2_zero()
    CALL test_norm2_empty()
    CALL test_norm2_small()
    CALL test_norm2_dim()

    CALL test_summary()

CONTAINS

    subroutine test_maxval_1d()
        type(FP8_E4M3), dimension(5) :: arr
        arr = [FP8_E4M3(1.0_real32), FP8_E4M3(5.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(-2.0_real32), FP8_E4M3(4.0_real32)]
        call check_fp8_e4m3_real64('maxval_1d', maxval(arr), 5.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_minval_1d()
        type(FP8_E4M3), dimension(5) :: arr
        arr = [FP8_E4M3(1.0_real32), FP8_E4M3(5.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(-2.0_real32), FP8_E4M3(4.0_real32)]
        call check_fp8_e4m3_real64('minval_1d', minval(arr), -2.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_maxval_2d()
        type(FP8_E4M3), dimension(2,3) :: arr
        arr = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), &
                       FP8_E4M3(3.0_real32), FP8_E4M3(10.0_real32), &
                       FP8_E4M3(5.0_real32), FP8_E4M3(-1.0_real32)], [2,3])
        call check_fp8_e4m3_real64('maxval_2d', maxval(arr), 10.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_minval_2d()
        type(FP8_E4M3), dimension(2,3) :: arr
        arr = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), &
                       FP8_E4M3(3.0_real32), FP8_E4M3(10.0_real32), &
                       FP8_E4M3(5.0_real32), FP8_E4M3(-1.0_real32)], [2,3])
        call check_fp8_e4m3_real64('minval_2d', minval(arr), -1.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_maxval_2d_dim()
        type(FP8_E4M3), dimension(2,3) :: arr
        type(FP8_E4M3), dimension(3) :: row_max
        type(FP8_E4M3), dimension(2) :: col_max
        arr = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), &
                       FP8_E4M3(3.0_real32), FP8_E4M3(10.0_real32), &
                       FP8_E4M3(5.0_real32), FP8_E4M3(-1.0_real32)], [2,3])
        row_max = maxval(arr, dim=1)
        call check_fp8_e4m3_real64('maxval_2d_dim1_c1', row_max(1), 2.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('maxval_2d_dim1_c2', row_max(2), 10.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('maxval_2d_dim1_c3', row_max(3), 5.0_real64, FP8_E4M3_TOL_TIGHT)

        col_max = maxval(arr, dim=2)
        call check_fp8_e4m3_real64('maxval_2d_dim2_r1', col_max(1), 5.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('maxval_2d_dim2_r2', col_max(2), 10.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_minval_2d_dim()
        type(FP8_E4M3), dimension(2,3) :: arr
        type(FP8_E4M3), dimension(3) :: row_min
        arr = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(-2.0_real32), &
                       FP8_E4M3(3.0_real32), FP8_E4M3(0.5_real32), &
                       FP8_E4M3(5.0_real32), FP8_E4M3(-1.0_real32)], [2,3])
        row_min = minval(arr, dim=1)
        call check_fp8_e4m3_real64('minval_2d_dim1_c1', row_min(1), -2.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('minval_2d_dim1_c2', row_min(2), 0.5_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('minval_2d_dim1_c3', row_min(3), -1.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_maxloc_1d()
        type(FP8_E4M3), dimension(5) :: arr
        arr = [FP8_E4M3(1.0_real32), FP8_E4M3(5.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(-2.0_real32), FP8_E4M3(4.0_real32)]
        call check_integer('maxloc_1d', maxloc(arr), 2)
    end subroutine

    subroutine test_minloc_1d()
        type(FP8_E4M3), dimension(5) :: arr
        arr = [FP8_E4M3(1.0_real32), FP8_E4M3(5.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(-2.0_real32), FP8_E4M3(4.0_real32)]
        call check_integer('minloc_1d', minloc(arr), 4)
    end subroutine

    subroutine test_maxval_1d_single()
        type(FP8_E4M3), dimension(1) :: arr
        arr = [FP8_E4M3(7.0_real32)]
        call check_fp8_e4m3_real64('maxval_1d_single', maxval(arr), 7.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_minval_1d_allneg()
        type(FP8_E4M3), dimension(3) :: arr
        arr = [FP8_E4M3(-1.0_real32), FP8_E4M3(-5.0_real32), FP8_E4M3(-3.0_real32)]
        call check_fp8_e4m3_real64('minval_1d_allneg', minval(arr), -5.0_real64, FP8_E4M3_TOL_TIGHT)
        call check_fp8_e4m3_real64('maxval_1d_allneg', maxval(arr), -1.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_1d()
        type(FP8_E4M3), dimension(4) :: arr
        arr = [FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32), FP8_E4M3(-12.0_real32), FP8_E4M3(1.0_real32)]
        call check_fp8_e4m3_real64('norm2_1d', norm2(arr), sqrt(170.0_real64), FP8_E4M3_TOL)
    end subroutine

    subroutine test_norm2_2d()
        type(FP8_E4M3), dimension(2,2) :: arr
        arr = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32)], [2,2])
        call check_fp8_e4m3_real64('norm2_2d', norm2(arr), sqrt(30.0_real64), FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_3d()
        type(FP8_E4M3), dimension(2,2,2) :: arr
        arr = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32), &
                       FP8_E4M3(5.0_real32), FP8_E4M3(6.0_real32), FP8_E4M3(7.0_real32), FP8_E4M3(8.0_real32)], [2,2,2])
        call check_fp8_e4m3_real64('norm2_3d', norm2(arr), sqrt(204.0_real64), FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_single()
        type(FP8_E4M3), dimension(1) :: arr
        arr = [FP8_E4M3(5.5_real32)]
        call check_fp8_e4m3_real64('norm2_single', norm2(arr), 5.5_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_zero()
        type(FP8_E4M3), dimension(5) :: arr
        arr = 0.0
        call check_fp8_e4m3_real64('norm2_zero', norm2(arr), 0.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_empty()
        type(FP8_E4M3), dimension(0) :: arr
        call check_fp8_e4m3_real64('norm2_empty', norm2(arr), 0.0_real64, FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_small()
        type(FP8_E4M3), dimension(2) :: arr
        type(FP8_E4M3) :: t
        t = tiny(FP8_E4M3(1.0))
        arr = [t, t]
        call check_fp8_e4m3_real64('norm2_small', norm2(arr), dble(t) * sqrt(2.0_real64), FP8_E4M3_TOL_TIGHT)
    end subroutine

    subroutine test_norm2_dim()
        type(FP8_E4M3), dimension(2) :: arr1
        type(FP8_E4M3), dimension(2,2) :: arr2
        type(FP8_E4M3), dimension(2,2,2) :: arr3
        type(FP8_E4M3), dimension(2) :: out2_d1, out2_d2
        type(FP8_E4M3), dimension(2,2) :: out3_d1, out3_d2, out3_d3

        arr1 = [FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32)]
        call check_fp8_e4m3_real64('norm2_dim_1d', norm2(arr1, dim=1), 5.0_real64, FP8_E4M3_TOL)

        arr2 = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32)], [2,2])
        out2_d1 = norm2(arr2, dim=1)
        call check_fp8_e4m3_real64('norm2_dim_2d_d1_c1', out2_d1(1), sqrt(5.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_2d_d1_c2', out2_d1(2), 5.0_real64, FP8_E4M3_TOL)

        out2_d2 = norm2(arr2, dim=2)
        call check_fp8_e4m3_real64('norm2_dim_2d_d2_r1', out2_d2(1), sqrt(10.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_2d_d2_r2', out2_d2(2), sqrt(20.0_real64), FP8_E4M3_TOL)

        arr3 = reshape([FP8_E4M3(1.0_real32), FP8_E4M3(2.0_real32), FP8_E4M3(3.0_real32), FP8_E4M3(4.0_real32), &
                       FP8_E4M3(5.0_real32), FP8_E4M3(6.0_real32), FP8_E4M3(7.0_real32), FP8_E4M3(8.0_real32)], [2,2,2])
        out3_d1 = norm2(arr3, dim=1)
        call check_fp8_e4m3_real64('norm2_dim_3d_d1_11', out3_d1(1,1), sqrt(5.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d1_21', out3_d1(2,1), 5.0_real64, FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d1_12', out3_d1(1,2), sqrt(61.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d1_22', out3_d1(2,2), sqrt(113.0_real64), FP8_E4M3_TOL)

        out3_d2 = norm2(arr3, dim=2)
        call check_fp8_e4m3_real64('norm2_dim_3d_d2_11', out3_d2(1,1), sqrt(10.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d2_21', out3_d2(2,1), sqrt(20.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d2_12', out3_d2(1,2), sqrt(74.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d2_22', out3_d2(2,2), 10.0_real64, FP8_E4M3_TOL)

        out3_d3 = norm2(arr3, dim=3)
        call check_fp8_e4m3_real64('norm2_dim_3d_d3_11', out3_d3(1,1), sqrt(26.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d3_21', out3_d3(2,1), sqrt(40.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d3_12', out3_d3(1,2), sqrt(58.0_real64), FP8_E4M3_TOL)
        call check_fp8_e4m3_real64('norm2_dim_3d_d3_22', out3_d3(2,2), sqrt(80.0_real64), FP8_E4M3_TOL)
    end subroutine

END PROGRAM test_fp8_e4m3_reductions
