!  SPDX-License-Identifier: LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program; if not, write to the Free Software Foundation,
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!



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

program test_level1_rotm
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
        type(DT), dimension(2) :: sx, sy
        type(DT) :: sa, sb, sc, sd
        type(DT), DIMENSION(5) :: sparam

        sx = [1.0, 2.0]
        sy = [3.0, 4.0]

        sa = 1.0
        sb = 1.0
        sc = 1.0
        sd = 3.0
        call rotmg(sa, sb, sc, sd, sparam)
        call rotm(n, sx, incx, sy, incy, sparam)

        call check_dt_real64("rotm_typical_sx1", 0_int64, sx(1), 3.33_real64, GENERIC_TOL)
        call check_dt_real64("rotm_typical_sx2", 0_int64, sx(2), 4.65_real64, GENERIC_TOL)
        call check_dt_real64("rotm_typical_sy1", 0_int64, sy(1), 0.0_real64, GENERIC_TOL)
        call check_dt_real64("rotm_typical_sy2", 0_int64, sy(2), -0.6666_real64, GENERIC_TOL)
    end subroutine

    subroutine test_n0()
        integer(int64) :: n = 0
        integer(int64) :: incx = 1, incy = 1
        type(DT), dimension(0) :: sx, sy, sparam

        call rotm(n, sx, incx, sy, incy, sparam)
    end subroutine

end program test_level1_rotm
