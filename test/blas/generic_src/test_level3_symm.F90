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

program test_level3_symm
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_combinations()
    call test_summary()

contains

    subroutine test_combinations()
        integer(int64) :: m, n
        integer(int64) :: lda, ldb, ldc
        type(DT) :: alpha, beta
        type(DT), allocatable :: a(:,:), b(:,:), c(:,:)
        character(len=1) :: side, uplo
        character(len=1) :: s_opts(2), u_opts(2)
        integer(int64) :: i, j, l_idx, j_idx, k_idx

        s_opts = ['L', 'R']
        u_opts = ['U', 'L']
        alpha = 1.0
        beta = 0.0

        do l_idx = 1, 2
            side = s_opts(l_idx)
            do j_idx = 1, 2
                uplo = u_opts(j_idx)

                if (side == 'L') then
                    m = 4
                    n = 3
                    lda = m
                    ldb = m
                    ldc = m
                else ! 'R'
                    m = 3
                    n = 4
                    lda = n
                    ldb = n
                    ldc = n
                end if

                allocate(a(lda, lda))
                allocate(b(ldb, n))
                allocate(c(ldc, n))

                ! Initialize symmetric matrix A
                do i = 1, lda
                    do k_idx = 1, lda
                        a(i,k_idx) = dble(i + k_idx)
                    end do
                end do
                do i = 1, ldb
                    do k_idx = 1, n
                        b(i,k_idx) = dble(i * k_idx)
                    end do
                end do
                c = 0.0

                call symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)

                ! Reference implementation
                do i = 1, m
                    do j = 1, n
                        block
                            real(real64) :: ref
                            ref = 0.0_real64
                            if (side == 'L') then
                                ! C = alpha*A*B + beta*C
                                do k_idx = 1, m
                                    ref = ref + dble(alpha) * dble(a(i,k_idx)) * dble(b(k_idx,j))
                                end do
                            else ! 'R'
                                ! C = alpha*B*A + beta*C
                                do k_idx = 1, n
                                    ref = ref + dble(alpha) * dble(b(i,k_idx)) * dble(a(k_idx,j))
                                end do
                            end if
                            call check_dt_real64("symm_combinations", int(i*100 + j, int64), c(i,j), ref, GENERIC_TOL)
                        end block
                    end do
                end do

                deallocate(a, b, c)
                end do
            end do
        end subroutine

end program test_level3_symm
