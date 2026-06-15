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

program test_level3_gemmtr
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_combinations()
    call test_summary()

contains

    subroutine test_combinations()
        integer(int64) :: n, k
        integer(int64) :: lda, ldb, ldc
        type(DT) :: alpha, beta
        type(DT), allocatable :: a(:,:), b(:,:), c(:,:)
        character(len=1) :: uplo, transa, transb
        character(len=1) :: u_opts(2), t_opts(3)
        integer(int64) :: i, j, l_idx, j_idx, t_idx, k_idx, js, je

        u_opts = ['U', 'L']
        t_opts = ['N', 'T', 'C']
        alpha = 1.0
        beta = 0.0

        do l_idx = 1, 2
            uplo = u_opts(l_idx)
            do j_idx = 1, 3
                transa = t_opts(j_idx)
                do t_idx = 1, 3
                    transb = t_opts(t_idx)

                    n = 4
                    k = 3
                    lda = n
                    ldb = n
                    ldc = n

                    allocate(a(lda, n))
                    allocate(b(ldb, n))
                    allocate(c(ldc, n))

                    do i = 1, size(a, 1)
                        do k_idx = 1, size(a, 2)
                            a(i,k_idx) = dble(i + k_idx)
                        end do
                    end do
                    do i = 1, size(b, 1)
                        do k_idx = 1, size(b, 2)
                            b(i,k_idx) = dble(i * k_idx)
                        end do
                    end do
                    c = 0.0

                    call gemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

                    ! Reference implementation
                    do i = 1, n
                        if (uplo.eq.'U') then
                            js = i
                            je = n
                        else
                            js = 1
                            je = i
                        end if

                        do j = js, je
                            block
                                real(real64) :: ref
                                ref = 0.0_real64
                                do k_idx = 1, k
                                    block
                                        real(real64) :: vala, valb
                                        if (transa == 'N') then
                                            vala = dble(a(i, k_idx))
                                        else if (transa == 'T') then
                                            vala = dble(a(k_idx, i))
                                        else ! 'C'
                                            vala = dble(a(k_idx, i))
                                        end if
                                        if (transb == 'N') then
                                            valb = dble(b(k_idx, j))
                                        else if (transb == 'T') then
                                            valb = dble(b(j, k_idx))
                                        else ! 'C'
                                            valb = dble(b(j, k_idx))
                                        end if
                                        ref = ref + dble(alpha) * vala * valb
                                    end block
                                end do
                                call check_dt_real64("gemmtr_combinations" // uplo // transa // transb, &
                                    & int(i*100 + j, int64), c(i,j), ref, GENERIC_TOL)
                            end block
                        end do
                    end do

                    deallocate(a, b, c)
                end do
            end do
        end do
    end subroutine

end program test_level3_gemmtr
