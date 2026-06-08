!  SPDX-License-Identifier LGPL-3.0-or-later
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

submodule (lpf_bf16) lpf_bf16_matmul
    use iso_fortran_env
    implicit none

    contains
        module pure function bf16_matmul_vm(a, b) result(c)
            type(bf16), intent(in) :: a(:)
            type(bf16), intent(in) :: b(:,:)
            type(bf16) :: c(size(b,2))

            integer(int64) ::  i,l
            integer(int64) ::  m, n, k
            type(bf16) :: temp

            m = size(a, 1)
            n = size(b, 1)
            k = size(b, 2)

            if ( m .ne. n ) error stop 'matmul: dimension of A and B not matching'

            do i = 1, k
                temp = 0.0
                do l = 1, n
                    temp = temp + a(l) * b(l,i)
                end do
                c(i) = temp
            end do

            return

        end function


        module pure function bf16_matmul_mv(a, b) result(c)
            type(bf16), intent(in) :: a(:,:)
            type(bf16), intent(in) :: b(:)
            type(bf16) :: c(size(a, 1))

            integer(int64) ::  i,j
            integer(int64) ::  m, n, k
            type(bf16) :: temp

            m = size(a, 1)
            n = size(a, 2)
            k = size(b, 1)

            if ( n .ne. k ) error stop 'matmul: dimension of A and B not matching'

            do i = 1, m
                temp = 0.0
                do j = 1, n
                    temp = temp + a(i, j) * b(j)
                end do
                c(i) = temp
            end do

            return

        end function

        module pure function bf16_matmul_mm(a, b) result(c)
            type(bf16), intent(in) :: a(:,:)
            type(bf16), intent(in) :: b(:,:)
            type(bf16) :: c(size(a,1),size(b,2))

            integer(int64) ::  i,j,l
            integer(int64) ::  m, n, k
            type(bf16) :: temp

            m = size(a, 1)
            n = size(a, 2)
            k = size(b, 1)

            if ( n .ne. k ) error stop 'matmul: dimension of A and B not matching'

            do i = 1, size(b, 2)
                do j = 1, m
                    temp = 0.0
                    do l = 1, k
                        temp = temp + a(j, l) * b(l,i)
                    end do
                    c(j,i) =  temp
                end do
            end do

            return

        end function


end submodule
