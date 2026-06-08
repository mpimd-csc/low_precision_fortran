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

submodule(lpf_lapack_fp16) lpf_lapack_potrfp_fp16

contains
    module subroutine potrfp(uplo, n, a, lda, ipiv, info)
        character, intent(in) :: uplo
        integer(lpf_default_int_kind), intent(in) :: n, lda
        integer(lpf_default_int_kind), intent(inout) :: info
        integer(lpf_default_int_kind), intent(inout) :: ipiv(*)
        type(fp16), intent(inout) :: a(lda, *)

        ! locals
        type(fp16) :: one, mone
        logical  ::          upper
        integer(lpf_default_int_kind)  ::          maxj, k, j
        type(fp16) :: max_val, tmp

        !     ..
        !     ..
        !
        !     test the input parameters.
        !

        one = 1.0
        mone = -1.0
        info = 0
        upper = lsame( uplo, 'u' )
        if( .not.upper .and. .not.lsame( uplo, 'l' ) ) then
            info = -1
        else if( n.lt.0 ) then
            info = -2
        else if( lda.lt.max( 1, n ) ) then
            info = -4
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'potrfp', -info )
            return
        end if

        !
        !     quick return if possible
        if( n.eq.0 ) return

        if ( upper ) then
            do k = 1, n
                max_val = a(k,k)
                maxj = k;
                do j = k+1, n
                    if (a(j,j) .gt. max_val) then
                        max_val = a(j,j);
                        maxj = j
                    end if
                end do

                if ( max_val .le. fp16(0.0)) then
                    info = k
                    return
                end if

                ! Swap if necessary
                ipiv(k) = maxj;
                if ( k .ne. maxj) then
                    call swap(k-1, a(1,k), 1_lpf_default_int_kind, a(1, maxj), 1_lpf_default_int_kind)
                    call swap(maxj-k-1, a(k,k+1), lda, a(k+1, maxj), 1_lpf_default_int_kind)
                    if (maxj.lt.n) call swap(n - maxj, a(k,maxj+1), lda, a(maxj,maxj+1), lda)
                    a(maxj, maxj) = a(k,k)
                    a(k,k) = max_val
                end if

                a (k,k) = sqrt(a(k,k))
                if (k .lt. n) then

                    call scal(n-k, one/a(k,k), a(k,k+1), lda)
                    call syr("u", n-k, mone, a(k,k+1), lda, a(k+1, k+1), lda)
                end if
            end do
        else
            !lower
            do k = 1, n
                max_val = a(k,k)
                maxj = k;
                do j = k+1, n
                    if (a(j,j) .gt. max_val) then
                        max_val = a(j,j);
                        maxj = j
                    end if
                end do

                if ( max_val .le. fp16(0.0)) then
                    info = k
                    return
                end if

                ! Swap if necessary
                ipiv(k) = maxj;
                if ( k .ne. maxj) then
                    call swap(k-1, a(k,1), lda, a(maxj, 1), lda)
                    call swap(maxj-k-1, a(k+1,k), 1_lpf_default_int_kind, a(maxj, k+1), lda)
                    if (maxj .lt. n) call swap(n - maxj, a(maxj+1,k), 1_lpf_default_int_kind, a(maxj+1,maxj), &
                        & 1_lpf_default_int_kind)
                    a(maxj, maxj) = a(k,k)
                    a(k,k) = max_val
                end if

                a (k,k) = sqrt(a(k,k))

                if (k .lt. n) then
                    call scal(n-k, one/a(k,k), a(k+1,k), 1_lpf_default_int_kind)
                    call syr("l", n-k, mone, a(k+1,k), 1_lpf_default_int_kind, a(k+1, k+1), lda)
                end if
            end do

        end if
    end subroutine

end submodule
