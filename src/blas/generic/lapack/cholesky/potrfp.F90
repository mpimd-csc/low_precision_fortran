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
#define TYPEMOD lpf_fp8_e4m3
submodule(lpf_lapack_fp8_e4m3) lpf_lapack_potrfp_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
#define TYPEMOD lpf_blas_fp8_e5m2
submodule(lpf_lapack_fp8_e5m2) lpf_lapack_potrfp_fp8_e5m2
#endif
#ifdef LPF_FP16
#define TYPEMOD lpf_fp16
submodule(lpf_lapack_fp16) lpf_lapack_potrfp_fp16
#endif
#ifdef LPF_BF16
#define TYPEMOD lpf_bf16
submodule(lpf_lapack_bf16) lpf_lapack_potrfp_bf16
#endif
    use TYPEMOD
    use iso_c_binding
    implicit none


contains
    module subroutine potrfp_64(uplo, n, a_, lda, ipiv_, info)
        character, intent(in) :: uplo
        integer(int64), intent(in) :: n, lda
        integer(int64), intent(inout) :: info
        integer(int64), target, intent(inout) :: ipiv_(..)
        type(DT), target, intent(inout) :: a_(..)

        ! locals
        type(DT), pointer :: a(:,:)
        integer(int64), pointer :: ipiv(:)
        type(c_ptr) :: ptr
        type(DT) :: one, mone
        logical  ::          upper
        integer(int64)  ::          maxj, k, j
        type(DT) :: max_val, tmp

        one = 1.0
        mone = -1.0
        info = 0

        ptr = c_loc(a_)
        call c_f_pointer(ptr, a, [lda, n])
        ptr = c_loc(ipiv_)
        call c_f_pointer(ptr, ipiv, [n])


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

                if ( max_val .le. DT(0.0)) then
                    info = k
                    return
                end if

                ! Swap if necessary
                ipiv(k) = maxj;
                if ( k .ne. maxj) then
                    call swap(k-1, a(1,k), 1_int64, a(1, maxj), 1_int64)
                    call swap(maxj-k-1, a(k,k+1), lda, a(k+1, maxj), 1_int64)
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

                if ( max_val .le. DT(0.0)) then
                    info = k
                    return
                end if

                ! Swap if necessary
                ipiv(k) = maxj;
                if ( k .ne. maxj) then
                    call swap(k-1, a(k,1), lda, a(maxj, 1), lda)
                    call swap(maxj-k-1, a(k+1,k), 1_int64, a(maxj, k+1), lda)
                    if (maxj .lt. n) call swap(n - maxj, a(maxj+1,k), 1_int64, a(maxj+1,maxj), &
                        & 1_int64)
                    a(maxj, maxj) = a(k,k)
                    a(k,k) = max_val
                end if

                a (k,k) = sqrt(a(k,k))

                if (k .lt. n) then
                    call scal(n-k, one/a(k,k), a(k+1,k), 1_int64)
                    call syr("l", n-k, mone, a(k+1,k), 1_int64, a(k+1, k+1), lda)
                end if
            end do

        end if
    end subroutine

    module subroutine potrfp_32(uplo, n, a_, lda, ipiv_, info)
        character, intent(in) :: uplo
        integer(int32), intent(in) :: n, lda
        integer(int32), intent(inout) :: info
        integer(int32), target, intent(inout) :: ipiv_(..)
        type(DT), target, intent(inout) :: a_(..)

        integer(int64) :: iinfo
        integer(int64), dimension(n) :: ipiv64
        type(c_ptr) :: ptr
        integer(int32), pointer :: ipiv(:)

        ptr = c_loc(ipiv_)
        call c_f_pointer(ptr, ipiv, [n])

        call potrfp_64(uplo, int(n, int64), a_, int(lda, int64), ipiv64, iinfo)
        info = int(iinfo, int32)

        ipiv ( 1:n ) = int(ipiv64(1:n), int32)

    end subroutine
end submodule
