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

#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_laset_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_laset_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding


contains

!> \brief \b LASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       REAL               ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LASET initializes an m-by-n matrix A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set; the strictly lower
!>                      triangular part of A is not changed.
!>          = 'L':      Lower triangular part is set; the strictly upper
!>                      triangular part of A is not changed.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>          The constant to which the offdiagonal elements are to be set.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>          The constant to which the diagonal elements are to be set.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On exit, the leading m-by-n submatrix of A is set as follows:
!>
!>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!>
!>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!> \ingroup OTHERauxiliary
!

    subroutine laset_impl( uplo, m, n, alpha, beta, a )
        character, intent(in)      ::   uplo
        integer(int64), intent(in) ::   m, n
        type(DT), intent(in)       ::   alpha, beta
        type(DT), intent(out)       ::   a( :, : )

        integer(int64) :: i, j

        if( uplo .eq. 'u' .or. uplo .eq. 'U') then
            ! Set the strictly upper triangular or trapezoidal part of the
            ! array to ALPHA.
            do j = 2, n
                do i = 1, min( j-1, m )
                    a( i, j ) = alpha
                end do
            end do
        else if( uplo .eq. 'l' .or. uplo .eq. 'L' ) then
            ! Set the strictly lower triangular or trapezoidal part of the
            ! array to ALPHA.
            do j = 1, min( m, n )
                do i = j + 1, m
                    a( i, j ) = alpha
                end do
            end do
        else
            ! Set the leading m-by-n submatrix to ALPHA.
            do j = 1, n
                do i = 1, m
                    a( i, j ) = alpha
                end do
            end do
        end if

        ! Set the first min(M,N) diagonal elements to BETA.
        do i = 1, min( m, n )
            a( i, i ) = beta
        end do
        return
    end subroutine

    module subroutine laset_32(uplo, m, n, alpha, beta, a, lda)
        character(c_char), intent(in) :: uplo
        integer(int32), intent(in) :: m, n, lda
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(out) :: a(..)

        type(c_ptr) :: ptr
        type(DT), pointer :: a_(:,:)
        integer(int64) :: asize(2)

        asize(1) = lda
        asize(2) = n


        ptr = c_loc(a)
        call c_f_pointer(ptr, a_, asize)

        call laset_impl(uplo, int(m, int64), int(n, int64), alpha, beta, a_)

    end subroutine

    module subroutine laset_single_32(uplo, m, n, alpha, beta, a, lda)
        character(c_char), intent(in) :: uplo
        integer(int32), intent(in) :: m, n, lda
        real(real32), intent(in) :: alpha, beta
        type(DT), target, intent(out) :: a(..)

        type(c_ptr) :: ptr
        type(DT), pointer :: a_(:,:)
        integer(int64) :: asize(2)

        asize(1) = lda
        asize(2) = n


        ptr = c_loc(a)
        call c_f_pointer(ptr, a_, asize)

        call laset_impl(uplo, int(m, int64), int(n, int64), DT(alpha), DT(beta), a_)

    end subroutine


    module subroutine laset_64(uplo, m, n, alpha, beta, a, lda)
        character(c_char), intent(in) :: uplo
        integer(int64), intent(in) :: m, n, lda
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(out) :: a(..)

        type(c_ptr) :: ptr
        type(DT), pointer :: a_(:,:)
        integer(int64) :: asize(2)


        asize(1) = lda
        asize(2) = n


        ptr = c_loc(a)
        call c_f_pointer(ptr, a_, asize)

        call laset_impl(uplo, int(m, int64), int(n, int64), alpha, beta, a_)

    end subroutine

    module subroutine laset_single_64(uplo, m, n, alpha, beta, a, lda)
        character(c_char), intent(in) :: uplo
        integer(int64), intent(in) :: m, n, lda
        real(real32), intent(in) :: alpha, beta
        type(DT), target, intent(out) :: a(..)

        type(c_ptr) :: ptr
        type(DT), pointer :: a_(:,:)
        integer(int64) :: asize(2)


        asize(1) = lda
        asize(2) = n


        ptr = c_loc(a)
        call c_f_pointer(ptr, a_, asize)

        call laset_impl(uplo, int(m, int64), int(n, int64), DT(alpha), DT(beta), a_)

    end subroutine



end submodule
