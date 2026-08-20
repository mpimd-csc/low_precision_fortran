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
submodule (lpf_blas_fp8_e5m2) lpf_blas_lacpy_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_lacpy_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding


contains


!> \brief \b LACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download LACPY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slacpy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slacpy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slacpy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE LACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
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
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!>          or trapezoid is accessed; if UPLO = 'L', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
    subroutine lacpy_impl( uplo, m, n, a, b )
        character, intent(in)      ::   uplo
        integer(int64), intent(in) ::   m, n
        type(DT), intent(in)       ::   a( :, : )
        type(DT), intent(inout)    ::   b( :, : )

        integer(int64) :: i, j

        if( uplo .eq. 'u' .or. uplo .eq. 'U') then
            do j = 1, n
                do i = 1, min( j, m )
                    b( i, j ) = a( i, j )
                end do
            end do
        else if( uplo .eq. 'l' .or. uplo .eq. 'L' ) then
            do j = 1, n
                do i = j, m
                    b( i, j ) = a( i, j )
                end do
            end do
        else
            do j = 1, n
                do i = 1, m
                    b( i, j ) = a( i, j )
                end do
            end do
        end if
        return
    end subroutine

    module subroutine lacpy_32(uplo, m, n, a, lda, b, ldb)
        character(c_char), intent(in) :: uplo
        integer(int32), intent(in) :: m, n, lda, ldb
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(inout) :: b(..)

        type(c_ptr) :: ptr
        type(DT), pointer :: a_(:,:), b_(:,:)
        integer(int64) :: asize(2), bsize(2)

        asize(1) = lda
        asize(2) = n

        bsize(1) = ldb
        bsize(2) = n

        ptr = c_loc(a)
        call c_f_pointer(ptr, a_, asize)

        ptr = c_loc(b)
        call c_f_pointer(ptr, b_, bsize)

        call lacpy_impl(uplo, int(m, int64), int(n, int64), a_, b_)

    end subroutine

    module subroutine lacpy_64(uplo, m, n, a, lda, b, ldb)
        character(c_char), intent(in) :: uplo
        integer(int64), intent(in) :: m, n, lda, ldb
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(inout) :: b(..)

        type(c_ptr) :: ptr
        type(DT), pointer :: a_(:,:), b_(:,:)
        integer(int64) :: asize(2), bsize(2)


        asize(1) = lda
        asize(2) = n

        bsize(1) = ldb
        bsize(2) = n

        ptr = c_loc(a)
        call c_f_pointer(ptr, a_, asize)

        ptr = c_loc(b)
        call c_f_pointer(ptr, b_, bsize)

        call lacpy_impl(uplo, int(m, int64), int(n, int64), a_, b_)

    end subroutine


end submodule
