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
submodule(lpf_lapack_fp8_e4m3) lpf_lapack_potrf_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
#define TYPEMOD lpf_blas_fp8_e5m2
submodule(lpf_lapack_fp8_e5m2) lpf_lapack_potrf_fp8_e5m2
#endif
#ifdef LPF_FP16
#define TYPEMOD lpf_fp16
submodule(lpf_lapack_fp16) lpf_lapack_potrf_fp16
#endif
#ifdef LPF_BF16
#define TYPEMOD lpf_bf16
submodule(lpf_lapack_bf16) lpf_lapack_potrf_bf16
#endif

contains
    !> \brief \b potrf
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE bpoTRF( UPLO, N, A, LDA, INFO )
    !
    !      ..Scalar Arguments ..
    !       CHARACTER          UPLO
    !       INTEGER            INFO, LDA, N
    !       ..
    !      ..Array Arguments ..
    !       type(DT)               A( LDA, * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> bpoTRF computes the Cholesky factorization of a real symmetric
    !> positive definite matrix A.
    !>
    !> The factorization has the form
    !>    A = U**T * U,  if UPLO = 'U', or
    !>    A = L  * L**T,  if UPLO = 'L',
    !> where U is an upper triangular matrix and L is lower triangular.
    !>
    !> This is the block version of the algorithm, calling Level 3 BLAS.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] UPLO
    !> \verbatim
    !>          UPLO is CHARACTER*1
    !>          = 'U':  Upper triangle of A is stored;
    !>          = 'L':  Lower triangle of A is stored.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The order of the matrix A.  N >= 0.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is type(DT) array, dimension (LDA,N)
    !>          On entry, the symmetric matrix A.  If UPLO = 'U', the leadin
    !>          N-by-N upper triangular part of A contains the upper
    !>          triangular part of the matrix A, and the strictly lower
    !>          triangular part of A is not referenced.  If UPLO = 'L', the
    !>          leading N-by-N lower triangular part of A contains the lower
    !>          triangular part of the matrix A, and the strictly upper
    !>          triangular part of A is not referenced.
    !>
    !>          On exit, if INFO = 0, the factor U or L from the Cholesky
    !>          factorization A = U**T*U or A = L*L**T.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(1,N).
    !> \endverbatim
    !>
    !> \param[out] INFO
    !> \verbatim
    !>          INFO is INTEGER
    !>          = 0:  successful exit
    !>          < 0:  if INFO = -i, the i-th argument had an illegal value
    !>          > 0:  if INFO = i, the leading principal minor of order i
    !>                is not positive, and the factorization could not be
    !>                completed.
    !> \endverbatim
    !
    !  =====================================================================
    module subroutine potrf_64( uplo, n, a_, lda, info )
        character, intent(in) ::          uplo
        integer(int64), intent(in) :: lda, n
        integer(int64), intent(inout) :: info
        type(DT), target, intent(inout) ::  a_(..)

        type(DT), pointer :: a(:,:)
        type(c_ptr) :: ptr
        type(DT)  ::              one
        logical  ::          upper
        integer(int64)  ::          j, jb, nb

        one = 1.0
        info = 0

        ptr = c_loc(a_)
        call c_f_pointer(ptr, a, [lda, n])


        upper = lsame( uplo, 'u' )
        if( .not.upper .and. .not.lsame( uplo, 'l' ) ) then
            info = -1
        else if( n.lt.0 ) then
            info = -2
        else if( lda.lt.max( 1, n ) ) then
            info = -4
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'potrf', -info )
            return
        end if
        !
        !     quick return if possible
        !
        if( n.eq.0 )                                                      &
            &   return
        !
        !     determine the block size for this environment.
        !
        ! @todo Better block size handling
        ! Lapack would use something like:
        ! nb = ilaenv( 1, 'potrf', uplo, n, -1, -1, -1 )
        nb = 192

        if( nb.le.1 .or. nb.ge.n ) then
            !
            !        use unblocked code.
            !
            call potrf2( uplo, n, a(1, 1), lda, info )
        else
            !
            !        use blocked code.
            !
            if( upper ) then
                !
                !           compute the cholesky factorization a = u**t*u.
                !
                do j = 1, n, nb
                    !
                    !              update and factorize the current diagonal block and test
                    !              for non-positive-definiteness.
                    !
                    jb = min( nb, n-j+1 )
                    call syrk( 'upper', 'transpose', jb, j-1, -one,         &
                        &                     a( 1, j ), lda, one, a( j, j ), lda )
                    call potrf2( 'upper', jb, a( j, j ), lda, info )
                    if( info.ne.0 ) then
                        info  = info + j - 1
                        exit
                    end if

                    if( j+jb.le.n ) then
                        !
                        !                 compute the current block row.
                        !
                        ! call gemm ('t', 'n', n, n, n, one, a, lda, a, lda, one, A, lda)
                        call gemm( 'transpose', 'no transpose', jb,          &
                            &                        n-j-jb+1,                                 &
                            &                        j-1, -one, a( 1, j ), lda, a( 1, j+jb ),  &
                            &                        lda, one, a( j, j+jb ), lda )
                        call trsm( 'left', 'upper', 'transpose',             &
                            &                        'non-unit',                               &
                            &                        jb, n-j-jb+1, one, a( j, j ), lda,        &
                            &                        a( j, j+jb ), lda )
                    end if
                end do
                !
            else
                !
                !           compute the cholesky factorization a = l*l**t.
                !
                do j = 1, n, nb
                    !
                    !              update and factorize the current diagonal block and test
                    !              for non-positive-definiteness.
                    !
                    jb = min( nb, n-j+1 )
                    call syrk( 'lower', 'no transpose', jb, j-1, -one,      &
                        &                     a( j, 1 ), lda, one, a( j, j ), lda )
                    call potrf2( 'lower', jb, a( j, j ), lda, info )
                    if( info.ne.0 ) then
                        info = info+j -1
                        exit
                    end if
                    if( j+jb.le.n ) then
                        !
                        !                 compute the current block column.
                        !
                        call gemm( 'no transpose', 'transpose', n-j-jb+1,    &
                            &                        jb,                                       &
                            &                        j-1, -one, a( j+jb, 1 ), lda, a( j, 1 ),  &
                            &                        lda, one, a( j+jb, j ), lda )
                        call trsm( 'right', 'lower', 'transpose',            &
                            &                        'non-unit',                               &
                            &                        n-j-jb+1, jb, one, a( j, j ), lda,        &
                            &                        a( j+jb, j ), lda )
                    end if
                end do
            end if
        end if

        return
        !
        !     end of potrf
        !
    end subroutine

    module subroutine potrf_32( uplo, n, a_, lda, info )
        character, intent(in) ::          uplo
        integer(int32), intent(in) :: lda, n
        integer(int32), intent(inout) :: info
        type(DT), target, intent(inout) ::  a_(..)

        integer(int64) :: iinfo
        call potrf_64(uplo, int(n, int64), a_, int(lda, int64), iinfo)
        info = int(iinfo , int32)
    end subroutine

end submodule
