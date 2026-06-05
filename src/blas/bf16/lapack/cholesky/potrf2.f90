! SPDX-License-Identifier: LGPL-3.0-or-later
submodule(lpf_lapack_bf16) lpf_lapack_potrf2_bf16
contains
    !> \brief \b POTRF2
    !
    !  Definition:
    !  ===========
    !
    !       RECURSIVE SUBROUTINE bpoTRF2( UPLO, N, A, LDA, INFO )
    !
    !       .. Scalar Arguments ..
    !       CHARACTER          UPLO
    !       INTEGER            INFO, LDA, N
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
    !> bpoTRF2 computes the Cholesky factorization of a real symmetric
    !> positive definite matrix A using the recursive algorithm.
    !>
    !> The factorization has the form
    !>    A = U**T * U,  if UPLO = 'U', or
    !>    A = L  * L**T,  if UPLO = 'L',
    !> where U is an upper triangular matrix and L is lower triangular.
    !>
    !> This is the recursive version of the algorithm. It divides
    !> the matrix into four submatrices:
    !>
    !>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
    !>    A = [ -----|----- ]  with n1 = n/2
    !>        [  A21 | A22  ]       n2 = n-n1
    !>
    !> The subroutine calls itself to factor A11. Update and scale A21
    !> or A12, update A22 then call itself to factor A22.
    !>
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
    !>          A is REAL array, dimension (LDA,N)
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
    module recursive subroutine potrf2( uplo, n, a, lda, info )
        !     .. scalar arguments ..
        character, intent(in) ::          uplo
        integer(lpf_default_int_kind), intent(in) ::  lda, n
        integer(lpf_default_int_kind), intent(inout) :: info
        !     ..
        !     .. array arguments ..
        type(bf16), intent(inout) ::               a( lda, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. parameters ..
        type(bf16)   ::            one, zero
        !     ..
        !     .. local scalars ..
        logical ::           upper
        integer(lpf_default_int_kind) ::           n1, n2, iinfo
        !     ..
        !     test the input parameters
        !

        one = 1.0
        zero = 0.0
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
            call lpf_blas_xerbla( 'potrf2', -info )
            return
        end if
        !
        !     quick return if possible
        !
        if( n.eq.0 )   then
            return
        end if
        !
        !     n=1 case
        !
        if( n.eq.1 ) then
            !
            !        test for non-positive-definiteness
            !
            if( a( 1, 1 ).le.zero.or.isnan( a( 1, 1 ) ) ) then
                info = 1
                return
            end if
            !
            !        factor
            !
            a( 1, 1 ) = sqrt( a( 1, 1 ) )
            !
            !     use recursive code
            !
        else
            n1 = n/2
            n2 = n-n1
            !
            !        factor a11
            !
            call potrf2( uplo, n1, a( 1, 1 ), lda, iinfo )
            if ( iinfo.ne.0 ) then
                info = iinfo
                return
            end if
            !
            !        compute the cholesky factorization a = u**t*u
            !
            if( upper ) then
                !
                !           update and scale a12
                !
                call trsm( 'l', 'u', 't', 'n', n1, n2, one, a( 1, 1 ), lda, a( 1, n1+1 ), lda )
                !
                !           update and factor a22
                !
                call syrk( uplo, 't', n2, n1, -one, a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
                call potrf2( uplo, n2, a( n1+1, n1+1 ), lda, iinfo )
                if ( iinfo.ne.0 ) then
                    info = iinfo + n1
                    return
                end if
                !
                !        compute the cholesky factorization a = l*l**t
                !
            else
                !
                !           update and scale a21
                !
                call trsm( 'r', 'l', 't', 'n', n2, n1, one, a( 1, 1 ), lda, a( n1+1, 1 ), lda )
                !
                !           update and factor a22
                !
                call syrk( uplo, 'n', n2, n1, -one, a( n1+1, 1 ), lda, one, a( n1+1, n1+1 ), lda )
                call potrf2( uplo, n2, a( n1+1, n1+1 ), lda, iinfo )
                if ( iinfo.ne.0 ) then
                    info = iinfo + n1
                    return
                end if
            end if
        end if
        return
        !
        !     end of bpotrf2
        !
    end subroutine
end submodule
