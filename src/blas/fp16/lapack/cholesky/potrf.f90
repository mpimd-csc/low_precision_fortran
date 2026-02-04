submodule(lpf_lapack_fp16) lpf_lapack_potrf_fp16
contains
    !> \brief \b hpoTRF
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download hpoTRF + dependencies
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&fil
    !> [TGZ]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&fil
    !> [ZIP]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&fil
    !> [TXT]</a>
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE hpoTRF( UPLO, N, A, LDA, INFO )
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
    !> hpoTRF computes the Cholesky factorization of a real symmetric
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
    !  Authors:
    !  ========
    !
    !> \author Univ. of Tennessee
    !> \author Univ. of California Berkeley
    !> \author Univ. of Colorado Denver
    !> \author NAG Ltd.
    !
    !> \ingroup potrf
    !
    !  =====================================================================
    module subroutine potrf( uplo, n, a, lda, info )
        character, intent(in) ::          uplo
        integer(lpf_default_int_kind), intent(in) :: lda, n
        integer(lpf_default_int_kind), intent(inout) :: info
        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) ::  a( lda, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. parameters ..
        type(fp16)  ::              one
        !     ..
        !     .. local scalars ..
        logical  ::          upper
        integer(lpf_default_int_kind)  ::          j, jb, nb
        !     ..
        !     ..
        !
        !     test the input parameters.
        !

        one = 1.0
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
        ! nb = ilaenv( 1, 'hpotrf', uplo, n, -1, -1, -1 )
        ! @todo Better block size handling
        nb = 192

        if( nb.le.1 .or. nb.ge.n ) then
            !
            !        use unblocked code.
            !
            call potrf2( uplo, n, a, lda, info )
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
                    call hsyrk( 'upper', 'transpose', jb, j-1, -one,         &
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
                        call hgemm( 'transpose', 'no transpose', jb,          &
                            &                        n-j-jb+1,                                 &
                            &                        j-1, -one, a( 1, j ), lda, a( 1, j+jb ),  &
                            &                        lda, one, a( j, j+jb ), lda )
                        call htrsm( 'left', 'upper', 'transpose',             &
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
                    call hsyrk( 'lower', 'no transpose', jb, j-1, -one,      &
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
                        call hgemm( 'no transpose', 'transpose', n-j-jb+1,    &
                            &                        jb,                                       &
                            &                        j-1, -one, a( j+jb, 1 ), lda, a( j, 1 ),  &
                            &                        lda, one, a( j+jb, j ), lda )
                        call htrsm( 'right', 'lower', 'transpose',            &
                            &                        'non-unit',                               &
                            &                        n-j-jb+1, jb, one, a( j, j ), lda,        &
                            &                        a( j+jb, j ), lda )
                    end if
                end do
            end if
        end if

        return
        !
        !     end of hpotrf
        !
    end subroutine
end submodule
