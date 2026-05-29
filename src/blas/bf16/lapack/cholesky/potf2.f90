submodule(lpf_lapack_bf16) lpf_lapack_potf2_bf16

contains
    !> \brief \b SPOTF2 computes the Cholesky factorization of a symmetric/B
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> SPOTF2 computes the Cholesky factorization of a real symmetric
    !> positive definite matrix A.
    !>
    !> The factorization has the form
    !>    A = U**T * U ,  if UPLO = 'U', or
    !>    A = L  * L**T,  if UPLO = 'L',
    !> where U is an upper triangular matrix and L is lower triangular.
    !>
    !> This is the unblocked version of the algorithm, calling Level 2 BLAS.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] UPLO
    !> \verbatim
    !>          UPLO is CHARACTER*1
    !>          Specifies whether the upper or lower triangular part of the
    !>          symmetric matrix A is stored.
    !>          = 'U':  Upper triangular
    !>          = 'L':  Lower triangular
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
    !>          n by n upper triangular part of A contains the upper
    !>          triangular part of the matrix A, and the strictly lower
    !>          triangular part of A is not referenced.  If UPLO = 'L', the
    !>          leading n by n lower triangular part of A contains the lower
    !>          triangular part of the matrix A, and the strictly upper
    !>          triangular part of A is not referenced.
    !>
    !>          On exit, if INFO = 0, the factor U or L from the Cholesky
    !>          factorization A = U**T *U  or A = L*L**T.
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
    !>          = 0: successful exit
    !>          < 0: if INFO = -k, the k-th argument had an illegal value
    !>          > 0: if INFO = k, the leading principal minor of order k
    !>               is not positive, and the factorization could not be
    !>               completed.
    !> \endverbatim
    !
    module subroutine potf2( uplo, n, a, lda, info )
        implicit none

        character, intent(in) :: uplo
        integer(lpf_default_int_kind), intent(in)   :: lda, n
        integer(lpf_default_int_kind), intent(inout) :: info
        !     ..
        !     .. array arguments ..
        type(bf16), intent(inout) :: a( lda, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. parameters ..
        type(bf16)               :: one, zero
        !     ..
        !     .. local scalars ..
        logical            :: upper
        integer(lpf_default_int_kind)            :: j
        type(bf16)         :: ajj
        !     ..
        !
        !     test the input parameters.
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
            call lpf_blas_xerbla( 'potf2', -info )
            return
        end if
        !
        !     quick return if possible
        !
        if( n.eq.0 ) then
            return
        end if
        !
        if( upper ) then
            !
            !        compute the cholesky factorization a = u**t *u.
            !
            do j = 1, n
                !
                !           compute u(j,j) and test for non-positive-definiteness.
                !
                ajj = a( j, j ) - dot( j-1, a( 1, j ), 1_lpf_default_int_kind, a( 1, j ), 1_lpf_default_int_kind )
                if( ajj.le.zero.or.isnan( ajj ) ) then
                    a( j, j ) = ajj
                    info = j
                    exit
                end if
                ajj = sqrt( ajj )
                a( j, j ) = ajj
                !
                !           compute elements j+1:n of row j.
                !
                if( j.lt.n ) then
                    call gemv( 'transpose', j-1, n-j, -one, a( 1, j+1 ), lda, a( 1, j ), 1_lpf_default_int_kind, &
                        & one, a( j, j+1 ), lda )
                    call scal( n-j, one / ajj, a( j, j+1 ), lda )
                end if
            end do
        else
            !
            !        compute the cholesky factorization a = l*l**t.
            !
            do j = 1, n
                !
                !           compute l(j,j) and test for non-positive-definiteness.
                !
                ajj = a( j, j ) - dot( j-1, a( j, 1 ), lda, a( j, 1 ), lda )
                if( ajj.le.zero.or.isnan( ajj ) ) then
                    a( j, j ) = ajj
                    info = j
                    exit
                end if
                ajj = sqrt( ajj )
                a( j, j ) = ajj
                !
                !           compute elements j+1:n of column j.
                !
                if( j.lt.n ) then
                    call gemv( 'no transpose', n-j, j-1, -one, a( j+1, 1 ), lda, a( j, 1 ), lda, one, a( j+1, j ), &
                        & 1_lpf_default_int_kind )
                    call scal( n-j, one / ajj, a( j+1, j ), 1_lpf_default_int_kind )
                end if
            end do
        end if
        !
        return
        !
        !     end of bpotf2
        !
    end subroutine
end submodule
