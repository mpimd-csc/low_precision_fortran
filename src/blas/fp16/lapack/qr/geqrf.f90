submodule(lpf_lapack_fp16) lpf_lapack_geqrf_fp16

contains


    !> \brief \b GEQRF
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> \htmlonly
    !> Download GEQRF + dependencies
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&fil
    !> [TGZ]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&fil
    !> [ZIP]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&fil
    !> [TXT]</a>
    !> \endhtmlonly
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE GEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER            INFO, LDA, LWORK, M, N
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)               A( LDA, * ), TAU( * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> GEQRF computes a QR factorization of a real M-by-N matrix A:
    !>
    !>    A = Q * ( R ),
    !>            ( 0 )
    !>
    !> where:
    !>
    !>    Q is a M-by-M orthogonal matrix;
    !>    R is an upper-triangular N-by-N matrix;
    !>    0 is a (M-N)-by-N zero matrix, if M > N.
    !>
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
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
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,N)
    !>          On entry, the M-by-N matrix A.
    !>          On exit, the elements on and above the diagonal of the array
    !>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
    !>          upper triangular if m >= n); the elements below the diagonal
    !>          with the array TAU, represent the orthogonal matrix Q as a
    !>          product of min(m,n) elementary reflectors (see Further
    !>          Details).
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] TAU
    !> \verbatim
    !>          TAU is TYPE(FP16) array, dimension (min(M,N))
    !>          The scalar factors of the elementary reflectors (see Further
    !>          Details).
    !> \endverbatim
    !>
    !> \param[out] WORK
    !> \verbatim
    !>          WORK is TYPE(FP16) array, dimension (MAX(1,LWORK))
    !>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    !> \endverbatim
    !>
    !> \param[in] LWORK
    !> \verbatim
    !>          LWORK is INTEGER
    !>          The dimension of the array WORK.
    !>          LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
    !>          For optimum performance LWORK >= N*NB, where NB is
    !>          the optimal blocksize.
    !>
    !>          If LWORK = -1, then a workspace query is assumed; the routin
    !>          only calculates the optimal size of the WORK array, returns
    !>          this value as the first entry of the WORK array, and no erro
    !>          message related to LWORK is issued by XERBLA.
    !> \endverbatim
    !>
    !> \param[out] INFO
    !> \verbatim
    !>          INFO is INTEGER
    !>          = 0:  successful exit
    !>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
    !> \ingroup geqrf
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The matrix Q is represented as a product of elementary reflectors
    !>
    !>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
    !>
    !>  Each H(i) has the form
    !>
    !>     H(i) = I - tau * v * v**T
    !>
    !>  where tau is a real scalar, and v is a real vector with
    !>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    !>  and tau in TAU(i).
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine geqrf( norm, m, n, a, lda, tau, work, lwork, info ) bind(c, name="hgeqrf_")

        !
        !  -- lapack computational routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        integer(lpf_default_int_kind), intent(inout) ::            info, lwork
        integer(lpf_default_int_kind), intent(in) :: lda, m, n
        character(len=*), intent(in)  :: norm
        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout)  :: a( lda, * ), tau( * ), work( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. local scalars ..
        logical            lquery
        integer(lpf_default_int_kind)  ::   i, ib, iinfo, iws, k, ldwork, lwkopt, nb,  nbmin, nx
        !     ..

        !
        !     test the input arguments
        !
        k = min( m, n )
        info = 0
        !
        !     @todo Better ILAENV nb handling
        !
        nb = 64

        lquery = ( lwork.eq.-1 )
        if( m.lt.0 ) then
            info = -1
        else if( n.lt.0 ) then
            info = -2
        else if( lda.lt.max( 1, m ) ) then
            info = -4
        else if( .not.lquery ) then
            if( lwork.le.0 .or. ( m.gt.0 .and. lwork.lt.max( 1, n ) ) )    &
                &      info = -7
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'geqrf', -info )
            return
        else if( lquery ) then
            if( k.eq.0 ) then
                lwkopt = 1
            else
                lwkopt = n*nb
            end if
            lwork = lwkopt
            return
        end if
        !
        !     quick return if possible
        !
        if( k.eq.0 ) then
            work( 1 ) = 1
            return
        end if
        !
        nbmin = 2
        nx = 0
        ldwork = n
        if( nb.gt.1 .and. nb.lt.k ) then
            !
            !        determine when to cross over from blocked to unblocked code.
            !

            ! @todo. Cross over point nx needs to obtianed from ilanv, now: nx = nb
            nx = nb
        end if
        !
        if( nb.lt.k .and. nx.lt.k ) then
            !
            !        use blocked code initially
            !
            do i = 1, k - nx, nb
                ib = min( k-i+1, nb )
                !
                !           compute the qr factorization of the current block
                !           a(i:m,i:i+ib-1)
                !
                call geqr2( norm, m-i+1, ib, a( i, i ), lda, tau( i ), work, iinfo )
                if( i+ib.le.n ) then
                    !
                    !              form the triangular factor of the block reflector
                    !              h = h(i) h(i+1) . . . h(i+ib-1)
                    !
                    call larft( 'forward', 'columnwise', m-i+1, ib, a( i, i ), lda, tau( i ), work, ldwork )
                    !
                    !              apply h**t to a(i:m,i+ib:n) from the left
                    !
                    call larfb( 'left', 'transpose', 'forward', 'columnwise', m-i+1, n-i-ib+1, ib,   &
                        &                      a( i, i ), lda, work, ldwork, a( i, i+ib ), lda, work( ib+1 ), ldwork )
                end if
            end do
        else
            i = 1
        end if
        !
        !     use unblocked code to factor the last or only block.
        !
        if( i.le.k ) then
            call geqr2( norm, m-i+1, n-i+1, a( i, i ), lda, tau( i ), work, iinfo )
        end if
        !
        return
        !
        !     end of geqrf
        !
        end
end submodule

