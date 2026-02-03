submodule(lpf_lapack_fp16) lpf_lapack_ormqrv_fp16

contains


    !> \brief \b SORMQR
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> \htmlonly
    !> Download SORMQR + dependencies
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
    !       SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
    !                          WORK, LWORK, INFO )
    !
    !       .. Scalar Arguments ..
    !       CHARACTER          SIDE, TRANS
    !       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)               A( LDA, * ), C( LDC, * ), TAU( * ),
    !      $                   WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> SORMQR overwrites the general real M-by-N matrix C with
    !>
    !>                 SIDE = 'L'     SIDE = 'R'
    !> TRANS = 'N':      Q * C          C * Q
    !> TRANS = 'T':      Q**T * C       C * Q**T
    !>
    !> where Q is a real orthogonal matrix defined as the product of k
    !> elementary reflectors
    !>
    !>       Q = H(1) H(2) . . . H(k)
    !>
    !> as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N
    !> if SIDE = 'R'.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] SIDE
    !> \verbatim
    !>          SIDE is CHARACTER*1
    !>          = 'L': apply Q or Q**T from the Left;
    !>          = 'R': apply Q or Q**T from the Right.
    !> \endverbatim
    !>
    !> \param[in] TRANS
    !> \verbatim
    !>          TRANS is CHARACTER*1
    !>          = 'N':  No transpose, apply Q;
    !>          = 'T':  Transpose, apply Q**T.
    !> \endverbatim
    !>
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix C. M >= 0.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix C. N >= 0.
    !> \endverbatim
    !>
    !> \param[in] K
    !> \verbatim
    !>          K is INTEGER
    !>          The number of elementary reflectors whose product defines
    !>          the matrix Q.
    !>          If SIDE = 'L', M >= K >= 0;
    !>          if SIDE = 'R', N >= K >= 0.
    !> \endverbatim
    !>
    !> \param[in] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,K)
    !>          The i-th column must contain the vector which defines the
    !>          elementary reflector H(i), for i = 1,2,...,k, as returned by
    !>          SGEQRF in the first k columns of its array argument A.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.
    !>          If SIDE = 'L', LDA >= max(1,M);
    !>          if SIDE = 'R', LDA >= max(1,N).
    !> \endverbatim
    !>
    !> \param[in] TAU
    !> \verbatim
    !>          TAU is TYPE(FP16) array, dimension (K)
    !>          TAU(i) must contain the scalar factor of the elementary
    !>          reflector H(i), as returned by SGEQRF.
    !> \endverbatim
    !>
    !> \param[in,out] C
    !> \verbatim
    !>          C is TYPE(FP16) array, dimension (LDC,N)
    !>          On entry, the M-by-N matrix C.
    !>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
    !> \endverbatim
    !>
    !> \param[in] LDC
    !> \verbatim
    !>          LDC is INTEGER
    !>          The leading dimension of the array C. LDC >= max(1,M).
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
    !>          If SIDE = 'L', LWORK >= max(1,N);
    !>          if SIDE = 'R', LWORK >= max(1,M).
    !>          For good performance, LWORK should generally be larger.
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
    !> \ingroup unmqr
    !
    !  =====================================================================
    module subroutine ormqrv( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
        !
        !  -- lapack computational routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        character, intent(in) ::          side, trans
        integer(lpf_default_int_kind), intent(inout) ::            info, lwork
        integer(lpf_default_int_kind), intent(in)    ::   lda, ldc, m, n, k
        !     ..
        !     .. array arguments ..
        type(fp16), intent(in)  ::  a( lda, * ), tau(*)
        type(fp16), intent(inout) :: c( ldc, * ), work( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. parameters ..
        integer(lpf_default_int_kind) :: nbmax, ldt, tsize

        !     ..
        !     .. local scalars ..
        logical ::           left, lquery, notran
        integer(lpf_default_int_kind) :: i, i1, i2, i3, ib, ic, iinfo, iwt, jc, ldwork, lwkopt, mi, nb, nbmin, ni, nq, nw
        !     ..
        !     .. external functions ..
        !     ..
        !     .. intrinsic functions ..
        intrinsic          max, min
        !     ..
        !     .. executable statements ..
        !
        !     test the input arguments
        !
        nbmax = 64
        nb = 64
        ldt = nbmax+1
        tsize = ldt*nbmax
        info = 0
        left = lsame( side, 'l' )
        notran = lsame( trans, 'n' )
        lquery = ( lwork.eq.-1 )
        !
        !     nq is the order of q and nw is the minimum dimension of work
        !
        if( left ) then
            nq = m
            nw = max( 1, n )
        else
            nq = n
            nw = max( 1, m )
        end if
        if( .not.left .and. .not.lsame( side, 'r' ) ) then
            info = -1
        else if( .not.notran .and. .not.lsame( trans, 't' ) ) then
            info = -2
        else if( m.lt.0 ) then
            info = -3
        else if( n.lt.0 ) then
            info = -4
        else if( k.lt.0 .or. k.gt.nq ) then
            info = -5
        else if( lda.lt.max( 1, nq ) ) then
            info = -7
        else if( ldc.lt.max( 1, m ) ) then
            info = -10
        else if( lwork.lt.nw .and. .not.lquery ) then
            info = -12
        end if
        !
        if( info.eq.0 ) then
            !
            !        compute the workspace requirements
            !
            lwkopt = nw*nb + tsize
            lwork = lwkopt
        end if
        !
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'ormqr', -info )
            return
        else if( lquery ) then
            return
        end if
        !
        !     quick return if possible
        !
        if( m.eq.0 .or. n.eq.0 .or. k.eq.0 ) then
            work( 1 ) = 1
            return
        end if
        !
        nbmin = 2
        ldwork = nw
        if( nb.gt.1 .and. nb.lt.k ) then
            if( lwork.lt.lwkopt ) then
                nb = (lwork-tsize) / ldwork
                nbmin = 2
            end if
        end if
        !
        if( nb.lt.nbmin .or. nb.ge.k ) then
            !
            !        use unblocked code
            !
            call orm2rv( side, trans, m, n, k, a, lda, tau, c, ldc,        &
                &                work,                                             &
                &                iinfo )
        else
            !
            !        use blocked code
            !
            iwt = 1 + nw*nb
            if( ( left .and. .not.notran ) .or.                            &
                &       ( .not.left .and. notran ) ) then
                i1 = 1
                i2 = k
                i3 = nb
            else
                i1 = ( ( k-1 ) / nb )*nb + 1
                i2 = 1
                i3 = -nb
            end if
            !
            if( left ) then
                ni = n
                jc = 1
            else
                mi = m
                ic = 1
            end if
            !
            do  i = i1, i2, i3
                ib = min( nb, k-i+1 )
                !
                !           form the triangular factor of the block reflector
                !           h = h(i) h(i+1) . . . h(i+ib-1)
                !
                call larftv( 'forward', 'columnwise', nq-i+1, ib, a( i, i ),                                           &
                    &                   lda, tau( i ), work( iwt ), ldt )
                if( left ) then
                    !
                    !              h or h**t is applied to c(i:m,1:n)
                    !
                    mi = m - i + 1
                    ic = i
                else
                    !
                    !              h or h**t is applied to c(1:m,i:n)
                    !
                    ni = n - i + 1
                    jc = i
                end if
                !
                !           apply h or h**t
                !
                call larfbv( side, trans, 'forward', 'columnwise', mi,  ni,                                            &
                    &                   ib, a( i, i ), lda, work( iwt ), ldt,          &
                    &                   c( ic, jc ), ldc, work, ldwork )
            end do
        end if
        work( 1 )  = lwork
        return
        !
        !     end of sormqr
        !
    end subroutine
end submodule
