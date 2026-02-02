submodule(lpf_lapack_fp16) lpf_lapack_gemqrt_fp16

contains


    !> \brief \b HGEMQRT
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download HGEMQRT + dependencies
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
    !       SUBROUTINE HGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
    !                          C, LDC, WORK, INFO )
    !
    !       .. Scalar Arguments ..
    !       CHARACTER SIDE, TRANS
    !       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
    !       ..
    !       .. Array Arguments ..
    !       REAL   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HGEMQRT overwrites the general real M-by-N matrix C with
    !>
    !>                 SIDE = 'L'     SIDE = 'R'
    !> TRANS = 'N':      Q C            C Q
    !> TRANS = 'T':   Q**T C            C Q**T
    !>
    !> where Q is a real orthogonal matrix defined as the product of K
    !> elementary reflectors:
    !>
    !>       Q = H(1) H(2) . . . H(K) = I - V T V**T
    !>
    !> generated using the compact WY representation as returned by SGEQRT.
    !>
    !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
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
    !> \param[in] NB
    !> \verbatim
    !>          NB is INTEGER
    !>          The block size used for the storage of T.  K >= NB >= 1.
    !>          This must be the same value of NB used to generate T
    !>          in SGEQRT.
    !> \endverbatim
    !>
    !> \param[in] V
    !> \verbatim
    !>          V is REAL array, dimension (LDV,K)
    !>          The i-th column must contain the vector which defines the
    !>          elementary reflector H(i), for i = 1,2,...,k, as returned by
    !>          SGEQRT in the first K columns of its array argument A.
    !> \endverbatim
    !>
    !> \param[in] LDV
    !> \verbatim
    !>          LDV is INTEGER
    !>          The leading dimension of the array V.
    !>          If SIDE = 'L', LDA >= max(1,M);
    !>          if SIDE = 'R', LDA >= max(1,N).
    !> \endverbatim
    !>
    !> \param[in] T
    !> \verbatim
    !>          T is REAL array, dimension (LDT,K)
    !>          The upper triangular factors of the block reflectors
    !>          as returned by SGEQRT, stored as a NB-by-N matrix.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.  LDT >= NB.
    !> \endverbatim
    !>
    !> \param[in,out] C
    !> \verbatim
    !>          C is REAL array, dimension (LDC,N)
    !>          On entry, the M-by-N matrix C.
    !>          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q.
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
    !>          WORK is REAL array. The dimension of WORK is
    !>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
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
    !> \ingroup gemqrt
    !
    !  =====================================================================
    module subroutine gemqrt( side, trans, m, n, k, nb, v, ldv, t, ldt,     &
            &                   c, ldc, work, info )

            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(in) :: k, ldv, ldc, m, n, nb, ldt
            integer(lpf_default_int_kind), intent(out) :: info
            type(fp16), intent(in) ::  v( ldv, * ), t(ldt, *)
            type(fp16), intent(inout) :: c( ldc, * ), work( * )

            !     .. local scalars ..
            logical ::  left, right, tran, notran
            integer(lpf_default_int_kind) ::  i, ib, ldwork, kf, q

            !
            info   = 0
            left   = lsame( side,  'l' )
            right  = lsame( side,  'r' )
            tran   = lsame( trans, 't' )
            notran = lsame( trans, 'n' )
            !
            if( left ) then
                ldwork = max( 1, n )
                q = m
            else if ( right ) then
                ldwork = max( 1, m )
                q = n
            end if
            if( .not.left .and. .not.right ) then
                info = -1
            else if( .not.tran .and. .not.notran ) then
                info = -2
            else if( m.lt.0 ) then
                info = -3
            else if( n.lt.0 ) then
                info = -4
            else if( k.lt.0 .or. k.gt.q ) then
                info = -5
            else if( nb.lt.1 .or. (nb.gt.k .and. k.gt.0)) then
                info = -6
            else if( ldv.lt.max( 1, q ) ) then
                info = -8
            else if( ldt.lt.nb ) then
                info = -10
            else if( ldc.lt.max( 1, m ) ) then
                info = -12
            end if
            !
            if( info.ne.0 ) then
                call lpf_blas_xerbla( 'gemqrt', -info )
                return
            end if
            !
            !     .. quick return if possible ..
            !
            if( m.eq.0 .or. n.eq.0 .or. k.eq.0 ) return
            !
            if( left .and. tran ) then
                !
                do i = 1, k, nb
                    ib = min( nb, k-i+1 )
                    call larfb( 'l', 't', 'f', 'c', m-i+1, n, ib,              &
                        &                   v( i, i ), ldv, t( 1, i ), ldt,                &
                        &                   c( i, 1 ), ldc, work, ldwork )
                end do
                !
            else if( right .and. notran ) then
                !
                do i = 1, k, nb
                    ib = min( nb, k-i+1 )
                    call larfb( 'r', 'n', 'f', 'c', m, n-i+1, ib,              &
                        &                   v( i, i ), ldv, t( 1, i ), ldt,                &
                        &                   c( 1, i ), ldc, work, ldwork )
                end do
                !
            else if( left .and. notran ) then
                !
                kf = ((k-1)/nb)*nb+1
                do i = kf, 1, -nb
                    ib = min( nb, k-i+1 )
                    call larfb( 'l', 'n', 'f', 'c', m-i+1, n, ib,              &
                        &                   v( i, i ), ldv, t( 1, i ), ldt,                &
                        &                   c( i, 1 ), ldc, work, ldwork )
                end do
                !
            else if( right .and. tran ) then
                !
                kf = ((k-1)/nb)*nb+1
                do i = kf, 1, -nb
                    ib = min( nb, k-i+1 )
                    call larfb( 'r', 't', 'f', 'c', m, n-i+1, ib,              &
                        &                   v( i, i ), ldv, t( 1, i ), ldt,                &
                        &                   c( 1, i ), ldc, work, ldwork )
                end do
                !
            end if
            !
            return
            !
            !     end of gemqrt
            !
        end subroutine

    end submodule
