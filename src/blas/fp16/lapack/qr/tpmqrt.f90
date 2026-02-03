submodule(lpf_lapack_fp16) lpf_lapack_tpmqrt_fp16

contains


    !> \brief \b STPMQRT
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download STPMQRT + dependencies
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
    !       SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
    !                           A, LDA, B, LDB, WORK, INFO )
    !
    !       .. Scalar Arguments ..
    !       CHARACTER SIDE, TRANS
    !       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
    !      $          WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> STPMQRT applies a type(fp16) orthogonal matrix Q obtained from a
    !> "triangular-pentagonal" type(fp16) block reflector H to a general
    !> type(fp16) matrix C, which consists of two blocks A and B.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] SIDE
    !> \verbatim
    !>          SIDE is CHARACTER*1
    !>          = 'L': apply Q or Q^T from the Left;
    !>          = 'R': apply Q or Q^T from the Right.
    !> \endverbatim
    !>
    !> \param[in] TRANS
    !> \verbatim
    !>          TRANS is CHARACTER*1
    !>          = 'N':  No transpose, apply Q;
    !>          = 'T':  Transpose, apply Q^T.
    !> \endverbatim
    !>
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix B. M >= 0.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix B. N >= 0.
    !> \endverbatim
    !>
    !> \param[in] K
    !> \verbatim
    !>          K is INTEGER
    !>          The number of elementary reflectors whose product defines
    !>          the matrix Q.
    !> \endverbatim
    !>
    !> \param[in] L
    !> \verbatim
    !>          L is INTEGER
    !>          The order of the trapezoidal part of V.
    !>          K >= L >= 0.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] NB
    !> \verbatim
    !>          NB is INTEGER
    !>          The block size used for the storage of T.  K >= NB >= 1.
    !>          This must be the same value of NB used to generate T
    !>          in CTPQRT.
    !> \endverbatim
    !>
    !> \param[in] V
    !> \verbatim
    !>          V is TYPE(FP16) array, dimension (LDV,K)
    !>          The i-th column must contain the vector which defines the
    !>          elementary reflector H(i), for i = 1,2,...,k, as returned by
    !>          CTPQRT in B.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDV
    !> \verbatim
    !>          LDV is INTEGER
    !>          The leading dimension of the array V.
    !>          If SIDE = 'L', LDV >= max(1,M);
    !>          if SIDE = 'R', LDV >= max(1,N).
    !> \endverbatim
    !>
    !> \param[in] T
    !> \verbatim
    !>          T is TYPE(FP16) array, dimension (LDT,K)
    !>          The upper triangular factors of the block reflectors
    !>          as returned by CTPQRT, stored as a NB-by-K matrix.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.  LDT >= NB.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension
    !>          (LDA,N) if SIDE = 'L' or
    !>          (LDA,K) if SIDE = 'R'
    !>          On entry, the K-by-N or M-by-K matrix A.
    !>          On exit, A is overwritten by the corresponding block of
    !>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.
    !>          If SIDE = 'L', LDC >= max(1,K);
    !>          If SIDE = 'R', LDC >= max(1,M).
    !> \endverbatim
    !>
    !> \param[in,out] B
    !> \verbatim
    !>          B is TYPE(FP16) array, dimension (LDB,N)
    !>          On entry, the M-by-N matrix B.
    !>          On exit, B is overwritten by the corresponding block of
    !>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDB
    !> \verbatim
    !>          LDB is INTEGER
    !>          The leading dimension of the array B.
    !>          LDB >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] WORK
    !> \verbatim
    !>          WORK is TYPE(FP16) array. The dimension of WORK is
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
    !> \ingroup tpmqrt
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The columns of the pentagonal matrix V contain the elementary reflec
    !>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
    !>  trapezoidal block V2:
    !>
    !>        V = [V1]
    !>            [V2].
    !>
    !>  The size of the trapezoidal block V2 is determined by the parameter
    !>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first
    !>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper trian
    !>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
    !>
    !>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-
    !>                      [B]
    !>
    !>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N
    !>
    !>  The type(fp16) orthogonal matrix Q is formed from V and T.
    !>
    !>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
    !>
    !>  If TRANS='T' and SIDE='L', C is on exit replaced with Q^T * C.
    !>
    !>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
    !>
    !>  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q^T.
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine tpmqrt( side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info )
        implicit none
        !
        !  -- lapack computational routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        character, intent(in) :: side, trans
        integer(lpf_default_int_kind), intent(inout) ::    info
        integer(lpf_default_int_kind), intent(in) :: k, ldv, lda, ldb, m, n, l, nb, ldt
        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) ::   v( ldv, * ), a( lda, * ), b( ldb, * ), t( ldt, * ), work( * )
        !     ..
        !
        !  =====================================================================
        !
        !     ..
        !     .. local scalars ..
        logical            left, right, tran, notran
        integer(lpf_default_int_kind) :: i, ib, mb, lb, kf, ldaq, ldvq
        !     ..
        !     ..
        !
        !     .. test the input arguments ..
        !
        info   = 0
        left   = lsame( side,  'l' )
        right  = lsame( side,  'r' )
        tran   = lsame( trans, 't' )
        notran = lsame( trans, 'n' )
        !
        if ( left ) then
            ldvq = max( 1, m )
            ldaq = max( 1, k )
        else if ( right ) then
            ldvq = max( 1, n )
            ldaq = max( 1, m )
        end if
        if( .not.left .and. .not.right ) then
            info = -1
        else if( .not.tran .and. .not.notran ) then
            info = -2
        else if( m.lt.0 ) then
            info = -3
        else if( n.lt.0 ) then
            info = -4
        else if( k.lt.0 ) then
            info = -5
        else if( l.lt.0 .or. l.gt.k ) then
            info = -6
        else if( nb.lt.1 .or. (nb.gt.k .and. k.gt.0) ) then
            info = -7
        else if( ldv.lt.ldvq ) then
            info = -9
        else if( ldt.lt.nb ) then
            info = -11
        else if( lda.lt.ldaq ) then
            info = -13
        else if( ldb.lt.max( 1, m ) ) then
            info = -15
        end if
        !
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'tpmqrt', -info )
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
                mb = min( m-l+i+ib-1, m )
                if( i.ge.l ) then
                    lb = 0
                else
                    lb = mb-m+l-i+1
                end if
                call tprfb( 'l', 't', 'f', 'c', mb, n, ib, lb,             &
                    &                   v( 1, i ), ldv, t( 1, i ), ldt,                &
                    &                   a( i, 1 ), lda, b, ldb, work, ib )
            end do
            !
        else if( right .and. notran ) then
            !
            do i = 1, k, nb
                ib = min( nb, k-i+1 )
                mb = min( n-l+i+ib-1, n )
                if( i.ge.l ) then
                    lb = 0
                else
                    lb = mb-n+l-i+1
                end if
                call tprfb( 'r', 'n', 'f', 'c', m, mb, ib, lb,             &
                    &                   v( 1, i ), ldv, t( 1, i ), ldt,                &
                    &                   a( 1, i ), lda, b, ldb, work, m )
            end do
            !
        else if( left .and. notran ) then
            !
            kf = ((k-1)/nb)*nb+1
            do i = kf, 1, -nb
                ib = min( nb, k-i+1 )
                mb = min( m-l+i+ib-1, m )
                if( i.ge.l ) then
                    lb = 0
                else
                    lb = mb-m+l-i+1
                end if
                call tprfb( 'l', 'n', 'f', 'c', mb, n, ib, lb,             &
                    &                   v( 1, i ), ldv, t( 1, i ), ldt,                &
                    &                   a( i, 1 ), lda, b, ldb, work, ib )
            end do
            !
        else if( right .and. tran ) then
            !
            kf = ((k-1)/nb)*nb+1
            do i = kf, 1, -nb
                ib = min( nb, k-i+1 )
                mb = min( n-l+i+ib-1, n )
                if( i.ge.l ) then
                    lb = 0
                else
                    lb = mb-n+l-i+1
                end if
                call tprfb( 'r', 't', 'f', 'c', m, mb, ib, lb,             &
                    &                   v( 1, i ), ldv, t( 1, i ), ldt,                &
                    &                   a( 1, i ), lda, b, ldb, work, m )
            end do
            !
        end if
        !
        return
        !
        !     end of stpmqrt
        !
    end subroutine
end submodule
