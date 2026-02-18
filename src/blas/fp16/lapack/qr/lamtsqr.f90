submodule(lpf_lapack_fp16) lpf_lapack_lamtsqr_fp16

contains
    !> \brief \b SLAMTSQR
    !
    !  Definition:
    !  ===========
    !
    !      SUBROUTINE SLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T,
    !     $                     LDT, C, LDC, WORK, LWORK, INFO )
    !
    !
    !     .. Scalar Arguments ..
    !      CHARACTER         SIDE, TRANS
    !      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
    !     ..
    !     .. Array Arguments ..
    !      DOUBLE        A( LDA, * ), WORK( * ), C(LDC, * ),
    !     $                  T( LDT, * )
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !>      SLAMTSQR overwrites the general type(fp16) M-by-N matrix C with
    !>
    !>
    !>                 SIDE = 'L'     SIDE = 'R'
    !> TRANS = 'N':      Q * C          C * Q
    !> TRANS = 'T':      Q**T * C       C * Q**T
    !>      where Q is a type(fp16) orthogonal matrix defined as the product
    !>      of blocked elementary reflectors computed by tall skinny
    !>      QR factorization (SLATSQR)
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
    !>          The number of rows of the matrix A.  M >=0.
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
    !>          the matrix Q. M >= K >= 0;
    !>
    !> \endverbatim
    !>
    !> \param[in] MB
    !> \verbatim
    !>          MB is INTEGER
    !>          The block size to be used in the blocked QR.
    !>          MB > N. (must be the same as SLATSQR)
    !> \endverbatim
    !>
    !> \param[in] NB
    !> \verbatim
    !>          NB is INTEGER
    !>          The column block size to be used in the blocked QR.
    !>          N >= NB >= 1.
    !> \endverbatim
    !>
    !> \param[in] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,K)
    !>          The i-th column must contain the vector which defines the
    !>          blockedelementary reflector H(i), for i = 1,2,...,k, as
    !>          returned by SLATSQR in the first k columns of
    !>          its array argument A.
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
    !> \param[in] T
    !> \verbatim
    !>          T is TYPE(FP16) array, dimension
    !>          ( N * Number of blocks(CEIL(M-K/MB-K)),
    !>          The blocked upper triangular block reflectors stored in comp
    !>          as a sequence of upper triangular blocks.  See below
    !>          for further details.
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
    !>          (workspace) TYPE(FP16) array, dimension (MAX(1,LWORK))
    !>          On exit, if INFO = 0, WORK(1) returns the minimal LWORK.
    !> \endverbatim
    !>
    !> \param[in] LWORK
    !> \verbatim
    !>          LWORK is INTEGER
    !>          The dimension of the array WORK.
    !>          If MIN(M,N,K) = 0, LWORK >= 1.
    !>          If SIDE = 'L', LWORK >= max(1,N*NB).
    !>          If SIDE = 'R', LWORK >= max(1,MB*NB).
    !>
    !>          If LWORK = -1, then a workspace query is assumed; the routin
    !>          only calculates the minimal size of the WORK array, returns
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
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !> Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transfo
    !> representing Q as a product of other orthogonal matrices
    !>   Q = Q(1) * Q(2) * . . . * Q(k)
    !> where each Q(i) zeros out subdiagonal entries of a block of MB rows o
    !>   Q(1) zeros out the subdiagonal entries of rows 1:MB of A
    !>   Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A
    !>   Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N]
    !>   . . .
    !>
    !> Q(1) is computed by GEQRT, which represents Q(1) by Householder vecto
    !> stored under the diagonal of rows 1:MB of A, and by upper triangular
    !> block reflectors, stored in array T(1:LDT,1:N).
    !> For more information see Further Details in GEQRT.
    !>
    !> Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Household
    !> stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper trian
    !> block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N).
    !> The last Q(k) may use fewer rows.
    !> For more information see Further Details in TPQRT.
    !>
    !> For more details of the overall algorithm, see the description of
    !> Sequential TSQR in Section 2.2 of [1].
    !>
    !> [1] “Communication-Optimal Parallel and Sequential QR and LU Factor
    !>     J. Demmel, L. Grigori, M. Hoemmen, J. Langou,
    !>     SIAM J. Sci. Comput, vol. 34, no. 1, 2012
    !> \endverbatim
    !>
    !> \ingroup lamtsqr
    !>
    !  =====================================================================
    module subroutine lamtsqr( side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info )
        implicit none
        !
        !  -- lapack computational routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        character, intent(in) :: side, trans
        integer(lpf_default_int_kind), intent(inout) :: info, lwork
        integer(lpf_default_int_kind), intent(in) :: lda, m, n, k, mb, nb, ldt, ldc
        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) ::  a( lda, * ), work( * ), c( ldc, * ), t( ldt, * )
        !     ..
        !
        ! =====================================================================
        !
        !     ..
        !     .. local scalars ..
        logical ::           left, right, tran, notran, lquery
        integer(lpf_default_int_kind)  ::  i, ii, kk, lw, ctr, q, minmnk, lwmin
        !     ..
        !     ..
        !
        !     test the input arguments
        !
        info = 0
        lquery  = ( lwork.eq.-1 )
        notran  = lsame( trans, 'n' )
        tran    = lsame( trans, 't' )
        left    = lsame( side, 'l' )
        right   = lsame( side, 'r' )
        if( left ) then
            lw = n * nb
            q = m
        else
            lw = mb * nb
            q = n
        end if
        !
        minmnk = min( m, n, k )
        if( minmnk.eq.0 ) then
            lwmin = 1
        else
            lwmin = max( 1, lw )
        end if
        !
        if( .not.left .and. .not.right ) then
            info = -1
        else if( .not.tran .and. .not.notran ) then
            info = -2
        else if( m.lt.k ) then
            info = -3
        else if( n.lt.0 ) then
            info = -4
        else if( k.lt.0 ) then
            info = -5
        else if( k.lt.nb .or. nb.lt.1 ) then
            info = -7
        else if( lda.lt.max( 1, q ) ) then
            info = -9
        else if( ldt.lt.max( 1, nb ) ) then
            info = -11
        else if( ldc.lt.max( 1, m ) ) then
            info = -13
        else if( lwork.lt.lwmin .and. (.not.lquery) ) then
            info = -15
        end if
        !
        if( info.eq.0 ) then
            lwork = lwmin
        end if
        !
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'lamtsqr', -info )
            return
        else if( lquery ) then
            return
        end if
        !
        !     quick return if possible
        !
        if( minmnk.eq.0 ) then
            return
        end if
        !
        !     determine the block size if it is tall skinny or short and wide
        !
        if((mb.le.k).or.(mb.ge.max(m,n,k))) then
            call gemqrt( side, trans, m, n, k, nb, a, lda,                 &
                &        t, ldt, c, ldc, work, info )
            return
        end if
        !
        if(left.and.notran) then
            !
            !         multiply q to the last block of c
            !
            kk = mod((m-k),(mb-k))
            ctr = (m-k)/(mb-k)
            if (kk.gt.0) then
                ii=m-kk+1
                call tpmqrt('l','n',kk , n, k, 0_lpf_default_int_kind, nb, a(ii,1), lda,         &
                    &       t(1,ctr*k+1),ldt , c(1,1), ldc,                            &
                    &       c(ii,1), ldc, work, info )
            else
                ii=m+1
            end if
            !
            do i=ii-(mb-k),mb+1,-(mb-k)
                !
                !         multiply q to the current block of c (i:i+mb,1:n)
                !
                ctr = ctr - 1
                call tpmqrt('l','n',mb-k , n, k, 0_lpf_default_int_kind,nb, a(i,1), lda,         &
                    &         t(1, ctr * k + 1), ldt, c(1,1), ldc,                     &
                    &         c(i,1), ldc, work, info )
                !
            end do
            !
            !         multiply q to the first block of c (1:mb,1:n)
            !
            call gemqrt('l','n',mb , n, k, nb, a(1,1), lda, t             &
                &            ,ldt ,c(1,1), ldc, work, info )
            !
        else if (left.and.tran) then
            !
            !         multiply q to the first block of c
            !
            kk = mod((m-k),(mb-k))
            ii=m-kk+1
            ctr = 1
            call gemqrt('l','t',mb , n, k, nb, a(1,1), lda, t             &
                &            ,ldt ,c(1,1), ldc, work, info )
            !
            do i=mb+1,ii-mb+k,(mb-k)
                !
                !         multiply q to the current block of c (i:i+mb,1:n)
                !
                call tpmqrt('l','t',mb-k , n, k, 0_lpf_default_int_kind,nb, a(i,1), lda,          &
                    &       t(1,ctr * k + 1),ldt, c(1,1), ldc,                         &
                    &       c(i,1), ldc, work, info )
                ctr = ctr + 1
                !
            end do
            if(ii.le.m) then
                !
                !         multiply q to the last block of c
                !
                call tpmqrt('l','t',kk , n, k, 0_lpf_default_int_kind,nb, a(ii,1), lda,           &
                    &      t(1, ctr * k + 1), ldt, c(1,1), ldc,                        &
                    &      c(ii,1), ldc, work, info )
                !
            end if
            !
        else if(right.and.tran) then
            !
            !         multiply q to the last block of c
            !
            kk = mod((n-k),(mb-k))
            ctr = (n-k)/(mb-k)
            if (kk.gt.0) then
                ii=n-kk+1
                call tpmqrt('r','t',m , kk, k, 0_lpf_default_int_kind, nb, a(ii,1), lda,        &
                    &        t(1, ctr * k + 1), ldt, c(1,1), ldc,                      &
                    &        c(1,ii), ldc, work, info )
            else
                ii=n+1
            end if
            !
            do i=ii-(mb-k),mb+1,-(mb-k)
                !
                !         multiply q to the current block of c (1:m,i:i+mb)
                !
                ctr = ctr - 1
                call tpmqrt('r','t',m , mb-k, k, 0_lpf_default_int_kind,nb, a(i,1), lda,        &
                    &          t(1, ctr * k + 1), ldt, c(1,1), ldc,                    &
                    &          c(1,i), ldc, work, info )
                !
            end do
            !
            !         multiply q to the first block of c (1:m,1:mb)
            !
            call gemqrt('r','t',m , mb, k, nb, a(1,1), lda, t            &
                &              ,ldt ,c(1,1), ldc, work, info )
            !
        else if (right.and.notran) then
            !
            !         multiply q to the first block of c
            !
            kk = mod((n-k),(mb-k))
            ii=n-kk+1
            ctr = 1
            call gemqrt('r','n', m, mb , k, nb, a(1,1), lda, t            &
                &              ,ldt ,c(1,1), ldc, work, info )
            !
            do i=mb+1,ii-mb+k,(mb-k)
                !
                !         multiply q to the current block of c (1:m,i:i+mb)
                !
                call tpmqrt('r','n', m, mb-k, k, 0_lpf_default_int_kind,nb, a(i,1), lda,          &
                    &         t(1, ctr * k + 1),ldt, c(1,1), ldc,                      &
                    &         c(1,i), ldc, work, info )
                ctr = ctr + 1
                !
            end do
            if(ii.le.n) then
                !
                !         multiply q to the last block of c
                !
                call tpmqrt('r','n', m, kk , k, 0_lpf_default_int_kind,nb, a(ii,1), lda,          &
                    &        t(1, ctr * k + 1),ldt, c(1,1), ldc,                       &
                    &        c(1,ii), ldc, work, info )
                !
            end if
            !
        end if
        !
        return
        !
        !     end of slamtsqr
        !
    end subroutine
end submodule
