submodule(lpf_lapack_fp16) lpf_lapack_latsqr_fp16

contains


    !> \brief \b SLATSQR
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE SLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK,
    !                           LWORK, INFO)
    !
    !       .. Scalar Arguments ..
    !       INTEGER           INFO, LDA, M, N, MB, NB, LDT, LWORK
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)              A( LDA, * ), T( LDT, * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> SLATSQR computes a blocked Tall-Skinny QR factorization of
    !> a type(fp16) M-by-N matrix A for M >= N:
    !>
    !>    A = Q * ( R ),
    !>            ( 0 )
    !>
    !> where:
    !>
    !>    Q is a M-by-M orthogonal matrix, stored on exit in an implicit
    !>    form in the elements below the diagonal of the array A and in
    !>    the elements of the array T;
    !>
    !>    R is an upper-triangular N-by-N matrix, stored on exit in
    !>    the elements on and above the diagonal of the array A.
    !>
    !>    0 is a (M-N)-by-N zero matrix, and is not stored.
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
    !>          The number of columns of the matrix A. M >= N >= 0.
    !> \endverbatim
    !>
    !> \param[in] MB
    !> \verbatim
    !>          MB is INTEGER
    !>          The row block size to be used in the blocked QR.
    !>          MB > N.
    !> \endverbatim
    !>
    !> \param[in] NB
    !> \verbatim
    !>          NB is INTEGER
    !>          The column block size to be used in the blocked QR.
    !>          N >= NB >= 1.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,N)
    !>          On entry, the M-by-N matrix A.
    !>          On exit, the elements on and above the diagonal
    !>          of the array contain the N-by-N upper triangular matrix R;
    !>          the elements below the diagonal represent Q by the columns
    !>          of blocked V (see Further Details).
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] T
    !> \verbatim
    !>          T is TYPE(FP16) array,
    !>          dimension (LDT, N * Number_of_row_blocks)
    !>          where Number_of_row_blocks = CEIL((M-N)/(MB-N))
    !>          The blocked upper triangular block reflectors stored in comp
    !>          as a sequence of upper triangular blocks.
    !>          See Further Details below.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.  LDT >= NB.
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
    !>          LWORK >= 1, if MIN(M,N) = 0, and LWORK >= NB*N, otherwise.
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
    !> \ingroup latsqr
    !>
    !  =====================================================================
    module subroutine latsqr( norm, m, n, mb, nb, a, lda, t, ldt, work, lwork, info )
        implicit none
        !     .. scalar arguments ..
        integer(lpf_default_int_kind), intent(inout) :: info, lwork
        integer(lpf_default_int_kind), intent(in) :: lda, m, n, mb, nb, ldt
        character(len=*), intent(in)  :: norm

        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) ::  a( lda, * ), work( * ), t( ldt, * )
        !     ..
        !
        !  =====================================================================
        !
        !     ..
        !     .. local scalars ..
        logical ::            lquery
        integer(lpf_default_int_kind)  ::  i, ii, kk, ctr, minmn, lwmin

        !     ..
        !     test the input arguments
        !
        info = 0
        !
        lquery = ( lwork.eq.-1 )
        !
        minmn = min( m, n )
        if( minmn.eq.0 ) then
            lwmin = 1
        else
            lwmin = n*nb
        end if
        !
        if( m.lt.0 ) then
            info = -1
        else if( n.lt.0 .or. m.lt.n ) then
            info = -2
        else if( mb.lt.1 ) then
            info = -3
        else if( nb.lt.1 .or. ( nb.gt.n .and. n.gt.0 ) ) then
            info = -4
        else if( lda.lt.max( 1, m ) ) then
            info = -6
        else if( ldt.lt.nb ) then
            info = -8
        else if( lwork.lt.lwmin .and. (.not.lquery) ) then
            info = -10
        end if
        !
        if( info.eq.0 )  then
            lwork  =  lwmin
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'latsqr', -info )
            return
        else if( lquery ) then
            return
        end if
        !
        !     quick return if possible
        !
        if( minmn.eq.0 ) then
            return
        end if
        !
        !     the qr decomposition
        !
        if( (mb.le.n) .or. (mb.ge.m) ) then
            call geqrt( norm, m, n, nb, a, lda, t, ldt, work, info )
            return
        end if
        kk = mod((m-n),(mb-n))
        ii = m-kk+1
        !
        !     compute the qr factorization of the first block a(1:mb,1:n)
        !
        call geqrt( norm, mb, n, nb, a(1,1), lda, t, ldt, work, info )
        !
        ctr = 1
        do i = mb+1, ii-mb+n, (mb-n)
            !
            !       compute the qr factorization of the current block a(i:i+mb-n,1:n
            !
            call tpqrt( norm, mb-n, n, 0_lpf_default_int_kind, nb, a(1,1), lda, a( i, 1 ), lda, t(1, ctr * n + 1), ldt, work, info)
            ctr = ctr + 1
        end do
        !
        !     compute the qr factorization of the last block a(ii:m,1:n)
        !
        if( ii.le.m ) then
            call tpqrt( norm,kk, n, 0_lpf_default_int_kind, nb, a(1,1), lda, a( ii, 1 ), lda, t(1, ctr * n + 1), ldt,   work, info)
        end if
        !
        return
        !
        !     end of slatsqr
        !
    end subroutine
end submodule
