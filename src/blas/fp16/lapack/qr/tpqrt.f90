submodule(lpf_lapack_fp16) lpf_lapack_tpqrt_fp16

contains


    !> \brief \b STPQRT
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download STPQRT + dependencies
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
    !       SUBROUTINE STPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,
    !                          INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16) A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> STPQRT computes a blocked QR factorization of a type(fp16)
    !> "triangular-pentagonal" matrix C, which is composed of a
    !> triangular block A and pentagonal block B, using the compact
    !> WY representation for Q.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix B.
    !>          M >= 0.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix B, and the order of the
    !>          triangular matrix A.
    !>          N >= 0.
    !> \endverbatim
    !>
    !> \param[in] L
    !> \verbatim
    !>          L is INTEGER
    !>          The number of rows of the upper trapezoidal part of B.
    !>          MIN(M,N) >= L >= 0.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] NB
    !> \verbatim
    !>          NB is INTEGER
    !>          The block size to be used in the blocked QR.  N >= NB >= 1.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,N)
    !>          On entry, the upper triangular N-by-N matrix A.
    !>          On exit, the elements on and above the diagonal of the array
    !>          contain the upper triangular matrix R.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(1,N).
    !> \endverbatim
    !>
    !> \param[in,out] B
    !> \verbatim
    !>          B is TYPE(FP16) array, dimension (LDB,N)
    !>          On entry, the pentagonal M-by-N matrix B.  The first M-L row
    !>          are rectangular, and the last L rows are upper trapezoidal.
    !>          On exit, B contains the pentagonal matrix V.  See Further De
    !> \endverbatim
    !>
    !> \param[in] LDB
    !> \verbatim
    !>          LDB is INTEGER
    !>          The leading dimension of the array B.  LDB >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] T
    !> \verbatim
    !>          T is TYPE(FP16) array, dimension (LDT,N)
    !>          The upper triangular block reflectors stored in compact form
    !>          as a sequence of upper triangular blocks.  See Further Detai
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
    !>          WORK is TYPE(FP16) array, dimension (NB*N)
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
    !> \ingroup tpqrt
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The input matrix C is a (N+M)-by-N matrix
    !>
    !>               C = [ A ]
    !>                   [ B ]
    !>
    !>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentag
    !>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a
    !>  upper trapezoidal matrix B2:
    !>
    !>               B = [ B1 ]  <- (M-L)-by-N rectangular
    !>                   [ B2 ]  <-     L-by-N upper trapezoidal.
    !>
    !>  The upper trapezoidal matrix B2 consists of the first L rows of a
    !>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
    !>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
    !>
    !>  The matrix W stores the elementary reflectors H(i) in the i-th colum
    !>  below the diagonal (of A) in the (N+M)-by-N input matrix C
    !>
    !>               C = [ A ]  <- upper triangular N-by-N
    !>                   [ B ]  <- M-by-N pentagonal
    !>
    !>  so that W can be represented as
    !>
    !>               W = [ I ]  <- identity, N-by-N
    !>                   [ V ]  <- M-by-N, same form as B.
    !>
    !>  Thus, all of information needed for W is contained on exit in B, whi
    !>  we call V above.  Note that V has the same form as B; that is,
    !>
    !>               V = [ V1 ] <- (M-L)-by-N rectangular
    !>                   [ V2 ] <-     L-by-N upper trapezoidal.
    !>
    !>  The columns of V represent the vectors which define the H(i)'s.
    !>
    !>  The number of blocks is B = ceiling(N/NB), where each
    !>  block is of order NB except for the last block, which is of order
    !>  IB = N - (B-1)*NB.  For each of the B blocks, a upper triangular blo
    !>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB
    !>  for the last block) T's are stored in the NB-by-N matrix T as
    !>
    !>               T = [T1 T2 ... TB].
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine tpqrt( norm, m, n, l, nb, a, lda, b, ldb, t, ldt, work, info )
        implicit none
        !
        !  -- lapack computational routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        integer(lpf_default_int_kind), intent(inout) ::  info
        integer(lpf_default_int_kind), intent(in) :: lda, ldb, ldt, n, m, l, nb
        character(len=*), intent(in)  :: norm

        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) :: a( lda, * ), b( ldb, * ), t( ldt, * ), work( * )
        !     ..
        !
        ! =====================================================================
        !
        !     ..
        !     .. local scalars ..
        integer(lpf_default_int_kind) ::  i, ib, lb, mb, iinfo
        !     ..
        !
        !     test the input arguments
        !
        info = 0
        if( m.lt.0 ) then
            info = -1
        else if( n.lt.0 ) then
            info = -2
        else if( l.lt.0 .or. (l.gt.min(m,n) .and. min(m,n).ge.0)) then
            info = -3
        else if( nb.lt.1 .or. (nb.gt.n .and. n.gt.0)) then
            info = -4
        else if( lda.lt.max( 1, n ) ) then
            info = -6
        else if( ldb.lt.max( 1, m ) ) then
            info = -8
        else if( ldt.lt.nb ) then
            info = -10
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'tpqrt', -info )
            return
        end if
        !
        !     quick return if possible
        !
        if( m.eq.0 .or. n.eq.0 ) return
        !
        do i = 1, n, nb
            !
            !     compute the qr factorization of the current block
            !
            ib = min( n-i+1, nb )
            mb = min( m-l+i+ib-1, m )
            if( i.ge.l ) then
                lb = 0
            else
                lb = mb-m+l-i+1
            end if
            !
            call tpqrt2( norm, mb, ib, lb, a(i,i), lda, b( 1, i ), ldb,         &
                &                 t(1, i ), ldt, iinfo )
            !
            !     update by applying h^h to b(:,i+ib:n) from the left
            !
            if( i+ib.le.n ) then
                call tprfb( 'l', 't', 'f', 'c', mb, n-i-ib+1, ib, lb,      &
                    &                    b( 1, i ), ldb, t( 1, i ), ldt,               &
                    &                    a( i, i+ib ), lda, b( 1, i+ib ), ldb,         &
                    &                    work, ib )
            end if
        end do
        return
        !
        !     end of stpqrt
        !
    end subroutine
end submodule
