submodule(lpf_lapack_fp16) lpf_lapack_tpqrt2_fp16

contains


    !> \brief \b STPQRT2 computes a QR factorization of a type(fp16) or complex "t
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download STPQRT2 + dependencies
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
    !       SUBROUTINE STPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER   INFO, LDA, LDB, LDT, N, M, L
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)   A( LDA, * ), B( LDB, * ), T( LDT, * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> STPQRT2 computes a QR factorization of a type(fp16) "triangular-pentagonal"
    !> matrix C, which is composed of a triangular block A and pentagonal bl
    !> using the compact WY representation for Q.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The total number of rows of the matrix B.
    !>          M >= 0.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix B, and the order of
    !>          the triangular matrix A.
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
    !>          The N-by-N upper triangular factor T of the block reflector.
    !>          See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.  LDT >= max(1,N)
    !> \endverbatim
    !>
    !> \param[out] INFO
    !> \verbatim
    !>          INFO is INTEGER
    !>          = 0: successful exit
    !>          < 0: if INFO = -i, the i-th argument had an illegal value
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
    !> \ingroup tpqrt2
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
    !>  The (M+N)-by-(M+N) block reflector H is then given by
    !>
    !>               H = I - W * T * W^H
    !>
    !>  where W^H is the conjugate transpose of W and T is the upper triangu
    !>  factor of the block reflector.
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine tpqrt2( norm, m, n, l, a, lda, b, ldb, t, ldt, info )
        implicit none
        !
        !  -- lapack computational routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        integer(lpf_default_int_kind), intent(inout) ::   info
        integer(lpf_default_int_kind), intent(in) :: lda, ldb, ldt, n, m, l
        character(len=*), intent(in)  :: norm


        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) ::  a( lda, * ), b( ldb, * ), t( ldt, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. parameters ..
        type(fp16)  ::  one, zero
        !     ..
        !     .. local scalars ..
        integer(lpf_default_int_kind) ::  i, j, p, mp, np
        type(fp16) ::  alpha
        !     ..
        !     test the input arguments
        !
        one = fp16(1.0)
        zero = fp16(0.0)
        info = 0
        if( m.lt.0 ) then
            info = -1
        else if( n.lt.0 ) then
            info = -2
        else if( l.lt.0 .or. l.gt.min(m,n) ) then
            info = -3
        else if( lda.lt.max( 1, n ) ) then
            info = -5
        else if( ldb.lt.max( 1, m ) ) then
            info = -7
        else if( ldt.lt.max( 1, n ) ) then
            info = -9
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'tpqrt2', -info )
            return
        end if
        !
        !     quick return if possible
        !
        if( n.eq.0 .or. m.eq.0 ) return
        !
        do i = 1, n
            !
            !        generate elementary reflector h(i) to annihilate b(:,i)
            !
            p = m-l+min( l, i )
            call larfg( norm, p+1, a( i, i ), b( 1, i ), 1_lpf_default_int_kind, t( i, 1 ) )
            if( i.lt.n ) then
                !
                !           w(1:n-i) := c(i:m,i+1:n)^h * c(i:m,i) [use w = t(:,n)]
                !
                do j = 1, n-i
                    t( j, n ) = (a( i, i+j ))
                end do
                call hgemv( 't', p, n-i, one, b( 1, i+1 ), ldb,             &
                    &                  b( 1, i ), 1_lpf_default_int_kind, one, t( 1, n ), 1_lpf_default_int_kind )
                !
                !           c(i:m,i+1:n) = c(i:m,i+1:n) + alpha*c(i:m,i)*w(1:n-1)^h
                !
                alpha = -(t( i, 1 ))
                do j = 1, n-i
                    a( i, i+j ) = a( i, i+j ) + alpha*(t( j, n ))
                end do
                call hger( p, n-i, alpha, b( 1, i ), 1_lpf_default_int_kind,                     &
                    &           t( 1, n ), 1_lpf_default_int_kind, b( 1, i+1 ), ldb )
            end if
        end do
        !
        do i = 2, n
            !
            !        t(1:i-1,i) := c(i:m,1:i-1)^h * (alpha * c(i:m,i))
            !
            alpha = -t( i, 1 )

            do j = 1, i-1
                t( j, i ) = zero
            end do
            p = min( i-1, l )
            mp = min( m-l+1, m )
            np = min( p+1, n )
            !
            !        triangular part of b2
            !
            do j = 1, p
                t( j, i ) = alpha*b( m-l+j, i )
            end do
            call htrmv( 'u', 't', 'n', p, b( mp, 1 ), ldb,                 &
                &               t( 1, i ), 1_lpf_default_int_kind )
            !
            !        rectangular part of b2
            !
            call hgemv( 't', l, i-1-p, alpha, b( mp, np ), ldb,            &
                &               b( mp, i ), 1_lpf_default_int_kind, zero, t( np, i ), 1_lpf_default_int_kind )
            !
            !        b1
            !
            call hgemv( 't', m-l, i-1, alpha, b, ldb, b( 1, i ), 1_lpf_default_int_kind,        &
                &               one, t( 1, i ), 1_lpf_default_int_kind )
            !
            !        t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
            !
            call htrmv( 'u', 'n', 'n', i-1, t, ldt, t( 1, i ), 1_lpf_default_int_kind )
            !
            !        t(i,i) = tau(i)
            !
            t( i, i ) = t( i, 1 )
            t( i, 1 ) = zero
        end do

        !
        !     end of stpqrt2
        !
    end subroutine
end submodule
