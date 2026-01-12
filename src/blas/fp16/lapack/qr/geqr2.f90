submodule(lpf_lapack_fp16) lpf_lapack_geqr2_fp16

contains

!> \brief \b SGEQR2 computes the QR factorization of a general rectangul
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQR2 + dependencies
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
!       SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEQR2 computes a QR factorization of a real m-by-n matrix A:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a m-by-m orthogonal matrix;
!>    R is an upper-triangular n-by-n matrix;
!>    0 is a (m-n)-by-n zero matrix, if m > n.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(m,n) by n upper trapezoidal matrix R (R is
!>          upper triangular if m >= n); the elements below the diagonal
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of elementary reflectors (see Further Details).
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
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
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
!> \ingroup geqr2
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
      module subroutine geqr2( norm, m, n, a, lda, tau, work, info ) bind(c, name="hgeqr2_")
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd
!
!     .. scalar arguments ..
      integer(lpf_default_int_kind), intent(in) ::  lda, m, n
      integer(lpf_default_int_kind), intent(inout) ::  info
!     ..
!     .. array arguments ..
      character(len=*), intent(in)  :: norm
      type(fp16), intent(inout)     :: a( lda, * ), tau( * ), work( * )
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
      integer(lpf_default_int_kind)  ::      i, k
      integer(lpf_default_int_kind), parameter :: ione = 1
      type(fp16) :: aii
!
!     test the input arguments
!
      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call lpf_blas_xerbla( 'hgeqr2', -info )
         return
      end if
!
      k = min( m, n )
!
      do i = 1, k
!
!        generate elementary reflector h(i) to annihilate a(i+1:m,i)
!
         call larfg( norm, m-i+1, a( i, i ), a( min( i+1, m ), i ), ione, tau( i ) )
         if( i.lt.n ) then
!
!           apply h(i) to a(i:m,i+1:n) from the left
!
            aii = a(i,i)
            a(i,i) = fp16(1.0)
            call larf( 'left', m-i+1, n-i, a( i, i ), ione, tau( i ), a( i, i+1 ), lda, work )
            a(i,i) = aii
         end if
      end do
      return
!
!     end of sgeqr2
!
      end

end submodule
