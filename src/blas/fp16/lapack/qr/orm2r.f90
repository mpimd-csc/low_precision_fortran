submodule(lpf_lapack_fp16) lpf_lapack_orm2r_fp16

contains

!> \brief \b SORM2R multiplies a general matrix by the orthogonal matrix
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORM2R + dependencies
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
!       SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SORM2R overwrites the general real m by n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**T if SIDE = 'R' and TRANS = 'T',
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by SGEQRF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left
!>          = 'R': apply Q or Q**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'T': apply Q**T (Transpose)
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
!>          A is REAL array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          SGEQRF in the first k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
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
!>          TAU is REAL array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by SGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
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
!>          WORK is REAL array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
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
!> \ingroup unm2r
!
!  =====================================================================
      module subroutine orm2r(  side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
!
!  -- lapack computational routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd
!
!     .. scalar arguments ..
      character, intent(in) :: side, trans
      integer(lpf_default_int_kind), intent(in)   :: k, lda, ldc, m, n
      integer(lpf_default_int_kind), intent(inout) :: info
!     ..
!     .. array arguments ..
      type(fp16), intent(in)     ::          a( lda, * ), tau( * )
      type(fp16), intent(inout)     ::          c( ldc, * ), work( * )
!     ..
!
!  =====================================================================
!
!     .. local scalars ..
      logical            :: left, notran
      integer(lpf_default_int_kind) :: i, i1, i2, i3, ic, jc, mi, ni, nq
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
      info = 0
      left = lsame( side, 'l' )
      notran = lsame( trans, 'n' )
!
!     nq is the order of q
!
      if( left ) then
         nq = m
      else
         nq = n
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
      end if
      if( info.ne.0 ) then
         call lpf_blas_xerbla( 'horm2r', -info )
         return
      end if
!
!     quick return if possible
!
      if( m.eq.0 .or. n.eq.0 .or. k.eq.0 )                              &
     &   return
!
      if( ( left .and. .not.notran ) .or. ( .not.left .and. notran ) )  &
     &     then
         i1 = 1
         i2 = k
         i3 = 1
      else
         i1 = k
         i2 = 1
         i3 = -1
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
      do i = i1, i2, i3
         if( left ) then
!
!           h(i) is applied to c(i:m,1:n)
!
            mi = m - i + 1
            ic = i
         else
!
!           h(i) is applied to c(1:m,i:n)
!
            ni = n - i + 1
            jc = i
         end if
!
!        apply h(i)
!
         call larf1f( side, mi, ni, a( i, i ), 1, tau( i ), c( ic,     &
     &               jc ), ldc, work )
      end do
      return
!
!     end of sorm2r
!
      END
end submodule
