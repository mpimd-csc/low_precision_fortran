submodule (lpf_lapack_fp16) lpf_lapack_larf1f_fp16
    use lpf_fp16
    use lpf_types
    implicit none


contains

!> \brief \b SLARF1F applies an elementary reflector to a general rectan
!              matrix assuming v(1) = 1.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARF1F + dependencies
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
!       SUBROUTINE SLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       type(fp16)               TAU
!       ..
!       .. Array Arguments ..
!       type(fp16)               C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARF1F applies a real elementary reflector H to a real m by n matrix
!> C, from either the left or the right. H is represented in the form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector assuming v(1) = 1.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is type(fp16) array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is type(fp16)
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is type(fp16) array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
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
!>          WORK is type(fp16) array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
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
!> \ingroup larf1f
!
!  =====================================================================
      module subroutine larf1f( side, m, n, v, incv, tau, c, ldc, work )
!
!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd
!
!     .. scalar arguments ..
      character, intent(in) :: side
      integer(lpf_default_int_kind), intent(in) ::incv, ldc, m, n
      type(fp16), intent(in) :: tau
!     ..
!     .. array arguments ..
      type(fp16), intent(inout) ::  c( ldc, * ), work( * )
      type(fp16), intent(in) ::  v( * )
!     ..
!
!  =====================================================================
!
!     .. parameters ..
      type(fp16) ::      one, zero
!     ..
!     .. local scalars ..
      logical   ::     applyleft
      integer(lpf_default_int_kind) :: i, lastv, lastc
      integer(lpf_default_int_kind), parameter :: ione = 1

      !
!     .. executable statements ..
!

      one  = 1.0
      zero = 0.0
      applyleft = lsame( side, 'l' )
      lastv = 1
      lastc = 0
      if( tau.ne.zero ) then
!     set up variables for scanning v.  lastv begins pointing to the end
!     of v up to v(1).
         if( applyleft ) then
            lastv = m
         else
            lastv = n
         end if
         if( incv.gt.0 ) then
            i = 1 + (lastv-1) * incv
         else
            i = 1
         end if
!     look for the last non-zero row in v.
         do while( lastv.gt.1 .and. v( i ).eq.zero )
            lastv = lastv - 1
            i = i - incv
         end do
         if( applyleft ) then
!     scan for the last non-zero column in c(1:lastv,:).
            lastc = last_column(lastv, n, c, ldc)
         else
!     scan for the last non-zero row in c(:,1:lastv).
            lastc = last_row(m, lastv, c, ldc)
         end if
      end if
      if( lastc.eq.0 ) then
         return
      end if
      if( applyleft ) then
!
!        form  h * c
!
         if( lastv.eq.1 ) then
!
!           c(1,1:lastc) := ( 1 - tau ) * c(1,1:lastc)
!
            call hscal( lastc, one - tau, c, ldc )
         else
!
!        w(1:lastc,1) := c(2:lastv,1:lastc)**t * v(2:lastv,1)
!
         call hgemv( 'transpose', lastv - 1, lastc, one, c( 2, 1 ),     &
     &               ldc, v( 1 + incv ), incv, zero, work, ione )
!
!        w(1:lastc,1) += v(1,1) * c(1,1:lastc)**t
!
         call haxpy( lastc, one, c, ldc, work, ione )
!
!        c(1, 1:lastc) += - tau * v(1,1) * w(1:lastc,1)**t
!
         call haxpy( lastc, -tau, work, ione, c, ldc )
!
!        c(2:lastv,1:lastc) += - tau * v(2:lastv,1) * w(1:lastc,1)**t
!
         call hger( lastv - 1, lastc, -tau, v( 1 + incv ), incv,        &
     &              work, ione, c( 2, 1 ), ldc )
            end if
      else
!
!        form  c * h
!
         if( lastv.eq.1 ) then
!
!           c(1:lastc,1) := ( 1 - tau ) * c(1:lastc,1)
!
            call hscal( lastc, one - tau, c, ione )
         else
!
!           w(1:lastc,1) := c(1:lastc,2:lastv) * v(2:lastv,1)
!
            call hgemv( 'no transpose', lastc, lastv - 1, one,          &
     &                  c( 1, 2 ), ldc, v( 1 + incv ), incv, zero,      &
     &                  work, ione )
!
!           w(1:lastc,1) += v(1,1) * c(1:lastc,1)
!
            call haxpy( lastc, one, c, ione, work, ione )
!
!           c(1:lastc,1) += - tau * v(1,1) * w(1:lastc,1)
!
            call haxpy( lastc, -tau, work, ione, c, ione )
!
!           c(1:lastc,2:lastv) += - tau * w(1:lastc,1) * v(2:lastv)**t
!
            call hger( lastc, lastv - 1, -tau, work, ione,                 &
     &                 v( 1 + incv ), incv, c( 1, 2 ), ldc )
         end if
      end if
      return
!
!     end of larf1f
!
      end subroutine

end submodule
