submodule (lpf_lapack_fp16) lpf_lapack_larf_fp16
    use lpf_fp16
    implicit none


contains
    !> \brief \b SLARF applies an elementary reflector to a general rectangular matrix.
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !        http://www.netlib.org/lapack/explore-html/
    !
    !> Download SLARF + dependencies
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarf.f">
    !> [TGZ]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarf.f">
    !> [ZIP]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarf.f">
    !> [TXT]</a>
    !
    !  Definition:
    !  ===========
    !
    !   SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
    !
    !   .. Scalar Arguments ..
    !   CHARACTER          SIDE
    !   INTEGER            INCV, LDC, M, N
    !   type(fp16)               TAU
    !   ..
    !   .. Array Arguments ..
    !   type(fp16)               C( LDC, * ), V( * ), WORK( * )
    !   ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> SLARF applies a real elementary reflector H to a real m by n matrix
    !> C, from either the left or the right. H is represented in the form
    !>
    !>       H = I - tau * v * v**T
    !>
    !> where tau is a real scalar and v is a real vector.
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
    !> \ingroup larf
    !
    !  =====================================================================
    module subroutine larf( side, m, n, v, incv, tau, c, ldc, work )
        implicit none
        !
        !  -- lapack auxiliary routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        !
        ! .. scalar arguments ..
        character, intent(in) :: side
        integer(lpf_default_int_kind), intent(in) :: incv, ldc, m, n
        type(fp16), intent(in) :: tau
        ! ..
        ! .. array arguments ..
        type(fp16), intent(in) :: v(*)
        type(fp16), intent(inout) :: c( ldc, * ), work( * )

        ! ..
        !
        !  =====================================================================
        !
        ! .. parameters ..
        type(fp16) :: one, zero
        ! ..
        ! .. local scalars ..
        logical :: applyleft
        integer(lpf_default_int_kind) :: i, lastv, lastc
        ! ..

        one  = fp16(1.0e+0)
        zero = fp16(0.0e+0)

        applyleft = lsame( side, 'l' )
        lastv = 0
        lastc = 0
        if( tau.ne.zero ) then
            ! set up variables for scanning v.  lastv begins pointing to the end
            ! of v.
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
            ! look for the last non-zero row in v.
            do while( lastv.gt.0 .and. v( i ).eq.zero )
                lastv = lastv - 1
                i = i - incv
            end do
            if( applyleft ) then
                ! scan for the last non-zero column in c(1:lastv,:).
                lastc = last_column(lastv, n, c, ldc)
            else
                ! scan for the last non-zero row in c(:,1:lastv).
                lastc = last_row(m, lastv, c, ldc)
            end if
        end if
        ! note that lastc.eq.0 renders the blas operations null; no special
        ! case is needed at this level.
        if( applyleft ) then
            !
            !    form  h * c
            !
            if( lastv.gt.0 ) then
                !
                !       w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
                !
                call hgemv( 'transpose', lastv, lastc, one, c, ldc, v, &
                    incv, &
                    zero, work, 1 )
                !
                !       c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**t
                !
                call hger( lastv, lastc, -tau, v, incv, work, 1, c, ldc )
            end if
        else
            !
            !    form  c * h
            !
            if( lastv.gt.0 ) then
                !
                !       w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                !
                call hgemv( 'no transpose', lastc, lastv, one, c, ldc, &
                    v, incv, zero, work, 1 )
                !
                !       c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**t
                !
                call hger( lastc, lastv, -tau, work, 1, v, incv, c, ldc )
            end if
        end if
        return
        !
        ! end of hlarf
        !
    end subroutine larf


end submodule
