submodule (lpf_lapack_fp16) lpf_lapack_larfgv1_fp16
    use lpf_fp16
    use lpf_types
    implicit none


contains

    !> \brief \b HLARFG generates an elementary reflector (Householder matrix).
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !        http://www.netlib.org/lapack/explore-html/
    !
    !> Download HLARFG + dependencies
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfg.f">
    !> [TGZ]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfg.f">
    !> [ZIP]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfg.f">
    !> [TXT]</a>
    !
    !  Definition:
    !  ===========
    !
    !   SUBROUTINE HLARFG( NORM, N, ALPHA, V1, X, INCX, TAU )
    !
    !   .. Scalar Arguments ..
    !   INTEGER            INCX, N
    !   TYPE(FP16)               ALPHA, V1, TAU
    !   ..
    !   .. Array Arguments ..
    !   CHARACTER          NORM(*)
    !   TYPE(FP16)               X( * )
    !   ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HLARFG generates a real elementary reflector H of order n, such
    !> that
    !>
    !>       H * ( alpha ) = ( beta ),   H**T * H = I.
    !>           (   x   )   (   0  )
    !>
    !> where alpha and beta are scalars, and x is an (n-1)-element real
    !> vector. H is represented in the form
    !>
    !>       H = I - tau  *( v1 ) ( v1 v**T ) ,
    !>                    ( v  )
    !>
    !> where tau, v1 are real scalars and v is a real (n-1)-element
    !> vector. The vector (v1 v) is normalized, such that || ( v1 v ) ||_2 = 1
    !> and tau = 2
    !>
    !> If the elements of x are all zero, then tau = 0 and H is taken to be
    !> the unit matrix.
    !>
    !> Otherwise  1 <= tau <= 2.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] NORM
    !> \verbatim
    !>          NORM is CHARACTER(*)
    !>          If NORM = "H", the 2-norms are computed using
    !>          an FP32 accumulator. Otherwise, the accumulator
    !>          is in low precision as well.
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The order of the elementary reflector.
    !> \endverbatim
    !>
    !> \param[in,out] ALPHA
    !> \verbatim
    !>          ALPHA is TYPE(FP16)
    !>          On entry, the value alpha.
    !>          On exit, it is overwritten with the value beta.
    !> \endverbatim
    !>
    !> \param[out] V1
    !> \verbatim
    !>          ALPHA is V1
    !>          On entry, the value alpha.
    !>          On exit, it is overwritten with the value beta.
    !> \endverbatim

    !>
    !> \param[in,out] X
    !> \verbatim
    !>          X is TYPE(FP16) array, dimension
    !>                         (1+(N-2)*abs(INCX))
    !>          On entry, the vector x.
    !>          On exit, it is overwritten with the vector v.
    !> \endverbatim
    !>
    !> \param[in] INCX
    !> \verbatim
    !>          INCX is INTEGER
    !>          The increment between elements of X. INCX > 0.
    !> \endverbatim
    !>
    !> \param[out] TAU
    !> \verbatim
    !>          TAU is TYPE(FP16)
    !>          The value tau.
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
    !> \ingroup larfg
    !
    !  =====================================================================
    module subroutine larfgv1( norm, n, alpha, v1,  x, incx, tau )
        integer(lpf_default_int_kind), intent(in) :: incx, n
        character,  intent(in)    :: norm(*)
        type(fp16), intent(inout) :: alpha
        type(fp16), intent(out)   :: v1
        type(fp16), intent(out)   :: tau
        type(fp16), intent(inout) :: x( * )

        ! .. parameters ..

        type(fp16) :: one, zero
        ! ..
        ! .. local scalars ..
        integer(lpf_default_int_kind) :: j, knt
        type(fp16) :: beta, rsafmn, safmin, xnorm
        type(fp16) :: s
        logical :: high_precision_norm

        one = 1.0
        zero = 0.0

        high_precision_norm = lsame(norm(1), "H")
        if( n.le.1 ) then
            tau = zero
            v1  = zero
            return
        end if

        if ( high_precision_norm ) then
            xnorm = hnrm2_fp32( n-1, x, incx )
        else
            xnorm = hnrm2( n-1, x, incx )
        end if

        if( xnorm.eq.zero ) then
            !
            !    h  =  i
            !
            tau = zero
            v1  = zero
        else
            !
            !    general case
            !
            if ( high_precision_norm ) then
                beta = -sign( lapy2_fp32( alpha, xnorm ), alpha )
            else
                beta = -sign( lapy2( alpha, xnorm ), alpha )
            end if
            ! if ( alpha**2.0 .gt. xnorm ) then
            !     beta = sqrt(xnorm) * sqrt(fp16(1.0) + alpha**2.0 / xnorm)
            ! else
            !     beta = alpha * sqrt(fp16(1.0) + xnorm / alpha**2.0)
            ! end if
            ! beta = -sign(beta, alpha)

            safmin = lamch( 's' ) / lamch( 'e' )
            knt = 0
            if( abs( beta ).lt.safmin ) then
                !
                !       xnorm, beta may be inaccurate; scale x and recompute them
                !
                rsafmn = one / safmin
                10    continue
                knt = knt + 1
                call hscal( n-1, rsafmn, x, incx )
                beta = beta*rsafmn
                alpha = alpha*rsafmn
                if( (abs( beta ).lt.safmin) .and. (knt .lt. 20) ) &
                    go to 10
                !
                !       new beta is at most 1, at least safmin
                !
                if (high_precision_norm) then
                    xnorm = hnrm2_fp32( n-1, x, incx )
                    beta = -sign( lapy2_fp32( alpha, xnorm ), alpha )
                else
                    xnorm = hnrm2( n-1, x, incx )
                    beta = -sign( lapy2( alpha, xnorm ), alpha )
                end if

            end if


            tau = 2
            if ( high_precision_norm ) then
                s = one / lapy2_fp32(xnorm, (alpha-beta))
            else
                s = one / lapy2(xnorm, (alpha-beta))
            end if

            v1  = (alpha-beta) * s
            call hscal( n-1, s , x, incx )
            !
            !    if alpha is subnormal, it may lose relative accuracy
            !
            do j = 1, knt
                beta = beta*safmin
            end do
            alpha = beta
        end if
        !
        return
        !
        ! end of hlarfg
        !
    end subroutine larfgv1


end submodule
