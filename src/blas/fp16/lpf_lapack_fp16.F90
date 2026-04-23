module lpf_lapack_fp16
    use lpf_fp16
    use lpf_blas_fp16
    use iso_fortran_env, only : real32, real64
    use iso_c_binding
    use lpf_types
    implicit none


    !
    ! Norms
    !
    interface
        module function lange( norm, m, n, a, lda, work ) result(nrm)
            character, intent(in) :: norm
            integer(lpf_default_int_kind), intent(in)   :: lda, m, n
            type(fp16), intent(in)      :: a( lda, * )
            type(fp16), intent(inout)   :: work( * )
            type(fp16) :: nrm
        end function
    end interface

    interface
        module subroutine lassq( n, x, incx, sca, sumsq )
            integer(lpf_default_int_kind), intent(in) ::   incx, n
            real(real32), intent(inout) :: sca, sumsq
            type(fp16), intent(in)  :: x( * )
        end subroutine
    end interface


    !
    ! Cholesky Decompostion
    !
    interface
        module subroutine potf2( uplo, n, a, lda, info )
            character, intent(in) :: uplo
            integer(lpf_default_int_kind), intent(in)   :: lda, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout) :: a( lda, * )
        end subroutine
    end interface

    interface
        module subroutine potrf( uplo, n, a, lda, info )
         character, intent(in) ::          uplo
         integer(lpf_default_int_kind), intent(in) :: lda, n
         integer(lpf_default_int_kind), intent(inout) :: info
         type(fp16), intent(inout) ::  a( lda, * )
        end subroutine
    end interface

    interface
        module recursive subroutine potrf2( uplo, n, a, lda, info )
            character, intent(in) ::          uplo
            integer(lpf_default_int_kind), intent(in) ::  lda, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout) ::               a( lda, * )
        end subroutine
    end interface

    interface
        module subroutine hpotrfp(uplo, n, a, lda, ipiv, info)
            character, intent(in) :: uplo
            integer(lpf_default_int_kind), intent(in) :: n, lda
            integer(lpf_default_int_kind), intent(inout) :: info
            integer(lpf_default_int_kind), intent(inout) :: ipiv(*)
            type(fp16), intent(inout) :: a(lda, *)
        end subroutine
    end interface

    !
    ! Householder related functions
    !
    interface
        module subroutine larf( side, m, n, v, incv, tau, c, ldc, work )
            implicit none
            character, intent(in) :: side
            integer(lpf_default_int_kind), intent(in) :: incv, ldc, m, n
            type(fp16), intent(in) :: tau
            type(fp16), intent(in) :: v(*)
            type(fp16), intent(inout) :: c( ldc, * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine larfg( norm, n, alpha, x, incx, tau )
            integer(lpf_default_int_kind), intent(in) :: incx, n
            character,  intent(in)    :: norm(*)
            type(fp16), intent(inout) :: alpha
            type(fp16), intent(out)   ::  tau
            type(fp16), intent(inout) :: x( * )
        end subroutine
    end interface

    interface
        module subroutine larfgv1( norm, n, alpha, v1, x, incx, tau )
            integer(lpf_default_int_kind), intent(in) :: incx, n
            character,  intent(in)    :: norm(*)
            type(fp16), intent(inout) :: alpha
            type(fp16), intent(out)   :: v1
            type(fp16), intent(out)   :: tau
            type(fp16), intent(inout) :: x( * )
        end subroutine
    end interface


    interface
        module subroutine larf1f( side, m, n, v, incv, tau, c, ldc, work )
            character, intent(in) :: side
            integer(lpf_default_int_kind), intent(in) ::incv, ldc, m, n
            type(fp16), intent(in) :: tau
            type(fp16), intent(inout) ::  c( ldc, * ), work( * )
            type(fp16), intent(in) ::  v( * )
        end subroutine
    end interface

    interface
        module subroutine larfb( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork ) &
                & bind(c, name="hlarfb_")
            character, intent(in) :: direct, side, storev, trans
            integer(lpf_default_int_kind), intent(in) :: k, ldc, ldt, ldv, ldwork, m, n
            type(fp16), intent(in) :: t( ldt, * ), v( ldv, * )
            type(fp16), intent(inout) :: c( ldc, * ), work( ldwork, * )
        end subroutine
    end interface
    interface
        module recursive subroutine larft( direct, storev, n, k, v, ldv, tau, t, ldt ) bind(c, name="hlarft_")
            character, intent(in) ::         direct, storev
            integer(lpf_default_int_kind), intent(in)   ::     k, ldt, ldv, n
            type(fp16), intent(out) ::   t( ldt, * )
            type(fp16), intent(in)  ::   tau( * ), v( ldv, * )
        end subroutine
    end interface

    interface
        module subroutine larfbv( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork ) &
                & bind(c, name="hlarfbv_")
            character, intent(in) :: direct, side, storev, trans
            integer(lpf_default_int_kind), intent(in) :: k, ldc, ldt, ldv, ldwork, m, n
            type(fp16), intent(in) :: t( ldt, * ), v( ldv, * )
            type(fp16), intent(inout) :: c( ldc, * ), work( ldwork, * )
        end subroutine
    end interface
    interface
        module recursive subroutine larftv( direct, storev, n, k, v, ldv, tau, t, ldt ) bind(c, name="hlarftv_")
            character, intent(in) ::         direct, storev
            integer(lpf_default_int_kind), intent(in)   ::     k, ldt, ldv, n
            type(fp16), intent(out) ::   t( ldt, * )
            type(fp16), intent(in)  ::   tau( * ), v( ldv, * )
        end subroutine
    end interface

    interface
        module subroutine tprfb( side, trans, direct, storev, m, n, k, l,       &
                &                   v, ldv, t, ldt, a, lda, b, ldb, work, ldwork )
            character, intent(in) ::  direct, side, storev, trans
            integer(lpf_default_int_kind), intent(in) ::  k, l, lda, ldb, ldt, ldv, ldwork, m, n
            type(fp16), intent(inout)  ::  a( lda, * ), b( ldb, * ), t( ldt, * ), v( ldv, * ), work( ldwork, * )
        end subroutine
    end interface

    !
    ! Householder QR Decompositions
    !
    interface
        module subroutine geqr2( norm, m, n, a, lda, tau, work, info ) bind(c, name="hgeqr2_")
            integer(lpf_default_int_kind), intent(in) ::  lda, m, n
            integer(lpf_default_int_kind), intent(inout) ::  info
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout)     :: a( lda, * ), tau( * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine geqr2v( norm, m, n, a, lda, diagr, tau, work, info ) bind(c, name="hgeqr2v_")
            integer(lpf_default_int_kind), intent(in) ::  lda, m, n
            integer(lpf_default_int_kind), intent(inout) ::  info
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout)     :: a( lda, * ), diagr(*), tau( * ), work( * )
        end subroutine
    end interface


    interface
        module subroutine orm2r(  side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(in)   :: k, lda, ldc, m, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(in)     ::          a( lda, * ), tau( * )
            type(fp16), intent(inout)     ::          c( ldc, * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine orm2rv(  side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(in)   :: k, lda, ldc, m, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(in)     ::          a( lda, * ), tau( * )
            type(fp16), intent(inout)     ::          c( ldc, * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine geqrf( norm, m, n, a, lda, tau, work, lwork, info ) bind(c, name="hgeqrf_")
            integer(lpf_default_int_kind),  intent(inout) ::            info, lwork
            integer(lpf_default_int_kind), intent(in) :: lda, m, n
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout)  :: a( lda, * ), tau( * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine geqrfv( norm, m, n, a, lda, diagr, tau, work, lwork, info )
            integer(lpf_default_int_kind),  intent(inout) ::            info, lwork
            integer(lpf_default_int_kind), intent(in) :: lda, m, n
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout)  :: a( lda, * ), diagr(*), tau( * ), work( * )
        end subroutine
    end interface


    interface
        module subroutine gemqrt( side, trans, m, n, k, nb, v, ldv, t, ldt,     &
            &                   c, ldc, work, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(in) :: k, ldv, ldc, m, n, nb, ldt
            integer(lpf_default_int_kind), intent(out) :: info
            type(fp16), intent(in) ::  v( ldv, * ), t(ldt, *)
            type(fp16), intent(inout) :: c( ldc, * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine gemqrtv( side, trans, m, n, k, nb, v, ldv, t, ldt,     &
            &                   c, ldc, work, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(in) :: k, ldv, ldc, m, n, nb, ldt
            integer(lpf_default_int_kind), intent(out) :: info
            type(fp16), intent(in) ::  v( ldv, * ), t(ldt, *)
            type(fp16), intent(inout) :: c( ldc, * ), work( * )
        end subroutine
    end interface


    interface
        module subroutine ormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
            character, intent(in) ::          side, trans
            integer(lpf_default_int_kind), intent(inout) ::            info, lwork
            integer(lpf_default_int_kind), intent(in)    ::   lda, ldc, m, n, k
            type(fp16), intent(in)  ::  a( lda, * ), tau(*)
            type(fp16), intent(inout) :: c( ldc, * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine ormqrv( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info )
            character, intent(in) ::          side, trans
            integer(lpf_default_int_kind), intent(inout) ::            info, lwork
            integer(lpf_default_int_kind), intent(in)    ::   lda, ldc, m, n, k
            type(fp16), intent(in)  ::  a( lda, * ), tau(*)
            type(fp16), intent(inout) :: c( ldc, * ), work( * )
        end subroutine
    end interface

    interface
        module subroutine geqrt( norm, m, n, nb, a, lda, t, ldt, work, info )
            integer(lpf_default_int_kind), intent(in)  :: lda, ldt, m, n, nb
            integer(lpf_default_int_kind), intent(out) :: info
            type(fp16), intent(inout) :: a( lda, * ), t( ldt, * ), work( * )
            character(len=*), intent(in)  :: norm
        end subroutine
    end interface
    interface
        module subroutine geqrt2( norm, m, n, a, lda, t, ldt, info )
            integer(lpf_default_int_kind), intent(in) :: lda, ldt, m, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout) :: a( lda, * ), t( ldt, * )
            character(len=*), intent(in)  :: norm
        end subroutine
    end interface

    interface
        module recursive subroutine geqrt3( norm, m, n, a, lda, t, ldt, info )
            integer(lpf_default_int_kind), intent(in) :: lda, m, n, ldt
            integer(lpf_default_int_kind), intent(out) :: info
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout) ::  a( lda, * ), t( ldt, * )
        end subroutine
    end interface

    interface
        module subroutine geqrtv( norm, m, n, nb, a, lda, diagr, t, ldt, work, info )
            integer(lpf_default_int_kind), intent(in)  :: lda, ldt, m, n, nb
            integer(lpf_default_int_kind), intent(out) :: info
            type(fp16), intent(inout) :: a( lda, * ), diagr(*), t( ldt, * ), work( * )
            character(len=*), intent(in)  :: norm
        end subroutine
    end interface
    interface
        module subroutine geqrt2v( norm, m, n, a, lda, diagr, t, ldt, info )
            integer(lpf_default_int_kind), intent(in) :: lda, ldt, m, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout) :: a( lda, * ), diagr(*), t( ldt, * )
            character(len=*), intent(in)  :: norm
        end subroutine
    end interface

    interface
        module recursive subroutine geqrt3v( norm, m, n, a, lda, diagr ,t, ldt, info )
            integer(lpf_default_int_kind), intent(in) :: lda, m, n, ldt
            integer(lpf_default_int_kind), intent(out) :: info
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout) ::  a( lda, * ), diagr(*), t( ldt, * )
        end subroutine
    end interface

    !
    ! Tall-Skinny QR
    !
    interface
        module subroutine latsqr( norm, m, n, mb, nb, a, lda, t, ldt, work, lwork, info )
            character(len=*), intent(in)  :: norm
            integer(lpf_default_int_kind), intent(inout) :: info, lwork
            integer(lpf_default_int_kind), intent(in) :: lda, m, n, mb, nb, ldt
            type(fp16), intent(inout) ::  a( lda, * ), work( * ), t( ldt, * )
        end subroutine
    end interface
    interface
        module subroutine tpqrt( norm, m, n, l, nb, a, lda, b, ldb, t, ldt, work, info )
            integer(lpf_default_int_kind), intent(inout) ::  info
            integer(lpf_default_int_kind), intent(in) :: lda, ldb, ldt, n, m, l, nb
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout) :: a( lda, * ), b( ldb, * ), t( ldt, * ), work( * )
        end subroutine
    end interface
    interface
        module subroutine tpqrt2( norm, m, n, l, a, lda, b, ldb, t, ldt, info )
            integer(lpf_default_int_kind), intent(inout) ::   info
            integer(lpf_default_int_kind), intent(in) :: lda, ldb, ldt, n, m, l
            character(len=*), intent(in)  :: norm
            type(fp16), intent(inout) ::  a( lda, * ), b( ldb, * ), t( ldt, * )
        end subroutine
    end interface
    interface
        module subroutine tpmqrt( side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(inout) ::    info
            integer(lpf_default_int_kind), intent(in) :: k, ldv, lda, ldb, m, n, l, nb, ldt
            type(fp16), intent(inout) ::   v( ldv, * ), a( lda, * ), b( ldb, * ), t( ldt, * ), work( * )
        end subroutine
    end interface
    interface
        module subroutine lamtsqr( side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(inout) :: info, lwork
            integer(lpf_default_int_kind), intent(in) :: lda, m, n, k, mb, nb, ldt, ldc
            type(fp16), intent(inout) ::  a( lda, * ), work( * ), c( ldc, * ), t( ldt, * )
        end subroutine
    end interface

    !
    ! Cholesky QR
    !
    interface
        module subroutine gecholqr(m, n, a, lda, r, ldr, info)
            integer(lpf_default_int_kind), intent(in)    :: m, n, lda, ldr
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16) , intent(inout)   :: a(lda, *), r(ldr, *)
        end subroutine
    end interface
    interface
        module subroutine gecholqr_shift(m, n, a, lda, r, ldr, work, info)
            implicit none
            ! arguments
            integer(lpf_default_int_kind), intent(in)   :: m, n, lda, ldr
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout) :: a(lda, *), r(ldr, *), work(*)
        end subroutine
    end interface

    contains
        function lsame(ca,cb) result(out)
        character, intent(in) :: ca,cb
        logical :: out
        intrinsic :: ichar
        integer :: inta,intb,zcode

        out = ca .eq. cb
        if (out) return

        zcode = ichar('z')

        inta = ichar(ca)
        intb = ichar(cb)

        if (zcode.eq.90 .or. zcode.eq.122) then
            if (inta.ge.97 .and. inta.le.122) inta = inta - 32
            if (intb.ge.97 .and. intb.le.122) intb = intb - 32
        else if (zcode.eq.233 .or. zcode.eq.169) then
            if (inta.ge.129 .and. inta.le.137 .or. inta.ge.145 .and. inta.le.153 .or. inta.ge.162 .and. inta.le.169) &
                & inta = inta + 64
            if (intb.ge.129 .and. intb.le.137 .or. intb.ge.145 .and. intb.le.153 .or. intb.ge.162 .and. intb.le.169) &
                & intb = intb + 64
        else if (zcode.eq.218 .or. zcode.eq.250) then
            if (inta.ge.225 .and. inta.le.250) inta = inta - 32
            if (intb.ge.225 .and. intb.le.250) intb = intb - 32
        end if
        out = inta .eq. intb
    end function

    integer(lpf_default_int_kind) function last_column( m, n, a, lda )
        integer(lpf_default_int_kind), intent(in) ::  m, n, lda
        type(fp16), intent(in) :: a( lda, * )


        type(fp16) ::  zero
        integer(lpf_default_int_kind) :: i, ilaslc

        zero = 0.0

        if( n.eq.0 ) then
            ilaslc = n
        else if( a(1, n).ne.zero .or. a(m, n).ne.zero ) then
            ilaslc = n
        else
            do ilaslc = n, 1, -1
                do i = 1, m
                    if( a(i, ilaslc).ne.zero ) return
                end do
            end do
        end if
        last_column = ilaslc
        return
    end function

    integer(lpf_default_int_kind) function last_row( m, n, a, lda )
        integer(lpf_default_int_kind), intent(in)  ::  m, n, lda
        type(fp16), intent(in) :: a( lda, * )

        type(fp16) :: zero
        integer(lpf_default_int_kind) :: i, j, ilaslr

        zero = 0.0


        if( m.eq.0 ) then
            ilaslr = m
        elseif( a(m, 1).ne.zero .or. a(m, n).ne.zero ) then
            ilaslr = m
        else
            ilaslr = 0
            do j = 1, n
                i=m
                do while((a(max(i,1),j).eq.zero).and.(i.ge.1))
                    i=i-1
                enddo
                ilaslr = max( ilaslr, i )
            end do
        end if
        last_row = ilaslr
        return
    end function

    function lamch(cmach) result(out)
        character, intent(in) :: cmach
        type(fp16) :: out
        ! ..
        !
        ! =====================================================================
        !
        ! .. parameters ..
        type(fp16) :: one, zero

        ! ..
        ! .. local scalars ..
        type(fp16) :: rnd, eps, sfmin, small, rmach
        integer :: d
        ! ..
        ! ..
        ! .. intrinsic functions ..
        intrinsic          digits, epsilon, huge, maxexponent, &
            minexponent, radix, tiny

        one = 1.0
        zero = 0.0

        rnd = one
        !
        if( one.eq.rnd ) then
            eps = epsilon(zero) * 0.5
        else
            eps = epsilon(zero)
        end if
        !
        if( lsame( cmach, 'e' ) ) then
            rmach = eps
        else if( lsame( cmach, 's' ) ) then
            sfmin = tiny(zero)
            small = one / huge(zero)
            if( small.ge.sfmin ) then
                !
                !       use small plus a bit, to avoid the possibility of rounding
                !       causing overflow when computing  1/sfmin.
                !
                sfmin = small*( one+eps )
            end if
            rmach = sfmin
        else if( lsame( cmach, 'b' ) ) then
            rmach = radix(zero)
        else if( lsame( cmach, 'p' ) ) then
            rmach = eps * fp16(radix(zero))
        else if( lsame( cmach, 'n' ) ) then
            d = digits(zero)
            rmach = fp16(d)
        else if( lsame( cmach, 'r' ) ) then
            rmach = rnd
        else if( lsame( cmach, 'm' ) ) then
            rmach = minexponent(zero)
        else if( lsame( cmach, 'u' ) ) then
            rmach = tiny(zero)
        else if( lsame( cmach, 'l' ) ) then
            rmach = maxexponent(zero)
        else if( lsame( cmach, 'o' ) ) then
            rmach = huge(zero)
        else
            rmach = zero
        end if
        !
        out = rmach
        return
    end function

    function lapy2( x, y ) result(out)
        type(fp16), intent(in) :: x, y
        type(fp16) :: out

        ! local parameters
        type(fp16) :: zero
        type(fp16) :: one
        ! ..
        ! .. local scalars ..
        type(fp16) :: w, xabs, yabs, z, hugeval
        logical :: x_is_nan, y_is_nan

        ! ..
        ! .. executable statements ..
        !
        zero = 0.0
        one = 1.0

        x_is_nan = isnan( x )
        y_is_nan = isnan( y )
        if ( x_is_nan ) out = x
        if ( y_is_nan ) out = y
        hugeval = huge(x)
        !
        if ( .not.( x_is_nan.or.y_is_nan ) ) then
           xabs = abs( x )
           yabs = abs( y )
           w = max( xabs, yabs )
           z = min( xabs, yabs )
           if( z.eq.zero .or. w.gt.hugeval ) then
              out = w
           else
              out = w*sqrt( one+( z / w )**2 )
           end if
        end if
        return

    end function

    function lapy2_fp32( x, y ) result(out)
        type(fp16), intent(in) :: x, y
        type(fp16) :: out

        ! local parameters
        real(real32), parameter :: zero = 0.0
        real(real32), parameter :: one  = 1.0

        ! ..
        ! .. local scalars ..
        real(real32) :: w, xabs, yabs, z, hugeval
        logical :: x_is_nan, y_is_nan
        real(real32) :: x32, y32

        ! ..
        ! .. executable statements ..
        !
        x32 = real(x)
        y32 = real(y)

        x_is_nan = isnan( x32 )
        y_is_nan = isnan( y32 )
        if ( x_is_nan ) out = fp16(x32)
        if ( y_is_nan ) out = fp16(y32)
        hugeval = huge(x)
        !
        if ( .not.( x_is_nan.or.y_is_nan ) ) then
           xabs = abs( x32 )
           yabs = abs( y32 )
           w = max( xabs, yabs )
           z = min( xabs, yabs )
           if( z.eq.zero .or. w.gt.hugeval ) then
              out = fp16(w)
           else
              out = fp16(w*sqrt( one+( z / w )**2 ))
           end if
        end if
        return

    end function



    function lapy3( x, y, z ) result(out)
        type(fp16), intent(in) :: x, y, z
        type(fp16) :: out
        type(fp16) :: zero
        type(fp16) :: w, xabs, yabs, zabs, hugeval

        zero = 0.0

        hugeval = huge(x)
        xabs = abs( x )
        yabs = abs( y )
        zabs = abs( z )
        w = max( xabs, yabs, zabs )

        if( w.eq.zero .or. w.gt.hugeval ) then
            ! w can be zero for max(0,nan,0)
            ! adding all three entries together will make sure
            ! nan will not disappear.
            out =  xabs + yabs + zabs
        else
            out = w*sqrt( ( xabs / w )**2+( yabs / w )**2+ &
                ( zabs / w )**2 )
        end if
        return
        !
        ! end of hlapy3
        !
    end function lapy3





end module
