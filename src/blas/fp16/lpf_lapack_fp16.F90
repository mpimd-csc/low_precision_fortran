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
        module subroutine larf1f( side, m, n, v, incv, tau, c, ldc, work )
            character, intent(in) :: side
            integer(lpf_default_int_kind), intent(in) ::incv, ldc, m, n
            type(fp16), intent(in) :: tau
            type(fp16), intent(inout) ::  c( ldc, * ), work( * )
            type(fp16), intent(in) ::  v( * )
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
        module subroutine orm2r(  side, trans, m, n, k, a, lda, tau, c, ldc, work, info )
            character, intent(in) :: side, trans
            integer(lpf_default_int_kind), intent(in)   :: k, lda, ldc, m, n
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(in)     ::          a( lda, * ), tau( * )
            type(fp16), intent(inout)     ::          c( ldc, * ), work( * )
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
