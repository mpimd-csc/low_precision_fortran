!  SPDX-License-Identifier: LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program; if not, write to the Free Software Foundation,
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

#ifdef LPF_FP8_E4M3
#define DT fp8_e4m3
#define MODNAME lpf_lapack_fp8_e4m3
#define MODNAME_BLAS lpf_blas_fp8_e4m3
#define TYPEMOD lpf_fp8_e4m3
#endif

#ifdef LPF_FP8_E5M2
#define DT fp8_e5m2
#define MODNAME lpf_lapack_fp8_e5m2
#define MODNAME_BLAS lpf_blas_fp8_e5m2
#define TYPEMOD lpf_fp8_e5m2
#endif

#ifdef LPF_FP16
#define DT fp16
#define MODNAME lpf_lapack_fp16
#define MODNAME_BLAS lpf_blas_fp16
#define TYPEMOD lpf_fp16
#endif

#ifdef LPF_BF16
#define DT bf16
#define MODNAME lpf_lapack_bf16
#define MODNAME_BLAS lpf_blas_bf16
#define TYPEMOD lpf_bf16
#endif


!
! We cannot use MODNAME here, since CMake (with GNU Make) does not detect
! the correct module names here. Either we have to enforce the user to switch
! ninja or we have to use be workaround from below.
!
#ifdef LPF_FP8_E4M3
module lpf_lapack_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
module lpf_lapack_fp8_e5m2
#endif
#ifdef LPF_FP16
module lpf_lapack_fp16
#endif
#ifdef LPF_BF16
module lpf_lapack_bf16
#endif
    use TYPEMOD
    use MODNAME_BLAS
    use iso_fortran_env, only: real32, real64, int32, int64
    use iso_c_binding
    use lpf_types
    implicit none

    !
    ! Norms
    !
    interface lange
        module function lange_64( norm, m, n, xa, lda, xwork ) result(nrm)
            character, intent(in) :: norm
            integer(int64), intent(in)   :: lda, m, n
            type(DT), target, intent(in)      :: xa(..)
            type(DT), target, intent(inout)   :: xwork(..)
            type(DT) :: nrm
        end function

        module function lange_32( norm, m, n, xa, lda, xwork ) result(nrm)
            character, intent(in) :: norm
            integer(int32), intent(in)   :: lda, m, n
            type(DT), target, intent(in)      :: xa(..)
            type(DT), target, intent(inout)   :: xwork(..)
            type(DT) :: nrm
        end function
    end interface

    interface lassq
        module subroutine lassq_64( n, xx, incx, sca, sumsq )
            integer(int64), intent(in) ::   incx, n
            real(real32), intent(inout) :: sca, sumsq
            type(DT), target, intent(in)  :: xx(..)
        end subroutine

        module subroutine lassq_32( n, xx, incx, sca, sumsq )
            integer(int32), intent(in) ::   incx, n
            real(real32), intent(inout) :: sca, sumsq
            type(DT), target, intent(in)  :: xx(..)
        end subroutine
    end interface


    !
    ! Cholesky Decompostion
    !
    interface potf2
        module subroutine potf2_64( uplo, n, a_, lda, info )
            character, intent(in) :: uplo
            integer(int64), intent(in)   :: lda, n
            integer(int64), intent(inout) :: info
            type(DT), target, intent(inout) :: a_(..)
        end subroutine
        module subroutine potf2_32( uplo, n, a_, lda, info )
            character, intent(in) :: uplo
            integer(int32), intent(in)   :: lda, n
            integer(int32), intent(inout) :: info
            type(DT), target, intent(inout) :: a_(..)
        end subroutine
    end interface

    interface potrf
        module subroutine potrf_32( uplo, n, a_, lda, info )
         character, intent(in) ::          uplo
         integer(int32), intent(in) :: lda, n
         integer(int32), intent(inout) :: info
         type(DT), target, intent(inout) ::  a_(..)
        end subroutine
        module subroutine potrf_64( uplo, n, a_, lda, info )
         character, intent(in) ::          uplo
         integer(int64), intent(in) :: lda, n
         integer(int64), intent(inout) :: info
         type(DT), target, intent(inout) ::  a_(..)
        end subroutine
    end interface

    interface potrf2
        module subroutine potrf2_32( uplo, n, a_, lda, info )
            character, intent(in) ::          uplo
            integer(int32), intent(in) ::  lda, n
            integer(int32), intent(inout) :: info
            type(DT), target, intent(inout) ::    a_(..)
        end subroutine
        module recursive subroutine potrf2_64( uplo, n, a_, lda, info )
            character, intent(in) ::          uplo
            integer(int64), intent(in) ::  lda, n
            integer(int64), intent(inout) :: info
            type(DT), target, intent(inout) ::    a_(..)
        end subroutine
    end interface

    interface potrfp
        module subroutine potrfp_32(uplo, n, a_, lda, ipiv_, info)
            character, intent(in) :: uplo
            integer(int32), intent(in) :: n, lda
            integer(int32), intent(inout) :: info
            integer(int32), target, intent(inout) :: ipiv_(..)
            type(DT), target, intent(inout) :: a_(..)
        end subroutine
        module subroutine potrfp_64(uplo, n, a_, lda, ipiv_, info)
            character, intent(in) :: uplo
            integer(int64), intent(in) :: n, lda
            integer(int64), intent(inout) :: info
            integer(int64), target, intent(inout) :: ipiv_(..)
            type(DT), target, intent(inout) :: a_(..)
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

    integer(int32) function last_column( m, n, a, lda )
        integer(int32), intent(in) ::  m, n, lda
        type(DT), intent(in) :: a( lda, * )


        type(DT) ::  zero
        integer(int32) :: i, ilaslc

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

    integer(int32) function last_row( m, n, a, lda )
        integer(int32), intent(in)  ::  m, n, lda
        type(DT), intent(in) :: a( lda, * )

        type(DT) :: zero
        integer(int32) :: i, j, ilaslr

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
        type(DT) :: out
        ! ..
        !
        ! =====================================================================
        !
        !..parameters ..
        type(DT) :: one, zero

        ! ..
        !..local scalars ..
        type(DT) :: rnd, eps, sfmin, small, rmach
        integer :: d
        ! ..
        ! ..
        !..intrinsic functions ..
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
            rmach = eps * DT(radix(zero))
        else if( lsame( cmach, 'n' ) ) then
            d = digits(zero)
            rmach = DT(d)
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
        type(DT), intent(in) :: x, y
        type(DT) :: out

        ! local parameters
        type(DT) :: zero
        type(DT) :: one
        ! ..
        !..local scalars ..
        type(DT) :: w, xabs, yabs, z, hugeval
        logical :: x_is_nan, y_is_nan

        ! ..
        !..executable statements ..
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
        type(DT), intent(in) :: x, y
        type(DT) :: out

        ! local parameters
        real(real32), parameter :: zero = 0.0
        real(real32), parameter :: one  = 1.0

        ! ..
        !..local scalars ..
        real(real32) :: w, xabs, yabs, z, hugeval
        logical :: x_is_nan, y_is_nan
        real(real32) :: x32, y32

        ! ..
        !..executable statements ..
        !
        x32 = real(x)
        y32 = real(y)

        x_is_nan = isnan( x32 )
        y_is_nan = isnan( y32 )
        if ( x_is_nan ) out = DT(x32)
        if ( y_is_nan ) out = DT(y32)
        hugeval = huge(x)
        !
        if ( .not.( x_is_nan.or.y_is_nan ) ) then
           xabs = abs( x32 )
           yabs = abs( y32 )
           w = max( xabs, yabs )
           z = min( xabs, yabs )
           if( z.eq.zero .or. w.gt.hugeval ) then
              out = DT(w)
           else
              out = DT(w*sqrt( one+( z / w )**2 ))
           end if
        end if
        return

    end function



    function lapy3( x, y, z ) result(out)
        type(DT), intent(in) :: x, y, z
        type(DT) :: out
        type(DT) :: zero
        type(DT) :: w, xabs, yabs, zabs, hugeval

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
        ! end of blapy3
        !
    end function lapy3

end module
