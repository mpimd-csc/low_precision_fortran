!> \brief \b SAXPBY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SAXPBY(N,SA,SX,INCX,SB,SY,INCY)
!
!       .. Scalar Arguments ..
!       REAL SA,SB
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SAXPBY constant times a vector plus constant times a vector.
!>
!>    Y = ALPHA * X + BETA * Y
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] SA
!> \verbatim
!>           SA is REAL
!>           On entry, SA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in] SB
!> \verbatim
!>           SB is REAL
!>           On entry, SB specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!> \author Martin Koehler, MPI Magdeburg
!
!> \ingroup axpby
!
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_axpby_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_axpby_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine axpby_impl(n,sa,sx,incx,sb,sy,incy)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sa, sb
        type(DT), intent(in) :: sx(*)
        type(DT), intent(inout) :: sy(*)

        integer(int64) :: i, ix, iy

        if (n.le.0) return

        if (sa.eq.0.0 .and. sb.ne.0.0) then
            if (incy.eq.1) then
                do i = 1, n
                    sy(i) = sb*sy(i)
                end do
            else
                iy = 1
                if (incy.lt.0) iy = (-n+1)*incy + 1
                do i = 1, n
                    sy(iy) = sb*sy(iy)
                    iy = iy + incy
                end do
            end if
            return
        end if

        if (incx.eq.1 .and. incy.eq.1) then
            do i = 1, n
                sy(i) = sb*sy(i) + sa*sx(i)
            end do
        else
            ix = 1
            iy = 1
            if (incx.lt.0) ix = (-n+1)*incx + 1
            if (incy.lt.0) iy = (-n+1)*incy + 1
            do i = 1, n
                sy(iy) = sb*sy(iy) + sa*sx(ix)
                ix = ix + incx
                iy = iy + incy
            end do
        end if
    end subroutine

    module subroutine axpby_64(n,sa,sx,incx,sb,sy,incy)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sa, sb
        type(DT), target, intent(in) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        call axpby_impl(n, sa, psx, incx, sb, psy, incy)
    end subroutine

    module subroutine axpby_32(n,sa,sx,incx,sb,sy,incy)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sa, sb
        type(DT), target, intent(in) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        call axpby_impl(int(n, int64), sa, psx, int(incx, int64), sb, psy, int(incy, int64))
    end subroutine

end submodule
