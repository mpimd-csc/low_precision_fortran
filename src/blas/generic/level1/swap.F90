! SPDX-License-Identifier: LGPL-3.0-or-later
!> \brief \b SSWAP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       type(DT) SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SSWAP interchanges two vectors.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] SX
!> \verbatim
!>          SX is type(DT) array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in,out] SY
!> \verbatim
!>          SY is type(DT) array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
!> \endverbatim
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_swap_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_swap_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine swap_impl(n,sx,incx,sy,incy)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(inout) :: sx(*)
        type(DT), intent(inout) :: sy(*)

        type(DT) :: stemp
        integer(int64) :: i, ix, iy, m, mp1

        if (n.le.0) return

        if (incx.eq.1 .and. incy.eq.1) then
            m = mod(n, 3)
            if (m.ne.0) then
                do i = 1, m
                    stemp = sx(i)
                    sx(i) = sy(i)
                    sy(i) = stemp
                end do
                if (n.lt.3) return
            end if
            mp1 = m + 1
            do i = mp1, n, 3
                stemp = sx(i)
                sx(i) = sy(i)
                sy(i) = stemp
                stemp = sx(i+1)
                sx(i+1) = sy(i+1)
                sy(i+1) = stemp
                stemp = sx(i+2)
                sx(i+2) = sy(i+2)
                sy(i+2) = stemp
            end do
        else
            ix = 1
            iy = 1
            if (incx.lt.0) ix = (-n+1)*incx + 1
            if (incy.lt.0) iy = (-n+1)*incy + 1
            do i = 1, n
                stemp = sx(ix)
                sx(ix) = sy(iy)
                sy(iy) = stemp
                ix = ix + incx
                iy = iy + incy
            end do
        end if
    end subroutine

    module subroutine swap_64(n,sx,incx,sy,incy)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
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

        call swap_impl(n, psx, incx, psy, incy)
    end subroutine

    module subroutine swap_32(n,sx,incx,sy,incy)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        integer(int32) :: lenx, leny
        integer(int32) :: lx(1), ly(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        call swap_impl(int(n, int64), psx, int(incx, int64), psy, int(incy, int64))
    end subroutine

end submodule
