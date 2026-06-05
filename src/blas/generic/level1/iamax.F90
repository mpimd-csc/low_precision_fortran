! SPDX-License-Identifier: LGPL-3.0-or-later
!> \brief \b ISAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ISAMAX(N,SX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       type(DT) SX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ISAMAX finds the index of the first element having maximum absolut
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
!> \param[in] SX
!> \verbatim
!>          SX is type(DT) array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
!> \endverbatim
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_samax_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_samax_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function iamax_impl(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: sx(*)
        integer(int64) :: out

        type(DT) :: smax
        integer(int64) :: i, ix

        out = 0
        if (n.lt.1 .or. incx.le.0) return
        out = 1
        if (n.eq.1) return

        if (incx.eq.1) then
            smax = abs(sx(1))
            do i = 2, n
                if (abs(sx(i)).gt.smax) then
                    out = i
                    smax = abs(sx(i))
                end if
            end do
        else
            ix = 1
            smax = abs(sx(1))
            ix = ix + incx
            do i = 2, n
                if (abs(sx(ix)).gt.smax) then
                    out = i
                    smax = abs(sx(ix))
                end if
                ix = ix + incx
            end do
        end if
    end function

    module function iamax_64(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        integer(int64) :: out

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int64) :: len
        integer(int64) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        out = iamax_impl(n, sxi, incx)
    end function

    module function iamax_32(n,sx,incx) result(out)
        implicit none
        integer(int32), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        integer(int32) :: out

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int32) :: len
        integer(int32) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        out = int(iamax_impl(int(n, int64), sxi, int(incx, int64)), int32)
    end function

end submodule
