! SPDX-License-Identifier: LGPL-3.0-or-later
!> \brief \b SDSDOT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       type(DT) FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       REAL SB
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
!>   Compute the inner product of two vectors with extended
!>   precision accumulation.
!>
!>   Returns S.P. result with dot product accumulated in D.P.
!>   SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
!>   where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!>   defined in a similar way using INCY.
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
!> \param[in] SB
!> \verbatim
!>          SB is REAL
!>          single precision scalar to be added to inner product
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!>          SX is type(DT) array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!>          single precision vector with N elements
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in] SY
!> \verbatim
!>          SY is type(DT) array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!>          single precision vector with N elements
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          storage spacing between elements of SY
!> \endverbatim
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>    REFERENCES
!>
!>    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!>    Krogh, Basic linear algebra subprograms for Fortran
!>    usage, Algorithm No. 539, Transactions on Mathematical
!>    Software 5, 3 (September 1979), pp. 308-323.
!>
!>    REVISION HISTORY  (YYMMDD)
!>
!>    791001  DATE WRITTEN
!>    890531  Changed all specific intrinsics to generic.  (WRB)
!>    890831  Modified array declarations. (WRB)
!>    890831  REVISION DATE from Version 3.2
!>    891214  Prologue converted to Version 4.0 format.  (BAB)
!>    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!>    920501  Reformatted the REFERENCES section.  (WRB)
!>    070118  Reformat to LAPACK coding style
!> \endverbatim
!>
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_sbdot_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_sbdot_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function sbdot_impl(n,sb,sx,incx,sy,incy) result(out)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sb
        type(DT), intent(in) :: sx(*)
        type(DT), intent(in) :: sy(*)
        type(DT) :: out

        real(real32) :: dsdot
        integer(int64) :: i, kx, ky, ns

        dsdot = sb
        if (n.le.0) then
            out = dsdot
            return
        end if

        if (incx.eq.incy .and. incx.gt.0) then
            ns = n*incx
            do i = 1, ns, incx
                dsdot = dsdot + real(sx(i)*sy(i))
            end do
        else
            kx = 1
            ky = 1
            if (incx.lt.0) kx = 1 + (1-n)*incx
            if (incy.lt.0) ky = 1 + (1-n)*incy
            do i = 1, n
                dsdot = dsdot + real(sx(kx)*sy(ky))
                kx = kx + incx
                ky = ky + incy
            end do
        end if
        out = dsdot
    end function

    module function sbdot_64(n,sb,sx,incx,sy,incy) result(out)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sb
        type(DT), target, intent(in) :: sx(..)
        type(DT), target, intent(in) :: sy(..)
        type(DT) :: out

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

        out = sbdot_impl(n, sb, psx, incx, psy, incy)
    end function

    module function sbdot_32(n,sb,sx,incx,sy,incy) result(out)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), intent(in) :: sb
        type(DT), target, intent(in) :: sx(..)
        type(DT), target, intent(in) :: sy(..)
        type(DT) :: out

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

        out = sbdot_impl(int(n, int64), sb, psx, int(incx, int64), psy, int(incy, int64))
    end function

end submodule
