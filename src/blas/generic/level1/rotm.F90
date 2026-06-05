! SPDX-License-Identifier: LGPL-3.0-or-later
!> \brief \b SROTM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       type(DT) SPARAM(5),SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
!>
!>    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
!>    (SX**T)
!>
!>    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
!>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
!>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
!>
!>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!>    H=(          )    (          )    (          )    (          )
!>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!>    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
!>
!>    IF SFLAG IS NOT ONE OF THE LISTED ABOVE, THE BEHAVIOR IS UNDEFINED
!>    NANS IN SFLAG MAY NOT PROPAGATE TO THE OUTPUT.
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
!>
!> \param[in] SPARAM
!> \verbatim
!>          SPARAM is type(DT) array, dimension (5)
!>     SPARAM(1)=SFLAG
!>     SPARAM(2)=SH11
!>     SPARAM(3)=SH21
!>     SPARAM(4)=SH12
!>     SPARAM(5)=SH22
!> \endverbatim
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_rotm_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rotm_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rotm_impl(n,sx,incx,sy,incy,sparams)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(*)
        type(DT), target, intent(inout) :: sy(*)
        type(DT), intent(in) :: sparams(*)

        type(DT) :: sflag, sh11, sh12, sh21, sh22, two, w, z, zero
        integer(int64) :: i, kx, ky, nsteps

        zero = 0.0
        two = 2.0
        sflag = sparams(1)
        if (n.le.0 .or. (sflag + two .eq. zero)) return

        if (incx.eq.incy .and. incx.gt.0) then
            nsteps = n*incx
            if (sflag.lt.zero) then
                sh11 = sparams(2)
                sh12 = sparams(4)
                sh21 = sparams(3)
                sh22 = sparams(5)
                do i = 1, nsteps, incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z*sh12
                    sy(i) = w*sh21 + z*sh22
                end do
            else if (sflag.eq.zero) then
                sh12 = sparams(4)
                sh21 = sparams(3)
                do i = 1, nsteps, incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w + z*sh12
                    sy(i) = w*sh21 + z
                end do
            else
                sh11 = sparams(2)
                sh22 = sparams(5)
                do i = 1, nsteps, incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z
                    sy(i) = -w + sh22*z
                end do
            end if
        else
            kx = 1
            ky = 1
            if (incx.lt.0) kx = 1 + (1-n)*incx
            if (incy.lt.0) ky = 1 + (1-n)*incy

            if (sflag.lt.zero) then
                sh11 = sparams(2)
                sh12 = sparams(4)
                sh21 = sparams(3)
                sh22 = sparams(5)
                do i = 1, n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w*sh11 + z*sh12
                    sy(ky) = w*sh21 + z*sh22
                    kx = kx + incx
                    ky = ky + incy
                end do
            else if (sflag.eq.zero) then
                sh12 = sparams(4)
                sh21 = sparams(3)
                do i = 1, n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w + z*sh12
                    sy(ky) = w*sh21 + z
                    kx = kx + incx
                    ky = ky + incy
                end do
            else
                sh11 = sparams(2)
                sh22 = sparams(5)
                do i = 1, n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w*sh11 + z
                    sy(ky) = -w + sh22*z
                    kx = kx + incx
                    ky = ky + incy
                end do
            end if
        end if
    end subroutine

    module subroutine rotm_64(n,sx,incx,sy,incy,sparam)
        implicit none
        integer(int64), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)
        type(DT), target, intent(in) :: sparam(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        type(DT), CONTIGUOUS, pointer :: psparams(:)
        integer(int64) :: lenx, leny, lenp
        integer(int64) :: lx(1), ly(1), lp(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        lenp = 5
        lp(1) = lenp
        ptr = c_loc(sparam)
        call c_f_pointer(ptr, psparams, lp)

        call rotm_impl(n, psx, incx, psy, incy, psparams)
    end subroutine

    module subroutine rotm_32(n,sx,incx,sy,incy,sparam)
        implicit none
        integer(int32), intent(in) :: incx, incy, n
        type(DT), target, intent(inout) :: sx(..)
        type(DT), target, intent(inout) :: sy(..)
        type(DT), target, intent(in) :: sparam(..)

        type(DT), CONTIGUOUS, pointer :: psx(:)
        type(DT), CONTIGUOUS, pointer :: psy(:)
        type(DT), CONTIGUOUS, pointer :: psparams(:)
        integer(int32) :: lenx, leny, lenp
        integer(int32) :: lx(1), ly(1), lp(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(sx)
        call c_f_pointer(ptr, psx, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(sy)
        call c_f_pointer(ptr, psy, ly)

        lenp = 5
        lp(1) = lenp
        ptr = c_loc(sparam)
        call c_f_pointer(ptr, psparams, lp)

        call rotm_impl(int(n, int64), psx, int(incx, int64), psy, int(incy, int64), psparams)
    end subroutine

end submodule
