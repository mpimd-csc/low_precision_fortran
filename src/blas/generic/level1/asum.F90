!> \brief \b SASUM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SASUM(N,SX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SASUM takes the sum of the absolute values.
!>    uses unrolled loops for increment equal to one.
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
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
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
!> \ingroup asum
!
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
submodule (lpf_blas_fp8_e5m2) lpf_blas_asum_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_asum_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding


contains


    pure function asum_impl(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: sx(*)
        type(DT) :: out

        type(DT) :: stemp
        integer(int64) :: i,nincx

        out = 0.0
        stemp = 0.0

        if (n.le.0 .or. incx.le.0) return

        if (incx .eq. 1) then
            do i = 1,n
                stemp = stemp + abs(sx(i))
            end do
            out = stemp
            return
        else
            nincx = n*incx
            do i = 1,nincx,incx
                stemp = stemp + abs(sx(i))
            end do
            out = stemp
            return
        end if
    end function

    module function asum_64(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in), target :: sx(..)
        type(DT) :: out

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int64) :: len
        integer(int64) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        out = asum_impl(n, sxi, incx)
    end function


    module function asum_32(n,sx,incx) result(out)
        implicit none
        integer(int32), intent(in) :: incx, n
        type(DT), intent(in), target :: sx(..)
        type(DT) :: out

        type(DT), CONTIGUOUS, pointer :: sxi(:)
        integer(int32) :: len
        integer(int32) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, sxi, lx)

        out = asum_impl(int(n, int64), sxi, int(incx, int64))
    end function

end submodule
