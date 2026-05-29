!> \brief \b SNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    SNRM2 := sqrt( x'*x ).
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
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup nrm2
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_nrm2_fp32_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_nrm2_fp32_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function nrm2_fp32_impl(n,x,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: x(*)
        type(DT) :: out

        integer(int64) :: ix
        real(real32) :: norm , scale, ssq, absxi, one, zero

        one  = 1.0
        zero = 0.0

        if (n.lt.1 .or. incx.lt.1) then
            norm = 0.0
        else if (n.eq.1) then
            norm = abs(x(1))
        else
            scale = 0.0
            ssq = 1.0
            !        the following loop is equivalent to this call to the lapack
            !        auxiliary routine:
            !        call slassq( n, x, incx, scale, ssq )
            !
            do ix = 1,1 + (n-1)*incx,incx
                if (x(ix).ne.zero) then
                    absxi = real(abs(x(ix)))
                    if (scale.lt.absxi) then
                        ssq = one + ssq* (scale/absxi)**2
                        scale = absxi
                    else
                        ssq = ssq + (absxi/scale)**2
                    end if
                end if
            end do
            norm = scale*sqrt(ssq)
        end if

        out = norm
        return
    end function

    module function nrm2_fp32_64(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        type(DT) :: out

        type(DT), CONTIGUOUS, pointer :: px(:)
        integer(int64) :: len
        integer(int64) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, px, lx)

        out = nrm2_fp32_impl(n, px, incx)
    end function

    module function nrm2_fp32_32(n,sx,incx) result(out)
        implicit none
        integer(int32), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        type(DT) :: out

        type(DT), CONTIGUOUS, pointer :: px(:)
        integer(int32) :: len
        integer(int32) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, px, lx)

        out = nrm2_fp32_impl(int(n, int64), px, int(incx, int64))
    end function

end submodule
