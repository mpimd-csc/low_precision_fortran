!> \brief \b SAXPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       REAL SA
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
!>    SAXPY constant times a vector plus a vector.
!>    uses unrolled loops for increments equal to one.
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
!> \param[in] SA
!> \verbatim
!>          SA is REAL
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
!>         storage spacing between elements of SX
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
!
!> \ingroup axpy
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_axpy_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
 submodule (lpf_blas_fp8_e4m3) lpf_blas_axpy_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding


    contains

        pure subroutine axpy_impl(n,sa,sx,incx,sy,incy)
            implicit none

            type(DT) , intent(in) :: sa
            integer(int64), intent(in) :: incx, incy, n
            type(DT), intent(in) :: sx(*)
            type(DT), intent(inout) :: sy(*)

            !     .. local scalars ..
            integer(int64) i,ix,iy
            !     ..
            if (n.le.0) return
            if (sa.eq.0.0) return
            if (incx.eq.1 .and. incy.eq.1) then
                do i = 1,n
                    sy(i) = sy(i) + sa*sx(i)
                end do
            else
                !        code for unequal increments or equal increments
                !          not equal to 1
                !
                ix = 1
                iy = 1
                if (incx.lt.0) ix = (-n+1)*incx + 1
                if (incy.lt.0) iy = (-n+1)*incy + 1
                do i = 1,n
                    sy(iy) = sy(iy) + sa*sx(ix)
                    ix = ix + incx
                    iy = iy + incy
                end do
            end if
            return
        end subroutine

        module subroutine axpy_64(n,sa,sx,incx,sy,incy)
            implicit none

            type(DT) , intent(in) :: sa
            integer(int64), intent(in) :: incx, incy, n
            type(DT), target, intent(in) :: sx(..)
            type(DT), target, intent(inout) :: sy(..)

            type(DT), CONTIGUOUS, pointer, dimension(:) :: psx
            type(DT), CONTIGUOUS, pointer, dimension(:) :: psy

            type(c_ptr) :: ptr
            integer(int64) :: ls(1)
            integer(int64) :: total_size

            ! Flatten the array sx
            total_size = size(sx)
            ls(1) = total_size
            ptr = c_loc(sx)
            call c_f_pointer(ptr, psx, ls)

            total_size = size(sy)
            ls(1) = total_size
            ptr = c_loc(sy)
            call c_f_pointer(ptr, psy, ls)


            call axpy_impl(n, sa, psx, incx, psy, incy)

        end subroutine

        module subroutine axpy_32(n,sa,sx,incx,sy,incy)
            implicit none

            type(DT) , intent(in) :: sa
            integer(int32), intent(in) :: incx, incy, n
            type(DT), target, intent(in) :: sx(..)
            type(DT), target, intent(inout) :: sy(..)

            type(DT), CONTIGUOUS, pointer, DIMENSION(:) :: psx
            type(DT), CONTIGUOUS, pointer, DIMENSION(:) :: psy

            integer(int64) :: cn, cincx, cincy
            integer(int64) :: total_size
            type(c_ptr) :: ptr
            integer(int64) :: ls(1)

            cn = int(n, int64)
            cincx = int(incx, int64)
            cincy = int(incy, int64)

            ! Flatten the array sx
            total_size = size(sx)
            ls(1) = total_size
            ptr = c_loc(sx)
            call c_f_pointer(ptr, psx, ls)

            total_size = size(sy)
            ls(1) = total_size
            ptr = c_loc(sy)
            call c_f_pointer(ptr, psy, ls)


            call axpy_impl(cn, sa, psx, cincx, psy, cincy)

        end subroutine


    end submodule
