!> \brief \b SGER
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGER   performs the rank 1 operation
!>
!>    A := alpha*x*y**T + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declare
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
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
!> \ingroup ger
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_ger_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_ger_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine ger_impl(m,n,alpha,x,incx,y,incy,a,lda)
        implicit none
        integer(int64), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: x(*)
        type(DT), intent(in) :: y(*)
        type(DT), intent(inout) :: a(lda,*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, j, jy, kx

        zero = 0.0
        if ((m.eq.0) .or. (n.eq.0) .or. (alpha.eq.zero)) return

        if (incy.gt.0) then
            jy = 1
        else
            jy = 1 - (n-1)*incy
        end if

        if (incx.eq.1) then
            do j = 1, n
                if (y(jy).ne.zero) then
                    temp = alpha*y(jy)
                    do i = 1, m
                        a(i,j) = a(i,j) + x(i)*temp
                    end do
                end if
                jy = jy + incy
            end do
        else
            if (incx.gt.0) then
                kx = 1
            else
                kx = 1 - (m-1)*incx
            end if
            do j = 1, n
                if (y(jy).ne.zero) then
                    temp = alpha*y(jy)
                    ix = kx
                    do i = 1, m
                        a(i,j) = a(i,j) + x(ix)*temp
                        ix = ix + incx
                    end do
                end if
                jy = jy + incy
            end do
        end if
    end subroutine

    module subroutine ger_64(m,n,alpha,x,incx,y,incy,a,lda)
        implicit none
        integer(int64), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: a(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (m - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        call ger_impl(m, n, alpha, px, incx, py, incy, pa, lda)
    end subroutine

    module subroutine ger_32(m,n,alpha,x,incx,y,incy,a,lda)
        implicit none
        integer(int32), intent(in) :: incx, incy, lda, m, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: a(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(m, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (int(n, int64) - 1) * abs(int(incy, int64)))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        call ger_impl(int(m, int64), int(n, int64), alpha, px, int(incx, int64), py, int(incy, int64), pa, int(lda, int64))
    end subroutine

end submodule
