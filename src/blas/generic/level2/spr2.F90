!> \brief \b SSPR2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER INCX,INCY,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       REAL AP(*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPR2   performs the symmetric rank 2 operation
!>
!>    A := alpha*x*y**T + alpha*y*x**T + A,
!>
!> where alpha is a real scalar, x and y are n element vectors and A is
!> an n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
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
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
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
!> \param[in,out] AP
!> \verbatim
!>          AP is REAL array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the upper triangular part of the
!>           updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that and AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the lower triangular part of the
!>           updated matrix.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup hpr
!
!> \par Further Details:
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
!>  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_spr2_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_spr2_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    subroutine spr2_impl(uplo, n, alpha, x, incx, y, incy, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: x(*)
        type(DT), intent(in) :: y(*)
        type(DT), intent(inout) :: ap(*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, iy, j, jx, jy, k, kk, kx, ky

        zero = 0.0
        if ((n.eq.0) .or. (alpha.eq.zero)) return

        if (incx.le.0) then
            kx = 1 - (n-1)*incx
        else if (incx.ne.1) then
            kx = 1
        end if

        if (incy.le.0) then
            ky = 1 - (n-1)*incy
        else if (incy.ne.1) then
            ky = 1
        end if

        kk = 1
        if (uplo == 'U' .or. uplo == 'u') then
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    k = kk;
                    do i = 1, j
                        temp = alpha*(x(i)*y(j) + y(i)*x(j))
                        ap(k) = ap(k) + temp
                        k = k + 1
                    end do
                    kk = kk + j
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    k = kk
                    ix = kx
                    iy = ky
                    do i = 1, j
                        ap(k) = ap(k) + alpha*(x(ix)*y(jy) + y(iy)*x(jx))
                        ix = ix + incx
                        k = k + 1
                    end do
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                end do
            end if
        else
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    k = kk
                    do i = j, n
                        ap(k) = ap(k) + alpha*(x(i)*y(j) + y(i)*x(j))
                        k = k + 1
                    end do
                    kk = kk + n - j + 1
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    k = kk
                    ix = jx
                    iy = jy
                    do i = j, n
                        ap(k) = ap(k) + alpha*(x(ix)*y(jy) + y(iy)*x(jx))
                        ix = ix + incx
                        k = k + 1
                    end do
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + n - j + 1
                end do
            end if
        end if
    end subroutine

    module subroutine spr2_64(uplo, n, alpha, x, incx, y, incy, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: ap(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pap(:)
        integer(int64) :: lenx, leny, lenap
        integer(int64) :: lx(1), ly(1), lap(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        lenap = (n * (n + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        call spr2_impl(uplo, n, alpha, px, incx, py, incy, pap)
    end subroutine

    module subroutine spr2_32(uplo, n, alpha, x, incx, y, incy, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(in) :: y(..)
        type(DT), target, intent(inout) :: ap(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pap(:)
        integer(int64) :: lenx, leny, lenap
        integer(int64) :: lx(1), ly(1), lap(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (int(n, int64) - 1) * abs(int(incy, int64)))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        lenap = (int(n, int64) * (int(n, int64) + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        call spr2_impl(uplo, int(n, int64), alpha, px, int(incx, int64), py, int(incy, int64), pap)
    end subroutine

end submodule
