!> \brief \b SSPMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
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
!> SSPMV  performs the matrix-vector operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric matrix, supplied in packed form.
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
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on.
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
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y. On exit, Y is overwritten by the updated
!>           vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
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
!> \ingroup hpmv
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M
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
submodule (lpf_blas_fp8_e5m2) lpf_blas_spmv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_spmv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine spmv_impl(uplo, n, alpha, ap, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: ap(*)
        type(DT), intent(in) :: x(*)
        type(DT), intent(inout) :: y(*)

        type(DT) :: temp1, temp2, one, zero
        integer(int64) :: i, ix, iy, j, jx, jy, k, kk, kx, ky

        one = 1.0
        zero = 0.0

        if ((n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return

        if (incx.gt.0) then
            kx = 1
        else
            kx = 1 - (n-1)*incx
        end if

        if (incy.gt.0) then
            ky = 1
        else
            ky = 1 - (n-1)*incy
        end if

        if (beta.ne.one) then
            if (incy.eq.1) then
                if (beta.eq.zero) then
                    do i = 1, n
                        y(i) = zero
                    end do
                else
                    do i = 1, n
                        y(i) = beta*y(i)
                    end do
                end if
            else
                iy = ky
                if (beta.eq.zero) then
                    do i = 1, n
                        y(iy) = zero
                        iy = iy + incy
                    end do
                else
                    do i = 1, n
                        y(iy) = beta*y(iy)
                        iy = iy + incy
                    end do
                end if
            end if
        end if

        if (alpha.eq.zero) return
        kk = 1

        if (uplo == 'U' .or. uplo == 'u') then
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    temp1 = alpha*x(j)
                    temp2 = zero
                    k = kk
                    do i = 1, j - 1
                        y(i) = y(i) + temp1*ap(k)
                        temp2 = temp2 + ap(k)*x(i)
                        k = k + 1
                    end do
                    y(j) = y(j) + temp1*ap(kk+j-1) + alpha*temp2
                    kk = kk + j
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    temp1 = alpha*x(jx)
                    temp2 = zero
                    ix = kx
                    iy = ky
                    do k = kk, kk + j - 2
                        y(iy) = y(iy) + temp1*ap(k)
                        temp2 = temp2 + ap(k)*x(ix)
                        ix = ix + incx
                        iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*ap(kk+j-1) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                end do
            end if
        else
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    temp1 = alpha*x(j)
                    temp2 = zero
                    y(j) = y(j) + temp1*ap(kk)
                    k = kk + 1
                    do i = j + 1, n
                        y(i) = y(i) + temp1*ap(k)
                        temp2 = temp2 + ap(k)*x(i)
                        k = k + 1
                    end do
                    y(j) = y(j) + alpha*temp2
                    kk = kk + (n-j+1)
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    temp1 = alpha*x(jx)
                    temp2 = zero
                    y(jy) = y(jy) + temp1*ap(kk)
                    ix = jx
                    iy = jy
                    do k = kk + 1, kk + n - j
                        ix = ix + incx
                        iy = iy + incy
                        y(iy) = y(iy) + temp1*ap(k)
                        temp2 = temp2 + ap(k)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + (n-j+1)
                end do
            end if
        end if
    end subroutine

    module subroutine spmv_64(uplo, n, alpha, ap, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: ap(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: pap(:)
        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        integer(int64) :: lenx, leny, lenap
        integer(int64) :: lx(1), ly(1), lap(1)
        type(c_ptr) :: ptr

        lenap = (n * (n + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (n - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        call spmv_impl(uplo, n, alpha, pap, px, incx, beta, py, incy)
    end subroutine

    module subroutine spmv_32(uplo, n, alpha, ap, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, incy, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: ap(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: pap(:)
        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        integer(int64) :: lenx, leny, lenap
        integer(int64) :: lx(1), ly(1), lap(1)
        type(c_ptr) :: ptr

        lenap = (int(n, int64) * (int(n, int64) + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (int(n, int64) - 1) * abs(int(incy, int64)))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        call spmv_impl(uplo, int(n, int64), alpha, pap, px, int(incx, int64), beta, py, int(incy, int64))
    end subroutine

end submodule
