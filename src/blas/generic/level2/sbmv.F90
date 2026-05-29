!> \brief \b SSBMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER INCX,INCY,K,LDA,N
!       CHARACTER UPLO
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
!> SSBMV  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric band matrix, with k super-diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the band matrix A is being supplied as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  being supplied.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  being supplied.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry, K specifies the number of super-diagonals of the
!>           matrix A. K must satisfy  0 .le. K.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension ( LDA, N )
!>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!>           by n part of the array A must contain the upper triangular
!>           band part of the symmetric matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row
!>           ( k + 1 ) of the array, the first super-diagonal starting a
!>           position 2 in row k, and so on. The top left k by k triangl
!>           of the array A is not referenced.
!>           The following program segment will transfer the upper
!>           triangular part of a symmetric band matrix from conventiona
!>           full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = K + 1 - J
!>                    DO 10, I = MAX( 1, J - K ), J
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!>
!>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!>           by n part of the array A must contain the lower triangular
!>           band part of the symmetric matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row 1 of
!>           the array, the first sub-diagonal starting at position 1 in
!>           row 2, and so on. The bottom right k by k triangle of the
!>           array A is not referenced.
!>           The following program segment will transfer the lower
!>           triangular part of a symmetric band matrix from conventiona
!>           full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = 1 - J
!>                    DO 10, I = J, MIN( N, J + K )
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declare
!>           in the calling (sub) program. LDA must be at least
!>           ( k + 1 ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the
!>           vector x.
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
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the
!>           vector y. On exit, Y is overwritten by the updated vector y
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
!> \ingroup hbmv
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
submodule (lpf_blas_fp8_e5m2) lpf_blas_sbmv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_sbmv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine sbmv_impl(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, k, lda, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(in) :: x(*)
        type(DT), intent(inout) :: y(*)

        type(DT) :: temp1, temp2, one, zero
        integer(int64) :: i, ix, iy, j, jx, jy, l, kplus1, kx, ky

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

        if (uplo == 'U' .or. uplo == 'u') then
            kplus1 = k + 1
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    temp1 = alpha*x(j)
                    temp2 = zero
                    l = kplus1 - j
                    do i = max(1, j-k), j - 1
                        y(i) = y(i) + temp1*a(l+i,j)
                        temp2 = temp2 + a(l+i,j)*x(i)
                    end do
                    y(j) = y(j) + temp1*a(kplus1,j) + alpha*temp2
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    temp1 = alpha*x(jx)
                    temp2 = zero
                    ix = kx
                    iy = ky
                    l = kplus1 - j
                    do i = max(1, j-k), j - 1
                        y(iy) = y(iy) + temp1*a(l+i,j)
                        temp2 = temp2 + a(l+i,j)*x(ix)
                        ix = ix + incx
                        iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*a(kplus1,j) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    if (j.gt.k) then
                        kx = kx + incx
                        ky = ky + incy
                    end if
                end do
            end if
        else
            if ((incx.eq.1) .and. (incy.eq.1)) then
                do j = 1, n
                    temp1 = alpha*x(j)
                    temp2 = zero
                    y(j) = y(j) + temp1*a(1,j)
                    l = 1 - j
                    do i = j + 1, min(n, j+k)
                        y(i) = y(i) + temp1*a(l+i,j)
                        temp2 = temp2 + a(l+i,j)*x(i)
                    end do
                    y(j) = y(j) + alpha*temp2
                end do
            else
                jx = kx
                jy = ky
                do j = 1, n
                    temp1 = alpha*x(jx)
                    temp2 = zero
                    y(jy) = y(jy) + temp1*a(1,j)
                    l = 1 - j
                    ix = jx
                    iy = jy
                    do i = j + 1, min(n, j+k)
                        ix = ix + incx
                        iy = iy + incy
                        y(iy) = y(iy) + temp1*a(l+i,j)
                        temp2 = temp2 + a(l+i,j)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                end do
            end if
        end if
    end subroutine

    module subroutine sbmv_64(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, k, lda, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
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

        call sbmv_impl(uplo, n, k, alpha, pa, lda, px, incx, beta, py, incy)
    end subroutine

    module subroutine sbmv_32(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, incy, k, lda, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        integer(int64) :: total_size
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
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

        call sbmv_impl(uplo, int(n, int64), int(k, int64), &
            & alpha, pa, int(lda, int64), px, int(incx, int64), beta, py, int(incy, int64))
    end subroutine

end submodule
