!> \brief \b SGBMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,              &
!        BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER INCX,INCY,KL,KU,LDA,M,N
!       CHARACTER TRANS
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
!> SGBMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n band matrix, with kl sub-diagonals and ku super-diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!> \endverbatim
!>
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           On entry, KL specifies the number of sub-diagonals of the
!>           matrix A. KL must satisfy  0 .le. KL.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           On entry, KU specifies the number of super-diagonals of the
!>           matrix A. KU must satisfy  0 .le. KU.
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
!>           Before entry, the leading ( kl + ku + 1 ) by n part of the
!>           array A must contain the matrix of coefficients, supplied
!>           column by column, with the leading diagonal of the matrix i
!>           row ( ku + 1 ) of the array, the first super-diagonal
!>           starting at position 2 in row ku, the first sub-diagonal
!>           starting at position 1 in row ( ku + 2 ), and so on.
!>           Elements in the array A that do not correspond to elements
!>           in the band matrix (such as the top left ku by ku triangle)
!>           are not referenced.
!>           The following program segment will transfer a band matrix
!>           from conventional full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    K = KU + 1 - J
!>                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
!>                       A( K + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declare
!>           in the calling (sub) program. LDA must be at least
!>           ( kl + ku + 1 ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
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
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry, the incremented array Y must contain the
!>           vector y. On exit, Y is overwritten by the updated vector y
!>           If either m or n is zero, then Y not referenced and the fun
!>           performs a quick return.
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
!> \ingroup gbmv
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
submodule (lpf_blas_fp8_e5m2) lpf_blas_gbmv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_gbmv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    subroutine gbmv_impl(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: trans
        integer(int64), intent(in) :: incx, incy, kl, ku, lda, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(in) :: x(*)
        type(DT), intent(inout) :: y(*)

        type(DT) :: temp, one, zero
        integer(int64) :: i, ix, iy, j, jx, jy, k, kup1, kx, ky, lenx, leny

        one = 1.0
        zero = 0.0

        if ((m.eq.0) .or. (n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return

        if (trans == 'N' .or. trans == 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        if (incx.gt.0) then
            kx = 1
        else
            kx = 1 - (lenx-1)*incx
        end if

        if (incy.gt.0) then
            ky = 1
        else
            ky = 1 - (leny-1)*incy
        end if

        if (beta.ne.one) then
            if (incy.eq.1) then
                if (beta.eq.zero) then
                    do i = 1,leny
                        y(i) = zero
                    end do
                else
                    do i = 1,leny
                        y(i) = beta*y(i)
                    end do
                end if
            else
                iy = ky
                if (beta.eq.zero) then
                    do i = 1,leny
                        y(iy) = zero
                        iy = iy + incy
                    end do
                else
                    do i = 1,leny
                        y(iy) = beta*y(iy)
                        iy = iy + incy
                    end do
                end if
            end if
        end if
        if (alpha.eq.zero) return
        kup1 = ku + 1

        if (trans .eq. 'n') then
            jx = kx
            if (incy.eq.1) then
                do j = 1,n
                    temp = alpha*x(jx)
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        y(i) = y(i) + temp*a(k+i,j)
                    end do
                    jx = jx + incx
                end do
            else
                do j = 1,n
                    temp = alpha*x(jx)
                    iy = ky
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        y(iy) = y(iy) + temp*a(k+i,j)
                        iy = iy + incy
                    end do
                    jx = jx + incx
                    if (j.gt.ku) ky = ky + incy
                end do
            end if
        else
            jy = ky
            if (incx.eq.1) then
                do j = 1,n
                    temp = zero
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        temp = temp + a(k+i,j)*x(i)
                    end do
                    y(jy) = y(jy) + alpha*temp
                    jy = jy + incy
                end do
            else
                do j = 1,n
                    temp = zero
                    ix = kx
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        temp = temp + a(k+i,j)*x(ix)
                        ix = ix + incx
                    end do

                    y(jy) = y(jy) + alpha*temp
                    jy = jy + incy
                    if (j.gt.ku) kx = kx + incx
                end do
            end if
        end if

    end subroutine

    module subroutine gbmv_64(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: trans
        integer(int64), intent(in) :: incx, incy, kl, ku, lda, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)

        integer(int64) :: total_size
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        type(c_ptr) :: ptr

        if (trans == 'N' .or. trans == 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)


        lenx = 1 + ( (lenx - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (leny - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        call gbmv_impl(trans, m, n, kl, ku, alpha, pa, lda, px, incx, beta, py, incy)
    end subroutine

    module subroutine gbmv_32(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        implicit none
        character, intent(in) :: trans
        integer(int32), intent(in) :: incx, incy, kl, ku, lda, m, n
        type(DT), intent(in) :: alpha, beta
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: y(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: py(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: total_size
        integer(int64) :: lenx, leny
        integer(int64) :: lx(1), ly(1), la(1)
        type(c_ptr) :: ptr

        if (trans == 'N' .or. trans == 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        lenx = 1 + ( (lenx - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        leny = 1 + ( (leny - 1) * abs(incy))
        ly(1) = leny
        ptr = c_loc(y)
        call c_f_pointer(ptr, py, ly)

        call gbmv_impl(trans, int(m, int64), int(n, int64), int(kl, int64), int(ku, int64), &
            & alpha, pa, int(lda, int64), px, int(incx, int64), beta, py, int(incy, int64))
    end subroutine

end submodule
