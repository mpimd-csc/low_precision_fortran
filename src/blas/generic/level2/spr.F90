!> \brief \b SSPR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER INCX,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       REAL AP(*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>L
!> SSPR    performs the symmetric rank 1 operation
!>
!>    A := alpha*x*x**T + A,
!>
!> where alpha is a real scalar, x is an n element vector and A is an
!> n by n symmetric matrix, supplied in packed form.
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
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
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
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup hpr
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
submodule (lpf_blas_fp8_e5m2) lpf_blas_spr_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_spr_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine spr_impl(uplo, n, alpha, x, incx, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: x(*)
        type(DT), intent(inout) :: ap(*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, j, jx, k, kk, kx

        zero = 0.0
        if ((n.eq.0) .or. (alpha.eq.zero)) return

        if (incx.le.0) then
            kx = 1 - (n-1)*incx
        else if (incx.ne.1) then
            kx = 1
        end if

        kk = 1
        if (uplo == 'U' .or. uplo == 'u') then
            if (incx.eq.1) then
                do j = 1, n
                    if (x(j).ne.zero) then
                        temp = alpha*x(j)
                        k = kk
                        do i = 1, j
                            ap(k) = ap(k) + x(i)*temp
                            k = k + 1
                        end do
                    end if
                    kk = kk + j
                end do
            else
                jx = kx
                do j = 1, n
                    if (x(jx).ne.zero) then
                        temp = alpha*x(jx)
                        ix = kx
                        do k = kk, kk + j - 1
                            ap(k) = ap(k) + x(ix)*temp
                            ix = ix + incx
                        end do
                    end if
                    jx = jx + incx
                    kk = kk + j
                end do
            end if
        else
            if (incx.eq.1) then
                do j = 1, n
                    if (x(j).ne.zero) then
                        temp = alpha*x(j)
                        k = kk
                        do i = j, n
                            ap(k) = ap(k) + x(i)*temp
                            k = k + 1
                        end do
                    end if
                    kk = kk + n - j + 1
                end do
            else
                jx = kx
                do j = 1, n
                    if (x(jx).ne.zero) then
                        temp = alpha*x(jx)
                        ix = jx
                        do k = kk, kk + n - j
                            ap(k) = ap(k) + x(ix)*temp
                            ix = ix + incx
                        end do
                    end if
                    jx = jx + incx
                    kk = kk + n - j + 1
                end do
            end if
        end if
    end subroutine

    module subroutine spr_64(uplo, n, alpha, x, incx, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: ap(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: pap(:)
        integer(int64) :: lenx, lenap
        integer(int64) :: lx(1), lap(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        lenap = (n * (n + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        call spr_impl(uplo, n, alpha, px, incx, pap)
    end subroutine

    module subroutine spr_32(uplo, n, alpha, x, incx, ap)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, n
        type(DT), intent(in) :: alpha
        type(DT), target, intent(in) :: x(..)
        type(DT), target, intent(inout) :: ap(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: pap(:)
        integer(int64) :: lenx, lenap
        integer(int64) :: lx(1), lap(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        lenap = (int(n, int64) * (int(n, int64) + 1)) / 2
        lap(1) = lenap
        ptr = c_loc(ap)
        call c_f_pointer(ptr, pap, lap)

        call spr_impl(uplo, int(n, int64), alpha, px, int(incx, int64), pap)
    end subroutine

end submodule
