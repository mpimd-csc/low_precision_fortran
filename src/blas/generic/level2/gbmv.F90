! SPDX-License-Identifier: LGPL-3.0-or-later
!
! \brief General Banded Matrix-Vector Multiplication (GBMV)
!
! This routine performs the operation:
! y := alpha * A * x + beta * y
! where A is a general banded matrix.
!
! \param[in] trans Character specifying the transpose. 'N' for no transpose, 'T' for transpose.
! \param[in] m Number of rows of matrix A.
! \param[in] n Number of columns of matrix A.
! \param[in] kl Number of lower diagonals.
! \param[in] ku Number of upper diagonals.
! \param[in] alpha Scalar multiplier for the matrix-vector product.
! \param[in] a The banded matrix A, stored in column-major order.
! \param[in] lda Leading dimension of matrix A.
! \param[in] x Vector X.
! \param[in] incx Increment for the elements of x.
! \param[in] beta Scalar multiplier for vector y.
! \param[in,out] y Vector Y.
! \param[in] incy Increment for the elements of y.

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
