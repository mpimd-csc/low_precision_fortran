! SPDX-License-Identifier: LGPL-3.0-or-later
!
! \brief Triangular Banded Matrix-Vector Multiplication (TBMV)
!
! This routine performs the operation:
! x := A * x
! where A is a triangular banded matrix.
!
! \param[in] uplo Character specifying the triangular part of the matrix A to be used: 'U' for upper, 'L' for lower.
! \param[in] trans Character specifying the transpose. 'N' for no transpose, 'T' for transpose.
! \param[in] diag Character specifying whether the diagonal of A is unit: 'U' for unit, 'N' for non-unit.
! \param[in] n Order of matrix A.
! \param[in] k Number of diagonals.
! \param[in] a The triangular banded matrix A.
! \param[in] lda Leading dimension of matrix A.
! \param[in,out] x Vector X.
! \param[in] incx Increment for the elements of x.

#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_tbmv_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_tbmv_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine tbmv_impl(uplo, trans, diag, n, k, a, lda, x, incx)
        implicit none
        character, intent(in) :: uplo, trans, diag
        integer(int64), intent(in) :: incx, k, lda, n
        type(DT), intent(in) :: a(lda,*)
        type(DT), intent(inout) :: x(*)

        type(DT) :: temp, zero
        integer(int64) :: i, ix, j, jx, kplus1, kx, l
        logical :: nounit

        zero = 0.0
        if (n == 0) return

        nounit = (diag == 'N' .or. diag == 'n')

        if (incx <= 0) then
            kx = 1 - (n-1)*incx
        else if (incx /= 1) then
            kx = 1
        end if

        if (trans == 'N' .or. trans == 'n') then
            if (uplo == 'U' .or. uplo == 'u') then
                kplus1 = k + 1
                if (incx == 1) then
                    do j = 1, n
                        if (x(j) /= zero) then
                            temp = x(j)
                            l = kplus1 - j
                            do i = max(1, j-k), j - 1
                                x(i) = x(i) + temp*a(l+i,j)
                            end do
                            if (nounit) x(j) = x(j)*a(kplus1,j)
                        end if
                    end do
                else
                    jx = kx
                    do j = 1, n
                        if (x(jx) /= zero) then
                            temp = x(jx)
                            ix = kx
                            l = kplus1 - j
                            do i = max(1, j-k), j - 1
                                x(ix) = x(ix) + temp*a(l+i,j)
                                ix = ix + incx
                            end do
                            if (nounit) x(jx) = x(jx)*a(kplus1,j)
                        end if
                        jx = jx + incx
                        if (j > k) kx = kx + incx
                    end do
                end if
            else
                if (incx == 1) then
                    do j = n, 1, -1
                        if (x(j) /= zero) then
                            temp = x(j)
                            l = 1 - j
                            do i = min(n, j+k), j + 1, -1
                                x(i) = x(i) + temp*a(l+i,j)
                            end do
                            if (nounit) x(j) = x(j)*a(1,j)
                        end if
                    end do
                else
                    kx = kx + (n-1)*incx
                    jx = kx
                    do j = n, 1, -1
                        if (x(jx) /= zero) then
                            temp = x(jx)
                            ix = kx
                            l = 1 - j
                            do i = min(n, j+k), j + 1, -1
                                x(ix) = x(ix) + temp*a(l+i,j)
                                ix = ix - incx
                            end do
                            if (nounit) x(jx) = x(jx)*a(1,j)
                        end if
                        jx = jx - incx
                        if ((n-j) >= k) kx = kx - incx
                    end do
                end if
            end if
        else
            if (uplo == 'U' .or. uplo == 'u') then
                kplus1 = k + 1
                if (incx == 1) then
                    do j = n, 1, -1
                        temp = x(j)
                        l = kplus1 - j
                        if (nounit) temp = temp*a(kplus1,j)
                        do i = j - 1, max(1, j-k), -1
                            temp = temp + a(l+i,j)*x(i)
                        end do
                        x(j) = temp
                    end do
                else
                    kx = kx + (n-1)*incx
                    jx = kx
                    do j = n, 1, -1
                        temp = x(jx)
                        kx = kx - incx
                        ix = kx
                        l = kplus1 - j
                        if (nounit) temp = temp*a(kplus1,j)
                        do i = j - 1, max(1, j-k), -1
                            temp = temp + a(l+i,j)*x(ix)
                            ix = ix - incx
                        end do
                        x(jx) = temp
                        jx = jx - incx
                    end do
                end if
            else
                if (incx == 1) then
                    do j = 1, n
                        temp = x(j)
                        l = 1 - j
                        if (nounit) temp = temp*a(1,j)
                        do i = j + 1, min(n, j+k)
                            temp = temp + a(l+i,j)*x(i)
                        end do
                        x(j) = temp
                    end do
                else
                    jx = kx
                    do j = 1, n
                        temp = x(jx)
                        kx = kx + incx
                        ix = kx
                        l = 1 - j
                        if (nounit) temp = temp*a(1,j)
                        do i = j + 1, min(n, j+k)
                            temp = temp + a(l+i,j)*x(ix)
                            ix = ix + incx
                        end do
                        x(jx) = temp
                        jx = jx + incx
                    end do
                end if
            end if
        end if
    end subroutine

    module subroutine tbmv_64(uplo, trans, diag, n, k, a, lda, x, incx)
        implicit none
        character, intent(in) :: uplo, trans, diag
        integer(int64), intent(in) :: incx, k, lda, n
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(inout) :: x(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx
        integer(int64) :: total_size
        integer(int64) :: lx(1), la(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (n - 1) * abs(incx))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        call tbmv_impl(uplo, trans, diag, n, k, pa, lda, px, incx)
    end subroutine

    module subroutine tbmv_32(uplo, trans, diag, n, k, a, lda, x, incx)
        implicit none
        character, intent(in) :: uplo, trans, diag
        integer(int32), intent(in) :: incx, k, lda, n
        type(DT), target, intent(in) :: a(..)
        type(DT), target, intent(inout) :: x(..)

        type(DT), CONTIGUOUS, pointer :: px(:)
        type(DT), CONTIGUOUS, pointer :: pa(:)
        integer(int64) :: lenx
        integer(int64) :: total_size
        integer(int64) :: lx(1), la(1)
        type(c_ptr) :: ptr

        lenx = 1 + ( (int(n, int64) - 1) * abs(int(incx, int64)))
        lx(1) = lenx
        ptr = c_loc(x)
        call c_f_pointer(ptr, px, lx)

        total_size = size(a)
        la(1) = total_size
        ptr = c_loc(a)
        call c_f_pointer(ptr, pa, la)

        call tbmv_impl(uplo, trans, diag, int(n, int64), int(k, int64), pa, int(lda, int64), px, int(incx, int64))
    end subroutine

end submodule
