! SPDX-License-Identifier: LGPL-3.0-or-later
!
! \brief Symmetric Rank-1 Update (SPR)
!
! This routine performs the operation:
! A := A + alpha * x * x^T
! where A is a symmetric matrix in packed storage.
!
! \param[in] uplo Character specifying the part of the matrix A to be used: 'U' for upper, 'L' for lower.
! \param[in] n Order of matrix A.
! \param[in] alpha Scalar multiplier.
! \param[in] x Vector X.
! \param[in] incx Increment for the elements of x.
! \param[in,out] ap The symmetric matrix A in packed storage.

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
