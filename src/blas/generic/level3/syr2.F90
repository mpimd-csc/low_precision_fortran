#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_syr2_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_syr2_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine syr2_impl(uplo, n, alpha, x, incx, y, incy, a, lda)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, lda, n
        type(DT), intent(in) :: alpha
        type(DT), intent(in) :: x(*)
        type(DT), intent(in) :: y(*)
        type(DT), intent(inout) :: a(lda,*)

        type(DT) :: temp1, temp2, zero
        integer(int64) :: i, ix, iy, j, jx, jy, kx, ky

        zero = 0.0
        if ((n == 0) .or. (alpha == zero)) return

        if ((incx /= 1) .or. (incy /= 1)) then
            if (incx > 0) then
                kx = 1
            else
                kx = 1 - (n-1)*incx
            end if
            if (incy > 0) then
                ky = 1
            else
                ky = 1 - (n-1)*incy
            end if
            jx = kx
            jy = ky
        end if

        if (uplo == 'U' .or. uplo == 'u') then
            if ((incx == 1) .and. (incy == 1)) then
                do j = 1, n
                    if ((x(j) /= zero) .or. (y(j) /= zero)) then
                        temp1 = alpha*y(j)
                        temp2 = alpha*x(j)
                        do i = 1, j
                            a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
                        end do
                    end if
                end do
            else
                do j = 1, n
                    if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
                        temp1 = alpha*y(jy)
                        temp2 = alpha*x(jx)
                        ix = kx
                        iy = ky
                        do i = 1, j
                            a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                            ix = ix + incx
                            iy = iy + incy
                        end do
                    end if
                    jx = jx + incx
                    jy = jy + incy
                end do
            end if
        else
            if ((incx == 1) .and. (incy == 1)) then
                do j = 1, n
                    if ((x(j) /= zero) .or. (y(j) /= zero)) then
                        temp1 = alpha*y(j)
                        temp2 = alpha*x(j)
                        do i = j, n
                            a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
                        end do
                    end if
                end do
            else
                do j = 1, n
                    if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
                        temp1 = alpha*y(jy)
                        temp2 = alpha*x(jx)
                        ix = jx
                        iy = jy
                        do i = j, n
                            a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                            ix = ix + incx
                            iy = iy + incy
                        end do
                    end if
                    jx = jx + incx
                    jy = jy + incy
                end do
            end if
        end if
    end subroutine

    module subroutine syr2_64(uplo, n, alpha, x, incx, y, incy, a, lda)
        implicit none
        character, intent(in) :: uplo
        integer(int64), intent(in) :: incx, incy, lda, n
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

        call syr2_impl(uplo, n, alpha, px, incx, py, incy, pa, lda)
    end subroutine

    module subroutine syr2_32(uplo, n, alpha, x, incx, y, incy, a, lda)
        implicit none
        character, intent(in) :: uplo
        integer(int32), intent(in) :: incx, incy, lda, n
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

        call syr2_impl(uplo, int(n, int64), alpha, px, int(incx, int64), py, int(incy, int64), pa, int(lda, int64))
    end subroutine

end submodule
