submodule (lpf_blas_bf16) lpf_blas_bf16_scale_diag
    use lpf_bf16
    use iso_fortran_env, only: real32, real64
    implicit none

contains

    module subroutine scale_diag_bf16(m, n, a, lda, dl, dr, info)
        integer, intent(in) :: m, n, lda
        integer, intent(inout) :: info
        type(bf16), intent(inout), dimension(lda, *) :: a
        type(bf16), intent(out), dimension(*) :: dl, dr

        integer :: k

        ! Check arguments
        info = 0
        if ( m .lt. 0 ) then
            info = -1
        else if ( n .lt. 0 ) then
            info = -2
        else if ( lda .lt. max(1, m)) then
            info = -4
        end if

        if ( info .ne. 0) then
            call lpf_blas_xerbla("scale_diag / scale_diag_bf16", info)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if

        do k = 1, m
            dl(k) = bf16(1.0)/maxval(abs(a(k,1:n)))
            a(k,1:n) = a(k,1:n)  * dl(k)
        end do

        do k = 1, n
            dr(k) = bf16(1.0)/maxval(abs(a(1:m,k)))
            a(1:m,k) = a(1:m,k) * dr(k)
        end do
    end subroutine

end submodule

