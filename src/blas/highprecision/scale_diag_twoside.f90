module lpf_blas_scale_diag
    use iso_fortran_env, only: real32, real64
    use lpf_types
    implicit none

    private

    public :: scale_diag

    interface scale_diag
        procedure :: scale_diag_fp32
        procedure :: scale_diag_fp64
    end interface

contains

    subroutine scale_diag_fp32(m, n, a, lda, dl, dr, info)
        use lpf_xerbla
        integer(lpf_default_int_kind), intent(in) :: m, n, lda
        integer(lpf_default_int_kind), intent(inout) :: info
        real(real32), intent(inout), dimension(lda, *) :: a
        real(real32), intent(out), dimension(*) :: dl, dr

        integer(lpf_default_int_kind) :: k

        external :: lpf_blas_xerbla

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
            call lpf_blas_xerbla("scale_diag / scale_diag_fp32", info)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if
        do k = 1, m
            dl(k) = 1.0/maxval(abs(a(k,1:n)))
            a(k,1:n) = a(k,1:n)  * dl(k)
        end do


        do k = 1, n
            dr(k) = 1.0/maxval(abs(a(1:m,k)))
            a(1:m,k) = a(1:m,k) * dr(k)
        end do

    end subroutine

    subroutine scale_diag_fp64(m, n, a, lda, dl, dr, info)
        integer(lpf_default_int_kind), intent(in) :: m, n, lda
        integer(lpf_default_int_kind), intent(inout) :: info
        real(real64), intent(inout), dimension(lda, *) :: a
        real(real64), intent(out), dimension(*) :: dl, dr

        integer(lpf_default_int_kind) :: k

        external :: lpf_blas_xerbla

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
            call lpf_blas_xerbla("scale_diag / scale_diag_fp64", info)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if
        do k = 1, m
            dl(k) = 1.0D0/maxval(abs(a(k,1:n)))
            a(k,1:n) = a(k,1:n)  * dl(k)
        end do


        do k = 1, n
            dr(k) = 1.0D0/maxval(abs(a(1:m,k)))
            a(1:m,k) = a(1:m,k) * dr(k)
        end do


    end subroutine

end module
