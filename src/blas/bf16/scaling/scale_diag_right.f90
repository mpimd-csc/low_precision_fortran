! SPDX-License-Identifier: LGPL-3.0-or-later
submodule (lpf_blas_bf16) lpf_blas_bf16_scale_diag_right
    use lpf_bf16
    use iso_fortran_env, only: real32, real64
    use lpf_types
    implicit none

contains

    module subroutine scale_diag_right_bf16(m, n, a, lda, dr, info)
        integer(lpf_default_int_kind), intent(in) :: m, n, lda
        integer(lpf_default_int_kind), intent(inout) :: info
        type(bf16), intent(inout), dimension(lda, *) :: a
        type(bf16), intent(out), dimension(*) :: dr

        integer(lpf_default_int_kind) :: k,l
        type(bf16) :: value, max_value
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
            call lpf_blas_xerbla("scale_diag_right / scale_diag_right_bf16", info)
            return
        end if

        ! Quick Return
        if ( m .eq. 0 .or. n.eq.0 ) then
            return
        end if

        do k = 1, n
            max_value = abs(a(1, k))
            do l = 1, m
                value = abs(a(l,k))
                if ( value .gt. max_value ) then
                    max_value = value
                end if
            end do
            dr(k) = bf16(1.0)/max_value
            a(1:m,k) = a(1:m,k) * dr(k)
        end do
    end subroutine

end submodule

