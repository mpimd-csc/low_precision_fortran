submodule(lpf_lapack_fp16) lpf_lapack_gecholqr_fp16

contains

    module subroutine gecholqr(m, n, a, lda, r, ldr, info)
        ! arguments
        integer(lpf_default_int_kind), intent(in)    :: m, n, lda, ldr
        integer(lpf_default_int_kind), intent(inout) :: info
        type(fp16) , intent(inout)   :: a(lda, *), r(ldr, *)

        ! locals
        type(fp16) :: one, zero
        integer(lpf_default_int_kind) :: iinfo

        one = 1.0
        zero = 0.0

        ! check arguments
        info = 0

        if ( m .lt. 0 ) then
            info = -1
        else if ( n .lt. 0 .or. n .gt. m) then
            info = -2
        else if ( lda .lt. max(1, m)) then
            info = -4
        else if ( ldr .lt. max(1, n)) then
            info = -6
        end if

        if ( info .ne. 0 ) then
            call lpf_blas_xerbla("gecholqr", -info)
            return
        end if

        ! quick return
        if ( m.eq.0 .or. n.eq.0 ) then
            return
        endif

        call hsyrk( "upper", "trans", n, m, one, a, lda, zero, r, ldr)
        call potrf("upper", n, r, ldr, iinfo)
        call htrsm("right", "upper", "notrans", "nounit", m, n, one, r, ldr, a, lda )

        if (iinfo .ne. 0 ) then
            info = 1
        end if
        return

    end subroutine
end submodule
