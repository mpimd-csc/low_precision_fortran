module lpf_blas_bf16
    use lpf_bf16
    use iso_fortran_env
    use iso_c_binding
    implicit none
#ifdef LPF_INTEGER8
    integer, parameter, private :: lpf_blas_int_t = c_int64_t
#else
    integer, parameter, private :: lpf_blas_int_t = c_int32_t
#endif
    interface
        !
        ! Level 1
        !
        pure function basum(N,SX,INCX) bind(C, name = "lpf_blas_basum_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16) :: out
        end function

        pure subroutine baxpy(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_baxpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine

        pure subroutine baxpby(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_baxpby_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine

        pure subroutine bcopy(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bcopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine

        pure function bdot(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bdot_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(*) :: sx, sy
            type(bf16) :: out
        end function

        pure function bnrm2(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16) :: out
        end function

        pure function bnrm2_fp32(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fp32_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16) :: out
        end function


        pure subroutine brot(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_brot_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sc, ss
            type(bf16), intent(inout), dimension(*) :: sx, sy
        end subroutine

        pure subroutine brotg(SA, SB, SC, SS) bind(C, name = "lpf_blas_brotg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(out) :: sc, ss
        end subroutine

        pure subroutine brotm(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_brotm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(*) :: sparam
            type(bf16), intent(inout), dimension(*) :: sx, sy
        end subroutine

        pure subroutine brotmg(SA, SB, SC, SD, SPARAM) bind(C, name = "lpf_blas_brotmg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            type(bf16), intent(in) :: sa, sb, sc, sd
            type(bf16), intent(out), dimension(*) :: sparam
        end subroutine

        pure subroutine bscal(N,SA,SX,INCX) bind(C, name = "lpf_blas_bscal_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(bf16), intent(in) :: sa
            type(bf16), intent(inout), dimension(*) :: sx
        end subroutine

        pure function hsbdot(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hsbdot_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sb
            type(bf16), intent(in), dimension(*) :: sx, sy
            type(bf16) :: out
        end function

        pure subroutine bswap(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bswap_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(inout), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine

        pure function ibamax(N,SX,INCX) bind(C, name = "lpf_blas_ibamax_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            integer(lpf_blas_int_t) :: out
        end function

        !
        ! Level 2
        !
        subroutine bgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: kl
            integer(lpf_blas_int_t), intent(in) :: ku
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine bgbmv

        subroutine bgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgemv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine bgemv

        subroutine bger(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_bger_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
            type(bf16), intent(inout) :: a(lda,*)
        end subroutine bger

        subroutine bsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine bsbmv

        subroutine bspmv(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bspmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: ap(*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine bspmv

        subroutine bspr(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_bspr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(inout) :: ap(*)
        end subroutine bspr

        subroutine bspr2(uplo,n,alpha,x,incx,y,incy,ap)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
            type(bf16), intent(inout) :: ap(*)
        end subroutine bspr2

        subroutine bsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsymv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine bsymv

        subroutine bsyr(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_bsyr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(inout) :: a(lda,*)
        end subroutine bsyr

        subroutine bsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(bf16), intent(in) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
            type(bf16), intent(inout) :: a(lda,*)
        end subroutine bsyr2

        subroutine btbmv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine btbmv

        subroutine btbsv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine btbsv

        subroutine btpmv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: ap(*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine btpmv

        subroutine btpsv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: ap(*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine btpsv

        subroutine btrmv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine btrmv

        subroutine btrsv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine btrsv

        !
        ! Level 3
        !
        subroutine bgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bgemm

        subroutine bgemm_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemm_fp32_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            use iso_fortran_env
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(ldc,*)
        end subroutine bgemm_32


        subroutine bsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bsymm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bsymm

        subroutine bsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) bind(C, name = "lpf_blas_bsyrk_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bsyrk

        subroutine bsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bsyr2k

        subroutine btrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrmm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: b(ldb,*)
        end subroutine btrmm

        subroutine btrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrsm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: b(ldb,*)
        end subroutine btrsm

        subroutine bgemmtr(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemmtr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bgemmtr

    end interface


    ! Auxillary Routines
    interface
        pure subroutine blacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_blacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(lda, *)
            type(bf16), intent(inout) :: b(ldb, *)
        end subroutine

        pure subroutine blaset(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_blaset_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda
            type(bf16), intent(in) :: alpha, beta
            type(bf16), intent(inout) :: a(lda, *)
        end subroutine


        pure subroutine b2slacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_b2slacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(lda, *)
            real(c_float), intent(in) :: b(ldb, *)
        end subroutine

        pure subroutine s2blacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2blacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(lda, *)
            type(bf16), intent(in) :: b(ldb, *)
        end subroutine

        pure subroutine b2scopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_b2scopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: x(*)
            real(c_float), intent(in) :: y(*)
        end subroutine

        pure subroutine s2bcopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2bcopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(*)
            type(bf16), intent(in) :: y(*)
        end subroutine

    end interface

end module
