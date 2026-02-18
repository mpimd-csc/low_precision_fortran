module lpf_blas_bf16
    use lpf_bf16
    use iso_fortran_env, only: real32, real64
    use iso_c_binding
    use lpf_types
    implicit none

    !
    ! Level 1
    !

    interface asum
                pure function basum(N,SX,INCX) bind(C, name = "lpf_blas_basum_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16) :: out
        end function
    end interface

    interface axpy
        pure subroutine baxpy(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_baxpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine
    end interface

    interface axpby
        pure subroutine baxpby(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_baxpby_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine
    end interface

    interface copy
        pure subroutine bcopy(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bcopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine
    end interface

    interface dot
        pure function bdot(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bdot_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(*) :: sx, sy
            type(bf16) :: out
        end function
    end interface

    interface nrm2
        pure function bnrm2(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16) :: out
        end function
    end interface

    interface nrm_fp32
        pure function bnrm2_fp32(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fp32_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            type(bf16) :: out
        end function
    end interface

    interface rot
        pure subroutine brot(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_brot_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sc, ss
            type(bf16), intent(inout), dimension(*) :: sx, sy
        end subroutine
    end interface

    interface rotg
        pure subroutine brotg(SA, SB, SC, SS) bind(C, name = "lpf_blas_brotg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(out) :: sc, ss
        end subroutine
    end interface

    interface rotm
        pure subroutine brotm(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_brotm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(*) :: sparam
            type(bf16), intent(inout), dimension(*) :: sx, sy
        end subroutine
    end interface

    interface rotmg
        pure subroutine brotmg(SA, SB, SC, SD, SPARAM) bind(C, name = "lpf_blas_brotmg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            type(bf16), intent(in) :: sa, sb, sc, sd
            type(bf16), intent(out), dimension(*) :: sparam
        end subroutine
    end interface

    interface scal
        pure subroutine bscal(N,SA,SX,INCX) bind(C, name = "lpf_blas_bscal_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(bf16), intent(in) :: sa
            type(bf16), intent(inout), dimension(*) :: sx
        end subroutine
    end interface

    interface dot
        pure function hsbdot(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hsbdot_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sb
            type(bf16), intent(in), dimension(*) :: sx, sy
            type(bf16) :: out
        end function
    end interface

    interface swap
        pure subroutine bswap(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bswap_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(inout), dimension(*) :: sx
            type(bf16), intent(inout), dimension(*) :: sy
        end subroutine
    end interface

    interface iamax
        pure function ibamax(N,SX,INCX) bind(C, name = "lpf_blas_ibamax_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(bf16), intent(in), dimension(*) :: sx
            integer(lpf_default_c_int_kind) :: out
        end function
    end interface

        !
        ! Level 2
        !

    interface gbmv
        subroutine bgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: kl
            integer(lpf_default_c_int_kind), intent(in) :: ku
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine bgbmv
    end interface

    interface gemv
        subroutine bgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgemv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine bgemv
    end interface

    interface ger
        subroutine bger(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_bger_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
            type(bf16), intent(inout) :: a(lda,*)
        end subroutine bger
    end interface

    interface sbmv
        subroutine bsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine bsbmv
    end interface

    interface spmv
        subroutine bspmv(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bspmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: ap(*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine bspmv
    end interface

    interface spr
        subroutine bspr(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_bspr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(inout) :: ap(*)
        end subroutine bspr
    end interface

    interface spr2
        subroutine bspr2(uplo,n,alpha,x,incx,y,incy,ap)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
            type(bf16), intent(inout) :: ap(*)
        end subroutine bspr2
    end interface

    interface symv
        subroutine bsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsymv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine bsymv
    end interface

    interface syr
        subroutine bsyr(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_bsyr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(inout) :: a(lda,*)
        end subroutine bsyr
    end interface

    interface syr2
        subroutine bsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(bf16), intent(in) :: y(*)
            integer(lpf_default_c_int_kind), intent(in) :: incy
            type(bf16), intent(inout) :: a(lda,*)
        end subroutine bsyr2
    end interface

    interface tbmv
        subroutine btbmv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine btbmv
    end interface

    interface tbsv
        subroutine btbsv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine btbsv
    end interface

    interface tpmv
        subroutine btpmv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: ap(*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine btpmv
    end interface

    interface tpsv
        subroutine btpsv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: ap(*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine btpsv
    end interface

    interface trmv
        subroutine btrmv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine btrmv
    end interface

    interface trsv
        subroutine btrsv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: x(*)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine btrsv
    end interface

        !
        ! Level 3
        !
    interface gemm
        subroutine bgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bgemm
    end interface

    interface gemm_fp32
        subroutine bgemm_fp32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemm_fp32_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(ldc,*)
        end subroutine bgemm_fp32
    end interface

    interface symm
        subroutine bsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bsymm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bsymm
    end interface

    interface syrk
        subroutine bsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) bind(C, name = "lpf_blas_bsyrk_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bsyrk
    end interface

    interface syrk2
        subroutine bsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bsyr2k
    end interface

    interface trmm
        subroutine btrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrmm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: b(ldb,*)
        end subroutine btrmm
    end interface

    interface trsm
        subroutine btrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrsm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(inout) :: b(ldb,*)
        end subroutine btrsm
    end interface

    interface gemmtr
        subroutine bgemmtr(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemmtr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(lda,*)
            type(bf16), intent(in) :: b(ldb,*)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(ldc,*)
        end subroutine bgemmtr

    end interface


    ! Auxillary Routines
    interface lacpy
        pure subroutine blacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_blacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(lda, *)
            type(bf16), intent(inout) :: b(ldb, *)
        end subroutine
    end interface

    interface laset
        pure subroutine blaset(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_blaset_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda
            type(bf16), intent(in) :: alpha, beta
            type(bf16), intent(inout) :: a(lda, *)
        end subroutine
    end interface

    interface lacpy
        pure subroutine b2slacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_b2slacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(lda, *)
            real(c_float), intent(in) :: b(ldb, *)
        end subroutine
    end interface

    interface lacpy
        pure subroutine s2blacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2blacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(lda, *)
            type(bf16), intent(in) :: b(ldb, *)
        end subroutine
    end interface

    interface copy
        pure subroutine b2scopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_b2scopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: x(*)
            real(c_float), intent(in) :: y(*)
        end subroutine
    end interface

    interface copy
        pure subroutine s2bcopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2bcopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(*)
            type(bf16), intent(in) :: y(*)
        end subroutine

    end interface

    !
    ! Submodules
    !
    interface scale_diag
        module subroutine scale_diag_bf16(m, n, a, lda, dl, dr, info)
            integer(lpf_default_int_kind), intent(in) :: m, n, lda
            integer(lpf_default_int_kind), intent(inout) :: info
            type(bf16), intent(inout), dimension(lda, *) :: a
            type(bf16), intent(out), dimension(*) :: dl, dr
        end subroutine
    end interface scale_diag

    interface scale_diag_right
        module subroutine scale_diag_right_bf16(m, n, a, lda, dr, info)
            integer(lpf_default_int_kind), intent(in) :: m, n, lda
            integer(lpf_default_int_kind), intent(inout) :: info
            type(bf16), intent(inout), dimension(lda, *) :: a
            type(bf16), intent(out), dimension(*) :: dr
        end subroutine
    end interface scale_diag_right


end module
