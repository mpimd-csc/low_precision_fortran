module lpf_blas_fp16
    use lpf_fp16
    use iso_fortran_env, only : real32, real64
    use iso_c_binding
    use lpf_types
    implicit none

    interface asum
        !
        ! Level 1
        !
        pure function hasum(N,SX,INCX) bind(C, name = "lpf_blas_hasum_fortran_dyn_rank") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
    end interface

    interface axpy
        pure subroutine haxpy(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_haxpy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface axpby
        pure subroutine haxpby(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_haxpby_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface copy
        pure subroutine hcopy(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hcopy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface dot
        pure function hdot(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hdot_fortran_dyn_rank") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sx, sy
            type(fp16) :: out
        end function
    end interface

    interface nrm2
        pure function hnrm2(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fortran_dyn_rank") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
    end interface

    interface nrm2_fp32
        pure function hnrm2_fp32(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fp32_fortran_dyn_rank") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
    end interface

    interface rot
        pure subroutine hrot(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_hrot_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sc, ss
            type(fp16), intent(inout), dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotg
        pure subroutine hrotg(SA, SB, SC, SS) bind(C, name = "lpf_blas_hrotg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(out) :: sc, ss
        end subroutine
    end interface

    interface rotm
        pure subroutine hrotm(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_hrotm_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sparam
            type(fp16), intent(inout), dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotmg
        pure subroutine hrotmg(SA, SB, SC, SD, SPARAM) bind(C, name = "lpf_blas_hrotmg_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            type(fp16), intent(in) :: sa, sb, sc, sd
            type(fp16), intent(out), dimension(..) :: sparam
        end subroutine
    end interface

    interface scal
        pure subroutine hscal(N,SA,SX,INCX) bind(C, name = "lpf_blas_hscal_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(fp16), intent(in) :: sa
            type(fp16), intent(inout), dimension(..) :: sx
        end subroutine
    end interface

    interface dot
        pure function hshdot(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hshdot_fortran_dyn_rank") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sb
            type(fp16), intent(in), dimension(..) :: sx, sy
            type(fp16) :: out
        end function
    end interface

    interface swap
        pure subroutine hswap(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hswap_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(inout), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface iamax
        pure function ihamax(N,SX,INCX) bind(C, name = "lpf_blas_ihamax_fortran_dyn_rank") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            integer(lpf_default_c_int_kind) :: out
        end function
    end interface

    !
    ! Level 2 - dyn_rank interfaces
    !
    interface gbmv
        subroutine hgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgbmv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: kl
            integer(lpf_default_c_int_kind), intent(in) :: ku
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine hgbmv
    end interface

    interface gemv
        subroutine hgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgemv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine hgemv
    end interface

    interface ger
        subroutine hger(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hger_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
            type(fp16), intent(inout) :: a(..)
        end subroutine hger
    end interface

    interface sbmv
        subroutine hsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsbmv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine hsbmv
    end interface

    interface spmv
        subroutine hspmv(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hspmv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine hspmv
    end interface

    interface spr
        subroutine hspr(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_hspr_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(inout) :: ap(..)
        end subroutine hspr
    end interface

    interface spr2
        subroutine hspr2(uplo,n,alpha,x,incx,y,incy,ap) bind(C, name = "lpf_blas_hspr2_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
            type(fp16), intent(inout) :: ap(..)
        end subroutine hspr2
    end interface

    interface symv
        subroutine hsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsymv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
        end subroutine hsymv
    end interface

    interface syr
        subroutine hsyr(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_hsyr_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(inout) :: a(..)
        end subroutine hsyr
    end interface

    interface syr2
        subroutine hsyr2(uplo,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hsyr2_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(lpf_default_c_int_kind), intent(in) :: incy
            type(fp16), intent(inout) :: a(..)
        end subroutine hsyr2
    end interface

    interface tbmv
        subroutine htbmv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbmv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine htbmv
    end interface

    interface tbsv
        subroutine htbsv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbsv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine htbsv
    end interface

    interface tpmv
        subroutine htpmv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpmv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(inout) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine htpmv
    end interface

    interface tpsv
        subroutine htpsv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpsv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(inout) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine htpsv
    end interface

    interface trmv
        subroutine htrmv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrmv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine htrmv
    end interface

    interface trsv
        subroutine htrsv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrsv_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(lpf_default_c_int_kind), intent(in) :: incx
        end subroutine htrsv
    end interface

    !
    ! Level 3
    !
    interface gemm
        subroutine hgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hgemm_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hgemm

    end interface


    interface gemm_fp32
        subroutine hgemm_fp32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hgemm_fp32_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(..)
        end subroutine hgemm_fp32
    end interface

    interface symm
        subroutine hsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsymm_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsymm
    end interface

    interface syrk
        subroutine hsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsyrk_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsyrk
    end interface

    interface syr2k
        subroutine hsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsyr2k_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsyr2k
    end interface

    interface trmm
        subroutine htrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) &
                bind(C, name = "lpf_blas_htrmm_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine htrmm
    end interface

    interface trsm
        subroutine htrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) &
                bind(C, name = "lpf_blas_htrsm_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: diag
            integer(lpf_default_c_int_kind), intent(in) :: m
            integer(lpf_default_c_int_kind), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine htrsm
    end interface

    interface gemmtr
        subroutine hgemmtr(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hgemmtr_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: ldc
            integer(lpf_default_c_int_kind), intent(in) :: ldb
            integer(lpf_default_c_int_kind), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(lpf_default_c_int_kind), intent(in) :: n
            integer(lpf_default_c_int_kind), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hgemmtr
    end interface


    !
    ! Auxillary Routines
    !
    interface lacpy
        pure subroutine hlacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_hlacpy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine hlacpy

         pure subroutine h2slacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_h2slacpy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            character(kind = c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(..)
            real(c_float), intent(in) :: b(..)
        end subroutine

        pure subroutine s2hlacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2hlacpy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
        end subroutine s2hlacpy
    end interface

    interface laset
        pure subroutine hlaset(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_hlaset_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            character(kind=c_char), intent(in) :: uplo
            integer(lpf_default_c_int_kind), intent(in) :: m, n, lda
            type(fp16), intent(in) :: alpha, beta
            type(fp16), intent(inout) :: a(..)
        end subroutine hlaset
    end interface

    interface copy
        pure subroutine h2scopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_h2scopy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: x(..)
            real(c_float), intent(in) :: y(..)
        end subroutine h2scopy

        pure subroutine s2hcopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2hcopy_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_default_c_int_kind

            integer(lpf_default_c_int_kind), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(..)
            type(fp16), intent(in) :: y(..)
        end subroutine s2hcopy
    end interface

    !
    ! Submodules
    !
    interface scale_diag
        module subroutine scale_diag_fp16(m, n, a, lda, dl, dr, info)
            integer(lpf_default_int_kind), intent(in) :: m, n, lda
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dl, dr
        end subroutine
    end interface scale_diag

    interface scale_diag_right
        module subroutine scale_diag_rightfp16(m, n, a, lda, dr, info)
            integer(lpf_default_int_kind), intent(in) :: m, n, lda
            integer(lpf_default_int_kind), intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dr
        end subroutine
    end interface scale_diag_right


end module
