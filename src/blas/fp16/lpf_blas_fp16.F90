module lpf_blas_fp16
    use lpf_fp16
    use iso_fortran_env, only : real32, real64
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
        pure function hasum(N,SX,INCX) bind(C, name = "lpf_blas_hasum_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(*) :: sx
            type(fp16) :: out
        end function

        pure subroutine haxpy(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_haxpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa
            type(fp16), intent(in), dimension(*) :: sx
            type(fp16), intent(inout), dimension(*) :: sy
        end subroutine

        pure subroutine haxpby(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_haxpby_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(in), dimension(*) :: sx
            type(fp16), intent(inout), dimension(*) :: sy
        end subroutine

        pure subroutine hcopy(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hcopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(*) :: sx
            type(fp16), intent(inout), dimension(*) :: sy
        end subroutine

        pure function hdot(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hdot_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(*) :: sx, sy
            type(fp16) :: out
        end function

        pure function hnrm2(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(*) :: sx
            type(fp16) :: out
        end function

        pure function hnrm2_fp32(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fp32_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(*) :: sx
            type(fp16) :: out
        end function


        pure subroutine hrot(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_hrot_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sc, ss
            type(fp16), intent(inout), dimension(*) :: sx, sy
        end subroutine

        pure subroutine hrotg(SA, SB, SC, SS) bind(C, name = "lpf_blas_hrotg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(out) :: sc, ss
        end subroutine

        pure subroutine hrotm(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_hrotm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(*) :: sparam
            type(fp16), intent(inout), dimension(*) :: sx, sy
        end subroutine

        pure subroutine hrotmg(SA, SB, SC, SD, SPARAM) bind(C, name = "lpf_blas_hrotmg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            type(fp16), intent(in) :: sa, sb, sc, sd
            type(fp16), intent(out), dimension(*) :: sparam
        end subroutine

        pure subroutine hscal(N,SA,SX,INCX) bind(C, name = "lpf_blas_hscal_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(fp16), intent(in) :: sa
            type(fp16), intent(inout), dimension(*) :: sx
        end subroutine

        pure function hshdot(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hshdot_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sb
            type(fp16), intent(in), dimension(*) :: sx, sy
            type(fp16) :: out
        end function

        pure subroutine hswap(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hswap_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(inout), dimension(*) :: sx
            type(fp16), intent(inout), dimension(*) :: sy
        end subroutine

        pure function ihamax(N,SX,INCX) bind(C, name = "lpf_blas_ihamax_fortran") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(*) :: sx
            integer(lpf_blas_int_t) :: out
        end function

        !
        ! Level 2
        !
        subroutine hgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: kl
            integer(lpf_blas_int_t), intent(in) :: ku
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine hgbmv

        subroutine hgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgemv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine hgemv

        subroutine hger(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hger_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
            type(fp16), intent(inout) :: a(lda,*)
        end subroutine hger

        subroutine hsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine hsbmv

        subroutine hspmv(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hspmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: ap(*)
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine hspmv

        subroutine hspr(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_hspr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(inout) :: ap(*)
        end subroutine hspr

        subroutine hspr2(uplo,n,alpha,x,incx,y,incy,ap)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
            type(fp16), intent(inout) :: ap(*)
        end subroutine hspr2

        subroutine hsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsymv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
        end subroutine hsymv

        subroutine hsyr(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_hsyr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(inout) :: a(lda,*)
        end subroutine hsyr

        subroutine hsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
            type(fp16), intent(in) :: y(*)
            integer(lpf_blas_int_t), intent(in) :: incy
            type(fp16), intent(inout) :: a(lda,*)
        end subroutine hsyr2

        subroutine htbmv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine htbmv

        subroutine htbsv(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine htbsv

        subroutine htpmv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: ap(*)
            type(fp16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine htpmv

        subroutine htpsv(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: ap(*)
            type(fp16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine htpsv

        subroutine htrmv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrmv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine htrmv

        subroutine htrsv(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrsv_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(inout) :: x(*)
            integer(lpf_blas_int_t), intent(in) :: incx
        end subroutine htrsv

        !
        ! Level 3
        !
        subroutine hgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hgemm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: b(ldb,*)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(ldc,*)
        end subroutine hgemm

        subroutine hgemm_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hgemm_fp32_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: b(ldb,*)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(ldc,*)
        end subroutine hgemm_32


        subroutine hsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hsymm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: b(ldb,*)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(ldc,*)
        end subroutine hsymm

        subroutine hsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) bind(C, name = "lpf_blas_hsyrk_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(ldc,*)
        end subroutine hsyrk

        subroutine hsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: trans
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: b(ldb,*)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(ldc,*)
        end subroutine hsyr2k

        subroutine htrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_htrmm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(inout) :: b(ldb,*)
        end subroutine htrmm

        subroutine htrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_htrsm_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: side
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: diag
            integer(lpf_blas_int_t), intent(in) :: m
            integer(lpf_blas_int_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(inout) :: b(ldb,*)
        end subroutine htrsm

        subroutine hgemmtr(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hgemmtr_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: ldc
            integer(lpf_blas_int_t), intent(in) :: ldb
            integer(lpf_blas_int_t), intent(in) :: lda
            character(c_char), dimension(*), intent(in) :: uplo
            character(c_char), dimension(*), intent(in) :: transa
            character(c_char), dimension(*), intent(in) :: transb
            integer(lpf_blas_int_t), intent(in) :: n
            integer(lpf_blas_int_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(lda,*)
            type(fp16), intent(in) :: b(ldb,*)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(ldc,*)
        end subroutine hgemmtr

    end interface


    ! Auxillary Routines
    interface
        pure subroutine hlacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_hlacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(lda, *)
            type(fp16), intent(inout) :: b(ldb, *)
        end subroutine

        pure subroutine hlaset(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_hlaset_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda
            type(fp16), intent(in) :: alpha, beta
            type(fp16), intent(inout) :: a(lda, *)
        end subroutine


        pure subroutine h2slacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_h2slacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(lda, *)
            real(c_float), intent(in) :: b(ldb, *)
        end subroutine

        pure subroutine s2hlacpy(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2hlacpy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            character(c_char), intent(in), dimension(*) :: uplo
            integer(lpf_blas_int_t), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(lda, *)
            type(fp16), intent(in) :: b(ldb, *)
        end subroutine

        pure subroutine h2scopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_h2scopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: x(*)
            real(c_float), intent(in) :: y(*)
        end subroutine

        pure subroutine s2hcopy(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2hcopy_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            import :: lpf_blas_int_t

            integer(lpf_blas_int_t), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(*)
            type(fp16), intent(in) :: y(*)
        end subroutine

    end interface

    !
    ! Submodules
    !
    interface scale_diag
        module subroutine scale_diag_fp16(m, n, a, lda, dl, dr, info)
            integer, intent(in) :: m, n, lda
            integer, intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dl, dr
        end subroutine
    end interface scale_diag



end module
