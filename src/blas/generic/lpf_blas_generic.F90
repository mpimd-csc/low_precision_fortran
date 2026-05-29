#ifdef LPF_FP8_E4M3
#define DT fp8_e4m3
#define MODNAME lpf_blas_fp8_e4m3
#define TYPEMOD lpf_fp8_e4m3
#endif

#ifdef LPF_FP8_E5M2
#define DT fp8_e5m2
#define MODNAME lpf_blas_fp8_e5m2
#define TYPEMOD lpf_fp8_e5m2
#endif

!
! We cannot use MODNAME here, since CMake (with GNU Make) does not detect
! the correct module names here. Either we have to enforce the user to switch
! ninja or we have to use be workaround from below.
!
#ifdef LPF_FP8_E4M3
module lpf_blas_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
module lpf_blas_fp8_e5m2
#endif
    use TYPEMOD
    use iso_fortran_env, only: real32, real64, int32, int64
    use iso_c_binding
    use lpf_types
    implicit none

    !
    ! Level 1
    !

    interface asum
        module function asum_32(N,SX,INCX) result(out)
            integer(int32), intent(in) :: n, incx
            type(DT), intent(in), dimension(..), target :: sx
            type(DT) :: out
        end function
        module function asum_64(N,SX,INCX) result(out)
            integer(int64), intent(in) :: n, incx
            type(DT), intent(in), dimension(..), target :: sx
            type(DT) :: out
        end function
    end interface

    interface axpy
        module subroutine axpy_32(N,SA,SX,INCX,SY,INCY)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: sa
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
        module subroutine axpy_64(N,SA,SX,INCX,SY,INCY)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: sa
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
    end interface

    interface axpby
        module subroutine axpby_32(N,SA,SX,INCX,SB,SY,INCY)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: sa, sb
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
        module subroutine axpby_64(N,SA,SX,INCX,SB,SY,INCY)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: sa, sb
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
    end interface

    interface copy
        module subroutine copy_32(N,SX,INCX,SY,INCY)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
        module subroutine copy_64(N,SX,INCX,SY,INCY)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
    end interface

    interface dot
        module function dot_32(N,SX,INCX,SY,INCY) result(out)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in), target, dimension(..) :: sx, sy
            type(DT) :: out
        end function
        module function dot_64(N,SX,INCX,SY,INCY) result(out)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in), target, dimension(..) :: sx, sy
            type(DT) :: out
        end function
    end interface

    interface nrm2
        module function nrm2_32(N,SX,INCX) result(out)
            integer(int32), intent(in) :: n, incx
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT) :: out
        end function
        module function nrm2_64(N,SX,INCX) result(out)
            integer(int64), intent(in) :: n, incx
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT) :: out
        end function
    end interface

    interface nrm2_fp32
        module function nrm2_fp32_32(N,SX,INCX) result(out)
            integer(int32), intent(in) :: n, incx
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT) :: out
        end function
        module function nrm2_fp32_64(N,SX,INCX) result(out)
            integer(int64), intent(in) :: n, incx
            type(DT), intent(in), target, dimension(..) :: sx
            type(DT) :: out
        end function
    end interface

    interface rot
        module subroutine rot_32(N,SX,INCX,SY,INCY, C, S)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: c, s
            type(DT), intent(inout), target, dimension(..) :: sx, sy
        end subroutine
        module subroutine rot_64(N,SX,INCX,SY,INCY, C, S)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: c, s
            type(DT), intent(inout), target, dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotg
        module subroutine rotg(SA, SB, SC, SS)
            type(DT), intent(inout) :: sa, sb
            type(DT), intent(out) :: sc, ss
        end subroutine
    end interface

    interface rotm
        module subroutine rotm_32(N,SX,INCX,SY,INCY, SPARAM)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in), target, dimension(..) :: sparam
            type(DT), intent(inout), target, dimension(..) :: sx, sy
        end subroutine
        module subroutine rotm_64(N,SX,INCX,SY,INCY, SPARAM)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in), target, dimension(..) :: sparam
            type(DT), intent(inout), target, dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotmg
        module subroutine rotmg(SA, SB, SC, SD, SPARAM)
            type(DT), intent(inout) :: sa, sb, sc, sd
            type(DT), intent(out), target, dimension(..) :: sparam
        end subroutine
    end interface

    interface scal
        module subroutine scal_32(N,SA,SX,INCX)
            integer(int32), intent(in) :: n, incx
            type(DT), intent(in) :: sa
            type(DT), intent(inout), target, dimension(..) :: sx
        end subroutine
        module subroutine scal_64(N,SA,SX,INCX)
            integer(int64), intent(in) :: n, incx
            type(DT), intent(in) :: sa
            type(DT), intent(inout), target, dimension(..) :: sx
        end subroutine
    end interface

    interface dot
        module function sbdot_32(N,SB,SX,INCX,SY,INCY) result(out)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: sb
            type(DT), intent(in), target, dimension(..) :: sx, sy
            type(DT) :: out
        end function
        module function sbdot_64(N,SB,SX,INCX,SY,INCY) result(out)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: sb
            type(DT), intent(in), target, dimension(..) :: sx, sy
            type(DT) :: out
        end function
    end interface

    interface swap
        module subroutine swap_32(N,SX,INCX,SY,INCY)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(inout), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
        module subroutine swap_64(N,SX,INCX,SY,INCY)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(inout), target, dimension(..) :: sx
            type(DT), intent(inout), target, dimension(..) :: sy
        end subroutine
    end interface

    interface iamax
        module function iamax_32(N,SX,INCX) result(out)
            integer(int32), intent(in) :: n, incx
            type(DT), intent(in), target, dimension(..) :: sx
            integer(int32) :: out
        end function
        module function iamax_64(N,SX,INCX) result(out)
            integer(int64), intent(in) :: n, incx
            type(DT), intent(in), target, dimension(..) :: sx
            integer(int64) :: out
        end function
    end interface

        !
        ! Level 2 - dyn_rank interfaces
        !

    interface gbmv
        module subroutine gbmv_32(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: trans
            integer(int32), intent(in) :: m, n, kl, ku, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
        module subroutine gbmv_64(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: trans
            integer(int64), intent(in) :: m, n, kl, ku, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
    end interface

    interface gemv
        module subroutine gemv_32(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: trans
            integer(int32), intent(in) :: m, n, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
        module subroutine gemv_64(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: trans
            integer(int64), intent(in) :: m, n, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
    end interface

    interface ger
        module subroutine ger_32(m,n,alpha,x,incx,y,incy,a,lda)
            integer(int32), intent(in) :: m, n, lda, incx, incy
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..), y(..)
            type(DT), target, intent(inout) :: a(..)
        end subroutine
        module subroutine ger_64(m,n,alpha,x,incx,y,incy,a,lda)
            integer(int64), intent(in) :: m, n, lda, incx, incy
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..), y(..)
            type(DT), target, intent(inout) :: a(..)
        end subroutine
    end interface

    interface sbmv
        module subroutine sbmv_32(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, k, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
        module subroutine sbmv_64(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, k, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
    end interface

    interface spmv
        module subroutine spmv_32(uplo,n,alpha,ap,x,incx,beta,y,incy)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: ap(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
        module subroutine spmv_64(uplo,n,alpha,ap,x,incx,beta,y,incy)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: ap(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
    end interface

    interface spr
        module subroutine spr_32(uplo,n,alpha,x,incx,ap)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, incx
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..)
            type(DT), target, intent(inout) :: ap(..)
        end subroutine
        module subroutine spr_64(uplo,n,alpha,x,incx,ap)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, incx
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..)
            type(DT), target, intent(inout) :: ap(..)
        end subroutine
    end interface

    interface spr2
        module subroutine spr2_32(uplo,n,alpha,x,incx,y,incy,ap)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..), y(..)
            type(DT), target, intent(inout) :: ap(..)
        end subroutine
        module subroutine spr2_64(uplo,n,alpha,x,incx,y,incy,ap)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..), y(..)
            type(DT), target, intent(inout) :: ap(..)
        end subroutine
    end interface

    interface symv
        module subroutine symv_32(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
        module subroutine symv_64(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, lda, incx, incy
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), x(..)
            type(DT), target, intent(inout) :: y(..)
        end subroutine
    end interface

    interface syr
        module subroutine syr_32(uplo,n,alpha,x,incx,a,lda)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, lda, incx
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..)
            type(DT), target, intent(inout) :: a(..)
        end subroutine
        module subroutine syr_64(uplo,n,alpha,x,incx,a,lda)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, lda, incx
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..)
            type(DT), target, intent(inout) :: a(..)
        end subroutine
    end interface

    interface syr2
        module subroutine syr2_32(uplo,n,alpha,x,incx,y,incy,a,lda)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: n, lda, incx, incy
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..), y(..)
            type(DT), target, intent(inout) :: a(..)
        end subroutine
        module subroutine syr2_64(uplo,n,alpha,x,incx,y,incy,a,lda)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: n, lda, incx, incy
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: x(..), y(..)
            type(DT), target, intent(inout) :: a(..)
        end subroutine
    end interface

    interface tbmv
        module subroutine tbmv_32(uplo,trans,diag,n,k,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, k, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
        module subroutine tbmv_64(uplo,trans,diag,n,k,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int64), intent(in) :: n, k, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
    end interface

    interface tbsv
        module subroutine tbsv_32(uplo,trans,diag,n,k,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, k, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
        module subroutine tbsv_64(uplo,trans,diag,n,k,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int64), intent(in) :: n, k, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
    end interface

    interface tpmv
        module subroutine tpmv_32(uplo,trans,diag,n,ap,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, incx
            type(DT), target, intent(in) :: ap(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
        module subroutine tpmv_64(uplo,trans,diag,n,ap,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int64), intent(in) :: n, incx
            type(DT), target, intent(in) :: ap(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
    end interface

    interface tpsv
        module subroutine tpsv_32(uplo,trans,diag,n,ap,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, incx
            type(DT), target, intent(in) :: ap(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
        module subroutine tpsv_64(uplo,trans,diag,n,ap,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int64), intent(in) :: n, incx
            type(DT), target, intent(in) :: ap(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
    end interface

    interface trmv
        module subroutine trmv_32(uplo,trans,diag,n,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
        module subroutine trmv_64(uplo,trans,diag,n,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int64), intent(in) :: n, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
    end interface

    interface trsv
        module subroutine trsv_32(uplo,trans,diag,n,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int32), intent(in) :: n, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
        module subroutine trsv_64(uplo,trans,diag,n,a,lda,x,incx)
            character(c_char), intent(in) :: uplo, trans, diag
            integer(int64), intent(in) :: n, lda, incx
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: x(..)
        end subroutine
    end interface

        !
        ! Level 3 - dyn_rank interfaces
        !
    interface gemm
        module subroutine gemm_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: transa, transb
            integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
        module subroutine gemm_64(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: transa, transb
            integer(int64), intent(in) :: m, n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
    end interface

    interface gemm_fp32
        module subroutine gemm_fp32_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: transa, transb
            integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: a(..), b(..)
            real(real32), intent(in) :: beta
            real(real32), target, intent(inout) :: c(..)
        end subroutine
        module subroutine gemm_fp32_64(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: transa, transb
            integer(int64), intent(in) :: m, n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: a(..), b(..)
            real(real32), intent(in) :: beta
            real(real32), target, intent(inout) :: c(..)
        end subroutine
    end interface

    interface symm
        module subroutine symm_32(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: side, uplo
            integer(int32), intent(in) :: m, n, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
        module subroutine symm_64(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: side, uplo
            integer(int64), intent(in) :: m, n, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
    end interface

    interface syrk
        module subroutine syrk_32(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
            character(c_char), intent(in) :: uplo, trans
            integer(int32), intent(in) :: n, k, lda, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
        module subroutine syrk_64(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
            character(c_char), intent(in) :: uplo, trans
            integer(int64), intent(in) :: n, k, lda, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
    end interface

    interface syr2k
        module subroutine syr2k_32(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: uplo, trans
            integer(int32), intent(in) :: n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
        module subroutine syr2k_64(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: uplo, trans
            integer(int64), intent(in) :: n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
    end interface

    interface trmm
        module subroutine trmm_32(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
            character(c_char), intent(in) :: side, uplo, transa, diag
            integer(int32), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: b(..)
        end subroutine
        module subroutine trmm_64(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
            character(c_char), intent(in) :: side, uplo, transa, diag
            integer(int64), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: b(..)
        end subroutine
    end interface

    interface trsm
        module subroutine trsm_32(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
            character(c_char), intent(in) :: side, uplo, transa, diag
            integer(int32), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: b(..)
        end subroutine
        module subroutine trsm_64(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
            character(c_char), intent(in) :: side, uplo, transa, diag
            integer(int64), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: alpha
            type(DT), target, intent(in) :: a(..)
            type(DT), target, intent(inout) :: b(..)
        end subroutine
    end interface

    interface gemmtr
        module subroutine gemmtr_32(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: uplo, transa, transb
            integer(int32), intent(in) :: n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
        module subroutine gemmtr_64(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
            character(c_char), intent(in) :: uplo, transa, transb
            integer(int64), intent(in) :: n, k, lda, ldb, ldc
            type(DT), intent(in) :: alpha, beta
            type(DT), target, intent(in) :: a(..), b(..)
            type(DT), target, intent(inout) :: c(..)
        end subroutine
    end interface


    ! Auxiliary Routines - dyn_rank interfaces
    interface lacpy
        module subroutine lacpy_32(uplo, m, n, a, lda, b, ldb)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: a(..)
            type(DT), intent(inout) :: b(..)
        end subroutine
        module subroutine lacpy_64(uplo, m, n, a, lda, b, ldb)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: a(..)
            type(DT), intent(inout) :: b(..)
        end subroutine
    end interface

    interface laset
        module subroutine laset_32(uplo, m, n, alpha, beta, a, lda)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: m, n, lda
            type(DT), intent(in) :: alpha, beta
            type(DT), intent(inout) :: a(..)
        end subroutine
        module subroutine laset_64(uplo, m, n, alpha, beta, a, lda)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: m, n, lda
            type(DT), intent(in) :: alpha, beta
            type(DT), intent(inout) :: a(..)
        end subroutine
    end interface

    interface lacpy
        module subroutine dt2slacpy_32(uplo, m, n, a, lda, b, ldb)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: a(..)
            real(real32), intent(in) :: b(..)
        end subroutine
        module subroutine dt2slacpy_64(uplo, m, n, a, lda, b, ldb)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: m, n, lda, ldb
            type(DT), intent(in) :: a(..)
            real(real32), intent(in) :: b(..)
        end subroutine
    end interface

    interface lacpy
        module subroutine s2dtlacpy_32(uplo, m, n, a, lda, b, ldb)
            character(c_char), intent(in) :: uplo
            integer(int32), intent(in) :: m, n, lda, ldb
            real(real32), intent(in) :: a(..)
            type(DT), intent(in) :: b(..)
        end subroutine
        module subroutine s2dtlacpy_64(uplo, m, n, a, lda, b, ldb)
            character(c_char), intent(in) :: uplo
            integer(int64), intent(in) :: m, n, lda, ldb
            real(real32), intent(in) :: a(..)
            type(DT), intent(in) :: b(..)
        end subroutine
    end interface

    interface copy
        module subroutine dt2scopy_32(n, x, incx, y, incy)
            integer(int32), intent(in) :: n, incx, incy
            type(DT), intent(in) :: x(..)
            real(real32), intent(in) :: y(..)
        end subroutine
        module subroutine dt2scopy_64(n, x, incx, y, incy)
            integer(int64), intent(in) :: n, incx, incy
            type(DT), intent(in) :: x(..)
            real(real32), intent(in) :: y(..)
        end subroutine
    end interface

    interface copy
        module subroutine s2dtcopy_32(n, x, incx, y, incy)
            integer(int32), intent(in) :: n, incx, incy
            real(real32), intent(in) :: x(..)
            type(DT), intent(in) :: y(..)
        end subroutine
        module subroutine s2dtcopy_64(n, x, incx, y, incy)
            integer(int64), intent(in) :: n, incx, incy
            real(real32), intent(in) :: x(..)
            type(DT), intent(in) :: y(..)
        end subroutine
    end interface

    !
    ! Submodules
    !
    interface scale_diag
        module subroutine scale_diag_32(m, n, a, lda, dl, dr, info)
            integer(int32), intent(in) :: m, n, lda
            integer(int32), intent(inout) :: info
            type(DT), intent(inout), dimension(lda, *) :: a
            type(DT), intent(out), dimension(*) :: dl, dr
        end subroutine
        module subroutine scale_diag_64(m, n, a, lda, dl, dr, info)
            integer(int64), intent(in) :: m, n, lda
            integer(int64), intent(inout) :: info
            type(DT), intent(inout), dimension(lda, *) :: a
            type(DT), intent(out), dimension(*) :: dl, dr
        end subroutine
    end interface scale_diag

    interface scale_diag_right
        module subroutine scale_diag_right_32(m, n, a, lda, dr, info)
            integer(int32), intent(in) :: m, n, lda
            integer(int32), intent(inout) :: info
            type(DT), intent(inout), dimension(lda, *) :: a
            type(DT), intent(out), dimension(*) :: dr
        end subroutine
        module subroutine scale_diag_right_64(m, n, a, lda, dr, info)
            integer(int64), intent(in) :: m, n, lda
            integer(int64), intent(inout) :: info
            type(DT), intent(inout), dimension(lda, *) :: a
            type(DT), intent(out), dimension(*) :: dr
        end subroutine
    end interface scale_diag_right


end module
