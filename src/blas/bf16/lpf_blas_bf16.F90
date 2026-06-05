! SPDX-License-Identifier: LGPL-3.0-or-later
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
        pure function basum_64(N,SX,INCX) bind(C, name = "lpf_blas_basum_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16) :: out
        end function
        pure function basum_32(N,SX,INCX) bind(C, name = "lpf_blas_basum_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16) :: out
        end function
    end interface

    interface axpy
        pure subroutine baxpy_64(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_baxpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine baxpy_32(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_baxpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface axpby
        pure subroutine baxpby_64(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_baxpby_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine baxpby_32(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_baxpby_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface copy
        pure subroutine bcopy_64(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bcopy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine bcopy_32(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bcopy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface dot
        pure function bdot_64(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bdot_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(..) :: sx, sy
            type(bf16) :: out
        end function
        pure function bdot_32(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bdot_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(..) :: sx, sy
            type(bf16) :: out
        end function
    end interface

    interface nrm2
        pure function bnrm2_64(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16) :: out
        end function
        pure function bnrm2_32(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16) :: out
        end function
    end interface

    interface nrm2_fp32
        pure function bnrm2_fp32_64(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fp32_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16) :: out
        end function
        pure function bnrm2_fp32_32(N,SX,INCX) bind(C, name = "lpf_blas_bnrm2_fp32_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            type(bf16) :: out
        end function
    end interface

    interface rot
        pure subroutine brot_64(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_brot_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sc, ss
            type(bf16), intent(inout), dimension(..) :: sx, sy
        end subroutine
        pure subroutine brot_32(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_brot_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sc, ss
            type(bf16), intent(inout), dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotg
        pure subroutine brotg(SA, SB, SC, SS) bind(C, name = "lpf_blas_brotg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            type(bf16), intent(in) :: sa, sb
            type(bf16), intent(out) :: sc, ss
        end subroutine
    end interface

    interface rotm
        pure subroutine brotm_64(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_brotm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(..) :: sparam
            type(bf16), intent(inout), dimension(..) :: sx, sy
        end subroutine
        pure subroutine brotm_32(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_brotm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in), dimension(..) :: sparam
            type(bf16), intent(inout), dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotmg
        pure subroutine brotmg(SA, SB, SC, SD, SPARAM) bind(C, name = "lpf_blas_brotmg_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            type(bf16), intent(in) :: sa, sb, sc, sd
            type(bf16), intent(out), dimension(..) :: sparam
        end subroutine
    end interface

    interface scal
        pure subroutine bscal_64(N,SA,SX,INCX) bind(C, name = "lpf_blas_bscal_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx
            type(bf16), intent(in) :: sa
            type(bf16), intent(inout), dimension(..) :: sx
        end subroutine
        pure subroutine bscal_32(N,SA,SX,INCX) bind(C, name = "lpf_blas_bscal_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx
            type(bf16), intent(in) :: sa
            type(bf16), intent(inout), dimension(..) :: sx
        end subroutine
    end interface

    interface dot
        pure function hsbdot_64(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bsbdot_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sb
            type(bf16), intent(in), dimension(..) :: sx, sy
            type(bf16) :: out
        end function
        pure function hsbdot_32(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bsbdot_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: sb
            type(bf16), intent(in), dimension(..) :: sx, sy
            type(bf16) :: out
        end function
    end interface

    interface swap
        pure subroutine bswap_64(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bswap_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(inout), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine bswap_32(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_bswap_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(inout), dimension(..) :: sx
            type(bf16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface iamax
        pure function ibamax_64(N,SX,INCX) bind(C, name = "lpf_blas_ibamax_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            integer(c_int64_t) :: out
        end function
        pure function ibamax_32(N,SX,INCX) bind(C, name = "lpf_blas_ibamax_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx
            type(bf16), intent(in), dimension(..) :: sx
            integer(c_int32_t) :: out
        end function
    end interface

        !
        ! Level 2 - dyn_rank interfaces
        !

    interface gbmv
        subroutine bgbmv_64(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgbmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: kl
            integer(c_int64_t), intent(in) :: ku
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine bgbmv_64
        subroutine bgbmv_32(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgbmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: kl
            integer(c_int32_t), intent(in) :: ku
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine bgbmv_32
    end interface

    interface gemv
        subroutine bgemv_64(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgemv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine bgemv_64
        subroutine bgemv_32(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bgemv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine bgemv_32
    end interface

    interface ger
        subroutine bger_64(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_bger_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: y(..)
            integer(c_int64_t), intent(in) :: incy
            type(bf16), intent(inout) :: a(..)
        end subroutine bger_64
        subroutine bger_32(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_bger_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: y(..)
            integer(c_int32_t), intent(in) :: incy
            type(bf16), intent(inout) :: a(..)
        end subroutine bger_32
    end interface

    interface sbmv
        subroutine bsbmv_64(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsbmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine bsbmv_64
        subroutine bsbmv_32(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsbmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine bsbmv_32
    end interface

    interface spmv
        subroutine bspmv_64(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bspmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: ap(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine bspmv_64
        subroutine bspmv_32(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bspmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: ap(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine bspmv_32
    end interface

    interface spr
        subroutine bspr_64(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_bspr_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(inout) :: ap(..)
        end subroutine bspr_64
        subroutine bspr_32(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_bspr_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(inout) :: ap(..)
        end subroutine bspr_32
    end interface

    interface spr2
        subroutine bspr2_64(uplo,n,alpha,x,incx,y,incy,ap) bind(C, name = "lpf_blas_bspr2_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: y(..)
            integer(c_int64_t), intent(in) :: incy
            type(bf16), intent(inout) :: ap(..)
        end subroutine bspr2_64
        subroutine bspr2_32(uplo,n,alpha,x,incx,y,incy,ap) bind(C, name = "lpf_blas_bspr2_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: y(..)
            integer(c_int32_t), intent(in) :: incy
            type(bf16), intent(inout) :: ap(..)
        end subroutine bspr2_32
    end interface

    interface symv
        subroutine bsymv_64(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsymv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine bsymv_64
        subroutine bsymv_32(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_bsymv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine bsymv_32
    end interface

    interface syr
        subroutine bsyr_64(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_bsyr_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(inout) :: a(..)
        end subroutine bsyr_64
        subroutine bsyr_32(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_bsyr_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(inout) :: a(..)
        end subroutine bsyr_32
    end interface

    interface syr2
        subroutine bsyr2_64(uplo,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_bsyr2_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(bf16), intent(in) :: y(..)
            integer(c_int64_t), intent(in) :: incy
            type(bf16), intent(inout) :: a(..)
        end subroutine bsyr2_64
        subroutine bsyr2_32(uplo,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_bsyr2_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(bf16), intent(in) :: y(..)
            integer(c_int32_t), intent(in) :: incy
            type(bf16), intent(inout) :: a(..)
        end subroutine bsyr2_32
    end interface

    interface tbmv
        subroutine btbmv_64(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine btbmv_64
        subroutine btbmv_32(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine btbmv_32
    end interface

    interface tbsv
        subroutine btbsv_64(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbsv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine btbsv_64
        subroutine btbsv_32(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_btbsv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine btbsv_32
    end interface

    interface tpmv
        subroutine btpmv_64(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: ap(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine btpmv_64
        subroutine btpmv_32(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: ap(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine btpmv_32
    end interface

    interface tpsv
        subroutine btpsv_64(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpsv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: ap(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine btpsv_64
        subroutine btpsv_32(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_btpsv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: ap(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine btpsv_32
    end interface

    interface trmv
        subroutine btrmv_64(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine btrmv_64
        subroutine btrmv_32(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine btrmv_32
    end interface

    interface trsv
        subroutine btrsv_64(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrsv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine btrsv_64
        subroutine btrsv_32(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_btrsv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine btrsv_32
    end interface

        !
        ! Level 3 - dyn_rank interfaces
        !
    interface gemm
        subroutine bgemm_64(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: transb
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bgemm_64
        subroutine bgemm_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bgemm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: transb
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bgemm_32
    end interface

    interface gemm_fp32
        subroutine bgemm_fp32_64(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
            bind(C, name = "lpf_blas_bgemm_fp32_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: transb
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            real(c_float), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(..)
        end subroutine bgemm_fp32_64
        subroutine bgemm_fp32_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
            bind(C, name = "lpf_blas_bgemm_fp32_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: transb
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            real(c_float), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(..)
        end subroutine bgemm_fp32_32
    end interface

    interface symm
        subroutine bsymm_64(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bsymm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: side
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bsymm_64
        subroutine bsymm_32(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bsymm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: side
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bsymm_32
    end interface

    interface syrk
        subroutine bsyrk_64(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) bind(C, name = "lpf_blas_bsyrk_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bsyrk_64
        subroutine bsyrk_32(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) bind(C, name = "lpf_blas_bsyrk_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bsyrk_32
    end interface

    interface syr2k
        subroutine bsyr2k_64(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bsyr2k_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bsyr2k_64
        subroutine bsyr2k_32(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_bsyr2k_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bsyr2k_32
    end interface

    interface trmm
        subroutine btrmm_64(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrmm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: side
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: b(..)
        end subroutine btrmm_64
        subroutine btrmm_32(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrmm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: side
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: b(..)
        end subroutine btrmm_32
    end interface

    interface trsm
        subroutine btrsm_64(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrsm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: side
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: b(..)
        end subroutine btrsm_64
        subroutine btrsm_32(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) bind(C, name = "lpf_blas_btrsm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: side
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: b(..)
        end subroutine btrsm_32
    end interface

    interface gemmtr
        subroutine bgemmtr_64(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                & bind(C, name = "lpf_blas_bgemmtr_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: transb
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bgemmtr_64
        subroutine bgemmtr_32(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                & bind(C, name = "lpf_blas_bgemmtr_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(c_char), intent(in) :: uplo
            character(c_char), intent(in) :: transa
            character(c_char), intent(in) :: transb
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(bf16), intent(in) :: alpha
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
            type(bf16), intent(in) :: beta
            type(bf16), intent(inout) :: c(..)
        end subroutine bgemmtr_32
    end interface


    ! Auxiliary Routines - dyn_rank interfaces
    interface lacpy
        pure subroutine blacpy_64(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_blacpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: b(..)
        end subroutine
        pure subroutine blacpy_32(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_blacpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(..)
            type(bf16), intent(inout) :: b(..)
        end subroutine
    end interface

    interface laset
        pure subroutine blaset_64(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_blaset_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda
            type(bf16), intent(in) :: alpha, beta
            type(bf16), intent(inout) :: a(..)
        end subroutine
        pure subroutine blaset_32(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_blaset_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda
            type(bf16), intent(in) :: alpha, beta
            type(bf16), intent(inout) :: a(..)
        end subroutine
    end interface

    interface lacpy_b2s
        pure subroutine b2slacpy_64(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_b2slacpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(..)
            real(c_float), intent(in) :: b(..)
        end subroutine
        pure subroutine b2slacpy_32(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_b2slacpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda, ldb
            type(bf16), intent(in) :: a(..)
            real(c_float), intent(in) :: b(..)
        end subroutine
    end interface

    interface lacpy_s2b
        pure subroutine s2blacpy_64(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2blacpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
        end subroutine
        pure subroutine s2blacpy_32(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2blacpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            character(c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(..)
            type(bf16), intent(in) :: b(..)
        end subroutine
    end interface

    interface copy_b2s
        pure subroutine b2scopy_64(n, x, incx, y, incy) bind(C, name = "lpf_blas_b2scopy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: x(..)
            real(c_float), intent(in) :: y(..)
        end subroutine
        pure subroutine b2scopy_32(n, x, incx, y, incy) bind(C, name = "lpf_blas_b2scopy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(bf16), intent(in) :: x(..)
            real(c_float), intent(in) :: y(..)
        end subroutine
    end interface

    interface copy_s2b
        pure subroutine s2bcopy_64(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2bcopy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int64_t), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(..)
            type(bf16), intent(in) :: y(..)
        end subroutine
        pure subroutine s2bcopy_32(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2bcopy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_bf16
            integer(c_int32_t), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(..)
            type(bf16), intent(in) :: y(..)
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
