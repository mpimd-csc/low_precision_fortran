!  SPDX-License-Identifier: LGPL-3.0-or-later
!
!  This file is part of LPF, a Low Precision helper for Fortran
!  Copyright (C) 2025 Martin Koehler
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with this program; if not, write to the Free Software Foundation,
!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

module lpf_blas_fp16
    use lpf_fp16
    use iso_fortran_env, only : real32, real64
    use iso_c_binding
    use lpf_types
    implicit none

    interface asum
        pure function hasum_64(N,SX,INCX) bind(C, name = "lpf_blas_hasum_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
        pure function hasum_32(N,SX,INCX) bind(C, name = "lpf_blas_hasum_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
    end interface

    interface axpy
        pure subroutine haxpy_64(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_haxpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine haxpy_32(N,SA,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_haxpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface axpby
        pure subroutine haxpby_64(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_haxpby_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine haxpby_32(N,SA,SX,INCX,SB,SY,INCY) bind(C, name = "lpf_blas_haxpby_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface copy
        pure subroutine hcopy_64(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hcopy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine hcopy_32(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hcopy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface dot
        pure function hdot_64(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hdot_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sx, sy
            type(fp16) :: out
        end function
        pure function hdot_32(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hdot_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sx, sy
            type(fp16) :: out
        end function
    end interface

    interface nrm2
        pure function hnrm2_64(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
        pure function hnrm2_32(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
    end interface

    interface nrm2_fp32
        pure function hnrm2_fp32_64(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fp32_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
        pure function hnrm2_fp32_32(N,SX,INCX) bind(C, name = "lpf_blas_hnrm2_fp32_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            type(fp16) :: out
        end function
    end interface

    interface rot
        pure subroutine hrot_64(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_hrot_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sc, ss
            type(fp16), intent(inout), dimension(..) :: sx, sy
        end subroutine
        pure subroutine hrot_32(N,SX,INCX,SY,INCY, SC, SS) bind(C, name = "lpf_blas_hrot_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sc, ss
            type(fp16), intent(inout), dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotg
        pure subroutine hrotg(SA, SB, SC, SS) bind(C, name = "lpf_blas_hrotg_fortran")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            type(fp16), intent(in) :: sa, sb
            type(fp16), intent(out) :: sc, ss
        end subroutine
    end interface

    interface rotm
        pure subroutine hrotm_64(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_hrotm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sparam
            type(fp16), intent(inout), dimension(..) :: sx, sy
        end subroutine
        pure subroutine hrotm_32(N,SX,INCX,SY,INCY, SPARAM) bind(C, name = "lpf_blas_hrotm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in), dimension(..) :: sparam
            type(fp16), intent(inout), dimension(..) :: sx, sy
        end subroutine
    end interface

    interface rotmg
        pure subroutine hrotmg(SA, SB, SC, SD, SPARAM) bind(C, name = "lpf_blas_hrotmg_fortran_dyn_rank")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            type(fp16), intent(in) :: sa, sb, sc, sd
            type(fp16), intent(out), dimension(..) :: sparam
        end subroutine
    end interface

    interface scal
        pure subroutine hscal_64(N,SA,SX,INCX) bind(C, name = "lpf_blas_hscal_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx
            type(fp16), intent(in) :: sa
            type(fp16), intent(inout), dimension(..) :: sx
        end subroutine
        pure subroutine hscal_32(N,SA,SX,INCX) bind(C, name = "lpf_blas_hscal_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx
            type(fp16), intent(in) :: sa
            type(fp16), intent(inout), dimension(..) :: sx
        end subroutine
    end interface

    interface dot
        pure function hshdot_64(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hshdot_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sb
            type(fp16), intent(in), dimension(..) :: sx, sy
            type(fp16) :: out
        end function
        pure function hshdot_32(N,SB,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hshdot_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: sb
            type(fp16), intent(in), dimension(..) :: sx, sy
            type(fp16) :: out
        end function
    end interface

    interface swap
        pure subroutine hswap_64(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hswap_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(inout), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
        pure subroutine hswap_32(N,SX,INCX,SY,INCY) bind(C, name = "lpf_blas_hswap_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(inout), dimension(..) :: sx
            type(fp16), intent(inout), dimension(..) :: sy
        end subroutine
    end interface

    interface iamax
        pure function ihamax_64(N,SX,INCX) bind(C, name = "lpf_blas_ihamax_fortran_dyn_rank_64") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            integer(c_int64_t) :: out
        end function
        pure function ihamax_32(N,SX,INCX) bind(C, name = "lpf_blas_ihamax_fortran_dyn_rank_32") result(out)
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: n, incx
            type(fp16), intent(in), dimension(..) :: sx
            integer(c_int32_t) :: out
        end function
    end interface

    !
    ! Level 2 - dyn_rank interfaces
    !
    interface gbmv
        subroutine hgbmv_64(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgbmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: kl
            integer(c_int64_t), intent(in) :: ku
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine hgbmv_64
        subroutine hgbmv_32(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgbmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: kl
            integer(c_int32_t), intent(in) :: ku
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine hgbmv_32
    end interface

    interface gemv
        subroutine hgemv_64(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgemv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine hgemv_64
        subroutine hgemv_32(trans,m,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hgemv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine hgemv_32
    end interface

    interface ger
        subroutine hger_64(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hger_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(c_int64_t), intent(in) :: incy
            type(fp16), intent(inout) :: a(..)
        end subroutine hger_64
        subroutine hger_32(m,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hger_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(c_int32_t), intent(in) :: incy
            type(fp16), intent(inout) :: a(..)
        end subroutine hger_32
    end interface

    interface sbmv
        subroutine hsbmv_64(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsbmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine hsbmv_64
        subroutine hsbmv_32(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsbmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine hsbmv_32
    end interface

    interface spmv
        subroutine hspmv_64(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hspmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine hspmv_64
        subroutine hspmv_32(uplo,n,alpha,ap,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hspmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine hspmv_32
    end interface

    interface spr
        subroutine hspr_64(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_hspr_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(inout) :: ap(..)
        end subroutine hspr_64
        subroutine hspr_32(uplo,n,alpha,x,incx,ap) bind(C, name = "lpf_blas_hspr_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(inout) :: ap(..)
        end subroutine hspr_32
    end interface

    interface spr2
        subroutine hspr2_64(uplo,n,alpha,x,incx,y,incy,ap) bind(C, name = "lpf_blas_hspr2_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(c_int64_t), intent(in) :: incy
            type(fp16), intent(inout) :: ap(..)
        end subroutine hspr2_64
        subroutine hspr2_32(uplo,n,alpha,x,incx,y,incy,ap) bind(C, name = "lpf_blas_hspr2_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(c_int32_t), intent(in) :: incy
            type(fp16), intent(inout) :: ap(..)
        end subroutine hspr2_32
    end interface

    interface symv
        subroutine hsymv_64(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsymv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int64_t), intent(in) :: incy
        end subroutine hsymv_64
        subroutine hsymv_32(uplo,n,alpha,a,lda,x,incx,beta,y,incy) bind(C, name = "lpf_blas_hsymv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: y(..)
            integer(c_int32_t), intent(in) :: incy
        end subroutine hsymv_32
    end interface

    interface syr
        subroutine hsyr_64(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_hsyr_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(inout) :: a(..)
        end subroutine hsyr_64
        subroutine hsyr_32(uplo,n,alpha,x,incx,a,lda) bind(C, name = "lpf_blas_hsyr_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(inout) :: a(..)
        end subroutine hsyr_32
    end interface

    interface syr2
        subroutine hsyr2_64(uplo,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hsyr2_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int64_t), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(c_int64_t), intent(in) :: incy
            type(fp16), intent(inout) :: a(..)
        end subroutine hsyr2_64
        subroutine hsyr2_32(uplo,n,alpha,x,incx,y,incy,a,lda) bind(C, name = "lpf_blas_hsyr2_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: x(..)
            integer(c_int32_t), intent(in) :: incx
            type(fp16), intent(in) :: y(..)
            integer(c_int32_t), intent(in) :: incy
            type(fp16), intent(inout) :: a(..)
        end subroutine hsyr2_32
    end interface

    interface tbmv
        subroutine htbmv_64(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine htbmv_64
        subroutine htbmv_32(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine htbmv_32
    end interface

    interface tbsv
        subroutine htbsv_64(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbsv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine htbsv_64
        subroutine htbsv_32(uplo,trans,diag,n,k,a,lda,x,incx) bind(C, name = "lpf_blas_htbsv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine htbsv_32
    end interface

    interface tpmv
        subroutine htpmv_64(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine htpmv_64
        subroutine htpmv_32(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine htpmv_32
    end interface

    interface tpsv
        subroutine htpsv_64(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpsv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine htpsv_64
        subroutine htpsv_32(uplo,trans,diag,n,ap,x,incx) bind(C, name = "lpf_blas_htpsv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: ap(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine htpsv_32
    end interface

    interface trmv
        subroutine htrmv_64(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrmv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine htrmv_64
        subroutine htrmv_32(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrmv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine htrmv_32
    end interface

    interface trsv
        subroutine htrsv_64(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrsv_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int64_t), intent(in) :: incx
        end subroutine htrsv_64
        subroutine htrsv_32(uplo,trans,diag,n,a,lda,x,incx) bind(C, name = "lpf_blas_htrsv_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: x(..)
            integer(c_int32_t), intent(in) :: incx
        end subroutine htrsv_32
    end interface

    !
    ! Level 3
    !
    interface gemm
        subroutine hgemm_64(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hgemm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hgemm_64
        subroutine hgemm_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(C, name = "lpf_blas_hgemm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hgemm_32
    end interface

    interface gemm_fp32
        subroutine hgemm_fp32_64(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hgemm_fp32_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            real(c_float), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(..)
        end subroutine hgemm_fp32_64
        subroutine hgemm_fp32_32(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hgemm_fp32_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            real(c_float), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            real(c_float), intent(in) :: beta
            real(c_float), intent(inout) :: c(..)
        end subroutine hgemm_fp32_32
    end interface

    interface symm
        subroutine hsymm_64(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsymm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsymm_64
        subroutine hsymm_32(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsymm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsymm_32
    end interface

    interface syrk
        subroutine hsyrk_64(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsyrk_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsyrk_64
        subroutine hsyrk_32(uplo,trans,n,k,alpha,a,lda,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsyrk_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsyrk_32
    end interface

    interface syr2k
        subroutine hsyr2k_64(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsyr2k_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsyr2k_64
        subroutine hsyr2k_32(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hsyr2k_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: trans
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hsyr2k_32
    end interface

    interface trmm
        subroutine htrmm_64(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) &
                bind(C, name = "lpf_blas_htrmm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine htrmm_64
        subroutine htrmm_32(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) &
                bind(C, name = "lpf_blas_htrmm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine htrmm_32
    end interface

    interface trsm
        subroutine htrsm_64(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) &
                bind(C, name = "lpf_blas_htrsm_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: diag
            integer(c_int64_t), intent(in) :: m
            integer(c_int64_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine htrsm_64
        subroutine htrsm_32(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb) &
                bind(C, name = "lpf_blas_htrsm_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: side
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: diag
            integer(c_int32_t), intent(in) :: m
            integer(c_int32_t), intent(in) :: n
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine htrsm_32
    end interface

    interface gemmtr
        subroutine hgemmtr_64(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hgemmtr_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: ldc
            integer(c_int64_t), intent(in) :: ldb
            integer(c_int64_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(c_int64_t), intent(in) :: n
            integer(c_int64_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hgemmtr_64
        subroutine hgemmtr_32(uplo,transa,transb,n,k,alpha,a,lda,b,ldb,beta,c,ldc) &
                bind(C, name = "lpf_blas_hgemmtr_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: ldc
            integer(c_int32_t), intent(in) :: ldb
            integer(c_int32_t), intent(in) :: lda
            character(kind=c_char), intent(in) :: uplo
            character(kind=c_char), intent(in) :: transa
            character(kind=c_char), intent(in) :: transb
            integer(c_int32_t), intent(in) :: n
            integer(c_int32_t), intent(in) :: k
            type(fp16), intent(in) :: alpha
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
            type(fp16), intent(in) :: beta
            type(fp16), intent(inout) :: c(..)
        end subroutine hgemmtr_32
    end interface


    !
    ! Auxillary Routines
    !
    interface lacpy
        pure subroutine hlacpy_64(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_hlacpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine hlacpy_64
        pure subroutine hlacpy_32(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_hlacpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(..)
            type(fp16), intent(inout) :: b(..)
        end subroutine hlacpy_32

         pure subroutine h2slacpy_64(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_h2slacpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind = c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(..)
            real(c_float), intent(in) :: b(..)
        end subroutine h2slacpy_64
        pure subroutine h2slacpy_32(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_h2slacpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind = c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda, ldb
            type(fp16), intent(in) :: a(..)
            real(c_float), intent(in) :: b(..)
        end subroutine h2slacpy_32

        pure subroutine s2hlacpy_64(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2hlacpy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
        end subroutine s2hlacpy_64
        pure subroutine s2hlacpy_32(uplo, m, n, a, lda, b, ldb) bind(C, name = "lpf_blas_s2hlacpy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda, ldb
            real(c_float), intent(in) :: a(..)
            type(fp16), intent(in) :: b(..)
        end subroutine s2hlacpy_32
    end interface

    interface laset
        pure subroutine hlaset_64(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_hlaset_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind=c_char), intent(in) :: uplo
            integer(c_int64_t), intent(in) :: m, n, lda
            type(fp16), intent(in) :: alpha, beta
            type(fp16), intent(inout) :: a(..)
        end subroutine hlaset_64
        pure subroutine hlaset_32(uplo, m, n, alpha, beta, a, lda) bind(C, name = "lpf_blas_hlaset_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            character(kind=c_char), intent(in) :: uplo
            integer(c_int32_t), intent(in) :: m, n, lda
            type(fp16), intent(in) :: alpha, beta
            type(fp16), intent(inout) :: a(..)
        end subroutine hlaset_32
    end interface

    interface copy
        pure subroutine h2scopy_64(n, x, incx, y, incy) bind(C, name = "lpf_blas_h2scopy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: x(..)
            real(c_float), intent(in) :: y(..)
        end subroutine h2scopy_64
        pure subroutine h2scopy_32(n, x, incx, y, incy) bind(C, name = "lpf_blas_h2scopy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: n, incx, incy
            type(fp16), intent(in) :: x(..)
            real(c_float), intent(in) :: y(..)
        end subroutine h2scopy_32

        pure subroutine s2hcopy_64(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2hcopy_fortran_dyn_rank_64")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int64_t), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(..)
            type(fp16), intent(in) :: y(..)
        end subroutine s2hcopy_64
        pure subroutine s2hcopy_32(n, x, incx, y, incy) bind(C, name = "lpf_blas_s2hcopy_fortran_dyn_rank_32")
            use, intrinsic :: iso_c_binding
            use lpf_fp16

            integer(c_int32_t), intent(in) :: n, incx, incy
            real(c_float), intent(in) :: x(..)
            type(fp16), intent(in) :: y(..)
        end subroutine s2hcopy_32
    end interface

    !
    ! Submodules
    !
    interface scale_diag
        module subroutine scale_diag_fp16_32(m, n, a, lda, dl, dr, info)
            integer(int32), intent(in) :: m, n, lda
            integer(int32), intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dl, dr
        end subroutine

        module subroutine scale_diag_fp16_64(m, n, a, lda, dl, dr, info)
            integer(int64), intent(in) :: m, n, lda
            integer(int64), intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dl, dr
        end subroutine
    end interface scale_diag

    interface scale_diag_right
        module subroutine scale_diag_right_fp16_32(m, n, a, lda, dr, info)
            integer(int32), intent(in) :: m, n, lda
            integer(int32), intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dr
        end subroutine
        module subroutine scale_diag_right_fp16_64(m, n, a, lda, dr, info)
            integer(int64), intent(in) :: m, n, lda
            integer(int64), intent(inout) :: info
            type(fp16), intent(inout), dimension(lda, *) :: a
            type(fp16), intent(out), dimension(*) :: dr
        end subroutine
    end interface scale_diag_right


end module
