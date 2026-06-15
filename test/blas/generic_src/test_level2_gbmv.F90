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


#ifdef LPF_FP8_E4M3
#define DT fp8_e4m3
#define BLASMOD lpf_blas_fp8_e4m3
#define TYPEMOD lpf_fp8_e4m3
#endif
#ifdef LPF_FP8_E5M2
#define DT fp8_e5m2
#define BLASMOD lpf_blas_fp8_e5m2
#define TYPEMOD lpf_fp8_e5m2
#endif
#ifdef LPF_FP16
#define DT fp16
#define BLASMOD lpf_blas_fp16
#define TYPEMOD lpf_fp16
#endif
#ifdef LPF_BF16
#define DT bf16
#define BLASMOD lpf_blas_bf16
#define TYPEMOD lpf_bf16
#endif

program test_level2_gbmv
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical('N')
    call test_typical('T')
    call test_edge()
    call test_stride()

    call test_summary()

contains

    subroutine test_typical(trans)
        character(len=1) :: trans
        integer(int64) :: m = 3, n = 3
        integer(int64) :: kl = 1, ku = 1
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(3, 3) :: a
        type(DT), dimension(3) :: x, y
        integer :: i
        real(real64) :: dalpha, dbeta
        real(real64), dimension(3, 3) :: da
        real(real64), dimension(3) :: dx, dy


        alpha = 2.0
        beta = 1.0
        a = 1.0
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]

        dalpha = alpha
        dbeta  = beta
        da = a
        dx = x
        dy = y


        call gbmv(trans, m, n, kl, ku, alpha, a, kl+ku+1, x, incx, beta, y, incy)
        call dgbmv('N', m, n, kl, ku, dalpha, da, kl+ku+1, dx, incx, dbeta, dy, incy)

        do i = 1, 3
            call check_dt_real64("gbmv_typical_" // trans, int(i,int64), y(i), dy(i), GENERIC_TOL)
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: m = 0, n = 3
        integer(int64) :: kl = 1, ku = 1
        integer(int64) :: incx = 1, incy = 1
        type(DT) :: alpha, beta
        type(DT), dimension(0, 3) :: a
        type(DT), dimension(3) :: x, y

        alpha = 1.0
        beta = 1.0
        x = [1.0, 2.0, 3.0]
        y = [1.0, 1.0, 1.0]

        call gbmv('N', m, n, kl, ku, alpha, a, kl+ku+1, x, incx, beta, y, incy)
    end subroutine

    subroutine test_stride()
        integer(int64) :: m = 3, n = 3
        integer(int64) :: kl = 1, ku = 1
        integer(int64) :: incx = 2, incy = 2
        type(DT) :: alpha, beta
        type(DT), dimension(3, 3) :: a
        type(DT), dimension(6) :: x, y
        real(real64) :: dalpha, dbeta
        real(real64), dimension(3, 3) :: da
        real(real64), dimension(6) :: dx, dy

        integer :: i

        alpha = DT(1.0)
        beta = DT(1.0)
        a = 1.0
        x = [1.0, 0.0, 2.0, 0.0, 1.0, 0.0]
        y = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

        dalpha = alpha
        dbeta  = beta
        da = a
        dx = x
        dy = y

        call gbmv('N', m, n, kl, ku, alpha, a, kl+ku+1, x, incx, beta, y, incy)
        call dgbmv('N', m, n, kl, ku, dalpha, da, kl+ku+1, dx, incx, dbeta, dy, incy)

        do i = 1, 3
          call check_dt_real64("gbmv_stride", int(i,int64), y(2*i-1), dy(2*i-1), GENERIC_TOL)
        end do
    end subroutine

    subroutine dgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx, beta,y,incy)
        implicit none
        real(real64) :: alpha,beta
        integer(int64) ::incx,incy,kl,ku,lda,m,n
        character :: trans
        real(real64) :: a(lda,*),x(*),y(*)

        real(real64) one,zero
        parameter (one=1.0d+0,zero=0.0d+0)

        real(real64) :: temp
        integer(int64) :: i,ix,iy,j,jx,jy,k,kup1,kx,ky,lenx,leny

        if ((m.eq.0) .or. (n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return

        if (trans .eq. 'n') then
            lenx = n
            leny = m
        else
            lenx = m
            leny = n
        end if

        if (incx.gt.0) then
            kx = 1
        else
            kx = 1 - (lenx-1)*incx
        end if
        if (incy.gt.0) then
            ky = 1
        else
            ky = 1 - (leny-1)*incy
        end if

        if (beta.ne.one) then
            if (incy.eq.1) then
                if (beta.eq.zero) then
                    do i = 1,leny
                        y(i) = zero
                    end do
                else
                    do i = 1,leny
                        y(i) = beta*y(i)
                    end do
                end if
            else
                iy = ky
                if (beta.eq.zero) then
                    do i = 1,leny
                        y(iy) = zero
                        iy = iy + incy
                    end do
                else
                    do i = 1,leny
                        y(iy) = beta*y(iy)
                        iy = iy + incy
                    end do
                end if
            end if
        end if
        if (alpha.eq.zero) return
        kup1 = ku + 1
        if (trans .eq. 'n') then
            jx = kx
            if (incy.eq.1) then
                do j = 1,n
                    temp = alpha*x(jx)
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        y(i) = y(i) + temp*a(k+i,j)
                    end do
                    jx = jx + incx
                end do
            else
                do j = 1,n
                    temp = alpha*x(jx)
                    iy = ky
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        y(iy) = y(iy) + temp*a(k+i,j)
                        iy = iy + incy
                    end do
                    jx = jx + incx
                    if (j.gt.ku) ky = ky + incy
                end do
            end if
        else
            jy = ky
            if (incx.eq.1) then
                do j = 1,n
                    temp = zero
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        temp = temp + a(k+i,j)*x(i)
                    end do
                    y(jy) = y(jy) + alpha*temp
                    jy = jy + incy
                end do
            else
                do j = 1,n
                    temp = zero
                    ix = kx
                    k = kup1 - j
                    do i = max(1,j-ku),min(m,j+kl)
                        temp = temp + a(k+i,j)*x(ix)
                        ix = ix + incx
                    end do
                    y(jy) = y(jy) + alpha*temp
                    jy = jy + incy
                    if (j.gt.ku) kx = kx + incx
                end do
            end if
        end if
        return
    end subroutine

end program test_level2_gbmv
