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

#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_rotmg_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rotmg_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rotmg_impl(sd1,sd2,sx1,sy1,sparams)
        implicit none
        type(DT), intent(inout) :: sd1, sd2, sx1
        type(DT), intent(in) :: sy1
        type(DT), intent(out) :: sparams(*)

        type(DT) :: gam, gamsq, one, rgamsq, sflag, sh11, sh12, sh21, sh22, sp1, sp2, sq1, sq2, stemp, su, two, zero

        zero = 0.0
        one = 1.0
        two = 2.0
        gam = 16.0
        gamsq = 256.0
        rgamsq = 3.90625E-3

        if (sd1.lt.zero) then
            sflag = -one
            sh11 = zero
            sh12 = zero
            sh21 = zero
            sh22 = zero
            sd1 = zero
            sd2 = zero
            sx1 = zero
        else
            sp2 = sd2*sy1
            if (sp2.eq.zero) then
                sflag = -two
                sparams(1) = sflag
                return
            end if
            sp1 = sd1*sx1
            sq2 = sp2*sy1
            sq1 = sp1*sx1

            if (abs(sq1).gt.abs(sq2)) then
                sh21 = -sy1/sx1
                sh12 = sp2/sp1
                su = one - sh12*sh21
                if (su.gt.zero) then
                    sflag = zero
                    sd1 = sd1/su
                    sd2 = sd2/su
                    sx1 = sx1*su
                else
                    sflag = -one
                    sh11 = zero
                    sh12 = zero
                    sh21 = zero
                    sh22 = zero
                    sd1 = zero
                    sd2 = zero
                    sx1 = zero
                end if
            else
                if (sq2.lt.zero) then
                    sflag = -one
                    sh11 = zero
                    sh12 = zero
                    sh21 = zero
                    sh22 = zero
                    sd1 = zero
                    sd2 = zero
                    sx1 = zero
                else
                    sflag = one
                    sh11 = sp1/sp2
                    sh22 = sx1/sy1
                    su = one + sh11*sh22
                    stemp = sd2/su
                    sd2 = sd1/su
                    sd1 = stemp
                    sx1 = sy1*su
                end if
            end if

            if (sd1.ne.zero) then
                do while ((sd1.le.rgamsq) .or. (sd1.ge.gamsq))
                    if (sflag.eq.zero) then
                        sh11 = one
                        sh22 = one
                        sflag = -one
                    else
                        sh21 = -one
                        sh12 = one
                        sflag = -one
                    end if
                    if (sd1.le.rgamsq) then
                        sd1 = sd1*gam**2
                        sx1 = sx1/gam
                        sh11 = sh11/gam
                        sh12 = sh12/gam
                    else
                        sd1 = sd1/gam**2
                        sx1 = sx1*gam
                        sh11 = sh11*gam
                        sh12 = sh12*gam
                    end if
                end do
            end if

            if (sd2.ne.zero) then
                do while ( (abs(sd2).le.rgamsq) .or. (abs(sd2).ge.gamsq) )
                    if (sflag.eq.zero) then
                        sh11 = one
                        sh22 = one
                        sflag = -one
                    else
                        sh21 = -one
                        sh12 = one
                        sflag = -one
                    end if
                    if (abs(sd2).le.rgamsq) then
                        sd2 = sd2*gam**2
                        sh21 = sh21/gam
                        sh22 = sh22/gam
                    else
                        sd2 = sd2/gam**2
                        sh21 = sh21*gam
                        sh22 = sh22*gam
                    end if
                end do
            end if
        end if

        if (sflag.lt.zero) then
            sparams(2) = sh11
            sparams(3) = sh21
            sparams(4) = sh12
            sparams(5) = sh22
        else if (sflag.eq.zero) then
            sparams(3) = sh21
            sparams(4) = sh12
        else
            sparams(2) = sh11
            sparams(5) = sh22
        end if
        sparams(1) = sflag
    end subroutine

    module subroutine rotmg(sa,sb,sc,sd,sparam)
        implicit none
        type(DT), intent(inout) :: sa, sb, sc
        type(DT), intent(inout) :: sd
        type(DT), target, intent(out) :: sparam(..)

        type(DT), CONTIGUOUS, pointer :: psparams(:)
        integer(int32) :: lenp
        integer(int32) :: lp(1)
        type(c_ptr) :: ptr

        lenp = 5
        lp(1) = lenp
        ptr = c_loc(sparam)
        call c_f_pointer(ptr, psparams, lp)

        call rotmg_impl(sa, sb, sc, sd, psparams)
    end subroutine

end submodule
