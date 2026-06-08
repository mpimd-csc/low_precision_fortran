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
submodule (lpf_blas_fp8_e5m2) lpf_blas_rotg_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rotg_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rotg_impl(a,b,c,s)
        implicit none
        type(DT), intent(inout) :: a, b
        type(DT), intent(out) :: c, s

        type(DT) :: anorm, bnorm, scl, sigma, r, z
        type(DT) :: zero, one, safmin, safmax

        zero = 0.0
        one = 1.0
        safmin = DT(radix(zero))**max(minexponent(zero)-1, 1-maxexponent(zero))
        safmax = DT(radix(zero))**max(1-minexponent(zero), maxexponent(zero)-1)

        anorm = abs(a)
        bnorm = abs(b)
        if( bnorm == zero ) then
            c = one
            s = zero
            b = zero
        else if( anorm == zero ) then
            c = zero
            s = one
            a = b
            b = one
        else
            scl = min( safmax, max( safmin, anorm, bnorm ) )
            if( anorm > bnorm ) then
                sigma = sign(one, a)
            else
                sigma = sign(one, b)
            end if
            r = sigma*( scl*sqrt((a/scl)**2 + (b/scl)**2) )
            c = a/r
            s = b/r

            if( anorm > bnorm ) then
                z = s
            else if( c /= zero ) then
                z = one/c
            else
                z = one
            end if
            a = r
            b = z
        end if
    end subroutine

    module subroutine rotg(sa,sb,sc,ss)
        implicit none
        type(DT), intent(inout) :: sa, sb
        type(DT), intent(out) :: sc, ss

        call rotg_impl(sa, sb, sc, ss)
    end subroutine

end submodule
