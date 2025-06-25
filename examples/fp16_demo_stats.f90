!  SPDX-License-Identifier LGPL-3.0-or-later
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

PROGRAM FP16_DEMO
    USE FP16_SUPPORT

    TYPE(FP16) :: x1, x2, x3, x4
    x1 = 1.0
    WRITE(*,*) "EPSILON FP16     = ", EPSILON(x1)
    WRITE(*,*) "HUGE FP16        = ", HUGE(x1)
    WRITE(*,*) "TINY FP16        = ", TINY(x1)
    WRITE(*,*) "MINEXPONENT FP16 = ", MINEXPONENT(x1)
    WRITE(*,*) "MAXEXPONENT FP16 = ", MAXEXPONENT(x1)
    WRITE(*,*) "PRECISION   FP16 = ", PRECISION(x1)
    WRITE(*,*) "RANGE       FP16 = ", RANGE(x1)
    WRITE(*,*) "RADIX       FP16 = ", RADIX(x1)


END PROGRAM FP16_DEMO
