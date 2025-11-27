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

PROGRAM BF16_DEMO
    USE LPF_BF16

    TYPE(BF16) :: x1, x2, x3, x4

    x1 = 4
    x2 = -3.14159
    x3 = 4.1D1
    x4 = x2 + x3



    WRITE(*,*) "Default Output X1 = ", X1
    WRITE(*,*) "Default Output X2 = ", X2
    WRITE(*,*) "Default Output X2 = ", X3
    WRITE(*,*) "x2 + x3", X4
    WRITE(*,*) "ABS(X2)", ABS(X2), ABS(-4.0)
    WRITE(*,*) "Exponent ", X1, " = ", exponent(X1)
    WRITE(*,*) "Exponent ", BF16(0.0), " = ", exponent(BF16(0.0))
    WRITE(*,*) "Exponent ", BF16(1.0), " = ", exponent(BF16(1.0))
    WRITE(*,*) "Exponent ", BF16(2.0), " = ", exponent(BF16(2.0))
    WRITE(*,*) "Exponent ", BF16(0.00006103515625D0), " = ", exponent(BF16(0.00006103515625D0))
    WRITE(*,*) "Exponent ", BF16(0.99951172D0), " = ", exponent(BF16(0.99951172D0))

    WRITE(*,'(A, DT)')      "Output with DT      X2 = ", X2
    WRITE(*,'(A, DT(1))')   "Output with DT(1)   X2 = ", X2
    WRITE(*,'(A, DT(6,3))') "Output with DT(6,3) X2 = ", X2
END PROGRAM BF16_DEMO
