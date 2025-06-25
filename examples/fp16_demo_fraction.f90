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

    TYPE(FP16) :: x1, x2, x3
    INTEGER:: e1, e2, e3
    TYPE(FP16) :: f1, f2, f3
    x1 = 4
    x2 = -0.5
    x3 = 2560

    e1 = exponent(x1)
    e2 = exponent(x2)
    e3 = exponent(x3)

    f1 = fraction(x1)
    f2 = fraction(x2)
    f3 = fraction(x3)

    WRITE(*,'(A, DT(12,7), A, I3, A, DT(10,7))') "Exponent ", X1, " = ", e1, " fraction = ", f1
    WRITE(*,'(A, DT(12,7), A, I3, A, DT(10,7))') "Exponent ", X2, " = ", e2, " fraction = ", f2
    WRITE(*,'(A, DT(12,7), A, I3, A, DT(10,7))') "Exponent ", X3, " = ", e3, " fraction = ", f3


END PROGRAM FP16_DEMO
