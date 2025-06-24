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
