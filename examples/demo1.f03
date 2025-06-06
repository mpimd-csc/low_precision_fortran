PROGRAM FP16_DEMO
    USE FP16_SUPPORT

    TYPE(FP16) :: x1, x2, x3, x4

    x1 = 4
    x2 = -3.14159
    x3 = 4.1D1
    x4 = x2 + x3



    WRITE(*,*) "Default Output X1 = ", X1
    WRITE(*,*) "Default Output X2 = ", X2
    WRITE(*,*) "Default Output X2 = ", X3
    WRITE(*,*) "x2 + x3", X4
    WRITE(*,*) "ABS(X2)", ABS(X2)

    WRITE(*,'(A, DT)')      "Output with DT      X2 = ", X2
    WRITE(*,'(A, DT(1))')   "Output with DT(1)   X2 = ", X2
    WRITE(*,'(A, DT(6,3))') "Output with DT(6,3) X2 = ", X2
END PROGRAM FP16_DEMO
