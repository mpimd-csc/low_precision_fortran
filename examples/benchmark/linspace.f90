SUBROUTINE test_linspace()
    USE FP16_SUPPORT
    IMPLICIT NONE

    INTEGER, PARAMETER :: N = 1024
    TYPE(FP16) :: x(N+1)
    TYPE(FP16) :: y(N+1)
    TYPE(FP16) :: a, b, h
    INTEGER :: K

    a = 0
    b = 1
    h = (b - a) / FP16(N)

    DO K = 1, N+1
        x(k) = a + FP16(K-1) * h
        y(k) = FP16(0.625) * x(k)*x(k) + FP16(0.125)*x(k)  + FP16(1.0)
    END DO

    ! WRITE(*,*) x, y
END SUBROUTINE


