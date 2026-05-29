program fp8
    use LPF_FP8_E5M2
    use lpf_blas_fp8_e5m2

    type(FP8_E5M2) :: x(4)
    type(FP8_E5M2) :: as

    x(1) = 1.0
    x(2) = 1.0
    x(3) = 2.0
    x(4) = 3.0

    as = asum(4, x, 1)

    write(*,* ) "as = ", as

end program
