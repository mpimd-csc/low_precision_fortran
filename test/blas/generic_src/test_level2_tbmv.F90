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

program test_level2_tbmv
    use BLASMOD
    use TYPEMOD
    use generic_test_utils
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    call test_typical('U', 'N', 'N')
    call test_typical('L', 'N', 'N')
    call test_typical('U', 'T', 'N')
    call test_typical('L', 'T', 'N')

    call test_typical('U', 'N', 'U')
    call test_typical('L', 'N', 'U')
    call test_typical('U', 'T', 'U')
    call test_typical('L', 'T', 'U')

    call test_stride('U', 'N', 'N')
    call test_stride('L', 'N', 'N')
    call test_stride('U', 'T', 'N')
    call test_stride('L', 'T', 'N')

    call test_stride('U', 'N', 'U')
    call test_stride('L', 'N', 'U')
    call test_stride('U', 'T', 'U')
    call test_stride('L', 'T', 'U')



    call test_edge()

    call test_summary()

contains

    subroutine test_typical(uplo, trans, diag)
        character(len=1) :: uplo, trans, diag
        integer(int64) :: n = 3
        integer(int64) :: k = 1
        integer(int64) :: incx = 1
        type(DT), dimension(2, 3) :: a
        type(DT), dimension(3,3) :: ao, bo
        type(DT), dimension(3) :: x, x_old
        integer(int64) :: i, j, m


        ao = 0.0

        do i = 1, 3
            ao(i,i) = 1.0
            if ( i .lt. 3) ao(i, i+1) = 0.5
        end do


        if (uplo.eq.'U') then
            do j = 1, n
                m = k + 1 - j
                do i = max(1, j-k), j
                    a(m+i, j) = ao(i,j)
                end do
            end do
        else
            bo = TRANSPOSE(ao)
            ao = bo
            do j = 1, n
                m = 1 - j
                do i = j, min(n, j+k)
                    a(m+i, j) = ao (i,j)
                end do
            end do
        end if

        x = [1.0, 2.0, 3.0]

        call gemv(trans, n, n, DT(1.0), ao, n, x, incx, DT(0.0), x_old, incx)
        call tbmv(uplo, trans, diag, n, k, a, k+1, x, incx)

        do i = 1, 3
            call check_dt_real64("tbmv_typical_" // uplo // trans // diag, int(i,int64), x(i), dble(x_old(i)), GENERIC_TOL)
        end do
    end subroutine

    subroutine test_edge()
        integer(int64) :: n = 0
        integer(int64) :: k = 1
        integer(int64) :: incx = 1
        type(DT), dimension(1, 0) :: a
        type(DT), dimension(0) :: x

        call tbmv('U', 'N', 'N', n, k, a, k+1, x, incx)
    end subroutine

    subroutine test_stride(uplo,trans,diag)
        character(len=1) :: uplo, trans, diag
        integer(int64) :: n = 3
        integer(int64) :: k = 1
        integer(int64) :: incx = 2
        type(DT), dimension(2, 3) :: a
        type(DT), dimension(3,3) :: ao, bo
        type(DT), dimension(5) :: x, x_old
        integer(int64)  :: i, j, m


        ao = 0.0

        do i = 1, 3
            ao(i,i) = 1.0
            if ( i .lt. 3) ao(i, i+1) = 0.5
        end do


        if (uplo.eq.'U') then
            do j = 1, n
                m = k + 1 - j
                do i = max(1, j-k), j
                    a(m+i, j) = ao(i,j)
                end do
            end do
        else
            bo = TRANSPOSE(ao)
            ao = bo
            do j = 1, n
                m = 1 - j
                do i = j, min(n, j+k)
                    a(m+i, j) = ao (i,j)
                end do
            end do
        end if

        x = [1.0, 0.0, 2.0, 0.0, 3.0]

        call gemv(trans, n, n, DT(1.0), ao, n, x, incx, DT(0.0), x_old, incx)
        call tbmv(uplo, trans, diag, n, k, a, k+1, x, incx)

        do i = 1, 5
            call check_dt_real64("tbmv_stride_" // uplo // trans // diag, int(i,int64), x(i), dble(x_old(i)), GENERIC_TOL)
        end do

    end subroutine

end program test_level2_tbmv
