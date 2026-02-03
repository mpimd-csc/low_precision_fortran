subroutine qrf2(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    k = min(m,n)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(work(max(m,n)))

    st = qr_stats_init(m,n, house_norm, "geqr2")

    A(1:m, 1:n) = srcA(1:m,1:n)


    tic = lpf_get_wtime()
    call geqr2(house_norm, m, n, A, m, tau, work, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, work)

end subroutine



subroutine qrf2v(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: diagr(:)
    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    k = min(m,n)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(diagr(k))
    allocate(work(max(m,n)))

    st = qr_stats_init(m,n, house_norm, "geqr2v")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call geqr2_v1(house_norm, m, n, A, m, diagr, tau, work, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)
    do l = 1, k
        R(l,l) = diagr(l)
    end do

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    call orm2rv("L", "N", m, m, k, A, m, tau, Q, m, work, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, diagr, work)

end subroutine

subroutine qrf(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm
    integer(lpf_default_int_kind) :: nb, lwork

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    nb = 64
    lwork =  nb * max(m,n) + (nb+1)**2

    k = min(m,n)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(work(lwork))

    st = qr_stats_init(m,n, house_norm, "geqrf")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call geqrf(house_norm, m, n, A, m, tau, work, lwork, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
    call ormqr("L", "N", m, m, k, A, m, tau, Q, m, work, lwork, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, work)

end subroutine

subroutine qrfv(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm
    integer(lpf_default_int_kind) :: nb, lwork

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: work(:)
    type(fp16), allocatable :: diagr(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    nb = 64
    lwork =  nb * max(m,n) + (nb+1)**2

    k = min(m,n)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(work(lwork))
    allocate(diagr(min(m,n)))

    st = qr_stats_init(m,n, house_norm, "geqrfv")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call geqrfv(house_norm, m, n, A, m, diagr, tau, work, lwork, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)
    do l = 1, k
        R(l,l) = diagr(l)
    end do

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    ! call orm2rv("L", "N", m, m, k, A, m, tau, Q, m, work, info )
    call ormqrv("L", "N", m, m, k, A, m, tau, Q, m, work, lwork, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, work, diagr )

end subroutine

subroutine qrt(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda, ldt
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm
    integer(lpf_default_int_kind) :: nb, lwork

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: T(:,:)

    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    nb = 64
    lwork =  nb * max(m,n) + (nb+1)**2
    ldt = nb
    k = min(m,n)
    nb = min(nb, k)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(T(ldt, k))
    allocate(work(lwork))

    st = qr_stats_init(m,n, house_norm, "geqrt")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call geqrt(house_norm, m, n, nb, A, m, T, ldt, work, info)
    toc = lpf_get_wtime()

    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
    call gemqrt("L", "N", m, m, k, nb, A, m, T, ldt, Q, m, work, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, T,  work)

end subroutine



subroutine qrtv(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda, ldt
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm
    integer(lpf_default_int_kind) :: nb, lwork

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: T(:,:)

    type(fp16), allocatable :: diagr(:)
    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    nb = 64
    lwork =  nb * max(m,n) + (nb+1)**2
    ldt = nb
    k = min(m,n)
    nb = min(nb, k)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(T(ldt, k))
    allocate(work(lwork))
    allocate(diagr(k))

    st = qr_stats_init(m,n, house_norm, "geqrtv")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call geqrtv(house_norm, m, n, nb, A, m, diagr,  T, ldt, work, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)
    do l = 1, k
        R(l,l) = diagr(l)
    end do


    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
    call gemqrtv("L", "N", m, m, k, nb, A, m, T, ldt, Q, m, work, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, T,  diagr,  work)

end subroutine

subroutine tsqr(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm
    integer(lpf_default_int_kind) :: mb, nb, lwork, num_blocks, ldt

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: T(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    nb = min(n,64)
    mb = 2048
    ldt = nb

    lwork =  nb * max(m,n) + (nb+1)**2
    num_blocks = ((m-n)/(mb-n))+1
    k = min(m,n)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(T(ldt, n * num_blocks))
    allocate(tau(k))
    allocate(work(lwork))

    st = qr_stats_init(m,n, house_norm, "tsqr")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call latsqr(house_norm, m, n, mb, nb, A, m, T, ldt, work, lwork, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
    ! call ormqr("L", "N", m, m, k, A, m, tau, Q, m, work, lwork, info )
    lwork =  nb * max(m,n) + (nb+1)**2

    call lamtsqr("L", "N", m, m, k, mb, nb, A, m, T, ldt, Q, m, work, lwork, info)

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, T,  work)

end subroutine




program qr_benchmark
    use qr_stat
    implicit none

    integer :: m, n, runs

    real (real32), allocatable :: A(:,:)
    character(len = 1024) :: arg
    type(qr_stats) :: st
    integer, dimension(4) :: seed =(/1,1,1,1/)
    external slarnv
    integer :: qr_type

    m = 100
    n = 100
    runs = 5

    if ( command_argument_count() .ge. 2) then
        call get_command_argument(1,  arg )
        read(arg, *) m
        call get_command_argument(2,  arg )
        read(arg, *) n
    end if
    if ( command_argument_count() .eq. 3) then
        call get_command_argument(3,  arg )
        read(arg, *) runs
    endif

    write(*,'("# m    = ", I8)') m
    write(*,'("# n    = ", I8)') n
    write(*,'("# runs = ", I8)') runs

    allocate(A(m,n))

    ! call RANDOM_NUMBER(A)
    call slarnv(2, seed, m * n, A)
    ! A = 1

    call print_qr_stats_header()

    do qr_type = 11, 14
        select case(qr_type)
            case (1)
                call qrf2 (m, n, 'L', A, m, st)
            case (2)
                call qrf2 (m, n, 'H', A, m, st)
            case (3)
                call qrf2v(m, n, 'L', A, m, st)
            case (4)
                call qrf2v(m, n, 'H', A, m, st)
            case (5)
                call qrf  (m, n, 'L', A, m, st)
            case (6)
                call qrf  (m, n, 'H', A, m, st)
            case (7)
                call qrfv (m, n, 'L', A, m, st)
            case (8)
                call qrfv (m, n, 'H', A, m, st)
            case (9)
                call qrt  (m, n, 'L', A, m, st)
            case (10)
                call qrt  (m, n, 'H', A, m, st)
            case (11)
                call qrtv (m, n, 'L', A, m, st)
            case (12)
                call qrtv (m, n, 'H', A, m, st)
            case (13)
                call tsqr (m, n, 'L', A, m, st)
            case (14)
                call tsqr (m, n, 'H', A, m, st)

        end select

        call print_qr_stats(st)

    end do

    deallocate(A)


end program qr_benchmark
