module qr_calls


contains
    subroutine qrf2(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        type(qr_stats) :: st
        real(real32) :: srcA(lda, *)
        character :: house_norm
        logical:: econ_q

        external slange
        real(real32) :: slange

        type(fp16), allocatable :: A(:,:)
        type(fp16), allocatable :: R(:,:)
        type(fp16), allocatable :: Q(:,:)
        type(fp16), allocatable :: QQ(:,:)
        type(fp16), allocatable :: tau(:)
        type(fp16), allocatable :: work(:)

        integer(lpf_default_int_kind) :: k, l, info, col_q, run
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        k = min(m,n)
        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if

        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(tau(k))
        allocate(work(max(m,n)))

        st = qr_stats_init(m,n, house_norm, "geqr2")

        A(1:m, 1:n) = srcA(1:m,1:n)


        toc = 0.0
        do run = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)

            tic = lpf_get_wtime()
            call geqr2(house_norm, m, n, A, m, tau, work, info)
            toc = toc + (lpf_get_wtime()-tic)
        end do


        st % walltime = toc / real(runs, kind = real64)
        st % info = info

        ! Get R
        R = 0.0
        call hlacpy("u", k, n, A, m , R, k)

        ! Get Q
        Q = 0.0
        QQ = 0.0
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        call orm2r("L", "N", m, col_q, k, A, m, tau, Q, m, work, info )

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(col_q))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, work)

    end subroutine


    subroutine qrf2v(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        type(qr_stats) :: st
        real(real32) :: srcA(lda, *)
        character :: house_norm
        logical:: econ_q

        external slange
        real(real32) :: slange

        type(fp16), allocatable :: A(:,:)
        type(fp16), allocatable :: R(:,:)
        type(fp16), allocatable :: Q(:,:)
        type(fp16), allocatable :: QQ(:,:)
        type(fp16), allocatable :: tau(:)
        type(fp16), allocatable :: work(:)
        type(fp16), allocatable :: diagr(:)

        integer(lpf_default_int_kind) :: k, l, info, col_q, run
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        k = min(m,n)
        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if

        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(tau(k))
        allocate(diagr(k))

        allocate(work(max(m,n)))

        st = qr_stats_init(m,n, house_norm, "geqr2v")

        A(1:m, 1:n) = srcA(1:m,1:n)


        toc = 0.0
        do run = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)

            tic = lpf_get_wtime()
            call geqr2v(house_norm, m, n, A, m, diagr, tau, work, info)
            toc = toc + (lpf_get_wtime()-tic)
        end do


        st % walltime = toc / real(runs, kind = real64)
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
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        call orm2rv("L", "N", m, col_q, k, A, m, tau, Q, m, work, info )

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(col_q))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, diagr, work)

    end subroutine

    subroutine qrf(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        logical :: econ_q
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

        integer(lpf_default_int_kind) :: k, l, info, col_q, run
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        nb = 64
        lwork =  nb * max(m,n) + (nb+1)**2

        k = min(m,n)
        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(tau(k))
        allocate(work(lwork))

        st = qr_stats_init(m,n, house_norm, "geqrf")


        toc = 0.0
        do run = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)
            tic = lpf_get_wtime()
            call geqrf(house_norm, m, n, A, m, tau, work, lwork, info)
            toc = toc  +(lpf_get_wtime() - tic)
        end do



        st % walltime = toc / real(runs, kind=real64)
        st % info = info

        ! Get R
        R = 0.0
        call hlacpy("u", k, n, A, m , R, k)

        ! Get Q
        Q = 0.0
        QQ = 0.0
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
        call ormqr("L", "N", m, col_q, k, A, m, tau, Q, m, work, lwork, info )

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, work)

    end subroutine


    subroutine qrfv(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        logical :: econ_q
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
        type(fp16), allocatable :: diagr(:)
        type(fp16), allocatable :: work(:)

        integer(lpf_default_int_kind) :: k, l, info, col_q, run
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        nb = 64
        lwork =  nb * max(m,n) + (nb+1)**2

        k = min(m,n)
        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(tau(k))
        allocate(diagr(k))
        allocate(work(lwork))

        st = qr_stats_init(m,n, house_norm, "geqrfv")


        toc = 0.0
        do run = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)
            tic = lpf_get_wtime()
            call geqrfv(house_norm, m, n, A, m, diagr, tau, work, lwork, info)
            toc = toc  +(lpf_get_wtime() - tic)
        end do



        st % walltime = toc / real(runs, kind=real64)
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
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
        call ormqrv("L", "N", m, col_q, k, A, m, tau, Q, m, work, lwork, info )

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, diagr, work)

    end subroutine

    subroutine qrt(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs, ldt
        type(qr_stats) :: st
        logical :: econ_q
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

        integer(lpf_default_int_kind) :: k, l, info, run, col_q
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        nb = 64
        lwork =  nb * max(m,n) + (nb+1)**2
        ldt = nb
        k = min(m,n)
        nb = min(nb, k)

        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(tau(k))
        allocate(T(ldt, k))
        allocate(work(lwork))

        st = qr_stats_init(m,n, house_norm, "geqrt")


        toc = 0.0
        do run = 1, runs

            A(1:m, 1:n) = srcA(1:m,1:n)
            tic = lpf_get_wtime()
            call geqrt(house_norm, m, n, nb, A, m, T, ldt, work, info)
            toc = toc + ( lpf_get_wtime() - tic)
        end do

        st % walltime = toc / real ( runs, kind = real64)
        st % info = info

        ! Get R
        R = 0.0
        call hlacpy("u", k, n, A, m , R, k)

        ! Get Q
        Q = 0.0
        QQ = 0.0
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        call gemqrt("L", "N", m, col_q, k, nb, A, m, T, ldt, Q, m, work, info )

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, T,  work)

    end subroutine

    subroutine qrtv(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs, ldt
        type(qr_stats) :: st
        logical :: econ_q
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
        type(fp16), allocatable :: diagr(:)
        type(fp16), allocatable :: T(:,:)

        type(fp16), allocatable :: work(:)

        integer(lpf_default_int_kind) :: k, l, info, run, col_q
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        nb = 64
        lwork =  nb * max(m,n) + (nb+1)**2
        ldt = nb
        k = min(m,n)
        nb = min(nb, k)

        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(tau(k))
        allocate(diagr(k))
        allocate(T(ldt, k))
        allocate(work(lwork))

        st = qr_stats_init(m,n, house_norm, "geqrtv")


        toc = 0.0
        do run = 1, runs

            A(1:m, 1:n) = srcA(1:m,1:n)
            tic = lpf_get_wtime()
            call geqrtv(house_norm, m, n, nb, A, m, diagr, T, ldt, work, info)
            toc = toc + ( lpf_get_wtime() - tic)
        end do

        st % walltime = toc / real ( runs, kind = real64)
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
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        call gemqrtv("L", "N", m, col_q, k, nb, A, m, T, ldt, Q, m, work, info )

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, diagr, T,  work)

    end subroutine

    subroutine tsqr(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        logical :: econ_q
        type(qr_stats) :: st
        real(real32) :: srcA(lda, *)
        character :: house_norm
        integer(lpf_default_int_kind) :: mb, nb, lwork, num_blocks, ldt, run, col_q

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

        if ( econ_q ) then
            col_q = k
        else
            col_q = m
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(T(ldt, n * num_blocks))
        allocate(tau(k))
        allocate(work(lwork))

        st = qr_stats_init(m,n, house_norm, "tsqr")


        toc = 0.0
        do run  = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)

            tic = lpf_get_wtime()
            call latsqr(house_norm, m, n, mb, nb, A, m, T, ldt, work, lwork, info)
            toc = toc + (lpf_get_wtime() - tic)
        end do



        st % walltime = toc / real(runs, kind=real64)
        st % info = info

        ! Get R
        R = 0.0
        call hlacpy("u", k, n, A, m , R, k)

        ! Get Q
        Q = 0.0
        QQ = 0.0
        do l = 1, col_q
            Q(l,l) = 1.0
            QQ(l,l) = 1.0
        end do
        ! call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )
        ! call ormqr("L", "N", m, m, k, A, m, tau, Q, m, work, lwork, info )
        lwork =  nb * max(m,n) + (nb+1)**2

        call lamtsqr("L", "N", m, col_q, k, mb, nb, A, m, T, ldt, Q, m, work, lwork, info)

        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, tau, T,  work)

    end subroutine

    subroutine cholqr(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        type(qr_stats) :: st
        real(real32) :: srcA(lda, *)
        character :: house_norm
        logical:: econ_q

        external slange
        real(real32) :: slange

        type(fp16), allocatable :: A(:,:)
        type(fp16), allocatable :: R(:,:)
        type(fp16), allocatable :: Q(:,:)
        type(fp16), allocatable :: QQ(:,:)

        integer(lpf_default_int_kind) :: k, l, info, col_q, run
        real(real64) :: tic, toc
        real(real32) :: nrm_A
        type(fp16) :: work(2)

        k = min(m,n)

        if ( econ_q ) then
            col_q = k
        else
            col_q = k
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))

        st = qr_stats_init(m,n, "X", "cholqr")

        R = 0.0
        toc = 0.0
        do run  = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)

            tic = lpf_get_wtime()
            call gecholqr(m, n, A, m, R, k, info)
            toc = toc + (lpf_get_wtime() - tic)
        end do



        st % walltime = toc / real(runs, kind=real64)
        st % info = info

        ! Get Q
        call hlacpy("all", m, n, A, m, Q, m)
        QQ = 0.0
        do l = 1, col_q
            QQ(l,l) = 1.0
        end do
        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ)

    end subroutine

    subroutine cholqr_shift(m, n, house_norm, srcA, lda, econ_q, runs, st)
        use qr_stat
        use lpf_fp16
        use lpf_lapack_fp16
        use lpf_types
        use iso_fortran_env, only: real32, real64
        implicit none

        integer(lpf_default_int_kind) :: m, n, lda, runs
        type(qr_stats) :: st
        real(real32) :: srcA(lda, *)
        character :: house_norm
        logical:: econ_q

        external slange
        real(real32) :: slange

        type(fp16), allocatable :: A(:,:)
        type(fp16), allocatable :: R(:,:)
        type(fp16), allocatable :: Q(:,:)
        type(fp16), allocatable :: QQ(:,:)
        type(fp16), allocatable :: work(:)

        integer(lpf_default_int_kind) :: k, l, info, col_q, run
        real(real64) :: tic, toc
        real(real32) :: nrm_A

        k = min(m,n)

        if ( econ_q ) then
            col_q = k
        else
            col_q = k
        end if


        allocate(A(m,n))
        allocate(R(k,n))
        allocate(Q(m,col_q))
        allocate(QQ(col_q,col_q))
        allocate(work(m*n))

        st = qr_stats_init(m,n, "X", "cholqr_shift")

        R = 0.0
        toc = 0.0
        do run  = 1, runs
            A(1:m, 1:n) = srcA(1:m,1:n)

            tic = lpf_get_wtime()
            call gecholqr_shift(m, n, A, m, R, k, work, info)
            toc = toc + (lpf_get_wtime() - tic)
        end do



        st % walltime = toc / real(runs, kind=real64)
        st % info = info

        ! Get Q
        call hlacpy("all", m, n, A, m, Q, m)
        QQ = 0.0
        do l = 1, col_q
            QQ(l,l) = 1.0
        end do
        call hgemm("T", "N", col_q, col_q, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, col_q)

        st % orth_err = lange("F", col_q, col_q, QQ, col_q, work) / sqrt(real(m))


        A(1:m, 1:n) = srcA(1:m,1:n)

        call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
        nrm_A = slange("F", m, n, srcA, lda, work)
        st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

        deallocate(A, R, Q, QQ, work)

    end subroutine


end module

