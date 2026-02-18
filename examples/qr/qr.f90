program qr_benchmark
    use qr_stat
    use qr_calls
    use lpf_types
    implicit none

    integer(lpf_default_int_kind) :: m, n, runs

    real (real32), allocatable :: A(:,:)
    character(len = 1024) :: arg
    type(qr_stats) :: st
    integer, dimension(4) :: seed =(/1,1,1,1/)
    external slarnv
    integer :: qr_type
    logical :: econ_q
    integer :: type_start

    m = 100
    n = 100
    runs = 1

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

    econ_q = .true.

    if ( m .le. 2000 ) then
        type_start = 1
    else
        type_start = 5
    end if

    do qr_type = type_start, 16
        select case(qr_type)
            case (1)
                call qrf2 (m, n, 'L', A, m, econ_q, runs, st)
            case (2)
                call qrf2 (m, n, 'H', A, m, econ_q, runs, st)
            case (3)
                call qrf2v(m, n, 'L', A, m, econ_q, runs, st)
            case (4)
                call qrf2v(m, n, 'H', A, m, econ_q, runs, st)
            case (5)
                call qrf  (m, n, 'L', A, m, econ_q, runs, st)
            case (6)
                call qrf  (m, n, 'H', A, m, econ_q, runs, st)
            case (7)
                call qrfv (m, n, 'L', A, m, econ_q, runs, st)
            case (8)
                call qrfv (m, n, 'H', A, m, econ_q, runs, st)
            case (9)
                call qrt  (m, n, 'L', A, m, econ_q, runs, st)
            case (10)
                call qrt  (m, n, 'H', A, m, econ_q, runs, st)
            case (11)
                call qrtv (m, n, 'L', A, m, econ_q, runs, st)
            case (12)
                call qrtv (m, n, 'H', A, m, econ_q, runs, st)
            case (13)
                call tsqr (m, n, 'L', A, m, econ_q, runs, st)
            case (14)
                call tsqr (m, n, 'H', A, m, econ_q, runs, st)
            case (15)
                call cholqr (m, n, 'H', A, m, econ_q, runs, st)
            case (16)
                call cholqr_shift (m, n, 'H', A, m, econ_q, runs, st)



        end select

        call print_qr_stats(st)

    end do

    deallocate(A)


end program qr_benchmark
