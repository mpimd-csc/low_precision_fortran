module qr_stat
    use iso_fortran_env, only: real32, real64
    use lpf_types
    implicit none

    type :: qr_stats
        character(len=1) :: house_norm
        character(len=15):: algorithm
        real(real32) :: orth_err
        real(real32) :: reconstruct_err
        real(real64) :: walltime
        integer(lpf_default_int_kind) :: info
        integer(lpf_default_int_kind) :: m, n
    end type

    interface qr_stats
        module procedure :: qr_stats_init
    end interface
contains

    function qr_stats_init(m ,n, house_norm, algorithm )
        integer(lpf_default_int_kind) :: m, n
        type(qr_stats) :: qr_stats_init
        character :: house_norm
        character(len=*) :: algorithm

        qr_stats_init % house_norm  = house_norm
        qr_stats_init % algorithm   = algorithm
        qr_stats_init % m = m
        qr_stats_init % n = n
        qr_stats_init % orth_err = 0.0
        qr_stats_init % reconstruct_err = 0.0
        qr_stats_init % walltime = 0.0
        return
    end function

    subroutine print_qr_stats( st )
        type(qr_stats) :: st

        write(*, '( I7, " | ", I7, " | ", A7, " | ", ES20.13, " | ", ES20.13, " | ", F20.5, " | ", A15)') st%m, st%n, &
            & st%house_norm, &
            & st%reconstruct_err, &
            & st%orth_err, st%walltime, st % algorithm
    end subroutine

    subroutine print_qr_stats_header()
     write(*, '( A7, " | ", A7, " | ", A7, " | ", A20, " | ", A20, " | ", A20, " | ", A20 )') "M", "N", "H-Norm", &
         & "||A-QR||_2 rel.", &
         & "||Q**TQ-I||_2", " Time [s]", ""
     write(*, '( A7, "-|-", A7, "-|-", A7, "-|-", A20, "-|-", A20, "-|-", A20, "-|-", A20)') "--------", "--------", "-------", &
         & "------------------------", &
         & "-----------------------", "-----------------------------", "-----------------------------"




    end subroutine

end module
