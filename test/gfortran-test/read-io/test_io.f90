module test_io

    TYPE :: MYTYPE
        REAL :: value
    END TYPE

    INTERFACE read(formatted)
        MODULE PROCEDURE read_formatted
    END INTERFACE
    PUBLIC :: read(formatted)

contains

    ! Formatted Input
    SUBROUTINE read_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        CLASS(MYTYPE), INTENT(INOUT) :: dtv
        INTEGER, INTENT(IN)      :: unit
        CHARACTER(*), INTENT(IN) :: iotype
        INTEGER, INTENT(IN)      :: v_list(:)
        INTEGER, INTENT(OUT)     :: iostat
        CHARACTER(*), INTENT(INOUT) :: iomsg

        REAL   :: tmp

        READ(unit, FMT = *, IOSTAT=iostat, IOMSG=iomsg) tmp
        IF (iostat == 0) dtv%value = tmp
    END SUBROUTINE read_formatted

end module

PROGRAM MAIN
    USE test_io

    INTEGER, PARAMETER  :: NIN = 5
    TYPE(MYTYPE)       :: V11, V12, V13
    ! REAL       :: V11, V12, V13
    INTEGER            :: V21, V22, V23

    READ(NIN, FMT = *) V11, V23, V12
    READ(NIN, FMT = *) V21, V13, V22
    WRITE(*,*) "V1 = ", V11%value, V12%value, V13%value
    ! WRITE(*,*) "V1 = ", V11, V12, V13
    WRITE(*,*) "V2 = ", V21, V22, V23
END PROGRAM MAIN
