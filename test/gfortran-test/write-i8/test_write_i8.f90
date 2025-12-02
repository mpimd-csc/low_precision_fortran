! file: dt_io_example.f90
module mymod
  implicit none
  public

  type :: my_type
    real :: value
  end type my_type

  interface write(formatted)
        module procedure  write_formatted
  end interface


contains

  subroutine write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)

        class(my_type), intent(in) :: dtv
        integer, intent(in) :: unit
        character(len=*), intent(in) :: iotype
        integer(kind = 4), intent(in), target :: v_list(:)
        integer, intent(out) :: iostat
        character(len=*), intent(inout) :: iomsg

 	! local
	real :: xx
    character(len=40) :: pfmt

	xx = dtv % value

    write(*,*) "v_list = | ", v_list , " |"
    write(unit, *) xx
  end subroutine write_formatted

end module mymod

program test_dt_io
  use mymod
  implicit none

  type(my_type) :: x

  x%value = 3.14159265

  write (*, *) x
  write (*, '(DT(8,3))') x
  write (*, '(DT)')  x
  write (*, '(DT(5))') x

end program test_dt_io

