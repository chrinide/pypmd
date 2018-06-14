subroutine loadwfn(filename)

  implicit none

  character(len=*), intent(in) :: filename

  call rdwfn(filename)

  call filtergto ()

end subroutine loadwfn
