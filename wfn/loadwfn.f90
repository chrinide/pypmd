subroutine loadwfn(filename)

  use mod_param, only: verbose
  implicit none

  character(len=*), intent(in) :: filename

  call rdwfn (filename)
  if (verbose) call infowfn ()
  call filtergto ()
  call isorthowfn ()

end subroutine loadwfn
