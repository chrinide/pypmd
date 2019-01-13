program bzip

  use, intrinsic :: iso_c_binding, only: c_char
  use bzip

  implicit none

  integer, parameter :: mline = 132
  character(len=mline) :: arg
  type(bzFile) :: bzptr
  integer :: length
  character(len=mline) :: buf

  integer :: i, j, k, l
  real(kind=8) :: valor

  call get_command_argument(1, arg)
  if (len_trim(arg) == 0) stop
  write (*,*) "Bzip2 filename : ", trim(arg)

  bzptr = bzOpenR(trim(arg))
  length = 0
  call bzReadline (bzptr, buf, length)
  write (*,'(a)') '1-RDM'
  do
    length = 0
    call bzReadline (bzptr, buf, length)
    if (buf(1:5) .eq. '2-RDM') exit
    read (buf,*) i, j, valor
    write (*,'(i3,1x,i3,1x,f24.16)') i, j, valor
  end do
  write (*,'(a)') '2-RDM'
  do while (bzptr%eof .eqv. .false.)
    length = 0
    call bzReadline (bzptr, buf, length)
    read (buf,*) i, j, k, l, valor
    write (*,'(i3,1x,i3,1x,i3,1x,i3,f24.16)') i, j, k, l, valor
  end do
  call bzCloseR (bzptr)

end program bzip
