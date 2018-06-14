subroutine optsparam(var,val)

  use iso_fortran_env, only: uout=>output_unit
  use mod_io, only: equal, faterr, ferror, string
  use mod_param, only: verbose, debug
  implicit none

  character(len=*), intent(in) :: var
  logical :: val

  if (equal(var,'verbose')) then
    verbose = val
    write (uout,'(1x,a)') string('# *** Verbose mode is enabled')
  else if (equal(var,'debug')) then
    debug = val
    verbose = val
    write (uout,'(1x,a)') string('# !! WARNING: DEBUG MODE ENABLED !!')
    write (uout,'(1x,a)') string('# !! Lot of info will be printed !!')
    write (uout,'(1x,a)') string('# !! WARNING: DEBUG MODE ENABLED !!')
    write (uout,*)
  else
    call ferror ('optsparam', 'unknown option', faterr)
  end if

end subroutine optsparam
