subroutine optssym(var,val)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: equal, faterr, ferror, string
  use mod_param, only: verbose
  use mod_sym, only: toldist
  implicit none

  character(len=*), intent(in) :: var
  real(kind=rp) :: val

  if (equal(var,'toldist')) then
    toldist = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable toldist changed to :'), toldist
    end if
  else
    call ferror ('optssym', 'unknown option '//string(var), faterr)
  end if

end subroutine optssym
