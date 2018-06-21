subroutine optsgeom(var,val)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: equal, faterr, ferror, string
  use mod_param, only: verbose
  use mod_geom, only: covx
  implicit none

  character(len=*), intent(in) :: var
  real(kind=rp) :: val

  if (equal(var,'covx')) then
    covx = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable covx changed to :'), covx
    end if
  else
    call ferror ('optsgeom', 'unknown option '//string(var), faterr)
  end if

end subroutine optsgeom
