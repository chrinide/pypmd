subroutine optswfn(var,val)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_wfn, only: cuttz, epsocc, rmaxatom, epsortho 
  use mod_io, only: equal, faterr, ferror, string
  use mod_param, only: verbose
  implicit none

  character(len=*), intent(in) :: var
  real(kind=rp) :: val

  if (equal(var,'cuttz')) then
    cuttz = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable cuttz changed to :'), cuttz
    end if
  else if (equal(var,'epsocc')) then
    epsocc = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsocc changed to :'), epsocc
    end if
  else if (equal(var,'epsortho')) then
    epsortho = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsortho changed to :'), epsortho
    end if
  else if (equal(var,'rmaxatom')) then
    rmaxatom = abs(val)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxatom changed to :'), rmaxatom
    end if
  else
    call ferror ('optswfn', 'unknown option '//string(var), faterr)
  end if

end subroutine optswfn
