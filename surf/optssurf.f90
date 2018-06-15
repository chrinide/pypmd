subroutine optssurf(var,rval,ival)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_io, only: equal, faterr, ferror, string
  use mod_surf, only: steeper, ntrial, rprimer, epsilon, &
                      epsiscp, rmaxsurf
  use mod_param, only: verbose
  implicit none

  character(len=*), intent(in) :: var
  real(kind=rp), optional :: rval
  integer(kind=ip), optional :: ival

  if (equal(var,'steeper')) then
    steeper = abs(ival)
    if (verbose) then
      write (uout,'(1x,a,1x,i0)') string('# *** Variable steeper changed to :'), steeper
    end if
  else if (equal(var,'ntrial')) then
    ntrial = abs(ival)
    if (verbose) then
      write (uout,'(1x,a,1x,i0)') string('# *** Variable ntrial changed to :'), ntrial
    end if
  else if (equal(var,'rmaxsurf')) then
    rmaxsurf = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxsurf changed to :'), rmaxsurf
    end if
  else if (equal(var,'epsilon')) then
    epsilon = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsilon changed to :'), epsilon
    end if
  else if (equal(var,'epsiscp')) then
    epsiscp = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsiscp changed to :'), epsiscp
    end if
  else if (equal(var,'rprimer')) then
    rprimer = abs(rval)
    if (verbose) then
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rprimer changed to :'), rprimer
    end if
  else
    call ferror ('optssurf', 'unknown option', faterr)
  end if

end subroutine optssurf
