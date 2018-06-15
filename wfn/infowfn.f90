subroutine infowfn()

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: ip
  use mod_io, only: string, mline
  use mod_wfn, only: ncent, charge, atnam, xyz, nvirtual, noccupied, &
                     occupied, epsocc, epsortho
  implicit none

  integer(kind=ip) :: i

  write (uout,'(1x,a,1x,i0)') string('# Number of centers :'), ncent
  atnam = adjustl(atnam)
  do i = 1,ncent
    write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') '#', i, atnam(i)(1:2), charge(i), xyz(i,:)
  end do
  write (uout,'(1x,a,1x,e13.6)') string('# Occupied eps ='), epsocc
  write (uout,'(1x,a,1x,i0)') string('# Number of occupied orbitals ='), noccupied
  write (uout,*) string('#'), occupied(:)
  write (uout,'(1x,a,1x,i0)') string('# Number of virtual orbitals ='), nvirtual
  write (uout,'(1x,a,1x,e13.6)') string('# epsortho = '), epsortho

end subroutine infowfn
