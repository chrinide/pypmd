! read in the atomic coordinates from the molecular
! database or write them back to the database from the local
! matrices.
subroutine rwcoord(init, mol)

  use iso_fortran_env, only: uout=>output_unit
  use mod_prec, only: rp, ip
  use mod_datatm, only: wgatm
  use mod_wfn, only: ncent, charge, xyz
  use mod_param, only: verbose
  use mod_sym, only: sdim, ax, ay, az, atzmol
  use mod_memory, only: alloc
  implicit none

  logical, intent(in) :: init, mol

  integer(kind=ip) :: i
  real(kind=rp) :: cm(3), wt, twt 

  cm = 0.0_rp
  wt = 0.0_rp

  if (init) then
    if (mol) then
      sdim = ncent
      call alloc ('mod_sym', 'ax', ax, sdim)
      call alloc ('mod_sym', 'ay', ay, sdim)
      call alloc ('mod_sym', 'az', az, sdim)
      call alloc ('mod_sym', 'atzmol', atzmol, sdim)
      atzmol(:) = (int(charge(:)))
      ax(:) = xyz(:,1)
      ay(:) = xyz(:,2)
      az(:) = xyz(:,3)
      if (verbose) write (uout,600)
      ! compute the center of mass and charge
      do i = 1,sdim
        twt = wgatm(int(charge(i)))
        wt = wt + twt
        cm(1) = cm(1) + twt*xyz(i,1)
        cm(2) = cm(2) + twt*xyz(i,2)
        cm(3) = cm(3) + twt*xyz(i,3)
      end do
      cm = cm / wt
      if (verbose) write (uout,605) cm
      ! move the origin to the center of mass and store locally the
      ax(:) = xyz(:,1) - cm(1) 
      ay(:) = xyz(:,2) - cm(2) 
      az(:) = xyz(:,3) - cm(3) 
    end if
  end if

 600 format (/                                                        &
     1x, '++SYMRWCOORD (R):'/                                         &
     1x, 'Copy coordinates from the molecular database to the local ',&
     1x, 'symmetry arrays.')
 605 format (                                                         &
     1x, 'Center of mass:            ', 3f12.6)                        

end subroutine rwcoord
