subroutine init_surf ()

  use mod_prec, only: rp, ip
  use mod_wfn, only: ncent
  use mod_surf, only: steeper, ntrial, rprimer, inuc, &
                      epsilon, epsiscp, nangleb, &
                      allocate_space_for_cp, rmaxsurf 
  implicit none

  steeper = 1_ip
  ntrial = 11_ip
  rprimer = 0.4_rp
  inuc = 0_ip
  epsilon = 1d-5 
  epsiscp = 0.08_rp
  nangleb(1) = 1
  nangleb(2) = 434
  nangleb(3) = 0
  nangleb(4) = 0
  rmaxsurf = 10.0_rp

  call allocate_space_for_cp (ncent)

end subroutine
