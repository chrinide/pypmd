subroutine init_geom()
                
  use mod_prec, only: rp
  use mod_geom, only: covx
  implicit none

  covx = 1.2_rp
  
end subroutine init_geom
