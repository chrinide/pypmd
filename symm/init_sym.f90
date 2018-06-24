subroutine init_sym()

  use mod_prec, only: rp, ip
  use mod_sym, only: toldist, toleqvm, tolsng, tolisint, &
                     tolnull, toleigen, toldirty, toltriplet, &
                     errsng, erreigen, inf_order, mol_linear, &
                     mol_planar
  implicit none
    
  ! default values  
  toldist = 3e-5_rp
  toleqvm = 2.0_rp*toldist
  tolsng = 1e-5_rp
  tolisint = 3e-5_rp
  tolnull = 3e-6_rp
  toleigen = 3e-5_rp
  toldirty = .false.
  toltriplet = 0.1_rp
  errsng = 0.0_rp
  erreigen = 0.0_rp
  inf_order = 16_ip
  mol_linear = .false.
  mol_planar = .false.


end subroutine init_sym
