subroutine init_sym()

  implicit none
    
  ! default values  
  TOLdist = 3E-5_rp
  TOLeqvm = 2_rp*TOLdist
  TOLsng = 1E-5_rp
  TOLisint = 3E-5_rp
  TOLnull = 3E-6_rp
  TOLeigen = 3E-5_rp
  TOLdirty = .false.
  TOLtriplet = 0.1_rp
  ERRsng = 0_rp
  ERReigen = 0_rp
  inf_order = 16
  mol_linear = .false.
  mol_planar = .false.

end subroutine init_sym
