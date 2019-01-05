module mod_basis
  
  use mod_prec, only: ip, rp
  implicit none
  public

  integer(kind=ip), parameter :: mgrp = 500_ip
  integer(kind=ip), parameter :: ngtoh = 21_ip
  integer(kind=ip), parameter, private :: maxtype = 56_ip

  ! nlm keeps the nlm values of x^n y^l z^m gaussian
  integer(kind=ip) :: nlm(maxtype,3)
  integer(kind=ip) :: nprims_
  integer(kind=ip) :: nmo_

  real(kind=rp), allocatable, dimension(:,:) :: coeff_
  real(kind=rp), allocatable, dimension(:) :: oexp_
  real(kind=rp), allocatable, dimension(:) :: occ_
  integer(kind=ip), allocatable, dimension(:) :: icen_
  integer(kind=ip), allocatable, dimension(:) :: ityp_

  integer(kind=ip), allocatable, dimension(:) :: npc_
  integer(kind=ip), allocatable, dimension(:) :: ngroup_
  integer(kind=ip), allocatable, dimension(:,:) :: nzexp_
  integer(kind=ip), allocatable, dimension(:,:,:) :: nuexp_
  real(kind=rp), allocatable, dimension(:,:) :: rcutte_

contains

  subroutine allocate_space_for_basis ()
  
    use mod_memory, only: alloc
    use mod_mole, only: ncent_
    implicit none

    call alloc ('mod_basis', 'coeff', coeff_, nmo_, nprims_)
    call alloc ('mod_basis', 'oexp', oexp_, nprims_)
    call alloc ('mod_basis', 'occ', occ_, nmo_)
    call alloc ('mod_basis', 'icen', icen_, nprims_)
    call alloc ('mod_basis', 'ityp', ityp_, nprims_)

    call alloc ('mod_basis', 'npc', npc_, ncent_)
    call alloc ('mod_basis', 'ngroup', ngroup_, ncent_)
    call alloc ('mod_basis', 'nzexp', nzexp_, ncent_, mgrp)
    call alloc ('mod_basis', 'nuexp', nuexp_, ncent_, ngtoh, mgrp)

    call alloc ('mod_basis', 'rcutte', rcutte_, ncent_, mgrp)
 
  end subroutine allocate_space_for_basis

  subroutine deallocate_space_for_basis ()

    use mod_memory, only: free
    implicit none

    call free ('mod_basis', 'coeff', coeff_)
    call free ('mod_basis', 'oexp', oexp_)
    call free ('mod_basis', 'occ', occ_)
    call free ('mod_basis', 'icen', icen_)
    call free ('mod_basis', 'ityp', ityp_)

    call free ('mod_basis', 'npc', npc_)
    call free ('mod_basis', 'ngroup', ngroup_)
    call free ('mod_basis', 'nzexp', nzexp_)
    call free ('mod_basis', 'nuexp', nuexp_)

    call free ('mod_basis', 'rcutte', rcutte_)

  end subroutine deallocate_space_for_basis


  subroutine init_basis()
                
    implicit none

    !.p's
    nlm(2,1)=1      !px
    nlm(2,2)=0      !px
    nlm(2,3)=0      !px
    nlm(3,1)=0      !py
    nlm(3,2)=1      !py
    nlm(3,3)=0      !py
    nlm(4,1)=0      !pz
    nlm(4,2)=0      !pz
    nlm(4,3)=1      !pz
    !.d's
    nlm(5,1)=2      !xx
    nlm(6,2)=2      !yy
    nlm(7,3)=2      !zz
    nlm(8,1)=1      !xy
    nlm(8,2)=1
    nlm(9,1)=1      !xz
    nlm(9,3)=1
    nlm(10,2)=1     !yz
    nlm(10,3)=1
    !.f's
    nlm(11,1)=3     !xxx
    nlm(12,2)=3     !yyy
    nlm(13,3)=3     !zzz    
    nlm(14,1)=2     !xxy
    nlm(14,2)=1
    nlm(15,1)=2     !xxz
    nlm(15,3)=1
    nlm(16,2)=2     !yyz
    nlm(16,3)=1 
    nlm(17,1)=1     !xyy
    nlm(17,2)=2 
    nlm(18,1)=1     !xzz
    nlm(18,3)=2     
    nlm(19,2)=1     !yzz
    nlm(19,3)=2
    nlm(20,1)=1     !xyz
    nlm(20,2)=1 
    nlm(20,3)=1 
    !.g's
    nlm(21,1)=4     !xxxx
    nlm(22,2)=4     !yyyy
    nlm(23,3)=4     !zzzz
    nlm(24,1)=3     !xxxy
    nlm(24,2)=1
    nlm(25,1)=3     !xxxz
    nlm(25,3)=1 
    nlm(26,1)=1     !xyyy
    nlm(26,2)=3 
    nlm(27,2)=3     !yyyz
    nlm(27,3)=1 
    nlm(28,1)=1     !xzzz 
    nlm(28,3)=3 
    nlm(29,2)=1     !yzzz
    nlm(29,3)=3
    nlm(30,1)=2     !xxyy
    nlm(30,2)=2 
    nlm(31,1)=2     !xxzz
    nlm(31,3)=2
    nlm(32,2)=2     !yyzz 
    nlm(32,3)=2 
    nlm(33,1)=2     !xxyz 
    nlm(33,2)=1
    nlm(33,3)=1
    nlm(34,1)=1     !xyyz
    nlm(34,2)=2 
    nlm(34,3)=1
    nlm(35,1)=1     !xyzz
    nlm(35,2)=1
    nlm(35,3)=2
    !.h's
    nlm(36,1)=0
    nlm(36,2)=0
    nlm(36,3)=5

    nlm(37,1)=0
    nlm(37,2)=1
    nlm(37,3)=4

    nlm(38,1)=0
    nlm(38,2)=2
    nlm(38,3)=3

    nlm(39,1)=0
    nlm(39,2)=3
    nlm(39,3)=2

    nlm(40,1)=0
    nlm(40,2)=4
    nlm(40,3)=1

    nlm(41,1)=0
    nlm(41,2)=5
    nlm(41,3)=0

    nlm(42,1)=1
    nlm(42,2)=0
    nlm(42,3)=4

    nlm(43,1)=1
    nlm(43,2)=1
    nlm(43,3)=3

    nlm(44,1)=1
    nlm(44,2)=2
    nlm(44,3)=2

    nlm(45,1)=1
    nlm(45,2)=3
    nlm(45,3)=1

    nlm(46,1)=1
    nlm(46,2)=4
    nlm(46,3)=0

    nlm(47,1)=2
    nlm(47,2)=0
    nlm(47,3)=3
  
    nlm(48,1)=2
    nlm(48,2)=1
    nlm(48,3)=2
  
    nlm(49,1)=2
    nlm(49,2)=2
    nlm(49,3)=1

    nlm(50,1)=2
    nlm(50,2)=3
    nlm(50,3)=0

    nlm(51,1)=3
    nlm(51,2)=0
    nlm(51,3)=2

    nlm(52,1)=3
    nlm(52,2)=1
    nlm(52,3)=1

    nlm(53,1)=3
    nlm(53,2)=2
    nlm(53,3)=0

    nlm(54,1)=4
    nlm(54,2)=0
    nlm(54,3)=1

    nlm(55,1)=4
    nlm(55,2)=1
    nlm(55,3)=0

    nlm(56,1)=5
    nlm(56,2)=0
    nlm(56,3)=0
  
  end subroutine init_basis

end module mod_basis
